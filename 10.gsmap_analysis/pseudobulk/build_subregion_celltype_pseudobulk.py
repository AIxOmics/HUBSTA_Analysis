"""Build sample x celltype x subregion pseudobulk profiles.

Input cell-level files are expected to contain ``celltype``, ``regions``, and
``region`` columns in ``obs``. The output is the single h5ad consumed by
``run_pseudobulk_gsmap.py`` after scVI/UMAP embedding is added.
"""

from __future__ import annotations

import argparse
import gc
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.celltypes import add_major_celltype_columns  # noqa: E402
from shared.paths import CELLBIN_H5AD_DIR, PSEUDOBULK_SUBREGION_CELLTYPE_H5AD  # noqa: E402


OBS_COLS_TO_AGGREGATE = ["dnbCount", "area", "n_counts", "celltype_main_frac"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=CELLBIN_H5AD_DIR)
    parser.add_argument(
        "--output-h5ad",
        type=Path,
        default=PSEUDOBULK_SUBREGION_CELLTYPE_H5AD,
    )
    parser.add_argument("--pattern", default="*.celltype.region.h5ad")
    parser.add_argument("--min-cells", type=int, default=1000)
    return parser.parse_args()


def common_gene_index(h5ad_files: list[Path]) -> list[str]:
    print("[1/3] Collecting common genes ...")
    var_sets = []
    for fpath in h5ad_files:
        adata = sc.read_h5ad(fpath, backed="r")
        var_sets.append(set(adata.var.index.astype(str)))
        del adata
    common_genes = sorted(var_sets[0].intersection(*var_sets[1:]))
    print(f"  Common genes across {len(h5ad_files)} samples: {len(common_genes)}")
    return common_genes


def build_pseudobulk_rows(
    h5ad_files: list[Path],
    common_genes: list[str],
) -> tuple[pd.DataFrame, sp.csr_matrix]:
    print("[2/3] Aggregating per sample ...")
    all_obs_rows: list[dict] = []
    all_x_rows: list[sp.csr_matrix] = []

    for fpath in h5ad_files:
        sample_id = fpath.name.split("_", 1)[0]
        sample_name = fpath.stem
        print(f"  {sample_name}", flush=True)

        adata = sc.read_h5ad(fpath)
        local_gene2idx = {g: i for i, g in enumerate(adata.var.index.astype(str))}
        local_col_idx = np.array([local_gene2idx[g] for g in common_genes], dtype=np.int32)

        x = adata.X
        if not sp.issparse(x):
            x = sp.csr_matrix(x)
        elif not isinstance(x, sp.csr_matrix):
            x = x.tocsr()
        x_common = x[:, local_col_idx]

        obs = adata.obs.copy()
        obs["_pair"] = obs["celltype"].astype(str) + "||" + obs["regions"].astype(str)
        obs["_pos"] = np.arange(len(obs))

        for pair, group in obs.groupby("_pair", observed=True):
            row_pos = group["_pos"].values
            n_cells = len(row_pos)
            x_sum = np.asarray(x_common[row_pos].sum(axis=0)).ravel()
            celltype, subregion = pair.split("||", 1)
            group_num = group[OBS_COLS_TO_AGGREGATE].apply(pd.to_numeric, errors="coerce")

            row = {
                "sample_id": sample_id,
                "sample_name": sample_name,
                "celltype": celltype,
                "regions": subregion,
                "region": group["region"].mode()[0],
                "n_cells": n_cells,
                "pair": f"{celltype}|{subregion}",
            }
            for col in OBS_COLS_TO_AGGREGATE:
                row[f"mean_{col}"] = group_num[col].mean()

            all_obs_rows.append(row)
            all_x_rows.append(sp.csr_matrix(x_sum))

        del adata, x, x_common, obs
        gc.collect()

    return pd.DataFrame(all_obs_rows), sp.vstack(all_x_rows, format="csr")


def main() -> None:
    args = parse_args()
    h5ad_files = sorted(args.input_dir.glob(args.pattern))
    if not h5ad_files:
        raise FileNotFoundError(f"No h5ad files matched {args.input_dir / args.pattern}")
    print(f"Found {len(h5ad_files)} input h5ad files.")

    common_genes = common_gene_index(h5ad_files)
    final_obs, final_x = build_pseudobulk_rows(h5ad_files, common_genes)

    print("[3/3] Filtering and writing output ...")
    keep = final_obs["n_cells"] >= args.min_cells
    print(f"  Total rows: {len(final_obs)}")
    print(f"  Kept rows with >= {args.min_cells} cells: {keep.sum()}")

    final_obs = final_obs.loc[keep].reset_index(drop=True)
    final_x = final_x[keep.values]
    final_obs = add_major_celltype_columns(final_obs)
    final_obs.index = (
        final_obs["sample_id"]
        + "__"
        + final_obs["celltype"].str.replace(" ", "_", regex=False)
        + "__"
        + final_obs["regions"]
        .str.replace("/", "-", regex=False)
        .str.replace(" ", "_", regex=False)
    )
    if final_obs.index.duplicated().any():
        final_obs.index = final_obs.index + "__" + (
            final_obs.groupby(level=0).cumcount().astype(str)
        )

    adata = ad.AnnData(
        X=final_x,
        obs=final_obs,
        var=pd.DataFrame(index=pd.Index(common_genes, name="gene_id")),
    )

    args.output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output_h5ad, compression="gzip")
    print(f"Saved: {args.output_h5ad}")
    print(f"Final shape: {adata.shape}")


if __name__ == "__main__":
    main()

