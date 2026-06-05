"""Collapse sample-level pseudobulk rows into weighted mean region-celltype rows."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.celltypes import add_major_celltype_columns  # noqa: E402
from shared.paths import (  # noqa: E402
    PSEUDOBULK_SUBREGION_CELLTYPE_H5AD,
    PSEUDOBULK_WEIGHTED_MEAN_H5AD,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, default=PSEUDOBULK_SUBREGION_CELLTYPE_H5AD)
    parser.add_argument("--output-h5ad", type=Path, default=PSEUDOBULK_WEIGHTED_MEAN_H5AD)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    pb = ad.read_h5ad(args.input_h5ad)
    obs = pb.obs.copy()
    x = pb.X if isinstance(pb.X, sp.csr_matrix) else sp.csr_matrix(pb.X)
    pair_key = obs["celltype"].astype(str) + "||" + obs["regions"].astype(str)
    n_cells_arr = obs["n_cells"].values.astype(np.float64)

    wm_obs_rows = []
    wm_x_rows = []
    for pair, group_idx in pair_key.groupby(pair_key, observed=True).groups.items():
        pos = obs.index.get_indexer(group_idx)
        weights = n_cells_arr[pos]
        x_sum = np.asarray(x[pos].sum(axis=0)).ravel()
        total_cells = weights.sum()
        celltype, subregion = pair.split("||", 1)
        group_obs = obs.iloc[pos]
        wm_obs_rows.append(
            {
                "celltype": celltype,
                "regions": subregion,
                "region": group_obs["region"].mode()[0],
                "pair": f"{celltype}|{subregion}",
                "n_samples": len(pos),
                "total_n_cells": int(total_cells),
                "mean_dnbCount": (group_obs["mean_dnbCount"] * weights).sum() / total_cells,
                "mean_area": (group_obs["mean_area"] * weights).sum() / total_cells,
                "mean_n_counts": (group_obs["mean_n_counts"] * weights).sum() / total_cells,
            }
        )
        wm_x_rows.append(x_sum / total_cells)

    wm_obs = add_major_celltype_columns(pd.DataFrame(wm_obs_rows))
    wm_obs.index = (
        wm_obs["celltype"].str.replace(" ", "_", regex=False)
        + "__"
        + wm_obs["regions"].str.replace("/", "-", regex=False).str.replace(" ", "_", regex=False)
    )
    adata = ad.AnnData(
        X=sp.csr_matrix(np.vstack(wm_x_rows).astype(np.float32)),
        obs=wm_obs,
        var=pb.var.copy(),
    )
    adata.write_h5ad(args.output_h5ad, compression="gzip")
    print(f"Saved: {args.output_h5ad}")
    print(f"Final shape: {adata.shape}")


if __name__ == "__main__":
    main()
