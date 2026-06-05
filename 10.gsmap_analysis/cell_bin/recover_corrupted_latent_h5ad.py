"""Recover a gsMap latent h5ad whose X dataset is corrupted.

The known historical case is sample ``39_A03589A1D4_WT202403310046`` in
``fetal_brain_cell_bin_gsmap2d/find_latent_representations``. The script reads
obs/obsm/var from the damaged file with h5py, reconstructs X from the original
cell-bin h5ad, and overwrites the target file.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import h5py
import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.celltypes import add_major_celltype_columns  # noqa: E402
from shared.paths import CELLBIN_MAJOR_CELLTYPE_DIR, WORK_DIR  # noqa: E402


SAMPLE = "39_A03589A1D4_WT202403310046"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target-h5ad",
        type=Path,
        default=(
            WORK_DIR
            / "fetal_brain_cell_bin_gsmap2d/find_latent_representations"
            / f"{SAMPLE}_add_latent.h5ad"
        ),
    )
    parser.add_argument(
        "--original-h5ad",
        type=Path,
        default=CELLBIN_MAJOR_CELLTYPE_DIR / f"{SAMPLE}.celltype.region.h5ad",
    )
    return parser.parse_args()


def read_string_vector(dataset) -> list[str]:
    return [x.decode() if isinstance(x, bytes) else str(x) for x in dataset[:]]


def read_obs_from_h5(path: Path, index: list[str]) -> pd.DataFrame:
    obs_dict = {}
    with h5py.File(path, "r") as handle:
        obs_group = handle["obs"]
        for key in obs_group.keys():
            if key == "_index":
                continue
            dataset = obs_group[key]
            if "categories" in dataset:
                categories = read_string_vector(dataset["categories"])
                obs_dict[key] = pd.Categorical.from_codes(
                    dataset["codes"][:],
                    categories=categories,
                )
            else:
                values = dataset[:]
                if values.dtype.kind in ("S", "O"):
                    values = [v.decode() if isinstance(v, bytes) else v for v in values]
                obs_dict[key] = values
    return pd.DataFrame(obs_dict, index=index)


def main() -> None:
    args = parse_args()
    with h5py.File(args.target_h5ad, "r") as handle:
        target_genes = read_string_vector(handle["var"]["HUMAN_GENE_SYM"])
        cell_index = read_string_vector(handle["obs"]["_index"])
        obsm = {key: handle["obsm"][key][:] for key in handle["obsm"].keys()}

    obs = add_major_celltype_columns(read_obs_from_h5(args.target_h5ad, cell_index))
    original = ad.read_h5ad(args.original_h5ad)
    missing = [gene for gene in target_genes if gene not in original.var_names]
    if missing:
        raise KeyError(f"{len(missing)} target genes missing from original h5ad: {missing[:5]}")

    original_subset = original[cell_index, target_genes].copy()
    recovered = ad.AnnData(
        X=original_subset.X,
        obs=obs,
        var=pd.DataFrame(index=pd.Index(target_genes, name="HUMAN_GENE_SYM")),
        obsm=obsm,
    )
    recovered.write_h5ad(args.target_h5ad, compression="gzip")
    print(f"Recovered: {args.target_h5ad}")


if __name__ == "__main__":
    main()
