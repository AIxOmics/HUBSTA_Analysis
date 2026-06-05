"""Add major cell type annotations to gsMap latent h5ad files in-place."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.celltypes import add_major_celltype_columns  # noqa: E402
from shared.paths import WORK_DIR  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target-dir",
        type=Path,
        default=WORK_DIR / "fetal_brain_cell_bin_gsmap2d/find_latent_representations",
    )
    parser.add_argument("--pattern", default="*.h5ad")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    h5ad_files = sorted(args.target_dir.glob(args.pattern))
    print(f"Found {len(h5ad_files)} h5ad files in {args.target_dir}")
    for fpath in h5ad_files:
        adata = ad.read_h5ad(fpath)
        adata.obs = add_major_celltype_columns(adata.obs)
        adata.write_h5ad(fpath, compression="gzip")
        print(f"Updated in-place: {fpath}")


if __name__ == "__main__":
    main()

