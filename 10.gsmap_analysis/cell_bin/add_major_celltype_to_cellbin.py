"""Add major cell type annotations to cell-bin h5ad files."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.celltypes import add_major_celltype_columns  # noqa: E402
from shared.paths import CELLBIN_H5AD_DIR, CELLBIN_MAJOR_CELLTYPE_DIR  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=CELLBIN_H5AD_DIR)
    parser.add_argument("--output-dir", type=Path, default=CELLBIN_MAJOR_CELLTYPE_DIR)
    parser.add_argument("--pattern", default="*.h5ad")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    h5ad_files = sorted(args.input_dir.glob(args.pattern))
    print(f"Found {len(h5ad_files)} h5ad files in {args.input_dir}")
    for fpath in h5ad_files:
        adata = ad.read_h5ad(fpath)
        adata.obs = add_major_celltype_columns(adata.obs)
        out_path = args.output_dir / fpath.name
        adata.write_h5ad(out_path, compression="gzip")
        print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()

