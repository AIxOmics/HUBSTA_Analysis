"""Plot heatmaps from gsMap Cauchy combination results."""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.paths import WORK_DIR  # noqa: E402


matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import seaborn as sns  # noqa: E402

warnings.filterwarnings("ignore")

ANNOTATION_TYPES = [
    "celltype",
    "major_celltype",
    "regions",
    "subregion_major_celltype_pair",
]
MEASURES = {
    "Cauchy": "neg_log10_cauchy",
    "Median": "neg_log10_median",
    "Q95": "q95",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=WORK_DIR / "fetal_brain_cell_bin_gsmap2d/cauchy_combination",
    )
    parser.add_argument("--output-dir", type=Path, default=None)
    return parser.parse_args()


def load_annotation(data_dir: Path, annotation_type: str) -> pd.DataFrame:
    files = sorted(data_dir.glob(f"*.{annotation_type}.cauchy.csv"))
    if not files:
        raise FileNotFoundError(f"No files found for annotation type: {annotation_type}")
    df = pd.concat((pd.read_csv(path) for path in files), ignore_index=True)
    df["neg_log10_cauchy"] = -np.log10(df["p_cauchy"].clip(lower=1e-300))
    df["neg_log10_median"] = -np.log10(df["p_median"].clip(lower=1e-300))
    df["q95"] = df["top_95_quantile"]
    return df


def draw_clustermap(mat: pd.DataFrame, title: str, out_path: Path, vmin=0, vmax=None) -> None:
    if mat.empty:
        return
    if vmax is None:
        vmax = np.nanpercentile(mat.values, 98)
    n_rows, n_cols = mat.shape
    figsize = (max(8, n_cols * 0.35 + 3), max(5, n_rows * 0.45 + 3))
    grid = sns.clustermap(
        mat.fillna(0),
        method="average",
        metric="euclidean",
        cmap="RdBu_r",
        vmin=vmin,
        vmax=vmax,
        figsize=figsize,
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"label": "-log10(p)", "shrink": 0.6},
        dendrogram_ratio=(0.12, 0.12),
    )
    grid.ax_heatmap.set_title(title, fontsize=12, fontweight="bold", pad=10)
    grid.ax_heatmap.set_xlabel("Annotation")
    grid.ax_heatmap.set_ylabel("Trait")
    plt.setp(grid.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=7)
    plt.setp(grid.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8)
    grid.figure.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(grid.figure)
    print(f"Saved: {out_path}")


def normalize_per_trait(mat: pd.DataFrame) -> pd.DataFrame:
    row_min = mat.min(axis=1)
    row_range = (mat.max(axis=1) - row_min).replace(0, 1)
    return mat.sub(row_min, axis=0).div(row_range, axis=0)


def draw_per_trait(df_trait: pd.DataFrame, trait: str, ann_type: str, out_path: Path) -> None:
    rows = df_trait.set_index("annotation")[
        ["neg_log10_cauchy", "neg_log10_median", "q95"]
    ].rename(
        columns={
            "neg_log10_cauchy": "Cauchy\n(-log10p)",
            "neg_log10_median": "Median\n(-log10p)",
            "q95": "Q95\n(-log10p)",
        }
    )
    if rows.shape[0] < 2:
        return
    grid = sns.clustermap(
        rows.fillna(0),
        method="average",
        metric="euclidean",
        col_cluster=False,
        row_cluster=True,
        cmap="RdBu_r",
        vmin=0,
        vmax=np.nanpercentile(rows.values, 98),
        figsize=(6, max(5, rows.shape[0] * 0.28 + 2)),
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"label": "-log10(p)", "shrink": 0.5},
        dendrogram_ratio=(0.15, 0.02),
    )
    grid.ax_heatmap.set_title(f"{trait} | {ann_type}", fontsize=11, fontweight="bold", pad=8)
    grid.ax_heatmap.set_xlabel("Measure")
    grid.ax_heatmap.set_ylabel("Annotation")
    grid.figure.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(grid.figure)
    print(f"Saved: {out_path}")


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir or args.data_dir / "heatmaps"
    per_trait_dir = output_dir / "per_trait"
    per_trait_dir.mkdir(parents=True, exist_ok=True)

    for ann_type in ANNOTATION_TYPES:
        df = load_annotation(args.data_dir, ann_type)
        for measure_label, measure_col in MEASURES.items():
            mat = df.pivot_table(
                index="trait",
                columns="annotation",
                values=measure_col,
                aggfunc="first",
            )
            draw_clustermap(
                mat,
                f"{measure_label} (-log10p) | {ann_type}",
                output_dir / f"all_traits_{ann_type}_{measure_label}.pdf",
            )
            draw_clustermap(
                normalize_per_trait(mat),
                f"{measure_label} trait-normalized | {ann_type}",
                output_dir / f"all_traits_{ann_type}_{measure_label}_traitnorm.pdf",
                vmin=0,
                vmax=1,
            )

        for trait in sorted(df["trait"].unique()):
            safe_trait = trait.replace(" ", "_")
            draw_per_trait(
                df[df["trait"] == trait].copy(),
                trait,
                ann_type,
                per_trait_dir / f"{safe_trait}_{ann_type}.pdf",
            )


if __name__ == "__main__":
    main()
