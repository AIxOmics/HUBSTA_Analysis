"""Plot CP/SP celltype-region Cauchy median p-values for ASD and IQ."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.paths import WORK_DIR  # noqa: E402


matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PROJECT_NAME = "20260319_cellbin_xgboost_gsMap2d_v1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=WORK_DIR / PROJECT_NAME / "cauchy_combination",
    )
    parser.add_argument("--project-name", default=PROJECT_NAME)
    parser.add_argument("--traits", nargs="+", default=["ASD", "IQ"])
    return parser.parse_args()


def plot_trait(ax, path: Path, trait: str) -> None:
    df = pd.read_csv(path)
    cp = df[df["annotation"].str.contains("CP", na=False)].copy()
    cp["celltype"] = cp["annotation"].str.split("|").str[0]
    cp["neg_log10_p_median"] = -np.log10(cp["p_median"])
    cp = cp.sort_values("neg_log10_p_median", ascending=False)

    reject, _, _, _ = multipletests(cp["p_median"].values, alpha=0.05, method="fdr_bh")
    fdr_threshold = cp["p_median"].values[reject].max() if reject.any() else None

    colors = plt.cm.tab10(np.linspace(0, 1, len(cp)))
    ax.bar(cp["celltype"], cp["neg_log10_p_median"], color=colors, edgecolor="white", linewidth=0.5)
    if fdr_threshold is not None:
        ax.axhline(-np.log10(fdr_threshold), color="red", linestyle="--", linewidth=1, label="FDR = 0.05")
    else:
        ax.text(
            0.5,
            0.95,
            "No FDR-significant cell types",
            transform=ax.transAxes,
            ha="center",
            va="top",
            color="red",
            fontsize=9,
        )
    ax.set_title(f"{trait} CP/SP cell types", fontsize=13, fontweight="bold")
    ax.set_xlabel("Cell type")
    ax.set_ylabel("-log10(median p-value)")
    ax.tick_params(axis="x", rotation=45)
    ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main() -> None:
    args = parse_args()
    fig, axes = plt.subplots(1, len(args.traits), figsize=(7 * len(args.traits), 5))
    if len(args.traits) == 1:
        axes = [axes]
    for ax, trait in zip(axes, args.traits):
        plot_trait(
            ax,
            args.work_dir / f"{args.project_name}_{trait}.celltype_region.cauchy.csv",
            trait,
        )
    plt.tight_layout()
    out = args.work_dir / f"{args.project_name}_CP_SP_barplot.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
