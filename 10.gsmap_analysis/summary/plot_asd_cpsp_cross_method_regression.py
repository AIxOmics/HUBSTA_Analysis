"""Compare ASD CP/SP gsMap signal with external pseudobulk LDSC signal."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.paths import WORK_DIR  # noqa: E402


matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import seaborn as sns  # noqa: E402

DEFAULT_CAUCHY = (
    WORK_DIR
    / "fetal_brain_cell_bin_gsmap2d/cauchy_combination/"
    "fetal_brain_cell_bin_gsmap2d_ASD.subregion_major_celltype_pair.cauchy.csv"
)
DEFAULT_PSEUDOBULK = Path(
    "/storage/yangjianLab/songliyang/SpatialData/gsMap3D_analysis/"
    "human_fetalbrain_pseudobulk/fetal_pesudo/spatial_ldsc/fetal_pesudo_ASD.csv"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cauchy-csv", type=Path, default=DEFAULT_CAUCHY)
    parser.add_argument("--pseudobulk-csv", type=Path, default=DEFAULT_PSEUDOBULK)
    parser.add_argument(
        "--measure",
        choices=["cauchy", "median"],
        default="cauchy",
        help="gsMap p-value column to compare.",
    )
    parser.add_argument("--output-prefix", default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    p_col = "p_cauchy" if args.measure == "cauchy" else "p_median"
    x_col = f"neg_log10_{args.measure}"

    cauchy = pd.read_csv(args.cauchy_csv)
    cauchy = cauchy[cauchy["annotation"].str.startswith("CP/SP|")].copy()
    cauchy["cell_type"] = cauchy["annotation"].str.split("|").str[1]
    cauchy[x_col] = -np.log10(cauchy[p_col])

    pseudo = pd.read_csv(args.pseudobulk_csv).rename(columns={"spot": "cell_type"})
    merged = cauchy[["cell_type", x_col]].merge(
        pseudo[["cell_type", "neg_log10_p"]],
        on="cell_type",
    )
    r, pval = stats.pearsonr(merged[x_col], merged["neg_log10_p"])

    fig, ax = plt.subplots(figsize=(5, 4.5))
    sns.regplot(
        data=merged,
        x=x_col,
        y="neg_log10_p",
        ax=ax,
        scatter_kws={"s": 60, "color": "#4C72B0", "edgecolors": "white", "linewidths": 0.5},
        line_kws={"color": "#C44E52", "linewidth": 1.5},
        ci=95,
    )
    for _, row in merged.iterrows():
        ax.text(row[x_col], row["neg_log10_p"], row["cell_type"], fontsize=7, ha="left", va="bottom")
    ax.text(
        0.05,
        0.95,
        f"Pearson r = {r:.3f}\np = {pval:.2e}",
        transform=ax.transAxes,
        fontsize=9,
        va="top",
        ha="left",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#cccccc", alpha=0.8),
    )
    ax.set_xlabel(f"-log10({args.measure} p) [gsMap CP/SP]")
    ax.set_ylabel("-log10(p) [Pseudobulk LDSC]")
    ax.set_title("ASD CP/SP cell types: gsMap vs pseudobulk LDSC")
    sns.despine()
    plt.tight_layout()

    prefix = args.output_prefix or f"ASD_CPSP_gsmap2d_{args.measure}_vs_pseudobulk_regression"
    out_pdf = Path(f"{prefix}.pdf")
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf.with_suffix(".png"), dpi=200, bbox_inches="tight")
    print(f"Saved: {out_pdf}")
    print(f"Pearson r = {r:.4f}, p = {pval:.4g}")


if __name__ == "__main__":
    main()
