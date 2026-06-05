"""Plot major cell type composition pies from pseudobulk h5ad obs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import scanpy as sc

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.paths import PSEUDOBULK_SUBREGION_CELLTYPE_H5AD  # noqa: E402


matplotlib.use("Agg")
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

CELLTYPES = [
    "Glu",
    "Gaba",
    "Dopaminergic",
    "NPC_RG",
    "NPC_precursor",
    "NPC_APC",
    "NPC_GCP",
    "NPC_Cajal",
    "NPC_Crest",
    "NPC_IPC",
    "OPC",
    "VLMC",
    "Immune",
    "Choroid",
    "unknown",
]
PALETTE = [
    "#E64B35",
    "#4DBBD5",
    "#F39B7F",
    "#3C5488",
    "#7E6148",
    "#B09C85",
    "#91D1C2",
    "#8491B4",
    "#DC0000",
    "#00A087",
    "#FFDC91",
    "#74C476",
    "#FD8D3C",
    "#9ECAE1",
    "#AAAAAA",
]
COLOR_MAP = dict(zip(CELLTYPES, PALETTE))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=PSEUDOBULK_SUBREGION_CELLTYPE_H5AD)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--focus-region", default="CP/SP")
    return parser.parse_args()


def region_counts(obs, region=None):
    if region is not None:
        obs = obs[obs["regions"] == region]
    counts = obs.groupby("major_celltype", observed=True)["n_cells"].sum()
    return counts[counts > 0].sort_values(ascending=False)


def plot_focus_region(obs, region: str, output_dir: Path) -> None:
    counts = region_counts(obs, region)
    fig, ax = plt.subplots(figsize=(7, 7))
    colors = [COLOR_MAP.get(ct, "#CCCCCC") for ct in counts.index]
    _, _, autotexts = ax.pie(
        counts.values,
        labels=None,
        colors=colors,
        autopct=lambda pct: f"{pct:.1f}%" if pct >= 2 else "",
        pctdistance=0.78,
        startangle=90,
        wedgeprops=dict(linewidth=0.8, edgecolor="white"),
        textprops=dict(fontsize=9),
    )
    for autotext in autotexts:
        autotext.set_fontsize(8)
        autotext.set_color("white")
        autotext.set_fontweight("bold")
    ax.legend(
        handles=[
            mpatches.Patch(color=COLOR_MAP.get(ct, "#CCCCCC"), label=f"{ct} ({v / counts.sum() * 100:.1f}%)")
            for ct, v in counts.items()
        ],
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=9,
        frameon=False,
        title="Cell Type",
    )
    safe_region = region.replace("/", "")
    ax.set_title(f"Cell Type Distribution - {region}", fontsize=14, fontweight="bold", pad=20)
    plt.tight_layout()
    fig.savefig(output_dir / f"celltype_pie_{safe_region}.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / f"celltype_pie_{safe_region}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_all_subregions(obs, output_dir: Path) -> None:
    subregions = [r for r in obs["regions"].unique() if r != "non-tissue"]
    totals = {r: obs.loc[obs["regions"] == r, "n_cells"].sum() for r in subregions}
    subregions = sorted(subregions, key=lambda r: totals[r], reverse=True)
    ncols = 4
    nrows = int(np.ceil(len(subregions) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.8, nrows * 3.8))
    axes = axes.flatten()

    for ax, region in zip(axes, subregions):
        counts = region_counts(obs, region)
        if counts.empty:
            ax.axis("off")
            continue
        total = counts.sum()
        main = counts[counts / total >= 0.02].copy()
        other = counts[counts / total < 0.02].sum()
        if other > 0:
            main.loc["Others"] = other
        colors = [COLOR_MAP.get(ct, "#CCCCCC") for ct in main.index]
        _, _, autotexts = ax.pie(
            main.values,
            labels=None,
            colors=colors,
            autopct=lambda pct: f"{pct:.0f}%" if pct >= 4 else "",
            pctdistance=0.75,
            startangle=90,
            wedgeprops=dict(linewidth=0.6, edgecolor="white"),
        )
        for autotext in autotexts:
            autotext.set_fontsize(7)
            autotext.set_color("white")
            autotext.set_fontweight("bold")
        ax.set_title(region, fontsize=10, fontweight="bold", pad=6)

    for ax in axes[len(subregions) :]:
        ax.axis("off")

    fig.legend(
        handles=[mpatches.Patch(color=COLOR_MAP.get(ct, "#CCCCCC"), label=ct) for ct in CELLTYPES]
        + [mpatches.Patch(color="#CCCCCC", label="Others")],
        loc="lower center",
        ncol=8,
        fontsize=8,
        frameon=False,
        bbox_to_anchor=(0.5, -0.02),
        title="Major Cell Type",
    )
    fig.suptitle("Major Cell Type Proportion Across Brain Subregions", fontsize=15, fontweight="bold", y=1.01)
    plt.tight_layout()
    fig.savefig(output_dir / "celltype_pie_subregions.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / "celltype_pie_subregions.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    obs = sc.read_h5ad(args.h5ad).obs.copy()
    plot_focus_region(obs, args.focus_region, args.output_dir)
    plot_all_subregions(obs, args.output_dir)
    print(f"Saved cell type pie plots to {args.output_dir}")


if __name__ == "__main__":
    main()
