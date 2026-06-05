"""Plot ASD pseudobulk LDSC enrichment with BH FDR annotation."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


matplotlib.use("Agg")
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

DEFAULT_INPUT = Path(
    "/storage/yangjianLab/songliyang/SpatialData/gsMap3D_analysis/"
    "human_fetalbrain_pseudobulk/fetal_pesudo/spatial_ldsc/fetal_pesudo_ASD.csv"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-prefix", default="ASD_pseudobulk_neglogp_barplot_fdr05")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input_csv).dropna(subset=["neg_log10_p", "spot", "p"])
    df = df.sort_values("neg_log10_p", ascending=False).reset_index(drop=True)
    reject, pvals_fdr, _, _ = multipletests(df["p"], alpha=0.05, method="fdr_bh")
    df["fdr"] = pvals_fdr
    df["sig_fdr"] = reject
    threshold = -np.log10(df.loc[reject, "p"].max()) if reject.any() else None

    colors = ["#d62728" if sig else "#aec7e8" for sig in df["sig_fdr"]]
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(df["spot"], df["neg_log10_p"], color=colors, edgecolor="white", linewidth=0.5)
    if threshold is not None:
        ax.axhline(threshold, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Cell Type")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title("ASD GWAS Enrichment by Cell Type")
    plt.xticks(rotation=45, ha="right", fontsize=10)
    handles = [
        mpatches.Patch(color="#d62728", label="FDR < 0.05"),
        mpatches.Patch(color="#aec7e8", label="Not significant"),
    ]
    if threshold is not None:
        handles.append(
            plt.Line2D([0], [0], color="black", linestyle="--", linewidth=1, label="FDR 0.05 threshold")
        )
    ax.legend(handles=handles, fontsize=9)
    plt.tight_layout()
    out_png = Path(f"{args.output_prefix}.png")
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".pdf"), dpi=300, bbox_inches="tight")
    print(f"Saved: {out_png}")
    print(df[["spot", "p", "fdr", "neg_log10_p", "sig_fdr"]].to_string(index=False))


if __name__ == "__main__":
    main()
