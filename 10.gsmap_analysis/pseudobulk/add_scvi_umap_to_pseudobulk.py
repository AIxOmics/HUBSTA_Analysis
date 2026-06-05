"""Train scVI on pseudobulk profiles and add X_scVI / X_umap in-place."""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import scanpy as sc
import scipy.sparse as sp

import torch._library.fake_impl as _fim

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.paths import PSEUDOBULK_SUBREGION_CELLTYPE_H5AD  # noqa: E402


matplotlib.use("Agg")
warnings.filterwarnings("ignore")

_ORIGINAL_REGISTER = _fim.FakeImplHolder.register


def _safe_register(self, func, source, lib=None, allow_override=False):
    try:
        return _ORIGINAL_REGISTER(
            self,
            func,
            source,
            lib=lib,
            allow_override=allow_override,
        )
    except RuntimeError:
        return None


_fim.FakeImplHolder.register = _safe_register

import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import scvi  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=PSEUDOBULK_SUBREGION_CELLTYPE_H5AD)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--n-top-genes", type=int, default=4000)
    parser.add_argument("--n-latent", type=int, default=30)
    parser.add_argument("--max-epochs", type=int, default=200)
    return parser.parse_args()


def make_palette(labels, cmap_name="tab20"):
    uniq = sorted(labels.unique())
    cmap = plt.get_cmap(cmap_name, len(uniq))
    return {value: cmap(i) for i, value in enumerate(uniq)}


def plot_umap(adata: ad.AnnData, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    umap = adata.obsm["X_umap"]
    obs = adata.obs

    palette_sample = make_palette(obs["sample_id"])
    fig, ax = plt.subplots(figsize=(10, 8))
    for sample_id, color in palette_sample.items():
        mask = obs["sample_id"] == sample_id
        ax.scatter(umap[mask, 0], umap[mask, 1], c=[color], s=6, alpha=0.7, linewidths=0)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("Pseudobulk UMAP by sample_id")
    ax.legend(
        handles=[mpatches.Patch(color=c, label=l) for l, c in palette_sample.items()],
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=6,
        ncol=2,
        frameon=False,
    )
    plt.tight_layout()
    fig.savefig(output_dir / "umap_by_sample.png", dpi=150)
    plt.close(fig)

    major_order = sorted(obs["major_celltype"].unique())
    cmap = plt.get_cmap("tab20", len(major_order))
    palette_major = {ct: cmap(i) for i, ct in enumerate(major_order)}

    fig, ax = plt.subplots(figsize=(10, 8))
    for celltype, color in palette_major.items():
        mask = obs["major_celltype"] == celltype
        ax.scatter(umap[mask, 0], umap[mask, 1], c=[color], s=8, alpha=0.8, linewidths=0)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("Pseudobulk UMAP by major_celltype")
    ax.legend(
        handles=[mpatches.Patch(color=c, label=l) for l, c in palette_major.items()],
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=8,
        frameon=False,
    )
    plt.tight_layout()
    fig.savefig(output_dir / "umap_by_major_celltype.png", dpi=150)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    scvi.settings.seed = 42
    scvi.settings.verbosity = 20

    adata = ad.read_h5ad(args.h5ad)
    x = adata.X
    if not isinstance(x, sp.csr_matrix):
        x = sp.csr_matrix(x)
    adata.layers["counts"] = x.astype(np.float32)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=args.n_top_genes,
        batch_key="sample_id",
        flavor="seurat_v3",
        layer="counts",
    )

    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="sample_id")
    model = scvi.model.SCVI(
        adata,
        n_latent=args.n_latent,
        n_layers=2,
        n_hidden=128,
        gene_likelihood="nb",
        dispersion="gene-batch",
    )
    model.train(
        max_epochs=args.max_epochs,
        early_stopping=True,
        early_stopping_patience=20,
        train_size=0.9,
        batch_size=256,
        plan_kwargs={"lr": 1e-3},
    )

    adata.obsm["X_scVI"] = model.get_latent_representation()
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15, n_pcs=None)
    sc.tl.umap(adata, min_dist=0.3)
    adata.write_h5ad(args.h5ad, compression="gzip")
    plot_umap(adata, args.output_dir)
    print(f"Updated: {args.h5ad}")


if __name__ == "__main__":
    main()
