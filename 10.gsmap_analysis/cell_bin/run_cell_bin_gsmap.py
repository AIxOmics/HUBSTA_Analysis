"""Run gsMap QuickMode on one cell-bin sample with multiple GWAS traits.

This reproduces the historical ``fetal_brain_cell_bin_gsmap2d`` run. The
current result metadata shows that this project used sample
``51_B03605C2E5_WT202406020126`` and started from an existing latent h5ad.
"""

from __future__ import annotations

import os
from collections import OrderedDict
from pathlib import Path

os.environ.setdefault("XLA_PYTHON_CLIENT_MEM_FRACTION", "0.4")

from gsMap.config import DatasetType, QuickModeConfig
from gsMap.pipeline.quick_mode import run_quick_mode

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.paths import (  # noqa: E402
    BIO_REPLICATE_SUMSTATS_CONFIG,
    GSMAP_RESOURCE_DIR,
    SNP_GENE_WEIGHT_ADATA,
    WORK_DIR,
    ensure_work_dirs,
)
from shared.script_utils import configure_logging  # noqa: E402


PROJECT_NAME = "fetal_brain_cell_bin_gsmap2d"
SAMPLE = "51_B03605C2E5_WT202406020126"


def build_config() -> QuickModeConfig:
    ensure_work_dirs(WORK_DIR)
    add_latent_h5ad = (
        WORK_DIR
        / PROJECT_NAME
        / "find_latent_representations"
        / f"{SAMPLE}_add_latent.h5ad"
    )
    return QuickModeConfig(
        workdir=WORK_DIR,
        project_name=PROJECT_NAME,
        sample_h5ad_dict=OrderedDict({SAMPLE: add_latent_h5ad}),
        sumstats_config_file=BIO_REPLICATE_SUMSTATS_CONFIG,
        w_ld_dir=GSMAP_RESOURCE_DIR / "LDSC_resource/weights_hm3_no_hla",
        snp_gene_weight_adata_path=SNP_GENE_WEIGHT_ADATA,
        annotation="subregion_major_celltype_pair",
        spatial_key="align_spatial_3d",
        data_layer="X",
        n_cell_training=100000,
        latent_representation_niche="emb_niche",
        latent_representation_cell="emb_cell",
        spatial_neighbors=300,
        adjacent_slice_spatial_neighbors=200,
        homogeneous_neighbors=21,
        n_adjacent_slices=4,
        dataset_type=DatasetType.SPATIAL_2D,
        rank_read_workers=30,
        compute_input_queue_size=5,
        mkscore_write_workers=5,
        mkscore_compute_workers=5,
        mkscore_batch_size=20,
        find_homogeneous_batch_size=200,
        memmap_tmp_dir="/data",
        high_quality_cell_qc=False,
        cauchy_annotations=[
            "regions",
            "celltype",
            "major_celltype",
            "subregion_major_celltype_pair",
        ],
        spots_per_chunk_quick_mode=50,
        ldsc_compute_workers=10,
        ldsc_read_workers=15,
        use_gpu=True,
    )


def main() -> None:
    configure_logging()
    run_quick_mode(build_config())


if __name__ == "__main__":
    main()

