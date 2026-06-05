"""Run gsMap QuickMode on selected xgboost-labelled cell-bin samples."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("XLA_PYTHON_CLIENT_MEM_FRACTION", "0.4")

from gsMap.config import DatasetType, QuickModeConfig
from gsMap.pipeline.quick_mode import run_quick_mode

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.paths import (  # noqa: E402
    GSMAP_RESOURCE_DIR,
    IQ_SUMSTATS,
    SNP_GENE_WEIGHT_ADATA,
    WORK_DIR,
    XGBOOST_CELLBIN_H5AD_DIR,
    ensure_work_dirs,
    sorted_h5ad_dict,
)
from shared.script_utils import configure_logging  # noqa: E402


PROJECT_NAME = "20260319_cellbin_xgboost_gsMap2d_v1"
SELECTED_SAMPLES = [
    "63_B03606C1E3_WT202403310061",
    "71_A03595A4D6_WT202403310063",
    "75_D03468D1E3_WT202403310066",
]


def build_config() -> QuickModeConfig:
    ensure_work_dirs(WORK_DIR)
    return QuickModeConfig(
        workdir=WORK_DIR,
        project_name=PROJECT_NAME,
        sample_h5ad_dict=sorted_h5ad_dict(
            XGBOOST_CELLBIN_H5AD_DIR,
            selected_samples=SELECTED_SAMPLES,
        ),
        trait_name="IQ",
        sumstats_file=IQ_SUMSTATS,
        w_ld_dir=GSMAP_RESOURCE_DIR / "LDSC_resource/weights_hm3_no_hla",
        snp_gene_weight_adata_path=SNP_GENE_WEIGHT_ADATA,
        annotation="celltype_region",
        spatial_key="align_spatial_2d",
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
        cauchy_annotations=["region_h2", "celltype"],
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

