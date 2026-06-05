"""Run gsMap QuickMode on fetal brain single-cell data."""

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
    ASD_SUMSTATS,
    GSMAP_RESOURCE_DIR,
    IQ_SUMSTATS,
    SINGLE_CELL_H5AD,
    SNP_GENE_WEIGHT_ADATA,
    WORK_DIR,
    ensure_work_dirs,
)
from shared.script_utils import configure_logging  # noqa: E402


PROJECT_NAME = "20260319_fetal_brain_SC_gsmap"


def build_config() -> QuickModeConfig:
    ensure_work_dirs(WORK_DIR)
    return QuickModeConfig(
        workdir=WORK_DIR,
        project_name=PROJECT_NAME,
        sample_h5ad_dict=OrderedDict({"Fetal_SC": SINGLE_CELL_H5AD}),
        sumstats_config_dict={
            "IQ": IQ_SUMSTATS,
            "ASD": ASD_SUMSTATS,
        },
        w_ld_dir=GSMAP_RESOURCE_DIR / "LDSC_resource/weights_hm3_no_hla",
        snp_gene_weight_adata_path=SNP_GENE_WEIGHT_ADATA,
        annotation="annotation_major",
        data_layer="X",
        n_cell_training=100000,
        latent_representation_cell="X_harmony",
        spatial_neighbors=300,
        adjacent_slice_spatial_neighbors=200,
        homogeneous_neighbors=21,
        n_adjacent_slices=4,
        dataset_type=DatasetType.SCRNA_SEQ,
        rank_read_workers=30,
        compute_input_queue_size=5,
        mkscore_write_workers=5,
        mkscore_compute_workers=5,
        mkscore_batch_size=20,
        find_homogeneous_batch_size=200,
        memmap_tmp_dir="/data",
        high_quality_cell_qc=False,
        cauchy_annotations=[
            "annotation_sub",
            "region",
            "library_meta",
            "week",
        ],
        spots_per_chunk_quick_mode=50,
        ldsc_compute_workers=10,
        ldsc_read_workers=15,
        use_gpu=True,
        start_step="latent2gene",
    )


def main() -> None:
    configure_logging()
    run_quick_mode(build_config())


if __name__ == "__main__":
    main()

