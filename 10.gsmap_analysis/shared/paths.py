"""Central path configuration for local gsMap fetal brain analyses.

The original analysis used absolute paths on the Yangjian lab storage system.
They are centralized here so workflow scripts stay readable and can be adapted
by changing one file or overriding command-line arguments.
"""

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = ANALYSIS_DIR.parent

WORK_DIR = Path(
    "/storage/yangjianLab/chenwenhao/01_Project/01_Research/202312_gsMap/"
    "experiment/20260107_gsmap_3d_report_all_collected/human_brain"
)

GSMAP_RESOURCE_DIR = Path(
    "/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/resource/"
    "gsMap_resource"
)
GSMAP_EXAMPLE_DATA_DIR = Path(
    "/storage/yangjianLab/chenwenhao/01_Project/01_Research/202312_gsMap/"
    "data/gsMap_dev_data/online_resource/gsMap_example_data"
)
SNP_GENE_WEIGHT_ADATA = Path(
    "/storage/yangjianLab/chenwenhao/01_Project/01_Research/202312_gsMap/"
    "experiment/20251125_gsMap_Omni/RNA/1000GP3_GRCh37/"
    "1000GP3_GRCh37_gencode_v46_protein_coding_ldscore_weights.h5ad"
)
BIO_REPLICATE_SUMSTATS_CONFIG = Path(
    "/storage/yangjianLab/chenwenhao/01_Project/01_Research/202312_gsMap/"
    "experiment/20241106_Human_Fetal_Brain_HongtaoYu/"
    "02_03_bio_replcate_sumstats.yaml"
)

IQ_SUMSTATS = GSMAP_EXAMPLE_DATA_DIR / "GWAS/IQ_NG_2018.sumstats.gz"
ASD_SUMSTATS = Path(
    "/storage/yangjianLab/songliyang/GWAS_trait/LDSC/"
    "PGC_ASD_2019_NatGenet.sumstats.gz"
)

BIN100_H5AD_DIR = Path(
    "/storage/yangjianLab/songliyang/SpatialData/Data/Embryo/Human/YHT/"
    "count_data"
)
CELLBIN_H5AD_DIR = Path(
    "/storage/yangjianLab/chenwenhao/data/ST/human/brain/"
    "2024_YHT_fetal_brain/v1/3d/cellbin"
)
CELLBIN_MAJOR_CELLTYPE_DIR = Path(
    "/storage/yangjianLab/chenwenhao/data/ST/human/brain/"
    "2024_YHT_fetal_brain/v1/3d/cellbin_major_cell_type"
)
XGBOOST_CELLBIN_H5AD_DIR = Path(
    "/storage/yangjianLab/chenwenhao/data/ST/human/brain/"
    "2024_YHT_fetal_brain/v1/3d/xgboost_cell_type/cell_bin_wt_xgboost_ct"
)
SINGLE_CELL_H5AD = Path(
    "/storage/yangjianLab/chenwenhao/data/ST/human/brain/"
    "2024_YHT_fetal_brain/v1/single_cell_data/"
    "fetal_brain_SC_adata_correct_meta.h5ad"
)

PSEUDOBULK_SUBREGION_CELLTYPE_H5AD = WORK_DIR / "pseudobulk_subregion_celltype.h5ad"
PSEUDOBULK_WEIGHTED_MEAN_H5AD = WORK_DIR / "pseudobulk_weighted_mean.h5ad"


def sample_name_from_h5ad(path: Path) -> str:
    """Return sample name before .celltype suffix when present."""
    return path.stem.split(".celltype")[0]


def sorted_h5ad_dict(
    h5ad_dir: Path,
    pattern: str = "*.h5ad",
    selected_samples: list[str] | None = None,
) -> OrderedDict[str, Path]:
    """Collect h5ad files as an OrderedDict sorted by numeric sample prefix."""
    selected = set(selected_samples) if selected_samples else None
    items: list[tuple[str, Path]] = []
    for path in h5ad_dir.glob(pattern):
        name = sample_name_from_h5ad(path)
        if selected is not None and name not in selected:
            continue
        items.append((name, path))

    def sort_key(item: tuple[str, Path]) -> tuple[int, str]:
        name, _ = item
        try:
            return int(name.split("_", 1)[0]), name
        except ValueError:
            return 10**9, name

    return OrderedDict(sorted(items, key=sort_key))


def ensure_work_dirs(workdir: Path = WORK_DIR) -> Path:
    """Create workdir and logs directory, then return the logs path."""
    workdir.mkdir(parents=True, exist_ok=True)
    log_dir = workdir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    return log_dir

