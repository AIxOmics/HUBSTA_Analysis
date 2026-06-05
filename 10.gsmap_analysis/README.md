# gsMap Analysis

Scripts for gsMap-based genetic enrichment analysis of the fetal brain
development atlas. The analysis is organized by input data type, with generated
AnnData files, gsMap results, logs, and figures excluded from the repository.

## Layout

```text
10.gsmap_analysis/
├── bin100_3d/       # bin100 whole-brain 3D gsMap run
├── cell_bin/        # cell-bin gsMap runs, annotation utilities, CP/SP plots
├── pseudobulk/      # pseudobulk construction, scVI embedding, gsMap run
├── shared/          # shared paths and celltype helpers
├── single_cell/     # fetal brain single-cell gsMap run
└── summary/         # final comparison and manuscript-style figures
```

## Run Scripts

| Topic | Script | Aim |
| --- | --- | --- |
| bin100 3D | `bin100_3d/run_bin100_3d_gsmap.py` | Whole-brain 3D spatial gsMap on bin100 data. |
| cell-bin | `cell_bin/run_cell_bin_gsmap.py` | One-sample cell-bin multi-trait run. |
| xgboost cell-bin | `cell_bin/run_xgboost_cell_bin_gsmap.py` | Selected xgboost-labelled cell-bin run. |
| single-cell | `single_cell/run_single_cell_gsmap.py` | scRNA gsMap using `X_harmony`. |
| pseudobulk | `pseudobulk/run_pseudobulk_gsmap.py` | Region-celltype pseudobulk gsMap using `X_scVI` / `X_umap`. |

## Required Inputs

Edit `shared/paths.py` if the lab storage paths differ on your machine. The
scripts expect these inputs:

- bin100 spatial h5ad files
- cell-bin h5ad files with `celltype`, `region`, and `regions` annotations
- fetal brain single-cell h5ad with `X_harmony`
- gsMap LDSC resources
- SNP-gene weight AnnData
- GWAS summary statistics and sumstats config YAML

To rebuild the pseudobulk input from cell-bin h5ad files:

```bash
python pseudobulk/build_subregion_celltype_pseudobulk.py
python pseudobulk/add_scvi_umap_to_pseudobulk.py
python pseudobulk/run_pseudobulk_gsmap.py
```

For cell-bin inputs that need major celltype annotations:

```bash
python cell_bin/add_major_celltype_to_cellbin.py
```
