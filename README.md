# Tumor–myeloid compartmentalization of polyamine metabolism in Ewing sarcoma

Analysis code for the manuscript:

> Cheng X, Pan L. *Single-cell analysis of two Ewing sarcoma cohorts reveals tumor–myeloid compartmentalization of polyamine metabolism.* (submitted)

The script `ES_analysis.R` reproduces all results, figures and supplementary tables of the manuscript from public GEO data.

## Data

All data are public (Gene Expression Omnibus).

| Dataset | Type | Files needed | Where to put them |
| --- | --- | --- | --- |
| [GSE243347](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243347) | scRNA-seq, CEL-seq2 (cohort A) | `GSE243347_RAW.tar` | anywhere; selected in a file dialog on first run |
| [GSE277083](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277083) | scRNA-seq, 10x Genomics (cohort B) | `GSE277083_cell_metadata_primary.txt.gz`, `GSE277083_counts_matrix_primary.txt.gz` | anywhere; selected in a file dialog on first run |
| [GSE17674](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17674), [GSE17618](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17618), [GSE63155](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63155), [GSE63156](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63156) | bulk microarrays with survival | `<GSE>_series_matrix.txt.gz` | `ES_project/01_data/bulk/` (downloaded automatically if missing) |

GSE17618 is read only to document its overlap with GSE17674; it is not analyzed.

## Requirements

R ≥ 4.6 with Seurat 5, harmony, Matrix, data.table, survival, ggplot2, ggrepel, patchwork, openxlsx (CRAN) and AnnotationDbi, org.Hs.eg.db, reactome.db, hgu133plus2.db, huex10sttranscriptcluster.db (Bioconductor). Missing packages are installed by the script. Exact versions used are listed in `sessionInfo.txt`.

## How to run

1. Create a working directory and set it in R, e.g. `setwd("path/to/workdir")`.
2. Run `source("ES_analysis.R")` (or open it in RStudio and run all).
3. On the first run, select the single-cell input files when prompted. Building the cohort A object takes about 15–20 minutes.

The script creates this structure inside the working directory:

```
ES_project/
├── 01_data/      Seurat objects, gene sets, bulk series matrix files (bulk/)
├── 02_results/   result tables, cached intermediate results, Supplementary_Tables.xlsx
└── 03_figures/   all figures as PDF, PNG and TIFF (300 dpi)
```

Time-consuming steps are cached as `.rds` files in `02_results/`; delete a file to recompute that step. Random steps (module scores, random gene sets, down-sampling) use fixed seeds, so results are reproducible.

## Script structure

| Part | Content | Manuscript |
| --- | --- | --- |
| 0 | Packages, paths, helper functions, 14 metabolic gene sets | Supplementary Table 1 |
| 1–2 | Cohort A and B processing and annotation | Methods 2.1–2.2, Fig. 1 |
| 3 | Patient-level module screen | Methods 2.3, Fig. 2 |
| 4 | Gene-level analysis of 40 candidate genes | Methods 2.4, Fig. 3 |
| 5 | Genome-wide differential expression and Reactome over-representation | Methods 2.4, Fig. 4 |
| 6 | Exploratory bulk survival analysis | Methods 2.6, Supplementary Figs. 1–2 |
| 7 | Sensitivity analyses | Methods 2.5, Supplementary Table 2 |
| 8 | Supplementary tables (Excel) | Supplementary Tables 1–5 |
| 9 | Exploratory analysis of the ENPP1 direction | Methods 2.7, Supplementary Fig. 3, Supplementary Table 6 |

## Contact

Lei Pan (corresponding author), Panlei21316@163.com

## License

MIT (see `LICENSE`).
