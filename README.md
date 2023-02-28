## Tumor eradication by hetIL-15 locoregional administration is related to intratumoral accumulation of a novel CD103intCD11b+dendritic cell population


This repo contains scripts for analysis of scRNAseq, bulkRNAseq and Nanostring data as described in Stellas & Karaliota et al. Cell Reports 2023

The original analysis was run using the NIH Integrated Data Analysis Platform (NIDAP), which is hosted behind NIH firewall and unfortunately cannot be accessed by external collaborators. The underlying code and datasets have been extracted from the platform and collected here to describe the logic parameters used for the analyses presented in the paper. The code has some NIDAP platform specific comments & variables, but encapsulates the core logic and analysis flow.

R package versions used are below:

- Seurat v3.1.5
- SingleR v1.0.0
- heatmap v1.0.12
- MSigDB v6.2
- GSVA v1.30.0

Each workbook folder contains input data and code related to the indicated figures below
- datasets are in the "nidap_downloads" subfolder
- rds objects are in the "rds_output" subfolder
- run_pipeline_R.R indicates the order the scripts are intended to be run

### Workbook_1 : Nanostring volcano plots
- Figure 2A, 2B, 2C, S2A, S2B, S2C


### Workbook_2 : Nanostring pathway analysis heatmaps & estimated immune cell composition
- Figure 2E, 2F, 2G, S2D, S2E, 3A


### Workbook_3 : Nanostring pathway analysis heatmaps
- Figure 2D


### Workbook_4 : Bulk RNAseq PCA
- Figure S5A


### Workbook_5 : Bulk RNAseq pathway analysis heatmaps
- Figure S5B, S5D


### Workbook_6 : scRNAseq QC, normalization, merging, clustering, celltype assignment, celltype marker establishment, pathway analysis, CD24 RNA/CITEseq overlay
- Figure 5A, 5B, 5C, 5D, S6A, S6B


### Workbook_7 : scRNAseq violin plots
- Figure S6C

Questions can be directed to kate.goldfarbmuren@nih.gov
