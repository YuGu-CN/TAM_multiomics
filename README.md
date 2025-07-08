# 🧬 TAM_multiomics

**Integrated single-cell multi-omics analysis of tumor-associated macrophage (TAM) subtypes in hepatocellular carcinoma (HCC)**

This repository contains reproducible R scripts and helper functions accompanying our manuscript:

> **"Integrated single-cell analysis dissects regulatory mechanisms underlying tumor-associated macrophage plasticity in hepatocellular carcinoma"**
>
> ---

## 📁 Repository structure

| File / Folder              | Description |
|---------------------------|-------------|
| `scATAC_preprocess.R`     | Preprocessing pipeline for scATAC-seq using ArchR (incl. quality control and LSI embedding) |
| `scRNA_preprocess.R`      | Preprocessing pipeline for scRNA-seq using Seurat |
| `scATAC_macro.R`          | TAM subset annotation, label transfer, marker gene scoring, and embedding visualization |
| `scRNA_macro.R`           | TAM subset identification in scRNA-seq and preparation for integration |
| `archr_helpers.R`         | Custom helper functions for ArchR-based workflows |
| `seurat_helpers.R`        | Custom helper functions for Seurat processing |
| `plotting_config.R`       | ggplot2 theme and color palette configuration |
| `matrix_helpers.R`        | Utility functions for manipulating ArchR matrices |
| `misc_helpers.R`          | Miscellaneous functions used across scripts |
| `GO_wrappers.R`           | Enrichment and visualization wrappers for GO/KEGG analysis |
| `Figure.R`                | All code for generating main and supplementary figures |
| `README.md`               | This file |

---

## 🔧 Dependencies

Please install the following packages in R (≥4.1.0):

```r
# Required
library(ArchR)
library(Seurat)
library(ComplexHeatmap)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
```

## 🧾 Citation

If you use this code, please cite our manuscript:
Gu Y, Zhu W, Zhang Z, *et al.*  
**Integrated single-cell analysis dissects regulatory mechanisms underlying tumor-associated macrophage plasticity in hepatocellular carcinoma.**  
genes, *Under review*, 2025.
