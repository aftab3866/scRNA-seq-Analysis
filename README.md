# scRNA-seq-Analysis


---

## Table of Contents

1. [Doublet Filtering (DoubletFinder)](#1-doublet-filtering-doubletfinder)
2. [Harmony Integration](#2-harmony-integration)
3. [Seurat CCA Integration](#3-seurat-cca-integration)
4. [Cluster Identity & Differential Expression](#4-cluster-identity--differential-expression)
5. [Pseudo-bulk DEG (DESeq2)](#5-pseudo-bulk-deg-deseq2)
6. [Standard Seurat Workflow (NSCLC 20k DTC)](#6-standard-seurat-workflow-nsclc-20k-dtc)
7. [Repository Structure](#repository-structure)
8. [Requirements](#requirements)

---

## 1. Doublet Filtering (DoubletFinder)

Detects and removes doublets from a 10x PBMC dataset using DoubletFinder.

| | |
|---|---|
| **Script** | `scripts/doublet_filtering.R` |
| **Input** | `data/raw_feature_bc_matrix/` (`matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`) |
| **Output** | `results/objects/pbmc_filtered.rds`, `pbmc_singlets.rds`; `results/plots/pk_sweep.png`, `umap_doublets.png` |
| **Packages** | Seurat, DoubletFinder, tidyverse, ggplot2 |

**Steps:** Load counts → QC filter → Normalize/PCA/UMAP → pK sweep →
estimate doublets → run DoubletFinder → keep singlets.

```r
setwd("single_cell_doublets")
source("scripts/doublet_filtering.R")
