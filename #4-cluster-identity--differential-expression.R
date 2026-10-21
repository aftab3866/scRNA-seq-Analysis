================================================================================
SCRIPT:     single_cell_DEG / Cluster Identity & Differential Expression Pipeline
REVIEWER:   Dr. Aftab
CONTACT:    If any problem is encountered, please visit: https://aftab3866.github.io/
PURPOSE:    Downstream Seurat analysis of the Harmony-integrated IFNB single-cell
            dataset. This script (1) identifies cluster identity via marker genes,
            (2) finds conserved markers across conditions, and (3) finds
            differentially expressed (DE) markers between stimulation conditions.
INPUT:      ../ifnb_harmony.rds  (Harmony-integrated Seurat object)
OUTPUT:     Marker tables (FindAllMarkers, FindConservedMarkers, FindMarkers)
            and diagnostic FeaturePlots / DimPlots.
================================================================================

SECTION GUIDE
--------------
1. Setup & Data Loading
2. Visualization of clusters and conditions
3. FindAllMarkers        -> markers for every cluster
4. FindConservedMarkers  -> markers conserved across conditions
5. Cluster annotation    -> rename / set identities
6. FindMarkers           -> DE markers between conditions
7. Conserved vs DE feature plots

NOTES FROM DR. AFTAB
---------------------
- Always confirm `DefaultAssay()` before running FindMarkers/FindConservedMarkers.
- `slot` must match the test: normalized data ("data") for Wilcoxon,
  raw counts ("counts") for DESeq2. Mismatch = invalid results.
- `min.cutoff = 'q10'` uses the 10th percentile as the lower color scale limit;
  useful to avoid one extreme cell washing out the FeaturePlot gradient.
- Annotation is iterative (automatic + manual). Verify each cell type forms a
  distinct cluster before trusting downstream DE calls.
- If any problem arises, refer to https://aftab3866.github.io/
================================================================================
# =============================================================================
# Script:        single_cell_DEG.R
# Description:   Cluster identity, conserved markers, and differential expression
#                analysis on the Harmony-integrated IFNB single-cell dataset.
# Author:        [Original author]
# Reviewed by:   Dr. Aftab
# Contact:       If any problem is encountered, please visit:
#                https://aftab3866.github.io/
# Input:         ../ifnb_harmony.rds  (Harmony-integrated Seurat object)
# Output:        Marker tables and diagnostic UMAP/Feature plots
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------
# setwd("~/Desktop/demo/single_cell_DEG/script")   # working directory (adjust as needed)

set.seed(1234)                                      # reproducibility

library(Seurat)                                     # single-cell toolkit
library(tidyverse)                                  # data wrangling & plotting

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
ifnb_harmony <- readRDS('../ifnb_harmony.rds')      # Harmony-integrated object
str(ifnb_harmony)                                   # inspect object structure
View(ifnb_harmony@meta.data)                        # inspect cell metadata

# -----------------------------------------------------------------------------
# 2. Visualize clusters and conditions on UMAP
# -----------------------------------------------------------------------------
clusters  <- DimPlot(ifnb_harmony, reduction = 'umap',
                     group.by = 'seurat_clusters', label = TRUE)
condition <- DimPlot(ifnb_harmony, reduction = 'umap',
                     group.by = 'stim')

condition | clusters                                # side-by-side comparison

# -----------------------------------------------------------------------------
# 3. FindAllMarkers: markers for every cluster
# -----------------------------------------------------------------------------
# test.use = 'DESeq2' requires slot = 'counts' (raw counts), NOT normalized data.
FindAllMarkers(ifnb_harmony,
               logfc.threshold = 0.25,              # min log2 fold-change
               min.pct         = 0.1,               # min fraction of cells expressing
               only.pos        = TRUE,              # keep only upregulated markers
               test.use        = 'DESeq2',
               slot            = 'counts')

# -----------------------------------------------------------------------------
# 4. FindConservedMarkers: markers conserved across conditions (stim)
# -----------------------------------------------------------------------------
# Notes:
#   - `slot` depends on the test used. Default is the "data" slot (normalized).
#   - DefaultAssay(ifnb_harmony) <- 'RNA'  # uncomment if assay needs switching
DefaultAssay(ifnb_harmony)                          # confirm active assay

markers_cluster3 <- FindConservedMarkers(ifnb_harmony,
                                         ident.1      = 3,        # cluster of interest
                                         grouping.var = 'stim')   # condition to split by

head(markers_cluster3)

# Visualize a top conserved feature
FeaturePlot(ifnb_harmony, features = c('FCGR3A'), min.cutoff = 'q10')

# Explanation of min.cutoff quantile behavior:
seq(1, 5)
SetQuantile('q50', seq(1, 5))                       # 50th percentile
SetQuantile('q10', seq(1, 5))                       # 10th percentile

# -----------------------------------------------------------------------------
# 5. Cluster annotation (iterative step)
# -----------------------------------------------------------------------------
Idents(ifnb_harmony)                                # current identities
ifnb_harmony <- RenameIdents(ifnb_harmony, `3` = 'CD16 Mono')   # label cluster 3

DimPlot(ifnb_harmony, reduction = 'umap', label = TRUE)

# Cells already carry annotations in metadata (sanity check)
View(ifnb_harmony@meta.data)

# NOTE: Setting cluster identities is iterative (automatic + manual).
#       Ensure each cell type forms a distinct cluster before DE analysis.

# Use the provided Seurat annotations as identities (also a sanity check)
Idents(ifnb_harmony) <- ifnb_harmony@meta.data$seurat_annotations
Idents(ifnb_harmony)

DimPlot(ifnb_harmony, reduction = 'umap', label = TRUE)

# -----------------------------------------------------------------------------
# 6. FindMarkers: DE markers between conditions
# -----------------------------------------------------------------------------
# Combine cell type + condition into a single identity for pairwise DE
ifnb_harmony$celltype.cnd <- paste0(ifnb_harmony$seurat_annotations, '_',
                                    ifnb_harmony$stim)
View(ifnb_harmony@meta.data)
Idents(ifnb_harmony) <- ifnb_harmony$celltype.cnd

DimPlot(ifnb_harmony, reduction = 'umap', label = TRUE)

# DE: STIM vs CTRL within CD16 monocytes
b.interferon.response <- FindMarkers(ifnb_harmony,
                                     ident.1 = 'CD16 Mono_STIM',
                                     ident.2 = 'CD16 Mono_CTRL')
head(b.interferon.response)

# -----------------------------------------------------------------------------
# 7. Conserved features vs. condition-specific DE features
# -----------------------------------------------------------------------------
head(markers_cluster3)

FeaturePlot(ifnb_harmony,
            features   = c('FCGR3A', 'AIF1', 'IFIT1'),
            split.by   = 'stim',                    # facet by condition
            min.cutoff = 'q10')                     # 10th percentile lower bound