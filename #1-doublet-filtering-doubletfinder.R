# ================================================================================
#   SCRIPT:     doublet_filtering.R
# LOCATION:   scripts/   (run from project root)
# REVIEWER:   Dr. Aftab
# CONTACT:    If any problem is encountered, please visit:
#   https://aftab3866.github.io/
#   PURPOSE:    Detect and classify doublets in a 10x Genomics PBMC dataset using
# DoubletFinder. Loads a raw count matrix, performs QC filtering,
# runs the standard Seurat preprocessing workflow, estimates expected
# doublets, runs DoubletFinder, and visualizes singlet/doublet
# classifications on the UMAP.
# INPUT:      data/raw_feature_bc_matrix/{matrix.mtx.gz, features.tsv.gz,
#   barcodes.tsv.gz}
# OUTPUT:     Seurat object with DoubletFinder classification column in metadata
# (DF.classifications_*), plus diagnostic plots (ElbowPlot, pK sweep),
# and singlet-only object ready for downstream analysis.
# ================================================================================
#   
#   PROJECT DIRECTORY LAYOUT (recommended)
# --------------------------------------
#   single_cell_doublets/
#   ├── data/
#   │   └── raw_feature_bc_matrix/
#   │       ├── matrix.mtx.gz
# │       ├── features.tsv.gz
# │       └── barcodes.tsv.gz
# ├── scripts/
#   │   └── doublet_filtering.R        <- this script
# ├── results/
#   │   ├── objects/
#   │   │   ├── pbmc_filtered.rds
# │   │   └── pbmc_singlets.rds
# │   └── plots/
#   │       ├── elbow_plot.png
# │       ├── pk_sweep.png
# │       └── umap_doublets.png
# └── README.md
# 
# NAMING CONVENTIONS
# ------------------
#   - Scripts:   snake_case.R            (e.g., doublet_filtering.R)
# - Objects:   snake_case.rds          (e.g., pbmc_filtered.rds)
# - Plots:     snake_case.png          (e.g., umap_doublets.png)
# - Directories: lowercase, no spaces (data/, scripts/, results/)
# 
# RENAMED FROM ORIGINAL
# ---------------------
#   OLD                                          NEW
# -------------------------------------------  ----------------------------------
#   ~/Desktop/demo/single_cell_doublets/scripts  scripts/            (relative)
# /Users/kr/Desktop/demo/single_cell_doublets  .                   (project root)
# /data/raw_feature_bc_matrix/...              data/raw_feature_bc_matrix/...
# pbmc.seurat                                  pbmc_raw
# pbmc.seurat.filtered                         pbmc_filtered
# sweep.res.list_pbmc                          sweep_list
# sweep.stats_pbmc                             sweep_stats
# bcmvn_pbmc                                   pk_metrics
# 
# NOTES FROM DR. AFTAB
# ---------------------
#   - The doublet rate (0.076 here) must match the 10x loading density for YOUR
# experiment — do NOT copy blindly. 10x documents ~0.8% per 1,000 cells.
# - Homotypic adjustment is important: doublets between identical cell types
# are biologically indistinguishable and should not be "removed."
# - The DF.classifications column name encodes pN, pK, and nExp
# (e.g., DF.classifications_0.25_0.21_691). If pK or nExp change, the column
# name WILL change — update any downstream references (DimPlot, table,
#                                                      subsetting) accordingly.
# - Recommended: keep doublets labeled, then subset to singlets before
# downstream clustering/DE.
# - If any problem arises, refer to https://aftab3866.github.io/
#   ================================================================================
# =============================================================================
# Script:        doublet_filtering.R
# Location:      scripts/   (run from project root)
# Description:   Doublet detection and removal in a 10x PBMC dataset using
#                DoubletFinder on top of a standard Seurat workflow.
# Author:        [Original author]
# Reviewed by:   Dr. Aftab
# Contact:       If any problem is encountered, please visit:
#                https://aftab3866.github.io/
# Input:         data/raw_feature_bc_matrix/ (matrix.mtx.gz, features.tsv.gz,
#                barcodes.tsv.gz)
# Output:        results/objects/pbmc_filtered.rds, pbmc_singlets.rds;
#                results/plots/{elbow_plot.png, pk_sweep.png, umap_doublets.png}
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup — paths, directories, libraries
# -----------------------------------------------------------------------------
# Run from the project root:
#   setwd("~/Desktop/demo/single_cell_doublets")   # or your project root

library(Seurat)          # single-cell toolkit
library(ggplot2)         # plotting
library(tidyverse)       # data wrangling & pipes
library(DoubletFinder)   # doublet detection

# Directory layout (relative to project root)
data_dir    <- "data/raw_feature_bc_matrix"   # raw 10x input
objects_dir <- "results/objects"              # .rds outputs
plots_dir   <- "results/plots"                # .png diagnostics
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,   recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Load count matrix (10x Matrix Market format)
# -----------------------------------------------------------------------------
cts <- ReadMtx(
  mtx      = file.path(data_dir, "matrix.mtx.gz"),
  features = file.path(data_dir, "features.tsv.gz"),
  cells    = file.path(data_dir, "barcodes.tsv.gz"))

cts[1:10, 1:10]                                  # sanity-check matrix dimensions

# -----------------------------------------------------------------------------
# 2. Build Seurat object
# -----------------------------------------------------------------------------
pbmc_raw <- CreateSeuratObject(counts = cts)
str(pbmc_raw)                                    # inspect slots/assays

# -----------------------------------------------------------------------------
# 3. QC & filtering
# -----------------------------------------------------------------------------
# Compute mitochondrial fraction (human gene symbol pattern '^MT-')
pbmc_raw$mitoPercent <- PercentageFeatureSet(pbmc_raw, pattern = '^MT-')

# Filter: retain cells with sufficient counts/features and low mito content
pbmc_filtered <- subset(pbmc_raw,
                        subset = nCount_RNA   > 800  &
                          nFeature_RNA > 500  &
                          mitoPercent  < 10)

pbmc_raw                                         # pre-filter cell counts
pbmc_filtered                                    # post-filter cell counts

# -----------------------------------------------------------------------------
# 4. Standard Seurat preprocessing
# -----------------------------------------------------------------------------
pbmc_filtered <- NormalizeData(object = pbmc_filtered)
pbmc_filtered <- FindVariableFeatures(object = pbmc_filtered)
pbmc_filtered <- ScaleData(object = pbmc_filtered)
pbmc_filtered <- RunPCA(object = pbmc_filtered)

ElbowPlot(pbmc_filtered)                         # guide selection of #PCs

pbmc_filtered <- FindNeighbors(object = pbmc_filtered,
                               dims = 1:20)      # match PCA dims
pbmc_filtered <- FindClusters(object = pbmc_filtered)
pbmc_filtered <- RunUMAP(object = pbmc_filtered,
                         dims = 1:20)

# -----------------------------------------------------------------------------
# 5. pK parameter sweep (no ground-truth mode)
# -----------------------------------------------------------------------------
# paramSweep_v3 explores pN/pK/pANN combinations; summarizeSweep aggregates
# BCmetrics; find.pK returns the table used to choose optimal pK.
sweep_list  <- paramSweep_v3(pbmc_filtered, PCs = 1:20, sct = FALSE)
sweep_stats <- summarizeSweep(sweep_list, GT = FALSE)
pk_metrics  <- find.pK(sweep_stats)

# Plot BCmetric vs pK to visually locate the maximum
pk_plot <- ggplot(pk_metrics, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(title = "DoubletFinder pK sweep", x = "pK", y = "BCmetric")
ggsave(file.path(plots_dir, "pk_sweep.png"), pk_plot, width = 6, height = 4)

# Select the pK that maximizes BCmetric (optimal doublet detection)
pK <- pk_metrics %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))          # coerce to numeric scalar

# -----------------------------------------------------------------------------
# 6. Homotypic doublet proportion estimate
# -----------------------------------------------------------------------------
# Homotypic = doublets formed between cells of the SAME cluster/type.
# These are biologically indistinguishable from singlets and should be
# removed from the expected-doublet count.
annotations    <- pbmc_filtered@meta.data$seurat_clusters
homotypic_prop <- modelHomotypic(annotations)

nExp_poi <- round(0.076 * nrow(pbmc_filtered@meta.data))
# ^ 7.6% assumed doublet rate — ADJUST to match your 10x loading density.

nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))  # homotypic-adjusted

# -----------------------------------------------------------------------------
# 7. Run DoubletFinder
# -----------------------------------------------------------------------------
# sct = FALSE because we used the standard LogNormalize workflow.
# reuse.pANN = FALSE -> compute pANN from scratch.
pbmc_filtered <- doubletFinder_v3(
  pbmc_filtered,
  PCs        = 1:20,
  pN         = 0.25,            # proportion of artificial doublets
  pK         = pK,              # optimal pK from the sweep
  nExp       = nExp_poi_adj,    # homotypic-adjusted expected doublets
  reuse.pANN = FALSE,
  sct        = FALSE)

# -----------------------------------------------------------------------------
# 8. Visualize doublet classification
# -----------------------------------------------------------------------------
# IMPORTANT: The metadata column name encodes pN_pK_nExp and WILL change
# whenever pK or nExp change. Verify the actual column name first:
df_col <- grep("^DF.classifications", colnames(pbmc_filtered@meta.data),
               value = TRUE)
message("DoubletFinder column: ", df_col)

umap_doublets <- DimPlot(pbmc_filtered, reduction = 'umap', group.by = df_col) +
  ggtitle("DoubletFinder classifications")
ggsave(file.path(plots_dir, "umap_doublets.png"),
       umap_doublets, width = 7, height = 5)

# Count singlets vs doublets
print(table(pbmc_filtered@meta.data[[df_col]]))

# -----------------------------------------------------------------------------
# 9. Save results
# -----------------------------------------------------------------------------
saveRDS(pbmc_filtered,
        file = file.path(objects_dir, "pbmc_filtered.rds"))

# Keep only singlets for downstream clustering / DE
pbmc_singlets <- subset(pbmc_filtered,
                        cells = colnames(pbmc_filtered)[
                          pbmc_filtered@meta.data[[df_col]] == "Singlet"])
saveRDS(pbmc_singlets,
        file = file.path(objects_dir, "pbmc_singlets.rds"))

message("Done. Singlets retained: ", ncol(pbmc_singlets), " / ",
        ncol(pbmc_filtered))
