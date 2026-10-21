# ================================================================================
#   SCRIPT:     harmony_integration.R
# LOCATION:   scripts/   (run from project root)
# AUTHOR:     [Original author]
# REVIEWER:   Dr. Aftab
# CONTACT:    If any problem is encountered, please visit:
#   https://aftab3866.github.io/
#   PURPOSE:    Integrate single-cell RNA-seq data across experimental conditions
# (stim: CTRL vs STIM) using Harmony to remove condition-driven batch
# effects. Compares UMAP before vs after integration and produces a
# Harmony-corrected Seurat object for downstream analysis.
# INPUT:      SeuratData dataset "ifnb" (installed via InstallData)
# OUTPUT:     results/objects/ifnb_harmony.rds  (Harmony-integrated Seurat object)
# results/plots/{umap_before.png, umap_after.png,
#   umap_before_vs_after.png}
# ================================================================================
#   
#   PROJECT DIRECTORY LAYOUT (recommended)
# --------------------------------------
#   single_cell_integrate/
#   ├── scripts/
#   │   └── harmony_integration.R      <- this script
# ├── results/
#   │   ├── objects/
#   │   │   └── ifnb_harmony.rds
# │   └── plots/
#   │       ├── umap_before.png
# │       ├── umap_after.png
# │       └── umap_before_vs_after.png
# └── README.md
# 
# NAMING CONVENTIONS
# ------------------
#   - Scripts:     snake_case.R          (e.g., harmony_integration.R)
# - Objects:     snake_case.rds        (e.g., ifnb_harmony.rds)
# - Plots:       snake_case.png        (e.g., umap_before_vs_after.png)
# - Directories: lowercase, no spaces  (data/, scripts/, results/)
# 

# -----------------------------------------------------------------------------
# 0. Setup — seed, paths, libraries
# -----------------------------------------------------------------------------
# Run from the project root:
#   setwd("~/Desktop/demo/single_cell_integrate")   # or your project root

set.seed(1234)                                   # reproducibility

library(harmony)         # Harmony integration
library(Seurat)          # single-cell toolkit
library(SeuratData)      # benchmark dataset distribution
library(tidyverse)       # data wrangling & pipes
library(ggplot2)         # plotting

# Output directories
objects_dir <- "results/objects"
plots_dir   <- "results/plots"
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,   recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Get data
# -----------------------------------------------------------------------------
AvailableData()                                  # list installable datasets
InstallData("ifnb")                              # install once; safe to re-run
LoadData("ifnb")                                 # loads object named `ifnb`
str(ifnb)                                        # inspect structure

# -----------------------------------------------------------------------------
# 2. QC & filtering
# -----------------------------------------------------------------------------
# Mitochondrial fraction (human gene symbols '^MT-')
ifnb$mito_percent <- PercentageFeatureSet(ifnb, pattern = '^MT-')
View(ifnb@meta.data)                             # explore QC metrics

ifnb                                             # pre-filter cell count
ifnb_filtered <- subset(ifnb,
                        subset = nCount_RNA    > 800 &
                          nFeature_RNA  > 200 &
                          mito_percent  < 5)

# -----------------------------------------------------------------------------
# 3. Standard Seurat workflow (before integration)
# -----------------------------------------------------------------------------
ifnb_filtered <- NormalizeData(ifnb_filtered)
ifnb_filtered <- FindVariableFeatures(ifnb_filtered)
ifnb_filtered <- ScaleData(ifnb_filtered)
ifnb_filtered <- RunPCA(ifnb_filtered)

ElbowPlot(ifnb_filtered)                         # guide #PC selection

ifnb_filtered <- RunUMAP(ifnb_filtered,
                         dims = 1:20, reduction = 'pca')

# UMAP BEFORE integration — expect condition-driven separation
umap_before <- DimPlot(ifnb_filtered, reduction = 'umap', group.by = 'stim')

# -----------------------------------------------------------------------------
# 4. Run Harmony
# -----------------------------------------------------------------------------
# group.by.vars = 'stim' tells Harmony which covariate (batch/condition) to
# correct for. Use plot_convergence = TRUE while developing to verify the
# objective stabilizes; disable for production for speed/cleanliness.
ifnb_harmony <- ifnb_filtered %>%
  RunHarmony(group.by.vars = 'stim', plot_convergence = FALSE)

# Confirm the 'harmony' reduction is present
ifnb_harmony@reductions

# Inspect Harmony embeddings (cells x harmony components)
harmony_embed <- Embeddings(ifnb_harmony, "harmony")
harmony_embed[1:10, 1:10]

# -----------------------------------------------------------------------------
# 5. UMAP + clustering on Harmony embeddings (NOT PCA)
# -----------------------------------------------------------------------------
# Critical: `reduction = 'harmony'` and `dims = 1:20` must MATCH the Harmony
# embedding, otherwise the workflow silently reverts to PCA and integration is
# effectively discarded.
ifnb_harmony <- ifnb_harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# UMAP AFTER integration — conditions should now be well-mixed
umap_after <- DimPlot(ifnb_harmony, reduction = 'umap', group.by = 'stim')

# -----------------------------------------------------------------------------
# 6. Compare & save
# -----------------------------------------------------------------------------
umap_before | umap_after                       # side-by-side comparison

ggsave(file.path(plots_dir, "umap_before.png"),
       umap_before, width = 6, height = 5)
ggsave(file.path(plots_dir, "umap_after.png"),
       umap_after,  width = 6, height = 5)
ggsave(file.path(plots_dir, "umap_before_vs_after.png"),
       umap_before | umap_after, width = 12, height = 5)

saveRDS(ifnb_harmony,
        file = file.path(objects_dir, "ifnb_harmony.rds"))

message("Done. Harmony object saved with ",
        ncol(ifnb_harmony), " cells.")