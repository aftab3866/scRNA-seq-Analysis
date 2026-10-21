# ================================================================================
#   SCRIPT:     standard_seurat_workflow.R
# LOCATION:   scripts/   (run from project root)
# REVIEWER:   Dr. Aftab
# CONTACT:    If any problem is encountered, please visit:
#   https://aftab3866.github.io/
#   PURPOSE:    End-to-end standard Seurat workflow for a 10x scRNA-seq dataset:
#   QC, filtering, normalization, variable feature selection, scaling,
# PCA, clustering, and UMAP. Demonstrates the canonical Seurat
# pipeline from raw counts to a clustered, UMAP-embedded object.
# INPUT:      data/20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5
# (10x Genomics 20k NSCLC DTCs, 7 donors, 3' v3.1)
# OUTPUT:     results/objects/nsclc_seurat.rds
#             results/plots/{qc_violin.png, feature_scatter.png,
#                            variable_features.png, pca_heatmap.png,
#                            elbow_plot.png, umap_res0.5.png, umap_res0.1.png}
# ================================================================================
# 
# DATASET NOTE
# ------------
# The original header comment says the data source is the "10k human PBMCs
# multiome" URL, but the file and object names clearly indicate the dataset is
# "20k NSCLC DTCs from 7 donors, 3' v3.1". The PBMC URL appears to be a
# copy-paste leftover and has been corrected in the header block. Confirm the
# actual source and record it in the README before publication.
# 
# PROJECT DIRECTORY LAYOUT (recommended)
# --------------------------------------
# single_cell_RNASeq/
# ├── data/
# │   └── 20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5
# ├── scripts/
# │   └── standard_seurat_workflow.R      <- this script
# ├── results/
# │   ├── objects/
# │   │   └── nsclc_seurat.rds
# │   └── plots/
# │       ├── qc_violin.png
# │       ├── feature_scatter.png
# │       ├── variable_features.png
# │       ├── pca_heatmap.png
# │       ├── elbow_plot.png
# │       ├── umap_res0.5.png
# │       └── umap_res0.1.png
# └── README.md
# 
# NAMING CONVENTIONS
# ------------------
# - Scripts:     snake_case.R          (e.g., standard_seurat_workflow.R)
# - Objects:     snake_case.rds        (e.g., nsclc_seurat.rds)
# - Plots:       snake_case.png        (e.g., umap_res0.5.png)
# - Directories: lowercase, no spaces  (data/, scripts/, results/)
# 
# RENAMED FROM ORIGINAL
# ---------------------
# OLD                          NEW
# ---------------------------  ----------------------------------
# ~/Desktop/demo/single_cell_RNASeq/scripts   scripts/  (relative)
# nsclc.sparse.m               nsclc_sparse
# nsclc.seurat.obj             nsclc_seurat
# all.genes                    all_genes
# top10                        top10_hvg
# plot1                        hvg_plot
# percent.mt                   percent_mt
# 
# NOTES FROM DR. AFTAB
# ---------------------
# - The dataset is 10x h5 with multiple modalities (Gene Expression + Antibody
#   Capture + Multiplexing). The `$`Gene Expression`` extraction correctly
#   isolates the RNA modality — do NOT pass the full sparse matrix to
#   CreateSeuratObject.
# - QC thresholds (nFeature 200-2500, mito < 5%) are dataset-specific. The
#   2500 upper bound is aggressive and can remove large, transcriptionally
#   active cells — verify on YOUR data and justify in methods.
# - Scaling ALL genes (`features = all.genes`) is memory-heavy (~30k x 40k
#   dense matrix). In practice, scale only HVGs unless you specifically need
#   all genes scaled (e.g., for heatmaps of specific markers).
# - `ElbowPlot` is subjective; consider `JackStraw` or `PCHeatmap` for a more
#   principled choice of PCs. `dims = 1:15` should match the chosen elbow.
# - Clustering with multiple resolutions (`resolution = c(0.1, 0.3, 0.5, 0.7, 1)`)
#   creates multiple `RNA_snn_res.*` columns — useful for exploration, but
#   pick ONE resolution for downstream analysis and document why.
# - Idents is set to `RNA_snn_res.0.1` but the DimPlot shown is for
#   `RNA_snn_res.0.5` — verify which resolution you actually intend to use
#   downstream (the UMAP call uses whatever Idents is set to).
# - If any problem arises, refer to https://aftab3866.github.io/
# ================================================================================
# -----------------------------------------------------------------------------
# 0. Setup — paths, directories, libraries
# -----------------------------------------------------------------------------
# Run from the project root:
#   setwd("~/Desktop/demo/single_cell_RNASeq")   # or your project root

library(Seurat)          # single-cell toolkit
library(tidyverse)       # data wrangling & ggplot2

# Output directories
objects_dir <- "results/objects"
plots_dir   <- "results/plots"
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,   recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Load the NSCLC dataset (10x h5)
# -----------------------------------------------------------------------------
nsclc_sparse <- Read10X_h5(
  filename = 'data/20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5')

str(nsclc_sparse)                                # multi-modal: inspect names

# Extract ONLY the Gene Expression modality (ignore Antibody Capture /
# Multiplexing modalities — they must not go into the RNA assay).
cts <- nsclc_sparse$`Gene Expression`

# -----------------------------------------------------------------------------
# 2. Initialize the Seurat object (raw, non-normalized counts)
# -----------------------------------------------------------------------------
nsclc_seurat <- CreateSeuratObject(counts       = cts,
                                   project      = "NSCLC",
                                   min.cells    = 3,
                                   min.features = 200)

str(nsclc_seurat)
nsclc_seurat
# ~29,552 features across ~42,081 cells (before QC)

# -----------------------------------------------------------------------------
# 3. QC
# -----------------------------------------------------------------------------
View(nsclc_seurat@meta.data)

# Mitochondrial fraction (human '^MT-' symbols)
nsclc_seurat[["percent_mt"]] <- PercentageFeatureSet(nsclc_seurat,
                                                     pattern = "^MT-")
View(nsclc_seurat@meta.data)

# QC violin plots
qc_violin <- VlnPlot(nsclc_seurat,
                     features = c("nFeature_RNA", "nCount_RNA", "percent_mt"),
                     ncol = 3)
ggsave(file.path(plots_dir, "qc_violin.png"),
       qc_violin, width = 9, height = 4)

# Count vs feature scatter with linear trend
feature_scatter <- FeatureScatter(nsclc_seurat,
                                  feature1 = "nCount_RNA",
                                  feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
ggsave(file.path(plots_dir, "feature_scatter.png"),
       feature_scatter, width = 6, height = 5)

# -----------------------------------------------------------------------------
# 4. Filtering
# -----------------------------------------------------------------------------
# NOTE: upper bound 2500 on nFeature is dataset-specific — adjust per dataset.
nsclc_seurat <- subset(nsclc_seurat,
                       subset = nFeature_RNA > 200 &
                         nFeature_RNA < 2500 &
                         percent_mt   < 5)

# -----------------------------------------------------------------------------
# 5. Normalization
# -----------------------------------------------------------------------------
# Default = LogNormalize with scale.factor = 1e4. The explicit form is shown
# for clarity; the two are equivalent.
# nsclc_seurat <- NormalizeData(nsclc_seurat,
#                               normalization.method = "LogNormalize",
#                               scale.factor = 10000)
nsclc_seurat <- NormalizeData(nsclc_seurat)
str(nsclc_seurat)

# -----------------------------------------------------------------------------
# 6. Highly variable features (HVGs)
# -----------------------------------------------------------------------------
nsclc_seurat <- FindVariableFeatures(nsclc_seurat,
                                     selection.method = "vst",
                                     nfeatures = 2000)

# 10 most variable genes
top10_hvg <- head(VariableFeatures(nsclc_seurat), 10)

# Plot HVGs with top-10 labels
hvg_plot <- VariableFeaturePlot(nsclc_seurat)
hvg_plot <- LabelPoints(plot = hvg_plot, points = top10_hvg, repel = TRUE)
ggsave(file.path(plots_dir, "variable_features.png"),
       hvg_plot, width = 8, height = 5)

# -----------------------------------------------------------------------------
# 7. Scaling
# -----------------------------------------------------------------------------
# NOTE: scaling ALL genes is memory-heavy. Downstream use (PCA) needs only
# HVGs; scale all only if you need all-gene heatmaps later.
all_genes <- rownames(nsclc_seurat)
nsclc_seurat <- ScaleData(nsclc_seurat, features = all_genes)

str(nsclc_seurat)

# -----------------------------------------------------------------------------
# 8. Linear dimensionality reduction (PCA)
# -----------------------------------------------------------------------------
nsclc_seurat <- RunPCA(nsclc_seurat,
                       features = VariableFeatures(object = nsclc_seurat))

# Inspect top loadings per PC
print(nsclc_seurat[["pca"]], dims = 1:5, nfeatures = 5)

# Heatmap of top genes in PC1
pca_heatmap <- DimHeatmap(nsclc_seurat, dims = 1, cells = 500, balanced = TRUE)

# Elbow plot to choose #PCs (subjective; consider JackStraw for rigor)
elbow_plot <- ElbowPlot(nsclc_seurat)
ggsave(file.path(plots_dir, "elbow_plot.png"),
       elbow_plot, width = 6, height = 4)

# -----------------------------------------------------------------------------
# 9. Clustering
# -----------------------------------------------------------------------------
# dims = 1:15 should match the elbow above
nsclc_seurat <- FindNeighbors(nsclc_seurat, dims = 1:15)

# Explore multiple resolutions (creates RNA_snn_res.0.1, .0.3, .0.5, .0.7, 1)
nsclc_seurat <- FindClusters(nsclc_seurat,
                             resolution = c(0.1, 0.3, 0.5, 0.7, 1))
View(nsclc_seurat@meta.data)

# Quick look at res = 0.5
umap_res05_preview <- DimPlot(nsclc_seurat,
                              group.by = "RNA_snn_res.0.5", label = TRUE)

# -----------------------------------------------------------------------------
# 10. Set identity (choose ONE resolution for downstream analysis)
# -----------------------------------------------------------------------------
Idents(nsclc_seurat)                             # show current
Idents(nsclc_seurat) <- "RNA_snn_res.0.1"        # choose resolution here
Idents(nsclc_seurat)

# -----------------------------------------------------------------------------
# 11. Non-linear dimensionality reduction (UMAP)
# -----------------------------------------------------------------------------
# UMAP is computed on the SNN graph; dims must match FindNeighbors.
nsclc_seurat <- RunUMAP(nsclc_seurat, dims = 1:15)

# UMAP colored by current Idents (res 0.1)
umap_res01 <- DimPlot(nsclc_seurat, reduction = "umap", label = TRUE)
ggsave(file.path(plots_dir, "umap_res0.1.png"),
       umap_res01, width = 7, height = 5)

# UMAP colored by res 0.5 (for comparison) — set Idents temporarily
Idents(nsclc_seurat) <- "RNA_snn_res.0.5"
umap_res05 <- DimPlot(nsclc_seurat, reduction = "umap", label = TRUE)
ggsave(file.path(plots_dir, "umap_res0.5.png"),
       umap_res05, width = 7, height = 5)

# Restore chosen identity
Idents(nsclc_seurat) <- "RNA_snn_res.0.1"

# -----------------------------------------------------------------------------
# 12. Save
# -----------------------------------------------------------------------------
saveRDS(nsclc_seurat, file = file.path(objects_dir, "nsclc_seurat.rds"))

message("Done. Object saved with ", ncol(nsclc_seurat), " cells and ",
        nlevels(Idents(nsclc_seurat)), " clusters (res 0.1).")