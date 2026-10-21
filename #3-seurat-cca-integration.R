# =============================================================================
# Script:        seurat_cca_integration.R
# Location:      scripts/   (run from project root)
# Description:   Integrate multiple scRNA-seq datasets (HB patients;
#                background / PDX / tumor) to correct for batch effects using
#                Seurat's CCA-based integration. Compares UMAPs before vs after.
# Reviewed by:   Dr. Aftab
# Contact:       If any problem is encountered, please visit:
#                https://aftab3866.github.io/
# Input:         data/<SAMPLE>_filtered_feature_bc_matrix/{matrix.mtx.gz,
#                features.tsv.gz, barcodes.tsv.gz}
# Output:        results/objects/{merged_filtered.rds, seurat_integrated.rds};
#                results/plots/{umap_before.png, umap_after.png,
#                               umap_before_vs_after.png}
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup — paths, directories, libraries
# -----------------------------------------------------------------------------
# Run from the project root:
#   setwd("~/Desktop/demo/single_cell_integrate")   # or your project root

library(Seurat)          # single-cell toolkit
library(ggplot2)         # plotting
library(tidyverse)       # data wrangling (separate, pipes)
library(gridExtra)       # grid.arrange for panel layouts

# Output directories
objects_dir <- "results/objects"
plots_dir   <- "results/plots"
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,   recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Load per-sample 10x matrices
# -----------------------------------------------------------------------------
# Each subdirectory of data/ must end with _filtered_feature_bc_matrix.
# NOTE (fragile pattern, left as-is for parity): this uses assign() into the
# GLOBAL environment, so the resulting objects must be referenced by their
# bare names in the merge() call below.
dirs <- list.dirs(path = 'data/', recursive = FALSE, full.names = FALSE)

for (x in dirs) {
  name <- gsub('_filtered_feature_bc_matrix', '', x)   # strip suffix -> sample name
  
  cts <- ReadMtx(
    mtx      = file.path('data', x, 'matrix.mtx.gz'),
    features = file.path('data', x, 'features.tsv.gz'),
    cells    = file.path('data', x, 'barcodes.tsv.gz'))
  
  # create Seurat object per sample
  assign(name, CreateSeuratObject(counts = cts))
}

# -----------------------------------------------------------------------------
# 2. Merge datasets
# -----------------------------------------------------------------------------
# add.cell.ids prefixes each cell barcode with its sample name to prevent
# collisions after merge. NOTE: ls()[3:9] is positional and fragile — it
# happens to work only if exactly the 7 sample objects are at positions 3-9
# in the workspace. A safer pattern is shown at the end of this script.
merged_seurat <- merge(
  HB17_background,
  y = c(HB17_PDX, HB17_tumor, HB30_PDX, HB30_tumor,
        HB53_background, HB53_tumor),
  add.cell.ids = ls()[3:9],
  project = 'HB')

merged_seurat

# -----------------------------------------------------------------------------
# 3. QC & filtering
# -----------------------------------------------------------------------------
View(merged_seurat@meta.data)

# Create a sample column from rownames (cell names)
merged_seurat$sample <- rownames(merged_seurat@meta.data)

# Split into Patient / Type / Barcode.
# Assumes cell names look like: <Patient>_<Type>_<Barcode>
merged_seurat@meta.data <- separate(merged_seurat@meta.data,
                                    col  = 'sample',
                                    into = c('Patient', 'Type', 'Barcode'),
                                    sep  = '_')

# Mitochondrial fraction (human '^MT-' symbols)
merged_seurat$mito_percent <- PercentageFeatureSet(merged_seurat, pattern = '^MT-')

# (explore QC here via VlnPlot / FeatureScatter as needed)

# Filter low-quality cells
merged_filtered <- subset(merged_seurat,
                          subset = nCount_RNA    > 800 &
                            nFeature_RNA  > 500 &
                            mito_percent  < 10)

merged_filtered
merged_seurat

# -----------------------------------------------------------------------------
# 4. Pre-integration workflow — expect batch effects
# -----------------------------------------------------------------------------
merged_filtered <- NormalizeData(object = merged_filtered)
merged_filtered <- FindVariableFeatures(object = merged_filtered)
merged_filtered <- ScaleData(object = merged_filtered)
merged_filtered <- RunPCA(object = merged_filtered)

ElbowPlot(merged_filtered)                       # guide #PC selection

merged_filtered <- FindNeighbors(object = merged_filtered, dims = 1:20)
merged_filtered <- FindClusters(object = merged_filtered)
merged_filtered <- RunUMAP(object = merged_filtered, dims = 1:20)

# BEFORE integration: expect strong Patient separation (batch effect)
plot_before_patient <- DimPlot(merged_filtered, reduction = 'umap',
                               group.by = 'Patient')
plot_before_type    <- DimPlot(merged_filtered, reduction = 'umap',
                               group.by = 'Type',
                               cols = c('red', 'green', 'blue'))

grid.arrange(plot_before_patient, plot_before_type, ncol = 2, nrow = 1)

# -----------------------------------------------------------------------------
# 5. Integration to correct for batch effects (CCA)
# -----------------------------------------------------------------------------
# Split into per-patient objects; each must be re-normalized + HVG-selected
# so that integration anchors are computed on comparable features.
obj_list <- SplitObject(merged_filtered, split.by = 'Patient')
for (i in seq_along(obj_list)) {
  obj_list[[i]] <- NormalizeData(object = obj_list[[i]])
  obj_list[[i]] <- FindVariableFeatures(object = obj_list[[i]])
}

# Select integration features (recurrent HVGs across datasets)
features <- SelectIntegrationFeatures(object.list = obj_list)

# Find integration anchors (CCA). For large datasets consider
# `reduction = 'rpca'` to reduce memory/time.
anchors <- FindIntegrationAnchors(object.list    = obj_list,
                                  anchor.features = features)

# Integrate data (creates an 'integrated' assay)
seurat_integrated <- IntegrateData(anchorset = anchors)

# -----------------------------------------------------------------------------
# 6. Scale, PCA, UMAP on the integrated assay, and visualize
# -----------------------------------------------------------------------------
seurat_integrated <- ScaleData(object = seurat_integrated)
seurat_integrated <- RunPCA(object = seurat_integrated)
# dims = 1:50 is typical for CCA integration (higher than PCA-only 1:20)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:50)

# AFTER integration: patient separation should be largely gone
plot_after_patient <- DimPlot(seurat_integrated, reduction = 'umap',
                              group.by = 'Patient')
plot_after_type    <- DimPlot(seurat_integrated, reduction = 'umap',
                              group.by = 'Type',
                              cols = c('red', 'green', 'blue'))

# -----------------------------------------------------------------------------
# 7. Compare & save
# -----------------------------------------------------------------------------
grid.arrange(plot_before_patient, plot_before_type,
             plot_after_patient,  plot_after_type,
             ncol = 2, nrow = 2)

ggsave(file.path(plots_dir, "umap_before.png"),
       grid.arrange(plot_before_patient, plot_before_type, ncol = 2, nrow = 1),
       width = 12, height = 5)
ggsave(file.path(plots_dir, "umap_after.png"),
       grid.arrange(plot_after_patient, plot_after_type, ncol = 2, nrow = 1),
       width = 12, height = 5)
ggsave(file.path(plots_dir, "umap_before_vs_after.png"),
       grid.arrange(plot_before_patient, plot_before_type,
                    plot_after_patient,  plot_after_type,
                    ncol = 2, nrow = 2),
       width = 12, height = 10)

saveRDS(merged_filtered,   file = file.path(objects_dir, "merged_filtered.rds"))
saveRDS(seurat_integrated, file = file.path(objects_dir, "seurat_integrated.rds"))

message("Done. Integration complete for ", ncol(seurat_integrated), " cells.")

# -----------------------------------------------------------------------------
# OPTIONAL REFACTOR — safer loading (recommended for Dr. Aftab's review)
# -----------------------------------------------------------------------------
# The assign()/ls() pattern above is fragile. A robust alternative:
#
# sample_names <- gsub('_filtered_feature_bc_matrix', '', dirs)
# obj_list_raw <- lapply(seq_along(dirs), function(i) {
#   x <- dirs[i]
#   cts <- ReadMtx(
#     mtx      = file.path('data', x, 'matrix.mtx.gz'),
#     features = file.path('data', x, 'features.tsv.gz'),
#     cells    = file.path('data', x, 'barcodes.tsv.gz'))
#   CreateSeuratObject(counts = cts)
# })
# names(obj_list_raw) <- sample_names
#
# merged_seurat <- merge(obj_list_raw[[1]],
#                        y = obj_list_raw[-1],
#                        add.cell.ids = sample_names,
#                        project = 'HB')
#
# This avoids polluting the global environment and removes the positional
# ls()[3:9] dependency.
# =============================================================================