# ================================================================================
#   SCRIPT:     pseudobulk_deg.R
# LOCATION:   scripts/   (run from project root)
# AUTHOR:     [Original author]
# REVIEWER:   Dr. Aftab
# CONTACT:    If any problem is encountered, please visit:
#   https://aftab3866.github.io/
#   PURPOSE:    Perform pseudo-bulk differential gene expression analysis on the
# Kang et al. (2018) IFN-beta PBMC dataset. Cells are aggregated by
# (cell type x sample) to form pseudo-bulk samples, then DESeq2 is
# run to compare Stimulated vs Control within a chosen cell type
# (default: B cells).
# INPUT:      ExperimentHub dataset EH2259 (Kang et al., 2018) via
# ExperimentHub::query(eh, "Kang")
# OUTPUT:     results/objects/dds_bcell.rds, res_bcell.rds, counts_bcell.rds,
# colData_bcell.rds; results/tables/res_bcell.csv
# ================================================================================
#   
#   PROJECT DIRECTORY LAYOUT (recommended)
# --------------------------------------
#   single_cell_DEG/
#   ├── scripts/
#   │   └── pseudobulk_deg.R           <- this script
# ├── results/
#   │   ├── objects/
#   │   │   ├── seu_filtered.rds
# │   │   ├── counts_bcell.rds
# │   │   ├── colData_bcell.rds
# │   │   ├── dds_bcell.rds
# │   │   └── res_bcell.rds
# │   ├── tables/
#   │   │   └── res_bcell.csv
# │   └── plots/
#   │       ├── umap_cell.png
# │       └── umap_condition.png
# └── README.md
# 
# NAMING CONVENTIONS
# ------------------
#   - Scripts:     snake_case.R          (e.g., pseudobulk_deg.R)
# - Objects:     snake_case.rds        (e.g., res_bcell.rds)
# - Tables:      snake_case.csv        (e.g., res_bcell.csv)
# - Plots:       snake_case.png        (e.g., umap_condition.png)
# - Directories: lowercase, no spaces  (data/, scripts/, results/)
# 
# RENAMED FROM ORIGINAL
# ---------------------
#   OLD                          NEW
# ---------------------------  ----------------------------------
#   ~/Desktop/demo/single_cell_DEG   .  (project root, relative)
# seu.obj                      seu_raw
# seu.filtered                 seu_filtered
# cts                          agg_counts
# cts.t                        agg_counts_t
# cts.split                    agg_counts_split
# cts.split.modified           agg_counts_split_mod
# counts_bcell                 counts_bcell      (kept)
# cell_plot / cond_plot        umap_cell / umap_condition
# mitoPercent                  mito_percent
# res                          res_bcell
# 
# NOTES FROM DR. AFTAB
# ---------------------
#   - The Kang dataset's counts assay is called "originalexp" (not "RNA") after
#   as.Seurat() conversion from a SingleCellExperiment. `DefaultAssay()` returns
#   "originalexp" — this MUST be referenced explicitly in AggregateExpression
#   and any downstream calls.
# - QC thresholds (nFeature 200-2500, nCount > 800, mito < 5%, singlet only)
#   are dataset-specific — adapt to YOUR data.
# - Pseudo-bulk aggregation is by (cell type x sample). This collapses cells of
#   the same type+donor+condition into one profile, yielding n = #donors per
#   condition (~6 per arm here). This is the correct unit of replication for
#   DESeq2 — do NOT run DESeq2 on single-cell counts directly.
# - `keep <- rowSums(counts(dds)) >= 10` removes genes with very low coverage
#   across pseudo-bulk samples. Adjust (e.g., >= 10 in ALL samples, or a
#   filterByExpr-based approach) if you need stricter filtering.
# - The default comparison assumes colData$condition is a factor with Control
#   as the reference level. Verify with `resultsNames(dds)` — the coefficient
#   name is "condition_Stimulated_vs_Control".
# - To analyze a different cell type, change `cell_type <- "B cells"` (must
#   match a value in the `cell` metadata column exactly).
# - If any problem arises, refer to https://aftab3866.github.io/
# ================================================================================

# -----------------------------------------------------------------------------
# 0. Setup — paths, directories, libraries
# -----------------------------------------------------------------------------
# Run from the project root:
#   setwd("~/Desktop/demo/single_cell_DEG")   # or your project root

library(ExperimentHub)   # fetch Kang dataset from Bioconductor
library(Seurat)          # single-cell toolkit
library(DESeq2)          # pseudo-bulk DE
library(tidyverse)       # data wrangling

# Output directories
objects_dir <- "results/objects"
tables_dir  <- "results/tables"
plots_dir   <- "results/plots"
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,   recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Get data (Kang et al. via ExperimentHub)
# -----------------------------------------------------------------------------
eh <- ExperimentHub()
query(eh, "Kang")                                # locate the Kang dataset

sce     <- eh[["EH2259"]]                        # SingleCellExperiment
seu_raw <- as.Seurat(sce, data = NULL)           # convert to Seurat
View(seu_raw@meta.data)

# -----------------------------------------------------------------------------
# 2. QC & filtering
# -----------------------------------------------------------------------------
# Mitochondrial fraction (human '^MT-' symbols)
seu_raw$mito_percent <- PercentageFeatureSet(seu_raw, pattern = '^MT-')
View(seu_raw@meta.data)

# Filter: gene count bounds, min UMI, low mito, singlets only.
# NOTE: assay columns are suffixed "_originalexp" because the Kang data uses
# an assay named "originalexp" (not "RNA").
seu_filtered <- subset(seu_raw,
                       subset = nFeature_originalexp > 200 &
                         nFeature_originalexp < 2500 &
                         nCount_originalexp   > 800 &
                         mito_percent         < 5   &
                         multiplets == 'singlet')

seu_raw                                          # pre-filter cell counts
seu_filtered                                     # post-filter cell counts

# -----------------------------------------------------------------------------
# 3. Standard Seurat workflow (for visualization only)
# -----------------------------------------------------------------------------
seu_filtered <- NormalizeData(seu_filtered)
seu_filtered <- FindVariableFeatures(seu_filtered)
seu_filtered <- ScaleData(seu_filtered)
seu_filtered <- RunPCA(seu_filtered)

ElbowPlot(seu_filtered)                          # guide #PC selection

seu_filtered <- RunUMAP(seu_filtered, dims = 1:20)

# Visualize cell type and condition
umap_cell      <- DimPlot(seu_filtered, reduction = 'umap',
                          group.by = 'cell', label = TRUE)
umap_condition <- DimPlot(seu_filtered, reduction = 'umap',
                          group.by = 'stim')

umap_cell | umap_condition

# -----------------------------------------------------------------------------
# 4. Pseudo-bulk workflow
# -----------------------------------------------------------------------------
# Step 1: aggregate raw counts by (cell type x sample).
# Each resulting column is a pseudo-bulk profile: sum of counts across all
# cells of a given type in a given sample.
View(seu_filtered@meta.data)
seu_filtered$samples <- paste0(seu_filtered$stim, seu_filtered$ind)

DefaultAssay(seu_filtered)                       # should be "originalexp"

agg_counts <- AggregateExpression(
  seu_filtered,
  group.by        = c("cell", "samples"),        # cell type x sample
  assays          = 'originalexp',
  slot            = "counts",
  return.seurat   = FALSE)

agg_counts <- agg_counts$originalexp             # extract the matrix

# Step 2: transpose so rows = genes, cols = pseudo-bulk samples
agg_counts_t <- t(agg_counts)
agg_counts_t <- as.data.frame(agg_counts_t)

# Step 3: row names look like "<cell type>_<sample>"; split by cell type
split_rows        <- gsub('_.*', '', rownames(agg_counts_t))
agg_counts_split  <- split.data.frame(agg_counts_t, f = factor(split_rows))

# Step 4: clean row names and transpose back -> genes x samples per cell type
agg_counts_split_mod <- lapply(agg_counts_split, function(x) {
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))   # strip cell-type prefix
  t(x)
})

# Example of the regex used above:
# gsub('.*_(.*)', '\\1', 'B cells_ctrl101')   # -> "ctrl101"

# -----------------------------------------------------------------------------
# 5. DESeq2 on B cells (default cell type)
# -----------------------------------------------------------------------------
cell_type <- "B cells"                           # change to analyze other types

# 1. Counts matrix (genes x samples for this cell type)
counts_bcell <- agg_counts_split_mod[[cell_type]]

# 2. Sample-level metadata
colData <- data.frame(samples = colnames(counts_bcell))

colData <- colData %>%
  mutate(condition = ifelse(grepl('stim', samples), 'Stimulated', 'Control')) %>%
  column_to_rownames(var = 'samples')

# (additional sample-level covariates can be merged here if needed)

# -----------------------------------------------------------------------------
# 6. Run DESeq2
# -----------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts_bcell,
                              colData   = colData,
                              design    = ~ condition)

# Filter out genes with very low coverage across pseudo-bulk samples
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

# Fit DESeq2 model
dds <- DESeq(dds)

# Confirm coefficient names (expect "condition_Stimulated_vs_Control")
resultsNames(dds)

# Extract results for Stimulated vs Control
res_bcell <- results(dds, name = "condition_Stimulated_vs_Control")
res_bcell

# -----------------------------------------------------------------------------
# 7. Save results
# -----------------------------------------------------------------------------
saveRDS(seu_filtered,    file = file.path(objects_dir, "seu_filtered.rds"))
saveRDS(counts_bcell,    file = file.path(objects_dir, "counts_bcell.rds"))
saveRDS(colData,         file = file.path(objects_dir, "colData_bcell.rds"))
saveRDS(dds,             file = file.path(objects_dir, "dds_bcell.rds"))
saveRDS(res_bcell,       file = file.path(objects_dir, "res_bcell.rds"))

write.csv(as.data.frame(res_bcell),
          file = file.path(tables_dir, "res_bcell.csv"))

ggsave(file.path(plots_dir, "umap_cell.png"),
       umap_cell, width = 7, height = 5)
ggsave(file.path(plots_dir, "umap_condition.png"),
       umap_condition, width = 7, height = 5)

message("Done. DESeq2 results for '", cell_type, "': ",
        sum(res_bcell$padj < 0.05, na.rm = TRUE), " DE genes (padj < 0.05).")









