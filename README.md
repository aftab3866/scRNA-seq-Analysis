# Single-Cell RNA-Seq Analysis Pipeline

A collection of six R scripts covering the standard scRNA-seq workflow: doublet
filtering, integration (Harmony & Seurat CCA), cluster identity & differential
expression, pseudo-bulk DESeq2, and the canonical Seurat pipeline.

**Reviewed by:** Dr. Aftab · **Contact / issues:** https://aftab3866.github.io/

---

## Graphical Abstract

---

## Scripts at a Glance

| # | Script | Purpose | Input → Output | Packages |
|---|--------|---------|----------------|----------|
| 1 | `single_cell_doublets/scripts/doublet_filtering.R` | Detect/remove doublets | 10x raw matrix → singlet-only `.rds` | Seurat, DoubletFinder, tidyverse, ggplot2 |
| 2 | `single_cell_integrate/scripts/harmony_integration.R` | Remove condition batch effect | `ifnb` (SeuratData) → `ifnb_harmony.rds` | Seurat, SeuratData, harmony, tidyverse, ggplot2 |
| 3 | `single_cell_integrate/scripts/seurat_cca_integration.R` | Remove patient batch effect | 7 samples → `seurat_integrated.rds` | Seurat, tidyverse, ggplot2, gridExtra |
| 4 | `single_cell_DEG/scripts/single_cell_DEG.R` | Cluster + condition markers | `ifnb_harmony.rds` → marker tables | Seurat, tidyverse |
| 5 | `single_cell_DEG/scripts/pseudobulk_deg.R` | Pseudo-bulk DESeq2 | Kang dataset → DE results `.csv` | ExperimentHub, Seurat, DESeq2, tidyverse |
| 6 | `single_cell_RNASeq/scripts/standard_seurat_workflow.R` | Canonical Seurat pipeline | 10x h5 → clustered `.rds` | Seurat, tidyverse |

---

## 1. Doublet Filtering (DoubletFinder)

Detects and removes doublets from a 10x PBMC dataset using DoubletFinder.

- **Script:** `scripts/doublet_filtering.R`
- **Input:** `data/raw_feature_bc_matrix/` (`matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`)
- **Output:** `results/objects/pbmc_filtered.rds`, `pbmc_singlets.rds`; `results/plots/pk_sweep.png`, `umap_doublets.png`
- **Steps:** Load counts → QC filter → Normalize/PCA/UMAP → pK sweep → estimate doublets → run DoubletFinder → keep singlets.
- **Run:** `setwd("single_cell_doublets"); source("scripts/doublet_filtering.R")`
- **Note:** Change the `0.076` doublet rate to match your 10x loading.

---

## 2. Harmony Integration

Removes condition (CTRL vs STIM) batch effects using Harmony.

- **Script:** `scripts/harmony_integration.R`
- **Input:** `ifnb` dataset from SeuratData
- **Output:** `results/objects/ifnb_harmony.rds`; `results/plots/umap_before.png`, `umap_after.png`, `umap_before_vs_after.png`
- **Steps:** Load ifnb → QC filter → PCA/UMAP → RunHarmony → re-run UMAP + clustering on `harmony` reduction → compare before/after.
- **Run:** `setwd("single_cell_integrate"); source("scripts/harmony_integration.R")`
- **Note:** UMAP and FindNeighbors must use `reduction = 'harmony'`.

---

## 3. Seurat CCA Integration

Removes batch effects across patients (HB17/HB30/HB53) using CCA.

- **Script:** `scripts/seurat_cca_integration.R`
- **Input:** `data/*_filtered_feature_bc_matrix/` (7 samples)
- **Output:** `results/objects/merged_filtered.rds`, `seurat_integrated.rds`; `results/plots/umap_before.png`, `umap_after.png`, `umap_before_vs_after.png`
- **Steps:** Load samples → merge → QC filter → UMAP (before) → SplitObject by Patient → FindIntegrationAnchors → IntegrateData → UMAP (after) → compare.
- **Run:** `setwd("single_cell_integrate"); source("scripts/seurat_cca_integration.R")`
- **Note:** CCA is memory-heavy; use `rpca` or Harmony for large data.

---

## 4. Cluster Identity & Differential Expression

Finds cluster markers, conserved markers, and condition DE on the Harmony-integrated IFNB object.

- **Script:** `scripts/single_cell_DEG.R`
- **Input:** `ifnb_harmony.rds`
- **Output:** Marker tables in the R session (FindAllMarkers, FindConservedMarkers, FindMarkers)
- **Steps:** UMAP clusters/conditions → FindAllMarkers → FindConservedMarkers (cluster 3, split by stim) → rename cluster 3 → FindMarkers (CD16 Mono STIM vs CTRL) → FeaturePlots.
- **Run:** `setwd("single_cell_DEG"); source("scripts/single_cell_DEG.R")`
- **Note:** `slot` must match the test — `"counts"` for DESeq2, `"data"` for Wilcoxon.

---

## 5. Pseudo-bulk DEG (DESeq2)

Pseudo-bulk DESeq2 analysis on the Kang IFN-β dataset.

- **Script:** `scripts/pseudobulk_deg.R`
- **Input:** Kang dataset via ExperimentHub (`EH2259`)
- **Output:** `results/objects/` (seu_filtered, counts_bcell, dds_bcell, res_bcell); `results/tables/res_bcell.csv`; `results/plots/umap_cell.png`, `umap_condition.png`
- **Steps:** Load Kang → QC filter → aggregate counts by (cell × sample) → split by cell type → DESeq2 on B cells (Stimulated vs Control).
- **Run:** `setwd("single_cell_DEG"); source("scripts/pseudobulk_deg.R")`
- **Note:** The assay is `originalexp`, not `RNA`. DESeq2 must run on pseudo-bulk counts, not single-cell counts.

---

## 6. Standard Seurat Workflow (NSCLC 20k DTC)

Canonical Seurat pipeline: raw counts → QC → clustering → UMAP.

- **Script:** `scripts/standard_seurat_workflow.R`
- **Input:** `data/20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5`
- **Output:** `results/objects/nsclc_seurat.rds`; `results/plots/` (QC, HVG, PCA, elbow, UMAP)
- **Steps:** Read10X_h5 → extract Gene Expression → QC filter → Normalize → HVG → Scale → PCA → Clusters → UMAP.
- **Run:** `setwd("single_cell_RNASeq"); source("scripts/standard_seurat_workflow.R")`
- **Note:** The h5 has multiple modalities — only use `` `Gene Expression` ``. Scaling all genes is memory-heavy; scale HVGs only if needed.

---

## Repository Structure

---

## Requirements

```r
# CRAN
install.packages(c("tidyverse", "ggplot2", "gridExtra", "harmony"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "ExperimentHub"))

# Seurat ecosystem
install.packages(c("Seurat", "SeuratData", "DoubletFinder"))

# SeuratData datasets (install once)
library(SeuratData); InstallData("ifnb")
