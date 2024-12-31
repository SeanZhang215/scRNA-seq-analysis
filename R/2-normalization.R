# =============================================================================
# Script: 2_normalization.R
# Description: Normalization, dimension reduction, and initial clustering of 
#              scRNA-seq data
# Author: Shaopeng Zhang
# =============================================================================
# Clean workspace
rm(list = ls(all.names = TRUE))
gc()
# Project Configuration ==================================================
PROJ_NAME <- ""              # Project identifier
PROJ_DESCRIPTION <- ""     # Project type/description
FULL_PROJ_NAME <- paste(PROJ_NAME, PROJ_DESCRIPTION, sep="_")

# Load Required Libraries ===============================================
suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(dplyr)
  library(patchwork)
  library(harmony)
})

# Global Settings =====================================================
set.seed(42)
options(
  max.print = .Machine$integer.max,
  scipen = 999,
  stringsAsFactors = FALSE,
  dplyr.summarise.inform = FALSE
)


# File Paths ========================================================
PROJECT_DIR <- file.path("path/to/projects", FULL_PROJ_NAME)
DATA_DIR <- file.path(PROJECT_DIR, "output")
OUTPUT_DIR <- file.path(PROJECT_DIR, "output")

# Load Data ========================================================
input_file <- file.path(DATA_DIR, paste0(FULL_PROJ_NAME, "_seurat_qc.qs"))
seurat_obj <- qread(input_file)

# Normalization ====================================================
message("Performing data normalization...")
seurat_obj <- seurat_obj %>% 
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(
    selection.method = "vst", 
    nfeatures = 2000, 
    verbose = FALSE
  ) %>% 
  ScaleData(verbose = FALSE)

# Dimension Reduction ==============================================
message("Running PCA...")
seurat_obj <- RunPCA(seurat_obj)

# PCA Quality Assessment
message("Generating PCA quality plots...")
ElbowPlot(seurat_obj)
DimPlot(seurat_obj, reduction = "pca", group.by = "Batch") +
  ggtitle("PCA before batch correction")

# Batch Correction with Harmony ===================================
message("Performing batch correction with Harmony...")
seurat_obj <- seurat_obj %>%
  RunHarmony(
    "Batch",              
    plot_convergence = TRUE,
    assay.use = "RNA"
  )

# Verify batch correction
DimPlot(seurat_obj, reduction = "harmony", group.by = "Batch") +
  ggtitle("Harmony after batch correction")

# Compare PCA and Harmony
p1 <- DimPlot(seurat_obj, reduction = "pca", group.by = "Batch") + 
  ggtitle("PCA")
p2 <- DimPlot(seurat_obj, reduction = "harmony", group.by = "Batch") + 
  ggtitle("Harmony")
p1 + p2

# Initial Clustering ==============================================
message("Performing initial clustering...")
seurat_obj <- seurat_obj %>%
  FindNeighbors(reduction = "harmony") %>%  
  FindClusters(resolution = 0.5)

# UMAP Generation ================================================
message("Generating UMAP...")
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:20) 

# Generate Multiple Resolution Clusterings ========================
message("Testing multiple clustering resolutions...")
resolutions <- c(0.05, 0.1, 0.3, 0.5, 0.75, 1, 1.25, 1.5)

for (res in resolutions) {
  message(sprintf("Testing resolution: %s", res))
  seurat_obj <- FindClusters(
    seurat_obj, 
    resolution = res, 
    algorithm = 1
  )
}

# Visualization =================================================
message("Generating visualization plots...")

# Multi-resolution clustering visualization
resolution_plots <- lapply(resolutions, function(res) {
  DimPlot(seurat_obj, 
          reduction = "umap", 
          group.by = paste0("RNA_snn_res.", res), 
          label = TRUE) & NoAxes()
})
print(wrap_plots(resolution_plots, ncol = 3))

# Basic UMAP visualizations
DimPlot(seurat_obj, reduction = "umap", group.by = "Treatment")
DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters")

# Compare batch effects in UMAP
DimPlot(seurat_obj, reduction = "umap", group.by = "Batch") +
  ggtitle("Batch distribution after correction")

# Save Results ==================================================
output_file <- file.path(OUTPUT_DIR, paste0(FULL_PROJ_NAME, "_normalized.qs"))
message("Saving normalized and batch-corrected Seurat object...")
qsave(seurat_obj, output_file)