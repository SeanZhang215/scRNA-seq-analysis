# =============================================================================
# Script: 02_quality_ctrl.R
# Description: Quality control and filtering of single-cell RNA-seq data
# Author: Shaopeng Zhang
# =============================================================================
# Clean workspace
rm(list = ls(all.names = TRUE))
gc()


# Project Configuration ==================================================
# Define project-specific variables
PROJ_NAME <- ""         
PROJ_DESCRIPTION <- ""    
FULL_PROJ_NAME <- paste(PROJ_NAME, PROJ_DESCRIPTION, sep="_")

# Load Required Libraries ===============================================
suppressPackageStartupMessages({
  library(Seurat)    
  library(qs)       
  library(dplyr)     
})

# Global Settings =====================================================
set.seed(42)

# Configure R options
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
input_file <- file.path(DATA_DIR, paste0(FULL_PROJ_NAME, "_seurat.qs"))
seurat_obj <- qread(input_file)

message("Loaded Seurat object with dimensions:")
print(dim(seurat_obj))

# Calculate QC Metrics =============================================
# Add mitochondrial percentage 
message("\nCalculating mitochondrial gene percentage...")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

# Generate Pre-filtering QC Plots =================================
message("\nGenerating pre-filtering QC plots...")
pre_qc_plot <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
)
pre_qc_plot
# Print pre-filtering statistics
message("\nPre-filtering cell counts:")
print(paste0("Number of cells before filtering: ", ncol(seurat_obj)))

# Apply QC Filtering =============================================
# Define QC thresholds
qc_filters <- list(
  min_features = 5000,    
  max_features = 16000,  
  min_counts = 1000000,  
  max_counts = 8000000,  
  max_mt_percent = 5     
)

# Apply filters
message("\nApplying QC filters...")
seurat_qc <- subset(
  seurat_obj,
  subset = nFeature_RNA > qc_filters$min_features &
    nFeature_RNA < qc_filters$max_features &
    nCount_RNA > qc_filters$min_counts &
    nCount_RNA < qc_filters$max_counts &
    percent.mt < qc_filters$max_mt_percent
)

# Print post-filtering statistics
message("Post-filtering cell counts:")
print(paste0("Number of cells after filtering: ", ncol(seurat_qc)))
print(paste0("Cells removed: ", ncol(seurat_obj) - ncol(seurat_qc)))
print(paste0("Percentage of cells kept: ", 
             round(ncol(seurat_qc)/ncol(seurat_obj)*100, 2), "%"))

# Generate Post-filtering QC Plots ================================
message("\nGenerating post-filtering QC plots...")
post_qc_plot <- VlnPlot(
  seurat_qc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
)
post_qc_plot
# Save Results ==================================================
# Create output filenames
output_qs <- file.path(OUTPUT_DIR, paste0(FULL_PROJ_NAME, "_seurat_qc.qs"))
output_rds <- file.path(OUTPUT_DIR, paste0(FULL_PROJ_NAME, "_seurat_qc.rds"))

# Save filtered Seurat object
message("\nSaving QC-filtered Seurat object...")
qsave(seurat_qc, output_qs)
saveRDS(seurat_qc, output_rds)
message("Quality control complete!")