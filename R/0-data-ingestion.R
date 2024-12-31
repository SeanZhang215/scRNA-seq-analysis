# =============================================================================
# Script: 0_data_ingestion.R
# Description: Data ingestion and initial Seurat object creation for scRNA-seq 
#              analysis
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

# Load Required Libraries =====================================================
suppressPackageStartupMessages({
  library(tidyverse)  
  library(Seurat)   
  library(qs)       
})

# Global Settings ==========================================================
# Set random seed 
set.seed(42)

# Configure R options
options(
  max.print = .Machine$integer.max, 
  scipen = 999,                      
  stringsAsFactors = FALSE,          
  dplyr.summarise.inform = FALSE    
)


# File Paths =============================================================
# Define project directories
PROJECT_DIR <- file.path("path/to/projects", FULL_PROJ_NAME)
DATA_DIR <- file.path(PROJECT_DIR, "data/processed")
OUTPUT_DIR <- file.path(PROJECT_DIR, "output")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Data Injestoin ===========================================================
# Read metadata and count data
metadata <- read.csv(
  file.path(DATA_DIR, "metadata.csv"), 
  check.names = FALSE
)

raw_counts <- read.csv(
  file.path(DATA_DIR, "raw_counts.csv"),
  check.names = FALSE
)

# Data Processing ======================================================
# Extract count matrix and clean gene names
count_cols <- grep("^(M-S|F-S)", colnames(raw_counts), value = TRUE)
counts_matrix <- raw_counts[, count_cols]
rownames(counts_matrix) <- gsub("_", "-", raw_counts$Geneid)

# Create Seurat Object =================================================
# Initialize Seurat object with counts
seurat_obj <- CreateSeuratObject(
  counts = counts_matrix,
  project = FULL_PROJ_NAME,
  min.cells = 0,  
  min.features = 0
)
# Ingestion Checks ===============================================
print("Dataset dimensions:")
print(dim(seurat_obj))

print("Preview of count matrix:")
print(head(seurat_obj@assays$RNA$counts)[,1:5])

# Add Metadata =======================================================
# Verify metadata cell names
rownames(metadata) <- metadata$Cell
if (!all(rownames(metadata) %in% colnames(seurat_obj))) {
  stop("Cell names in metadata don't match Seurat object")
}

# Add metadata to Seurat object
seurat_obj <- AddMetaData(seurat_obj, metadata)

# Save Output =======================================================
# Save in both QS and RDS formats
output_qs <- file.path(OUTPUT_DIR, paste0(FULL_PROJ_NAME, "_seurat.qs"))
output_rds <- file.path(OUTPUT_DIR, paste0(FULL_PROJ_NAME, "_seurat.rds"))
message("\nSaving Seurat object...")
qsave(seurat_obj, output_qs)
saveRDS(seurat_obj, output_rds)
message("Data ingestion complete!")
