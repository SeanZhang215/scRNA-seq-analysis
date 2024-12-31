# =============================================================================
# Script: 3_cell_annotation.R
# Description: Cell type identification and annotation for neurons
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
  library(ggplot2)
  library(COSG)
})

# Global Settings ====================================================
set.seed(42)
options(
  max.print = .Machine$integer.max,
  scipen = 999,
  stringsAsFactors = FALSE,
  dplyr.summarise.inform = FALSE
)

# File Paths =======================================================
PROJECT_DIR <- file.path("path/to/projects", FULL_PROJ_NAME)
DATA_DIR <- file.path(PROJECT_DIR, "output")
OUTPUT_DIR <- file.path(PROJECT_DIR, "output")

# Load Data =======================================================
seurat_obj <- qread(file.path(DATA_DIR, paste0(FULL_PROJ_NAME, "_normalized.qs")))

# Check Glial Contamination ========================================
message("Checking for glial contamination...")
glial_markers <- c("Atp1b2", "Fabp7", "Sostdc1", "Timp3")

# Plot glial markers
DotPlot(seurat_obj, 
        features = glial_markers,
        group.by = "RNA_snn_res.1",
        scale = TRUE) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Remove Glial Cells ==============================================
message("Removing glial cells...")
# Define glial clusters based on marker expression
glial_clusters <- c(6)
if(length(glial_clusters) > 0) {
  seurat_obj <- subset(seurat_obj, 
                       idents = glial_clusters, 
                       invert = TRUE)
  
  # Rerun dimensional reduction on only-neurons data
  seurat_obj <- seurat_obj %>%
    RunPCA() %>%
    RunHarmony("Batch", plot_convergence = TRUE) %>%
    FindNeighbors(reduction = "harmony") %>%
    FindClusters(resolution = 1) %>%
    RunUMAP(reduction = "harmony", dims = 1:20)
}

# Neuronal Subtype Analysis =======================================
message("Analyzing neuronal subtypes...")
neuronal_subtypes <- list(
  mNP = c("P2rx3", "Mrgprd", "Gfra2"),
  mNFa = c("Necab2", "Nefh"),
  mNFb = c("Fam19a1", "Nefh"),
  mPEPa = c("Smr2", "Calca"),
  mPEPb = c("Trpa1", "Calca"),
  pNF = c("Spp1", "Nefh"),
  pPEP = c("Th", "Calca")
)

# Plot neuronal markers
all_markers <- unique(unlist(neuronal_subtypes))
print(DotPlot(seurat_obj, 
              features = all_markers,
              group.by = "RNA_snn_res.1",
              scale = TRUE) +
        coord_flip() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))

# Find Cluster Markers ============================================
# Set marker analysis parameters
marker_params <- list(
  find_all_markers = FALSE,    
  use_cosg = TRUE,         
  min_pct = 0.25,           
  logfc_threshold = 0.25,     
  n_genes_cosg = 200    
)

# Initialize list to store marker results
marker_results <- list()

# Method 1: FindAllMarkers
if(marker_params$find_all_markers) {
  message("Running FindAllMarkers...")
  marker_results$findall <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = marker_params$min_pct,
    logfc.threshold = marker_params$logfc_threshold,
    test.use = "wilcox"
  )
  
  top_markers_findall <- marker_results$findall %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
  
  message("\nTop markers from FindAllMarkers:")
  print(head(top_markers_findall))
}

# Method 2: COSG
if(marker_params$use_cosg) {
  message("\nRunning COSG analysis...")
  marker_results$cosg <- COSG::cosg(
    seurat_obj,
    groups = 'all',
    assay = 'RNA',
    slot = 'data',
    mu = 1,
    expressed_pct = marker_params$min_pct,
    remove_lowly_expressed = TRUE,
    n_genes_user = marker_params$n_genes_cosg
  )
  
  # Display top markers from COSG
  markers_df_cosg <- as.data.frame(marker_results$cosg$names)
  message("\nTop markers from COSG:")
  print(head(markers_df_cosg, 10))
}

# Compare results if both methods used
if(marker_params$find_all_markers && marker_params$use_cosg) {
  message("\nComparing top markers between methods...")
  for(cluster in unique(marker_results$findall$cluster)) {
    message(sprintf("\nCluster %s top markers:", cluster))
    findall_top <- marker_results$findall %>%
      filter(cluster == !!cluster) %>%
      top_n(5, avg_log2FC) %>%
      pull(gene)
    
    cosg_top <- markers_df_cosg[[as.character(cluster)]][1:5]
    common_genes <- intersect(findall_top, cosg_top)
    message("FindAllMarkers: ", paste(findall_top, collapse=", "))
    message("COSG: ", paste(cosg_top, collapse=", "))
    message("Common: ", paste(common_genes, collapse=", "))
  }
}
# Assign Cell Types ==============================================
message("Assigning cell types...")
# Initialize cell type assignments
n_clusters <- length(unique(seurat_obj@active.ident))
celltype <- data.frame(
  ClusterID = 0:(n_clusters-1),
  celltype = 'Unknown'
)

# Define cell type assignments based on marker expression
neuron_assignments <- list(
  "mNP" = c(0, 2),     # Clusters showing mNP markers
  "mNFa" = c(7),       # Clusters showing mNFa markers
  "mNFb" = c(3),       # Clusters showing mNFb markers
  "mPEPa" = c(5),      # Clusters showing mPEPa markers
  "mPEPb" = c(6),      # Clusters showing mPEPb markers
  "pNF" = c(4),        # Clusters showing pNF markers
  "pPEP" = c(1)        # Clusters showing pPEP markers
)

# Assign cell types
for(celltype_name in names(neuron_assignments)) {
  clusters <- neuron_assignments[[celltype_name]]
  celltype$celltype[celltype$ClusterID %in% clusters] <- celltype_name
}

# Add cell type annotations to Seurat object
seurat_obj$celltype <- "Unknown"
for(i in 1:nrow(celltype)){
  seurat_obj$celltype[seurat_obj@active.ident == celltype$ClusterID[i]] <- 
    celltype$celltype[i]
}

# Final Visualization ===========================================
# UMAP by cell type
DimPlot(seurat_obj, group.by = "celltype", label = TRUE) + 
  ggtitle("Cell Types")

# UMAP split by treatment
DimPlot(seurat_obj, 
        group.by = "celltype", 
        split.by = "Treatment", 
        label = TRUE)

FeaturePlot(seurat_obj, 
            features = unique(unlist(neuronal_subtypes)[1:4]),
            ncol = 2)

# Save Results ==================================================
output_file <- file.path(OUTPUT_DIR, 
                         paste0(FULL_PROJ_NAME, "_annotated.qs"))
message("Saving annotated Seurat object...")
qsave(seurat_obj, output_file)