# =============================================================================
# Script: 4_enrichment_analysis.R
# Description: Enrichment analysis
# Author: Shaopeng Zhang
# =============================================================================
# Clean workspace
rm(list = ls(all.names = TRUE))
gc()
# Project Configuration ==================================================
PROJ_NAME <- ""         
PROJ_DESCRIPTION <- "" 
FULL_PROJ_NAME <- paste(PROJ_NAME, PROJ_DESCRIPTION, sep="_")

# Load Required Libraries ===============================================
suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(dplyr)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(fgsea)
  library(msigdbr)
  library(org.Mm.eg.db) 
  
})

# Global Settings =====================================================
set.seed(42)
options(
  max.print = .Machine$integer.max,
  scipen = 999,
  stringsAsFactors = FALSE,
  dplyr.summarise.inform = FALSE
)

# Set up parallel processing
future::plan("multisession", workers = 4)
options(future.globals.maxSize = 10000 * 1024^3) # 10 GB

# File Paths ========================================================
PROJECT_DIR <- file.path("path/to/projects", FULL_PROJ_NAME)
DATA_DIR <- file.path(PROJECT_DIR, "output")
OUTPUT_DIR <- file.path(PROJECT_DIR, "output")

# Load Data ========================================================
seurat_obj <- qread(file.path(DATA_DIR, paste0(FULL_PROJ_NAME, "_annotated.qs")))

# Find Differentially Expressed Genes ================================
message("Finding differentially expressed genes...")
Idents(seurat_obj) <- "celltype"
de_results <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# Filter significant DEGs
sig_genes <- de_results %>%
  filter(avg_log2FC >= 0.8, p_val_adj < 0.05)

# Create gene lists by cluster
gene_lists <- split(sig_genes$gene, sig_genes$cluster)

# GO Enrichment Analysis ============================================
message("Performing GO enrichment analysis...")
go_results <- lapply(gene_lists, function(genes) {
  enrichGO(
    gene = genes,
    OrgDb = "org.Mm.eg.db",
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = 0.05
  )
})

# Visualize GO results
plot_go <- function(go_result, n_terms = 5) {
  if(nrow(go_result@result) > 0) {
    barplot(go_result, showCategory = n_terms, font.size = 10) +
      theme(axis.text.y = element_text(size = 8))
  }
}

go_plots <- lapply(names(go_results), function(cluster) {
  plot_go(go_results[[cluster]]) + ggtitle(cluster)
})
wrap_plots(go_plots, ncol = 3)

# KEGG Pathway Analysis ============================================
message("Performing KEGG pathway analysis...")
# Convert gene symbols to ENTREZ IDs
gene_lists_entrez <- lapply(gene_lists, function(genes) {
  bitr(genes, fromType = "SYMBOL", 
       toType = "ENTREZID", 
       OrgDb = "org.Mm.eg.db")$ENTREZID
})

# Run KEGG enrichment
kegg_results <- lapply(gene_lists_entrez, function(genes) {
  enrichKEGG(
    gene = genes,
    organism = "mmu",
    pvalueCutoff = 0.05
  )
})

# Create comparison plot for KEGG results
kegg_compare <- compareCluster(
  geneClusters = gene_lists_entrez,
  fun = "enrichKEGG",
  organism = "mmu"
)
dotplot(kegg_compare, showCategory = 5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# GSEA Analysis ===================================================
message("Performing GSEA analysis...")

# Prepare ranked gene lists
prepare_ranked_list <- function(de_results, cluster_id) {
  cluster_genes <- de_results %>%
    filter(cluster == cluster_id) %>%
    group_by(gene) %>%
    summarize(avg_log2FC = max(avg_log2FC)) %>%
    arrange(desc(avg_log2FC))
  
  gene_list <- setNames(cluster_genes$avg_log2FC, cluster_genes$gene)
  return(gene_list)
}

# Get hallmark gene sets
hallmark_sets <- msigdbr(
  species = "Mus musculus",  # Use "Homo sapiens" for human
  category = "H"
) %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Run GSEA for each cluster
run_gsea_analysis <- function(ranked_genes, gene_sets) {
  tryCatch({
    fgsea(
      pathways = gene_sets,
      stats = ranked_genes,
      minSize = 15,
      maxSize = 500,
      nperm = 1000
    )
  }, error = function(e) {
    message("Error in GSEA analysis: ", e$message)
    return(NULL)
  })
}

gsea_results <- list()
for(cluster in unique(de_results$cluster)) {
  ranked_genes <- prepare_ranked_list(de_results, cluster)
  gsea_results[[cluster]] <- run_gsea_analysis(ranked_genes, hallmark_sets)
}

# Create GSEA visualization
message("Creating GSEA visualization...")
create_gsea_plot <- function(gsea_results) {
  # Combine results from all clusters
  gsea_df <- do.call(rbind, lapply(names(gsea_results), function(cluster) {
    if (!is.null(gsea_results[[cluster]])) {
      res <- as.data.frame(gsea_results[[cluster]])
      res$cluster <- cluster
      return(res)
    }
  }))
  
  if (!is.null(gsea_df)) {
    ggplot(gsea_df, 
           aes(x = cluster, y = pathway, 
               size = -log10(padj), 
               color = NES)) +
      geom_point() +
      scale_color_gradient2(
        low = "blue", 
        mid = "white", 
        high = "red", 
        midpoint = 0
      ) +
      labs(
        title = "GSEA Results",
        x = "Cell Type",
        y = "Pathway",
        size = "-log10(FDR)",
        color = "NES"
      )
  }
}

gsea_plot <- create_gsea_plot(gsea_results)
gsea_plot
# Save Results ==================================================
message("Saving enrichment results...")
enrichment_results <- list(
  DEGs = de_results,
  GO = go_results,
  KEGG = kegg_results,
  GSEA = gsea_results
)

output_file <- file.path(OUTPUT_DIR, paste0(FULL_PROJ_NAME, "_enrichment.qs"))
qsave(enrichment_results, output_file)
