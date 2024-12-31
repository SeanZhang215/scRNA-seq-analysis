# Single-cell RNA-seq Analysis Pipeline

This repository contains a comprehensive pipeline for analyzing single-cell RNA sequencing data, implemented in both R (Seurat) and Python (Scanpy). 

## Pipeline Overview

1. **Data Ingestion**
   - Raw count matrix loading
   - Metadata integration
   - Initial quality metrics calculation

2. **Quality Control**
   - Cell filtering based on:
     - Gene count thresholds
     - Mitochondrial content
   - Quality metric visualization
   - Outlier removal

3. **Normalization & Dimensionality Reduction**
   - Data normalization
   - Feature selection
   - PCA computation
   - Batch correction using Harmony
   - UMAP visualization
   - Leiden clustering at multiple resolutions

4. **Cell Type Annotation**
   - Marker-based cell type identification
   - Differential expression analysis
   - Automated cell type assignment
   - Contamination removal
   - Cluster annotation visualization

5. **Enrichment Analysis**
   - Differential expression analysis
   - GO term enrichment
   - KEGG pathway analysis
   - Gene set enrichment analysis (GSEA)

## Implementation

### R Version (Using Seurat)
- Located in `R/` directory
- Based on Seurat workflow
- Uses Harmony for batch correction
- Implements COSG for marker detection

### Python Version (Using Scanpy)
- Located in `python/` directory
- Implemented in Jupyter notebooks
- Based on Scanpy framework
- Parallel implementation of the Seurat workflow

## Repository Structure
```
scRNA-seq-DRG/              
├── R/
│   ├── 0_data_ingestion.R       # Data loading and Seurat object creation
│   ├── 1_quality_ctrl.R         # Quality control and filtering
│   ├── 2_normalization.R        # Normalization and batch correction
│   ├── 3_cell_annotation.R      # Cell type identification
│   └── 4_enrichment_analysis.R  # GO, KEGG, and GSEA analysis
├── python/
│   ├── 0_data_ingestion.ipynb       # Data loading and AnnData creation
│   ├── 1_quality_control.ipynb      # Quality control and filtering
│   ├── 2_normalization.ipynb        # Normalization and batch correction
│   ├── 3_cell_annotation.ipynb      # Cell type identification
│   └── 4_enrichment_analysis.ipynb  # GO, KEGG, and GSEA analysis
└── README.md    
```

## Dependencies

### R
```r
- Seurat
- tidyverse
- harmony
- COSG
- clusterProfiler
- enrichplot
```

### Python
```python
- scanpy
- numpy
- pandas
- matplotlib
- seaborn
- harmonypy
- scipy
- gseapy
```

## Usage

Each script/notebook is numbered according to the analysis workflow. Execute them in order:

1. `0_data_ingestion`
2. `1_quality_control`
3. `2_normalization`
4. `3_cell_annotation`
5. `4_enrichment_analysis`

## Input Data Requirements

- Raw count matrix (genes × cells)
- Metadata file containing:
  - Cell IDs
  - Batch information
  - Treatment conditions
  - Other experimental metadata

## Output

The pipeline generates:
- Filtered and normalized expression data
- Quality control metrics and visualizations
- Dimensionality reduction plots
- Cell type annotations
- Enrichment analysis results

## Notes

- QC parameters should be adjusted based on your dataset
- Cell type markers are specific to neuronal subtypes in DRG
- Batch correction can be skipped if not needed
