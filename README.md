# Identification of Efficacy Driving Markers of Immunotherapy in Melanoma

This repository contains R scripts for analyzing single-cell RNA sequencing (scRNA-seq) data to identify cell subtypes, marker genes, and perform functional enrichment analysis for responder (R) and non-responder (NR) patient groups.

## File Overview

### Data Processing and Analysis Scripts

1. **Data_Preprocessing.R**
   - Processes raw scRNA-seq count data
   - Performs quality control, filtering, and normalization
   - Prepares metadata for downstream analysis

2. **Classification_and_Dimensionality_Reduction.R**
   - Normalizes and scales the data
   - Performs dimensionality reduction (PCA, t-SNE)
   - Classifies cells into major cell types

3. **B_cells_subsetting.R**
   - Extracts and analyzes B cell subsets
   - Performs consensus clustering to identify B cell subtypes
   - Identifies subtype-specific marker genes

4. **CD4_subsetting.R**
   - Processes CD4+ T cell data
   - Performs clustering and subtype identification
   - Analyzes CD4+ T cell marker genes

5. **CD8_cell_subsetting.R**
   - Focuses on CD8+ T cell re-clustering
   - Identifies CD8+ T cell subtypes and their markers
   - Prepares data for survival analysis

6. **T_cell_description.R**
   - Performs functional enrichment analysis (GO and KEGG) for T cells
   - Generates visualizations including chord diagrams and bar charts
   - Analyzes non-responder (NR) specific T cell populations

7. **T_cell_subsetting.R**
   - Subsets and clusters T cells (including γδ T cells)
   - Identifies T cell subtype markers
   - Performs correlation analysis between subtypes

8. **Treg_subsetting.R**
   - Processes regulatory T cells (Tregs)
   - Performs consensus clustering for Treg subtypes
   - Extracts Treg-specific marker genes

9. **pDCs_subsetting.R**
   - Analyzes plasmacytoid dendritic cells (pDCs)
   - Performs clustering and identifies pDC subtypes
   - Finds pDC subtype-specific marker genes

10. **marker_genes_RvsNR.R**
    - Compares marker genes between responder (R) and non-responder (NR) groups
    - Performs intersection analysis of marker genes
    - Generates visualizations for marker gene comparisons
    - Performs co-expression analysis for responder genes

11. **randomWalk2.R**
    - Implements random walk algorithm for PPI network analysis
    - Scores gene importance in protein-protein interaction networks
    - Performs statistical validation of results

12. **survival_analysis.R**
    - Prepares survival data from clinical metadata
    - Performs survival analysis for both responder and non-responder marker genes
    - Generates Kaplan-Meier survival curves

## Key Features

- **Cell Type Classification**: Identifies major immune cell types including B cells, CD4+ T cells, CD8+ T cells, T cells, Tregs, and pDCs
- **Subtype Identification**: Uses consensus clustering (ConsensusClusterPlus) to identify cell subtypes within each major cell type
- **Marker Gene Detection**: Finds differentially expressed genes (DEGs) for each cell type and subtype using Seurat's FindAllMarkers
- **Functional Enrichment**: Performs GO (Gene Ontology) and KEGG pathway enrichment analysis
- **Network Analysis**: Uses random walk algorithms on PPI networks to prioritize important genes
- **Survival Analysis**: Associates gene expression with patient survival outcomes
- **Visualization**: Generates various plots including violin plots, feature plots, heatmaps, bar charts, and survival curves

## Dependencies

The following R packages are required:
- Seurat
- dplyr
- ggplot2
- ConsensusClusterPlus
- clusterProfiler
- org.Hs.eg.db
- GOplot
- survival
- survminer
- Hmisc
- tidyr
- ggthemes
- stringr
- Matrix
- magrittr

## Usage Notes

1. The scripts are designed to be run sequentially
2. Input data files (count matrices, metadata) should be placed in the working directory
3. Output files (CSV, PNG, PDF, RData) will be generated in the same directory
4. Ensure all required packages are installed before running the scripts
5. Some scripts depend on output files from previous scripts

## Visualization Tools

- Intervene web tool: https://asntech.shinyapps.io/intervene/ (for intersection analysis visualization)
- Cytoscape (for network visualization)
- STRING (for PPI network construction)
