# Reproducible Analysis Archive for "A Diagnostic Analysis for Gaussian Mixture Models in Single-Cell Population Analysis"

**Author**: Quanjin (Terry) Su  
**Group**: Purdom Group, UC Berkeley  
**Date**: `r Sys.Date()`

---

## 1. Project Overview

This archive contains all the necessary code, data, and scripts to reproduce the analysis and figures for the summer research report titled, "A Diagnostic Analysis for Gaussian Mixture Models in Single-Cell Population Analysis."

The project investigates the properties of Gaussian Mixture Models (GMMs) as applied to single-cell RNA sequencing data within the GloScope framework. Using an HIV-PBMC dataset, this analysis explores the "matching problem"—the challenge of consistently mapping statistical GMM components to biological cell types across different samples.

The main findings are:
1.  GloScope's global divergence metric is robust to changes in GMM parameterization.
2.  Individual GMM components are flexible statistical constructs that prioritize modeling the overall data density rather than acting as consistent biological units, making a post-hoc decomposition of compositional vs. expression differences infeasible within the current framework.

This archive is self-contained and fully portable.

## 2. How to Run the Analysis

The entire analysis workflow is orchestrated by a single R Markdown file that follows the exact structure and logic of the original research paper. To reproduce all figures from the report, please follow these steps:

**Prerequisites:**
*   **R**: Version 4.0 or higher.
*   **RStudio**: Recommended for easily knitting the R Markdown document.
*   **Required R Packages**: Before running, ensure you have the necessary packages installed. You can install them by running the following command in an R console:
    ```r
    install.packages(c("here", "dplyr", "ggplot2", "readr", "Seurat", "pheatmap", "RColorBrewer", "viridis", "gridExtra", "GloScope", "knitr"))
    ```

**Steps to Reproduce:**

1.  **Open the R Markdown File**: Launch RStudio and open the `analysis_notebook.Rmd` file located in the root of this archive.
2.  **Knit the Document**: Click the "Knit" button in the RStudio interface.

This will execute the complete analysis workflow following the paper's structure:
- **Figure 1**: GloScope Divergence Heatmap (baseline results)
- **Figure 2**: MDS Plot of GloScope Divergences (baseline results)  
- **Table 1**: GMM Parameter Selection under Default Settings (k=1-9)
- **Table 2**: GMM Parameter Selection under Expanded Settings (k=6-20)
- **Figure 3**: BIC Curves for Model Selection Across Expanded k-Range
- **Figure 4**: Comparison of MDS Plots Before and After GMM Refitting
- **Figure 5**: PCA Visualization of GMM Components and Cell Types
- **Figure 6**: Heatmap of Distances Between GMM Component Centers

All figures will be displayed directly within the HTML output document, providing an integrated, reproducible analysis report.

## 3. Directory Structure

The archive is organized as follows:

```
Summer_Report_Archive/
│
├── analysis_notebook.Rmd     # The main R Markdown file to run the entire analysis.
│
├── README.md                 # This explanatory file.
│
├── code/                     # Contains all R scripts.
│   ├── HIV_PBMC/             # Scripts for individual analysis components.
│   │   ├── generate_bic_plot.r                    # BIC curve analysis (Figure 3)
│   │   ├── generate_6months_pca_plot.r            # PCA visualization (Figure 5)
│   │   ├── generate_6month_heatmap.r              # Distance heatmap (Figure 6)
│   │   ├── generate_mds_plot_for_default_k_range.r # 1-9k MDS analysis
│   │   └── generate_mds_plot_for_expanded_k_range.r # 6-20k MDS analysis
│   └── utils/                # Reusable utility functions.
│       ├── util_visualize_gloscope.r              # MDS and heatmap utilities
│       ├── util_bic_visualization.r               # BIC analysis utilities
│       ├── util_unified_condition_heatmaps.r      # Distance heatmap utilities
│       └── util_multi_patient_visualization.r     # PCA visualization utilities
│
├── data/                     # Contains all necessary input data.
│   ├── processed_data/       # Pre-processed Seurat object.
│   └── gloscope_results/     # GMM parameters and distance matrices from GloScope.
│       ├── 1-9k/             # GMM parameters for samples (k=1-9 range).
│       ├── 6-20k/            # GMM parameters for samples (k=6-20 range).
│       ├── defaultGMM_1-9k_distMat.rds    # Distance matrix (baseline)
│       ├── defaultGMM_6-20k_distMat.rds   # Distance matrix (expanded)
│       ├── gmm_fitted_models_1-9k.rds   # Fitted GMM models (Default)
│       └── gmm_fitted_models_6-20k.rds    # Fitted GMM models (Expand)
│
└── final_report/             # Output directory (created during knitting).
    └── (analysis outputs)    # All generated plots and data files.
```

## 4. Notes on Data

*   The data included in this archive is a subset of the full project's results, containing only the files necessary to reproduce the figures in the final report.
*   The `6-20k` identifier in file names corresponds to GMMs fitted with an expanded k-range (6-20), which was found to be optimal.
*   The `1-9k` identifier corresponds to GMMs fitted with the default k-range (1-9), used for comparison. The data for "6 Months" samples from this run is included for completeness.

---
