# TCGA Cancer Microbiome Analysis

Welcome to the **TCGA Cancer Microbiome Analysis** repository! This project provides a comprehensive pipeline for analyzing microbiome data from The Cancer Genome Atlas (TCGA), focusing on data integration, preprocessing, normalization, clustering, differential abundance analysis, survival analysis, and more. The analysis is implemented in R and documented in an R Markdown file, with results rendered for easy viewing.

## Overview

This repository contains an R Markdown document (`TCGA_Microbiome_Analysis.Rmd`) that processes TCGA microbiome data for a specified cancer type (e.g., TCGA-STAD). The pipeline:
- **Loads Data**: Integrates OTU counts, metadata, clinical subtypes, immune profiles, and coverage statistics.
- **Preprocesses**: Cleans and standardizes the data.
- **Normalizes**: Applies coverage-based normalization.
- **Clusters**: Groups samples into "High" and "Low" categories based on richness and abundance.
- **Analyzes Composition**: Uses ANCOM to identify differentially abundant species, visualized with volcano and bar plots.
- **Assesses Survival**: Merges with clinical data for survival analysis using Cox models and survival plots.
- **Modular Design**: Employs helper functions for flexibility across cancer types.

The rendered output is available as `TCGA_Microbiome_Analysis.md`, which includes all text, code, and plots.

## Files

- **`TCGA_Microbiome_Analysis.Rmd`**: The source R Markdown file containing the full analysis pipeline.
- **`TCGA_Microbiome_Analysis.md`**: The knitted Markdown file with rendered results, viewable directly on GitHub.
- **`TCGA_Microbiome_Analysis_files/`**: Folder containing plot images (e.g., `.png` files) referenced in the `.md` file.
- - **`Input_data/`**: Folder contains all the original data used for the analysis

## Viewing the Analysis

To view the rendered analysis with text, code, and plots:
- Click on [`TCGA_Microbiome_Analysis.md`](TCGA_Microbiome_Analysis.md) in this repository.
- GitHub will display the Markdown file, embedding all plots stored in the `TCGA_Microbiome_Analysis_files/` folder.

## Running the Analysis Locally

### Prerequisites
- **R**: Version 4.0 or higher.
- **RStudio**: Recommended for knitting and Git integration.
- **R Packages**: Install the required libraries listed in the `.Rmd` file:
  ```R
  install.packages(c("survminer", "survival", "ggplot2", "readxl", "data.table", "gplots", "dplyr", "gridExtra", "forcats", "flexsurv", "ciTools", "lemon", "vegan", "MicrobiotaProcess", "patchwork", "stringr", "compositions", "phyloseq", "ConQuR", "doParallel", "DESeq2", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "ggbeeswarm", "nlme"))
