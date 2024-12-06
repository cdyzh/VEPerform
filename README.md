# VEPerform

**VEPerform** is an interactive web application built with **Shiny** to analyze and visualize the performance of Variant Effect Predictors (VEPs) using Balanced Precision-Recall Curves (BPRCs).

## Features

- **Precision-Recall Curve Generation**: 
  - Visualize the performance of VEP scores across multiple thresholds using BRPCs.
  - Compare multiple predictors simultaneously.
  
- **Log Likelihood Ratio (LLR) Calculation (In Development)**:
  - Calculate LLR of pathogenicity using kernel density estimation (or Gaussian distributions) using the same set of variants.
  - Generate LLR plots for multiple VEPs.

- **Custom Data Upload**:
  - Upload your own CSV file with any number of VEPs for customized analysis.

- **Variant Filtering and Selection**:
  - Interactively filter and deselect variants based on ClinVar annotation quality.
  - Support for overriding Variants of Uncertain Significance (VUS) and Conflicting classification.

- **Export Options**:
  - Download plot PNG and PDF with metadata.
  - Download variants used as a CSV for further analysis.

## Usage

### 1. Data Input
- Use the **Basic** tab to generate BPRCs for a specific gene.
- Use the **Advanced** tab to upload your own reference sets or to retrieve variant annotations from OpenCRAVAT.

### 2. Generating PRC and LLR Plots
- Select one or more predictors to analyze.

### 3. Export Results
- Download plots and associated data for further analysis or reporting.

## Installation
- This is currently a Web App at https://cdyzh.shinyapps.io/veperform/. To run your own, clone this repository and extract full.csv from its folder. The package "yogiroc" is installed using remotes::install_github("cdyzh/yogiroc").