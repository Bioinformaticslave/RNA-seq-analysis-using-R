# RNA-seq Analysis Pipeline Using R

A modular RNA-seq analysis workflow implemented in R for differential gene expression and functional enrichment analysis.

## Overview

This repository provides a streamlined RNA-seq analysis pipeline covering data preparation, differential expression testing, visualization and pathway enrichment. The workflow is based on standard Bioconductor tools and has been adapted from teaching material developed for the Computational Omics course at Imperial College London.

## Features

* Data preparation and quality control
* Differential expression analysis using DESeq2
* Principal Component Analysis (PCA)
* Volcano plot generation
* Gene Set Enrichment Analysis (GSEA)
* Gene Ontology (GO) over-representation analysis
* Modular script structure for easy customization

## Workflow

```text
Raw count matrix
        ↓
Data preparation
        ↓
DESeq2 differential expression analysis
        ↓
PCA visualization
        ↓
Volcano plot visualization
        ↓
GSEA pathway enrichment
        ↓
GO enrichment analysis
```

## Repository Structure

| Script            | Purpose                              |
| ----------------- | ------------------------------------ |
| 00_config.R       | Configuration and parameter settings |
| 01_prepare_data.R | Import and prepare count data        |
| 02_run_deseq2.R   | Differential expression analysis     |
| 03_plot_pca.R     | PCA visualization                    |
| 04_plot_volcano.R | Volcano plot generation              |
| 05_run_gsea.R     | Gene Set Enrichment Analysis         |
| 06_run_go_ora.R   | Gene Ontology enrichment analysis    |

## Requirements

* R (≥ 4.0)
* DESeq2
* ggplot2
* clusterProfiler
* enrichplot
* fgsea
* dplyr

## Author

Developed as part of computational biology and bioinformatics training in DoLS, Imperial College London, with modifications to improve flexibility and reproducibility of RNA-seq analysis workflows.
