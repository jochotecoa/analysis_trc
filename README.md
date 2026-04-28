# TRC Analysis Project

## Description

This project contains a series of R scripts for the analysis of transcriptomics (TRC) and proteomics data. The analyses include differential expression, correlation analysis, and pathway analysis to investigate the relationship between different omics layers.

## Key Improvements

- **Centralized Utilities:** Shared functions are now consolidated in `utils.R`.
- **Decluttered Workspace:** Older versions of scripts have been moved to the `archive/` directory.
- **Improved Maintainability:** Main scripts have been refactored to use `utils.R` and include descriptive headers.

## Project Structure

- `utils.R`: Centralized utility functions used across the project.
- `archive/`: Historical versions of scripts.
- `best_cases/`: Scripts for analyzing best-case scenarios or specific subsets of data.
- `comparing_higuest_correlations_tpm_trc/`: Scripts for comparing high-correlation data between TPM and TRC normalized data.
- `correlation_transcr_prot/`: Scripts for analyzing the correlation between transcript and protein expression.
- `CPDB/`: Scripts related to the ConsensusPathDB (CPDB) analysis.
- `DEGs/`: Scripts for differential gene expression (DGE) analysis using tools like DESeq2, edgeR, and limma.
- `differentially_expressed_proteins/`: Scripts for analyzing differentially expressed proteins (DEPs).
- `GOrilla/`: Scripts for Gene Ontology (GO) enrichment analysis using GOrilla.
- `integrated_timeline/`: Scripts for creating integrated timelines of different omics data.
- `miRNA_effect/`: Scripts to analyze the effect of miRNAs on gene expression.
- `p.values_multiomics/`: Scripts for multi-omics p-value analysis.
- `protein/`: Scripts for processing proteomics data.
- `quality_control/`: Scripts for quality control of miRNA and RNA-Seq data.
- `time-shift_analysis/`: Scripts for time-shift analysis between different omics layers.
- `trc_protein_comparison_severaltranscripts/`: Scripts for comparing TRC and protein data for several transcripts.

## Key Scripts

- `analysis_trc.Rproj`: RStudio project file.
- `utils.R`: The main library of helper functions.
- `trc_protein_comparison_severaltranscripts/trc_protein_comparison_severaltranscriptsv33.R`: Latest multi-transcript comparison script.
- `time-shift_analysis/time-shift_analysisv8.R`: Latest time-shift analysis script.

## How to Run

1.  Open the `analysis_trc.Rproj` file in RStudio.
2.  The `utils.R` file will be sourced by the main scripts. Ensure you have the necessary permissions for package installation.
3.  Run the individual scripts in their respective directories.
