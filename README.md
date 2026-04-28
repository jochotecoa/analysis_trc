# TRC Analysis Project

## Description

This project contains a series of R scripts for the analysis of transcriptomics (TRC) and proteomics data. The analyses include differential expression, correlation analysis, and pathway analysis to investigate the relationship between different omics layers.

## Project Structure

The project is organized into the following main directories:

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
- `checkSeveralOmicsTrTFile.R`: Script to check omics data.
- `correlation_trc_vs_tpm.R`: Script for correlating TRC and TPM data.
- `functions.R`, `functions_JOA.R`: Contain custom functions used in the analysis.
- `heatmap_RNA.R`: Script for generating heatmaps.
- `mirna_effect_on_TPM_TRC_Protein.R`: Script to analyze miRNA effects.

## Dependencies

This project is based on R and RStudio. It requires several Bioconductor and CRAN packages, including:

- `DESeq2`
- `edgeR`
- `limma`
- `ggplot2`
- `pheatmap`
- and others.

To ensure reproducibility, it is recommended to use a version of R and Bioconductor consistent with the analysis period.

## How to Run

1.  Open the `analysis_trc.Rproj` file in RStudio.
2.  Install the required packages.
3.  Run the individual scripts as needed. The scripts are generally organized by analysis type, so you can run the scripts in the corresponding directory for a specific analysis.

**Note:** This README was generated based on the project's file structure. The descriptions are inferred and may require further detail from the original authors.
