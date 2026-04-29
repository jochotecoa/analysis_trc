# TRC Analysis Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![Conda](https://img.shields.io/badge/conda-environment-green.svg)](https://docs.conda.io/en/latest/)

## Overview

This repository contains a professional bioinformatics pipeline for the analysis and integration of Transcriptomics (TRC), regular RNA-Seq (TPM), and Proteomics data. The pipeline encompasses differential expression analysis, correlation analyses between multi-omics layers, and time-shift/half-life analysis, primarily applied to time-course experimental data.

## Project Structure

A clean, modular architecture allows for reproducible runs:

```
analysis_trc/
├── config.R             # Global paths and experimental parameters
├── utils.R              # Common helper functions
├── environment.yml      # Conda environment definition for reproducibility
├── README.md            # Project documentation
├── LICENSE              # MIT License
├── 01_DESeq2_RNA_Seq.R            # DE analysis (DESeq2)
├── 02_Limma_DEG_Analysis.R        # DE analysis (Limma)
├── 03_Correlation_TPM_vs_Protein.R # Transcript vs Protein baseline (TPM)
├── 04_Correlation_TRC_vs_Protein.R # Translation vs Protein baseline (TRC)
├── 05_Compare_Correlations.R      # Comparison between correlations
├── 06_Correlation_With_Shifts.R   # Dynamic time-shift enabled correlations
├── 07_Time_Shift_Analysis.R       # Half-life and time-shift distributions
├── data/                          # [Ignored] Raw and processed datasets
├── results/                       # [Ignored] Analysis outputs (TSV, RDA)
├── figures/                       # [Ignored] Output figures (PDF, PNG)
└── archive/                       # Historical script versions
```

## Installation & Setup

We highly recommend using Conda to manage your dependencies to ensure reproducibility.

1. **Clone the repository:**
   ```bash
   git clone https://github.com/jochotecoa/analysis_trc.git
   cd analysis_trc
   ```

2. **Create the Conda environment:**
   We have provided an `environment.yml` to automatically install R, Bioconductor packages, and CRAN dependencies.
   ```bash
   conda env create -f environment.yml
   ```

3. **Activate the environment:**
   ```bash
   conda activate analysis_trc
   ```

## Configuration

Before running the pipeline, update `config.R` with the correct paths for your local machine or high-performance computing (HPC) environment. You need to adjust variables like:
- `BASE_ANALYSIS_PATH`
- `BASE_NGS_PATH`
- `PROTEOMICS_BASE_PATH`
- `SALMON_INPUT_PATH`
- `TRC_OUTPUT_PATH`

## Pipeline Execution

The scripts are designed to be run sequentially from `01` to `07`.

1. **Differential Expression Analysis**
   - Run `01_DESeq2_RNA_Seq.R` and `02_Limma_DEG_Analysis.R` to compute differentially expressed genes/transcripts.
2. **Correlation Profiling**
   - Run `03` and `04` to measure basic correlation coefficients between RNA layers and Proteome layers.
3. **Advanced Time-Shift Modeling**
   - Run `05`, `06`, and `07` to analyze how long it takes for a transcriptomic shift to reflect at the proteomic level, integrating half-life calculations.

*(You can run them in RStudio via `analysis_trc.Rproj` or from the terminal using `Rscript 01_DESeq2_RNA_Seq.R`, etc.)*

## Output

All generated tables, lists of DEGs/DEPs, and intermediate `.Rda` files are saved into a specified `results/` folder (configured in `config.R`). All visualizations and plots are placed into `figures/`. *(Note: These folders are ignored by git to keep the repository lightweight).*