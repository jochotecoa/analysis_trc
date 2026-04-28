#' Configuration file for TRC Analysis Project
#' Centralizes base paths and common variables

# Base paths for data storage
BASE_ANALYSIS_PATH = '/share/analysis/hecatos/juantxo/'
BASE_NGS_PATH = '/ngs-data/data/hecatos/'
BASE_CEFIC_PATH = '/ngs-data-2/data/CEFIC/'

# Specific data directories
PROTEOMICS_BASE_PATH = file.path(BASE_NGS_PATH, 'Cardiac')
SCORE_OUTPUT_PATH = file.path(BASE_ANALYSIS_PATH, 'Score/output/Output_Run_mrna_SEPT2019/')
TABLEOMICS_PATH = file.path(BASE_ANALYSIS_PATH, 'tableomics/')

# Common experimental parameters
DEFAULT_COMPOUND = 'Con_UNTR'
DEFAULT_COMP = 'UNTR'
DEFAULT_TIMEPOINTS = c('002', '008', '024', '072')
DEFAULT_TRIPLICATES = 1:3
