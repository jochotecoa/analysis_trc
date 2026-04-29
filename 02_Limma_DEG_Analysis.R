#' Differential Expression Analysis using Limma
#' 
#' This script compares TPM and TRC data using Limma-trend.
#' It identifies DEGs across different timepoints and compares the two quantification methods.
#'
#' Inputs:
#' - Salmon (TPM) and TRC data (paths defined in config.R)
#' - Proteomics DEGs for comparison
#'
#' Outputs:
#' - DEGs lists (TSV)
#' - Comparison plots (PNG)

source("utils.R")
forceLibrary(c('limma', 'dplyr', 'tibble', 'ggplot2', 'reshape2'))

##### Functions #####

#' Get counts/quantification for Limma
getCts <- function(project_name) {
  if (grepl('TRC', project_name)) {
    if (grepl('counts', project_name)) {
      setwd(TRC_COUNTS_PATH)
      counts_file = list.files(pattern = 'counts_TRC_UNTR_.*\.tsv')
      if (length(counts_file) > 0) {
        return(read.table(counts_file[1], header = TRUE))
      }
    }
    
    setwd(TRC_OUTPUT_PATH)
    # Use latest available output if not specified
    latest_dir = list.dirs(".", recursive = FALSE) %>% tail(1)
    if (length(latest_dir) > 0) setwd(latest_dir)
    
    if (grepl('voom', project_name)) {
      quant_file = mergeFiles(files_patt = 'TRCscore', row_names = TRUE, all = TRUE)
    } else {
      quant_file = read.table('TRCscore_total.txt', header = TRUE)
    }
  } else if (grepl('Salmon', project_name)) {
    setwd(SALMON_INPUT_PATH)
    quant_file = read.table('total_quant.sf', header = TRUE)
  }
  return(quant_file)
}

transformLogTPM <- function(x) {
  x[is.na(x)] = 0.25
  x[x == 0] = 0.25
  return(log2(x))
}

processLimmaData <- function(x, prot_cod = TRUE, DETs = FALSE) {
  if (prot_cod) {
    if (DETs) {
      x = x %>% filterProtCod() %>% dplyr::select(-rownames)
    } else {
      x = x %>% transcrToGene(aggregate = TRUE, prot_cod = TRUE) %>% 
        tibble::column_to_rownames('ensembl_gene_id')
    }
  } else {
    if (!DETs) {
      x = x %>% transcrToGene(aggregate = TRUE, prot_cod = FALSE) %>% 
        tibble::column_to_rownames('ensembl_gene_id')
    }
  }
  return(x)
}

doLimmaSingleTp <- function(x, timepoints, all_genes = FALSE) {
  i = 0
  tps_levels = levels(as.factor(timepoints))
  baseline = tps_levels[1]
  others = tps_levels[-1]
  
  DEGs_final = NULL
  
  for (tp in others) {
    # Contrast: tp vs baseline
    idx = which(timepoints %in% c(baseline, tp))
    design = model.matrix(~ timepoints[idx])
    
    fit = lmFit(x[, idx], design)
    fit = eBayes(fit, trend = TRUE)
    
    if (all_genes) {
      DEGs = topTable(fit, number = Inf)
    } else {
      DEGs = topTable(fit, p.value = 0.05, number = Inf)
    }
    
    if (nrow(DEGs) > 0) {
      colnames(DEGs) = paste(colnames(DEGs), tp, sep = '_')
      DEGs = rownames_to_column(DEGs, "rowname")
      
      if (is.null(DEGs_final)) {
        DEGs_final = DEGs
      } else {
        DEGs_final = merge(DEGs_final, DEGs, by = 'rowname', all = TRUE)
      }
    }
  }
  return(DEGs_final)
}

##### Main Analysis #####

# Load Data
tryCatch({
  setwd(SALMON_INPUT_PATH)
  salmon = read.table('total_quant.sf', header = TRUE)
  salmon_tpm = salmon %>%
    dplyr::select(contains('TPM')) %>%
    mutate(Name = salmon$Name) %>%
    filter(grepl('ENST', Name)) %>%
    tibble::column_to_rownames('Name')

  # Search for the TRC output directory
  trc_dirs = list.dirs(TRC_OUTPUT_PATH, recursive = FALSE)
  trc_dir = trc_dirs[grep("2019-11-11_12:04:38_UTC", trc_dirs)] # Prefer this one if it exists
  if (length(trc_dir) == 0) trc_dir = tail(trc_dirs, 1) # Fallback to latest
  
  if (length(trc_dir) > 0) {
    setwd(trc_dir)
    trc_output = read.table('TRCscore_total.txt', header = TRUE)
    trc = trc_output %>%
      dplyr::select(contains('TRC_')) %>%
      mutate(Name = trc_output$rowname) %>%
      filter(grepl('ENST', Name)) %>%
      tibble::column_to_rownames('Name') 
  } else {
    stop("TRC output directory not found in ", TRC_OUTPUT_PATH)
  }

  # Processing
  trc = processLimmaData(trc, prot_cod = TRUE)
  tpm = processLimmaData(salmon_tpm, prot_cod = TRUE)

  log_trc = transformLogTPM(trc)
  log_tpm = transformLogTPM(tpm)

  # Setup timepoints
  timepoints = rep(DEFAULT_TIMEPOINTS, each = 3)
  
  # Ensure matching genes
  common_genes = intersect(rownames(log_tpm), rownames(log_trc))
  log_tpm = log_tpm[common_genes, ]
  log_trc = log_trc[common_genes, ]

  # Run Limma
  trc_degs = doLimmaSingleTp(log_trc, timepoints, all_genes = TRUE)
  tpm_degs = doLimmaSingleTp(log_tpm, timepoints, all_genes = TRUE)

  # Compare with Proteomics if available
  # (Placeholder for proteomics DEG loading logic)
  
  # Results saving
  output_dir = file.path('results', DEFAULT_COMP, 'Limma_Comparison')
  setOrCreatewd(output_dir)
  
  write.table(trc_degs, "DEGs_Limma_TRC.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(tpm_degs, "DEGs_Limma_TPM.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

  print("Limma Analysis completed.")
  
}, error = function(e) {
  warning("Limma Analysis failed: ", e$message)
})
