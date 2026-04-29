#' Compare TPM and TRC Correlation Results
#' 
#' This script compares the correlation results obtained from TPM and TRC analyses
#' to identify which transcripts show better correlation with protein expression.
#'
#' Inputs:
#' - Correlation results from 03_Correlation_TPM_vs_Protein.R
#' - Correlation results from 04_Correlation_TRC_vs_Protein.R

source("utils.R")

##### Variables #####
comp = if(exists("DEFAULT_COMP")) DEFAULT_COMP else 'UNTR'
results_dir = 'results'

tpm_file = file.path(results_dir, comp, 'targetRNA_TPM', 'correlation_results_between_targetRNA_TPM_and_protein.tsv')
trc_file = file.path(results_dir, comp, 'TRC', 'correlation_results_between_TRC_and_protein.tsv')

##### Analyses #####

if (file.exists(tpm_file) && file.exists(trc_file)) {
  corrs.tpm = read.table(file = tpm_file, header = T, sep = '\t')
  corrs.trc = read.table(file = trc_file, header = T, sep = '\t')

  # Ensure common proteins are compared
  common_prots = intersect(corrs.tpm$protlist, corrs.trc$protlist)
  corrs.tpm = corrs.tpm[corrs.tpm$protlist %in% common_prots, ]
  corrs.trc = corrs.trc[corrs.trc$protlist %in% common_prots, ]

  limit = quantile(corrs.tpm$corlist, probs = 0.9, na.rm = T)

  tpm.best = corrs.tpm[corrs.tpm$corlist >= limit, ]
  trc.best = corrs.trc[corrs.trc$corlist >= limit, ]
  common.best = sum(tpm.best$protlist %in% trc.best$protlist)

  print(paste("Proteins in top 10% correlation (TPM):", nrow(tpm.best)))
  print(paste("Proteins in top 10% correlation (TRC):", nrow(trc.best)))
  print(paste("Common proteins in top 10%:", common.best))

  # Plot comparison
  png(file.path(results_dir, comp, 'comparison_TPM_vs_TRC.png'))
  plot(corrs.tpm$corlist, corrs.trc$corlist, 
       xlab = 'TPM Correlation', ylab = 'TRC Correlation',
       main = 'TPM vs TRC Correlation with Protein Expression',
       pch = 19, col = rgb(0.1, 0.1, 0.8, 0.3))
  abline(0, 1, col = 'red', lty = 2)
  dev.off()

} else {
  warning("Correlation files not found. Please run scripts 03 and 04 first.")
}