################################## VARIABLES ##################################
kMain.dir = '/share/analysis/hecatos/juantxo/score_protein_analysis/'
kSpecific.dir = 'UNTR/shift-FALSE/targetRNA_TPM/correlation_results/'
comp = 'UNTR'
shift = 'shift-FALSE'
input = 'correlation_results'
file.name = 'minimumexpressedsamples-12.tsv'
################################## ANALYSES ###################################
setwd(paste(kMain.dir, comp, shift, 'targetRNA_TPM', input, sep = '/'))
corrs.tpm = read.table(
  file = file.name, 
  header = T, 
  sep = '\t', 
  stringsAsFactors = F)
corrs.tpm = corrs.tpm[complete.cases(corrs.tpm), ]

setwd(paste(kMain.dir, comp, shift, 'TRC', input, sep = '/'))
corrs.trc = read.table(
  file = file.name, 
  header = T, 
  sep = '\t',
  stringsAsFactors = F)
corrs.trc = corrs.trc[complete.cases(corrs.trc), ]


limit = quantile(corrs.tpm$corlist, probs = 0.9, na.rm = T)

tpm.best = corrs.tpm[corrs.tpm$corlist >= limit, ]
trc.best = corrs.trc[corrs.trc$corlist >= limit, ]
common.best = intersect(tpm.best, trc.best)
