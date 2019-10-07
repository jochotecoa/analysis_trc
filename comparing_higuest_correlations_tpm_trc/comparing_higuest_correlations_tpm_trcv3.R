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
common.best = intersect(tpm.best$protlist, trc.best$protlist)
trc.prot.diff = setdiff(trc.best$protlist, tpm.best$protlist)
trc.all.diff = trc.best[trc.prot.diff %in% trc.best$protlist, ]
setwd("/share/analysis/hecatos/juantxo/score_protein_analysis")
setwd("./UNTR/shift-FALSE/TRC")
write.table(x = trc.all.diff, 
            file = 'most_correlated_proteins_trc.tsv', 
            sep = '\t')
