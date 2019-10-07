get.values <- function(table, row.multiplier) {
  rows.index = seq(from = row.multiplier, to = nrow(table), by = 2)
  half.table = table[rows.index,]
  return(half.table)
}

setwd("/share/analysis/hecatos/juantxo/score_protein_analysis/UNTR/shift-FALSE/")

tpmprotein.table = read.table(file = 'UNTR_dataframe_containing_targetRNA_TPM_protein_1onseveral_until072_minimumexpressedsamples3_shift-FALSE_timeps-shifted-002.tsv', header = T, sep = '\t')
protein.table = get.values(table = tpmprotein.table, row.multiplier = 1)
tpm.table = get.values(table = tpmprotein.table, row.multiplier = 2)

nas.prot = is.na(protein.table$UNTR_The_002_1)
prot.test = protein.table$UNTR_The_002_1 <= mean(protein.table$UNTR_The_002_1, na.rm = T)
nas.prot = 
tpm.test = tpm.table$UNTR_The_002_1 <= mean(tpm.table$UNTR_The_002_1, na.rm = T)
all.test = prot.test + tpm.test
all.test  = all.test == 2
lm.transcprot = lm(formula = protein.table$UNTR_The_002_1[all.test] ~ tpm.table$UNTR_The_002_1[all.test])

setwd('/share/analysis/hecatos/juantxo/score_protein_analysis/plots/scatterplots/')
png(filename = 'plot_TPM_Protein_UNTRThe0021.png')

plot(x = tpm.table$UNTR_The_002_1, y = protein.table$UNTR_The_002_1, xlim = c(0, mean(tpm.table$UNTR_The_002_1)), ylim = c(0, mean(protein.table$UNTR_The_002_1, na.rm = T)))
dev.off()
lmplot2(lm.transcprot)
