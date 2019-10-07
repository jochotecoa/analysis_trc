get.values <- function(table, row.multiplier) {
  rows.index = seq(from = row.multiplier, to = nrow(table), by = 2)
  half.table = table[rows.index,]
  return(half.table)
}

setwd("/share/analysis/hecatos/juantxo/score_protein_analysis")

tpmprotein.table = read.table(file = 'UNTR_dataframe_containing_targetRNA_TPM_protein_1onseveral_until072_minimumexpressedsamples3_shift-FALSE_timeps-shifted-002.tsv', header = T, sep = '\t')
protein.table = get.values(table = tpmprotein.table, row.multiplier = 1)
tpm.table = get.values(table = tpmprotein.table, row.multiplier = 2)

plot(x = tpm.table$UNTR_The_002_1, y = protein.table$UNTR_The_002_1, xlim = c(0, mean(tpm.table$UNTR_The_002_1)), ylim = c(0, mean(protein.table$UNTR_The_002_1, na.rm = T)))
