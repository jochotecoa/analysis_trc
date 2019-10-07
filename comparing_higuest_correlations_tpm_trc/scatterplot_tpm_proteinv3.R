get.values <- function(table, row.multiplier) {
  rows.index = seq(from = row.multiplier, to = nrow(table), by = 2)
  half.table = table[rows.index,]
  return(half.table)
}
filter.data <- function(tpm, prot, outliers.rm = NULL, nonexpressed.rm = T) {
  num.cond = nas.prot = zeros.tpm = prot.test = tpm.test = 0
  if (nonexpressed.rm) {
    nas.prot = !is.na(prot)
    zeros.tpm = tpm != 0
    num.cond = num.cond + 2
  }
  if (!is.null(outliers.rm)) {
    prot.summ = summary(prot)
    prot.IQR = prot.summ[['3rd Qu.']] - prot.summ[['1st Qu.']]
    tpm.summ = summary(tpm)
    tpm.IQR = tpm.summ[['3rd Qu.']] - tpm.summ[['1st Qu.']]
    
    if (outliers.rm=='mild') {
      prot.test = prot >= (prot.summ[['1st Qu.']] - 1.5*prot.IQR) & prot <= (prot.summ[['3rd Qu.']] + 1.5*prot.IQR)
      prot.test[is.na(prot.test)] = FALSE
      tpm.test = tpm >= (tpm.summ[['1st Qu.']] - 1.5*tpm.IQR) & tpm <= (tpm.summ[['3rd Qu.']] + 1.5*tpm.IQR)
      tpm.test[is.na(tpm.test)] = FALSE
    }
    if (outliers.rm=='extreme') {
      prot.test = prot >= (prot.summ[['1st Qu.']] - 3*prot.IQR) & prot <= (prot.summ[['3rd Qu.']] + 3*prot.IQR)
      prot.test[is.na(prot.test)] = FALSE
      tpm.test = tpm >= (tpm.summ[['1st Qu.']] - 3*tpm.IQR) & tpm <= (tpm.summ[['3rd Qu.']] + 3*tpm.IQR)
      tpm.test[is.na(tpm.test)] = FALSE
    }
    num.cond = num.cond + 2
  }
  
  conditions = nas.prot + zeros.tpm + prot.test + tpm.test
  conditions = conditions == num.cond
  
  tpm.filt = tpm[conditions]
  prot.filt = prot[conditions]
  
  results = list('rna.tpm' = tpm.filt, 'protein' = prot.filt)
  return(results)
}
png.plot.line <- function(filename, x, y) {
  png(filename = filename)
  plot(x = x, y = y)
  abline(lm.transcprot)
  dev.off()
}
classify.translation.group <- function(df.data) {
  vector.translationgroups = c()
  for (row in 1:nrow(df.data)) {
    high.transl = medium.transl = low.transl = FALSE
    high.transl = (1413*df.data[row,1] + 13018.019) < df.data[row,2]
    medium.transl = ((1413*df.data[row,1]) > df.data[row,2]) & ((1413*df.data[row,1] - 90279.371) < df.data[row,2])
    low.transl = (1413*df.data[row,1] - 90279.371) > df.data[row,2]
    
    if (high.transl) {
      group = 'high'
    }
    if (medium.transl) {
      group = 'medium'
    }
    if (low.transl) {
      group = 'low'
    }
    vector.translationgroups = c(vector.translationgroups, group)
  }
  df.new = cbind.data.frame(df.data, vector.translationgroups)
  colnames(df.new)[3] = 'Translation group'
  return(df.new)
}


intercept.fixed = T

setwd("/share/analysis/hecatos/juantxo/score_protein_analysis/UNTR/shift-FALSE/")
tpmprotein.table = read.table(file = 'UNTR_dataframe_containing_targetRNA_TPM_protein_1onseveral_until072_minimumexpressedsamples3_shift-FALSE_timeps-shifted-002.tsv', header = T, sep = '\t')
protein.table = get.values(table = tpmprotein.table, row.multiplier = 1)
tpm.table = get.values(table = tpmprotein.table, row.multiplier = 2)

tpm.all = stack(tpm.table)[,1]
prot.all = stack(protein.table)[,1]

filtered.data = filter.data(tpm = tpm.all, prot = prot.all, outliers.rm = 'extreme', nonexpressed.rm = T)
if (intercept.fixed) {
  lm.transcprot = lm(formula = I(filtered.data$protein - 0) ~ 0 + filtered.data$rna.tpm)
} else {
  lm.transcprot = lm(formula = filtered.data$protein ~ filtered.data$rna.tpm)
}
# setwd('/share/analysis/hecatos/juantxo/score_protein_analysis/plots/scatterplots/')
# png.plot.line(filename = 'plot_TPM_Protein_UNTRThe0021_extreme-outliers_intercept-fixed.png', x = filtered.data$rna.tpm, y = filtered.data$protein)

df.data = as.data.frame(filtered.data)
df.data = classify.translation.group(df.data = df.data)

