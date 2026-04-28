
# Author: Juan Ochoteco Asensio

# This program plots expression (TPM) of RNA against the expression of their 
# translatable protein, in addition to several analyses that I did based on that
# plot

getValues <- function(table, row.multiplier) {
  rows.index = seq(from = row.multiplier, to = nrow(table), by = 2)
  half.table = table[rows.index,]
  return(half.table)
}
filterData <- function(tpm, prot.df, prot.col, outliers.rm = NULL, nonexpressed.rm = T) {
  prot = prot.df[, prot.col]
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
  prot.filt = prot.df[conditions, ]
  
  results = list('rna.tpm' = tpm.filt, 'protein' = prot.filt)
  return(results)
}
pngPlotLine <- function(filename, x, y) {
  png(filename = filename)
  plot(x = x, y = y)
  abline(lm.transcprot)
  dev.off()
}
classifyTranslationGroup <- function(df.data, lm, fractions = 3) {
  vector.translationgroups = c()
  b = as.numeric(lm[[1]])
  min.a = -(b*(max(df.data[,1], na.rm = T)))
  max.a = max(df.data[,2], na.rm = T)
  fraction = (max.a - min.a)/fractions
  a.high.limit = fraction
  a.low.limit = -fraction
  for (row in 1:nrow(df.data)) {
    high.transl = medium.transl = low.transl = FALSE
    y = df.data[row,2]
    x = df.data[row,1]
    high.transl = y >= (b*x + a.high.limit)
    medium.transl = (y < (b*x + a.high.limit)) & (y > (b*x + a.low.limit))
    low.transl = y <= (b*x + a.low.limit)
    
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
  colnames(df.new)[ncol(df.new)] = 'translation.group'
  return(df.new)
}
plotTranslGroups <- function(x, y, lm, fractions) {
  plot(x, y)
  b = as.numeric(lm[[1]])
  min.a = -(b*(max(x, na.rm = T)))
  max.a = max(y, na.rm = T)
  fraction = (max.a - min.a)/fractions
  a.high.limit = fraction
  a.low.limit = -fraction
  abline(lm.transcprot, col = 'blue')
  abline(a = min.a, b = b, col = 'red')
  abline(a = max.a, b = b, col = 'red')
  abline(a = a.high.limit, b = b, col = 'green')
  abline(a = a.low.limit, b = b, col = 'green')
  
}
findRownames <- function(df.pattern, col.pattern, df.x) {
  df.pattern[, 'prot.ID'] = NA
  for (x in 1:nrow(df.pattern)) {
    sample.x = grep(pattern = df.pattern[x, col.pattern], x = df.x)
    row.x = grep(pattern = df.pattern[x, col.pattern], x = df.x[,sample.x])
    if (length(row.x) > 1) {
      stop()
    }
    df.pattern[x, 'prot.ID'] = rownames(df.x)[row.x]
  }
}


kIntercept.fixed = T

setwd("/share/analysis/hecatos/juantxo/score_protein_analysis/UNTR/shift-FALSE/")
tpmprotein.table = read.table(
  file = 'UNTR_dataframe_containing_targetRNA_TPM_protein_1onseveral_until072_minimumexpressedsamples3_shift-FALSE_timeps-shifted-002.tsv', 
  header = T, 
  sep = '\t')
protein.table = getValues(table = tpmprotein.table, row.multiplier = 1)
tpm.table = getValues(table = tpmprotein.table, row.multiplier = 2)

tpm.all = stack(tpm.table)[,1]
prot.all = stack(protein.table)[,1]
prot.all = as.data.frame(prot.all)
prot.all[, 'prot.ID'] = rep(x = rownames(protein.table),
                            times = ncol(protein.table))
colnames(prot.all)[1] = 'prot.expr'

filtered.data = filterData(tpm = tpm.all, 
                           prot.df = prot.all, prot.col = 'prot.expr', 
                           outliers.rm = 'extreme', 
                           nonexpressed.rm = T)
if (kIntercept.fixed) {
  lm.transcprot = lm(formula = I(filtered.data$protein$prot.expr - 0) ~ 0 + filtered.data$rna.tpm)
} else {
  lm.transcprot = lm(formula = filtered.data$protein$prot.expr ~ filtered.data$rna.tpm)
}
# setwd('/share/analysis/hecatos/juantxo/score_protein_analysis/plots/scatterplots/')
# pngPlotLine(
#   filename = 'plot_TPM_Protein_UNTRThe0021_extreme-outliers_intercept-fixed.png', 
#   x = filtered.data$rna.tpm, 
#   y = filtered.data$protein)

df.data = as.data.frame(filtered.data)
df.data.classified = classifyTranslationGroup(
  df.data = df.data, 
  lm = lm.transcprot, fractions = 80)
df.data.rownames = findRownames(df.pattern = df.data.classified,
                                col.pattern = 2,
                                df.x = protein.table)


a = dput(df.data, file = 'data.txt')
lm.transcprot = lm(formula = I(a$protein.prot.expr - 0) ~ 0 + a$rna.tpm)
lm.transcprot

a = dput(df.data, file = 'data.txt')
lm.transcprot = lm(formula =protein.prot.expr ~ 0 + rna.tpm, data=a)
lm.transcprot
