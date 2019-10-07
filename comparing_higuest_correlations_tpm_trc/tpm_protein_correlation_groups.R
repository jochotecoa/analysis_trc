################################## FUNCTIONS ##################################

getValues <- function(table, row.multiplier) {
  rows.index = seq(from = row.multiplier, to = nrow(table), by = 2)
  half.table = table[rows.index,]
  return(half.table)
}
filterUnexpressed <- function(x, table = x) {
  rows.non.expressed = complete.cases(table)
  table.filtered = x[rows.non.expressed, ]
}
corPerRow <- function(x, y, names = rownames(x)) {
  cor.list = list()
  for (i in 1:nrow(x)) {
    x.numeric = as.numeric(x[i, ])
    y.numeric = as.numeric(y[i, ])
    if (identical(x.numeric, rep(0, 12))) {
      cor.row = 0
    } else {
      cor.row = cor(x = x.numeric, y = y.numeric)
    }

    if (length(cor.row) > 1) {
      stop('Comparing columns instead of rows')
    }
    cor.list = c(cor.list, cor.row)
  }
  names(cor.list) = names
  return(cor.list)
}
getTopPc <- function(x, percentage) {
  per.one = percentage/100
  inv.per.one = 1 - per.one
  limit = quantile(as.numeric(x), probs = inv.per.one)
  x.fraction = x[x >= limit]
}

################################## VARIABLES ##################################
kMain.dir = '/share/analysis/hecatos/juantxo/score_protein_analysis/'
kSpecific.dir = 'UNTR/shift-FALSE/targetRNA_TPM/data_tables/'
################################## ANALYSES ###################################
setwd(paste0(kMain.dir, kSpecific.dir))
tpmprotein.table = read.table(
  file = 'minimumexpressedsamples3.tsv', 
  header = T, 
  sep = '\t')
protein.table = getValues(table = tpmprotein.table, row.multiplier = 1)
tpm.table = getValues(table = tpmprotein.table, row.multiplier = 2)
protein.table.expressed = filterUnexpressed(x = protein.table, 
                                            table = c(tpm.table, protein.table))
tpm.table.expressed = filterUnexpressed(x = tpm.table, 
                                        table = c(tpm.table, protein.table))
correlation.values = corPerRow(x = tpm.table.expressed, y = protein.table.expressed)

top10pc.correlations = getTopPc(x = correlation.values, percentage = 10)
corrs.top10.stacked = stack(top10pc.correlations)
corrs.top10.df = data.frame(corrs.top10.stacked[, 1], 
                            row.names = corrs.top10.stacked[, 2])
