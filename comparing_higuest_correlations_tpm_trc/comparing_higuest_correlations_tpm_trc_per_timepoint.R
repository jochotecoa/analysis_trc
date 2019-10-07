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
    if (sd(x.numeric) == 0) {
      cor.row = 0
    } else {
      cor.row = cor(x = x.numeric, y = y.numeric)
    }
    if (length(cor.row) > 1) {
      stop('Comparing columns instead of rows')
    }
    cor.list = c(cor.list, cor.row)
  }
  cor.num = as.numeric(cor.list)
  cor.df = data.frame(cors = cor.num, row.names = names)
  return(cor.df)
}
filterByLimit <- function(x, limit) {
  filtered = x[x > limit, ]
  names = rownames(x)[x[, 1] > limit]
}

################################## VARIABLES ##################################
kMain.dir = '/share/analysis/hecatos/juantxo/score_protein_analysis/'
kSpecific.dir = 'UNTR/shift-FALSE/targetRNA_TPM/data_tables/'
file.name = 'minimumexpressedsamples3.tsv'
################################## ANALYSES ###################################
setwd(paste0(kMain.dir, kSpecific.dir))
tpmprotein.table = read.table(
  file = file.name, 
  header = T, 
  sep = '\t')

rownames(tpmprotein.table) = gsub('TRC', 'TPM', rownames(tpmprotein.table))

protein.table = getValues(table = tpmprotein.table, row.multiplier = 1)
tpm.table = getValues(table = tpmprotein.table, row.multiplier = 2)
protein.table.expressed = filterUnexpressed(x = protein.table, 
                                            table = c(tpm.table, protein.table))
tpm.table.expressed = filterUnexpressed(x = tpm.table, 
                                        table = c(tpm.table, protein.table))

correlation.values.002 = corPerRow(x = tpm.table.expressed[, 1:3], 
                               y = protein.table.expressed[, 1:3])
correlation.values.008 = corPerRow(x = tpm.table.expressed[, 4:6], 
                                   y = protein.table.expressed[, 4:6])
correlation.values.024 = corPerRow(x = tpm.table.expressed[, 5:9], 
                                  y = protein.table.expressed[, 5:9])
correlation.values.072 = corPerRow(x = tpm.table.expressed[, 10:12], 
                                   y = protein.table.expressed[, 10:12])


limit = quantile(correlation.values.002, probs = 0.9, na.rm = T)

tpm.002.best = filterByLimit(x = correlation.values.002, limit = limit)
tpm.008.best = filterByLimit(x = correlation.values.008, limit = limit)
tpm.024.best = filterByLimit(x = correlation.values.024, limit = limit)
tpm.072.best = filterByLimit(x = correlation.values.072, limit = limit)

Reduce(f = intersect, 
       x = list(tpm.002.best, tpm.008.best, tpm.024.best, tpm.072.best))
