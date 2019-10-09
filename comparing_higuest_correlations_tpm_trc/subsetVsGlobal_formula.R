###################

forceLibrary <- function(list.of.packages) {
  list.of.packages <- c("pbmcapply", "foreach", "doParallel", "biomaRt")
  new.packages.log = !(list.of.packages %in% installed.packages()[,"Package"])
  new.packages <- list.of.packages[new.packages.log]
  if (length(new.packages)) install.packages(new.packages)
  new.packages.log = !(list.of.packages %in% installed.packages()[,"Package"])
  new.packages <- list.of.packages[new.packages.log]
  if (length(new.packages)) {
    setRepositories(graphics = F, ind = 1:8)
    install.packages(new.packages)
  }
  lapply(list.of.packages, library, character.only = T)
}
greppend <- function(pattern, 
                     x,
                     output = x, 
                     output.col = colnames(output)) {
  pb = progressBar(min = 0, max = length(pattern))
  new.table = data.frame()
  for (variable in 1:length(pattern)) {
    pos = grep(pattern = pattern[variable], x = x)
    new.table = rbind.data.frame(new.table, output[pos, output.col])
    setTxtProgressBar(pb, value = variable)
  }
  close(pb)
  return(new.table)
}
plot.freq <- function(x, 
                      y = NULL, 
                      outliers.rm = F, 
                      x.lab = '', y.lab = '', fill = '', 
                      x.name = '', y.name = '', angle = 0, nbreaks = 10,
                      ...) {
  naToZero <- function(x) {
    x[is.na(x)] = 0
    return(x)
  }
  x = naToZero(x)
  if (outliers.rm) {
    x.old.min = min(x)
    x.old.max = max(x)
    Q1 = quantile(x, .25, na.rm = T)
    Q3 = quantile(x, .75, na.rm = T)
    IQ = Q3 - Q1
    lower.outer.fence = Q1 - 3*IQ
    upper.outer.fence = Q3 + 3*IQ
    x = x[x > lower.outer.fence & x < upper.outer.fence]
  }
  x.range = max(x) - min(x)
  x.breaks = seq(from = min(x), to = max(x), by = (x.range/nbreaks))
  if (outliers.rm) {
    x.breaks = c(x.old.min, x.breaks[2:nbreaks], x.old.max)
  }
  x.cut = cut(x, breaks = x.breaks, include.lowest = T)
  x.table = table(x.cut) / length(x) * 100
  if (!is.null(y)) {
    y = naToZero(y)
    y.cut = cut(y, breaks = x.breaks, include.lowest = T)
    y.table = table(y.cut) / length(y) * 100
    xy.table = rbind(x.table, y.table)
    xy.df = as.data.frame(xy.table)
    f = merge(stack(xy.df), stack(as.data.frame(t(xy.df))), by = 'values')
    if (x.name != '') {
      f[, 3] = gsub(pattern = 'x.table', replacement = x.name, x = f[, 3])
    }
    if (y.name != '') {
      f[, 3] = gsub(pattern = 'y.table', replacement = y.name, x = f[, 3])
    }
    ggplot(f, aes(x = f[, 2], y = f[, 1], fill = f[, 3])) +
      geom_bar(stat = "identity", position = "dodge") +
      theme(axis.text.x = element_text(angle = angle)) +
      labs(fill = fill, x = x.lab, y = y.lab)
    
  } else {
    x.df = as.data.frame(x.table)
    x.df$names = rownames(x.df)
    ggplot(x.df, aes(x = x.cut, y = x.df[, 2])) +
      geom_bar(stat = 'identity') +
      theme(axis.text.x = element_text(angle = angle)) +
      labs(x = x.lab, y = y.lab)
  }
}
rmStrings <- function(x) {
  if (is.null(dim(x))) {
    y = as.numeric(x)
    names(y) = names(x)
    return(y)
  } else {
    y = apply(X = x, MARGIN = 2, FUN = as.numeric)
    rownames(y) = rownames(x)
    return(y)
  }
}
naToZero <- function(x) {
  x[is.na(x)] = 0
  return(x)
}

################################
forceLibrary(list.of.packages = c('pbmcapply', 'ggplot2', 'biomaRt'))

################################
setwd('/share/analysis/hecatos/juantxo/Score/output/Output_RunJUN2019_MV')
trc_unt_002_1 = read.table('UNTR_002_1_TRCscore.txt', stringsAsFactors = F)

setwd('/share/analysis/hecatos/juantxo/Score/output/Output_Run_mrna_SEPT2019/')
setwd('V3/output/UNTR/TRCscore/')
files = list.files()
untr_trc.files = subset(files, grepl('UNTR', files))

# Fuse all tables by summing rowsxcolumns in common, and appending the rest
pb = progressBar(min = 0, max = length(untr_trc.files))
for (f in untr_trc.files) {
  f.table = read.table(file = f, stringsAsFactors = F)
  f.table = rmStrings(f.table)
  f.table = naToZero(f.table)
  if (f == untr_trc.files[1]) {
    i = 1
    all.sum = f.table
  } else {
    rows.n = rownames(f.table)[!(rownames(f.table) %in% rownames(all.sum))]
    cols.n = colnames(f.table)[!(colnames(f.table) %in% colnames(all.sum))]
    pb2 = progressBar(min = 0, max = length(rows.n))
    o = 1
    for (r in rows.n) { # Add rows that don't exist in the fused DF
      all.sum = rbind(all.sum, 0)
      rownames(all.sum)[length(rownames(all.sum))] = r
      setTxtProgressBar(pb = pb2, value = o, title = 'Adding new rows...')
      o = o + 1
    }
    close(pb2)
    for (c in cols.n) { # Add cols that don't exist in the fused DF
      all.sum = cbind(all.sum, 0)
      colnames(all.sum)[length(colnames(all.sum))] = c
    }
    cols = colnames(f.table)
    for (r in rownames(f.table)) { # Sum values of rows & cols in common
      all.sum[r, cols] = all.sum[r, cols] + f.table[r, cols]
    }
  }
  setTxtProgressBar(pb = pb, value = i)
  i = i + 1
}
close(pb)
setwd('/share/analysis/hecatos/juantxo/score_protein_analysis/UNTR')
write.table('all_sum.tsv')

setwd('/share/analysis/hecatos/juantxo/score_protein_analysis/UNTR/')
all.sum = read.table(file = 'all_sum.tsv')
all.avg = all.sum / length(untr_trc.files)


setwd("/share/analysis/hecatos/juantxo/score_protein_analysis/UNTR/shift-FALSE")
proteins.table = read.table('./TRC/most_correlated_proteins_trc.tsv')

mart.human = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                     dataset = 'hsapiens_gene_ensembl') 
# View(searchAttributes(mart = mart.human))
transcr.table = getBM(attributes = c('ensembl_transcript_id', 
                                          'uniprot_gn_id'),
                           filters = 'uniprot_gn_id',
                           values = list(proteins.table$protlist), 
                           mart = mart.human)
form_subs.table = greppend(pattern = transcr.table$ensembl_transcript_id, 
                           x = rownames(all.avg), 
                           output = all.avg)

plot.freq(x = all.avg$TRC,
          outliers.rm = T, 
          x.lab = 'Sum of TPM (miRNA)', 
          y.lab = '% of transcripts', 
          fill = 'Groups', 
          x.name = 'All transcripts', 
          y.name = 'TRC best proteins', angle = 90, nbreaks = 100)

x = naToZero(trc_unt_002_1$miRNA_TPM_avg)
x = x[x > 0]

setwd('../V3/')
new_trc = read.table('UNTR_002_1_TRCscore.txt', stringsAsFactors = F)
new_trc_pos = new_trc[!is.na(new_trc$int_miRNA_sum), ]
table(new_trc_pos$TRC < 0) / length(new_trc_pos$TRC)

new_trc$miRNA_CPM_mean = new_trc$sum_miRNA_TPM / new_trc$X.miRNAs
summary(new_trc$miRNA_CPM_mean)

new_trc_pos = new_trc[!is.na(new_trc$MITcount_sum), ]
summary(new_trc_pos$TRC)

new_trc$MITscore_mean = as.numeric(new_trc$sum_MITscore) / new_trc$X.miRNAs
new_trc$MITcount_mean = as.numeric(new_trc$count) / new_trc$X.miRNAs

all.avg.pos$MITcount_sum[all.avg.pos$MITcount_sum == 0] = NA
all.avg.pos = all.avg[!is.na(all.avg.pos$MITcount_sum), ]

trc_unt_002_1$MIT_avg = trc_unt_002_1$MITcount_sum / trc_unt_002_1$sum_miRNA_TPM
new_trc$MIT_avg = as.numeric(new_trc$sum_MITscore) / new_trc$X.miRNAs
