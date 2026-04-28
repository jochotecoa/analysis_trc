# setRepositories()
# install.packages('biomaRt')
library('biomaRt')

compare <- function(dataset1, dataset2, column) {
  difference = summary(dataset2[,column]) - summary(dataset1[,column])
  return(abs(difference))
}

freq.dist <- function(column, type = 'sd', from = '', to = '', rot_angle = 45, main ='', xlab = '', ylab = '') {
  sd = sd(column, na.rm = T)
  mn = mean(column, na.rm = T)
  mxm = max(column, na.rm = T)
  mnm = min(column, na.rm = T)
  if (type == 'sd') {
    bye = (2*sd + mn - (mn - 2*sd))/10
    brks = seq(from = mnm, to = mxm, by = bye)
  }
  if (type == 'minmax') {
    bye = (mxm - mnm)/10
    brks = seq(from = mnm, to = mxm, by = bye)
  }
  if (type == 'custom') {
    bye = (to - from)/10
    brks = seq(from = from, to = to, by = bye)
  }
  
  if (type == 'sd') {
    if (mnm < (mn - 2*sd)) {
      brks = c(mnm, brks)
    }
    if (mxm > (mn + 2*sd)) {
      brks = c(brks, mxm)
    }
  }
  ct = cut(column, breaks = brks)
  tbl = table(ct)
  plt <- barplot(tbl, col='steelblue', xaxt="n", main = main, xlab = xlab, ylab = ylab)
  text(plt, par("usr")[3], labels = names(tbl), srt = rot_angle, adj = c(1.1,1.1), xpd = TRUE, cex=0.6) 
  return(tbl)
}

setwd('/share/analysis/hecatos/juantxo/score_protein_analysis/')
datafile = read.table('UNTR_correlation_results_TRC_protein_1onseveral_until072_minimumexpressedsamples-12_shift-TRUE_timeps-shifted-002.tsv', header = T, stringsAsFactors = F)
shifts = as.numeric(names(table(datafile$shiftlist)))

for (s in shifts) {
  assign(x = paste0('timeshift',s), datafile$protlist[grepl(s, datafile$shiftlist)])
}
timeshift.analyzed = timeshift0
name.timeshift.analyzed = 'Time-shift 0'

halflives = read.csv('nature10098-s5_halflife_proteins.csv', header = T, sep = ';', stringsAsFactors = F)

marts = listMarts()
mart.human = useMart(biomart = marts$biomart[1])
datasets = listDatasets(mart = mart.human)
mart.human = useMart(biomart = marts$biomart[1], dataset = datasets$dataset[73])
filt = searchFilters(mart = mart.human, pattern = 'uniprot')
attribute = searchAttributes(mart = mart.human, pattern = 'Mouse protein')
ensembl.prot.ids.mouse = getBM(attributes = attribute$name,
                               filters = "uniprot_gn_id",
                               values = list(timeshift.analyzed), mart = mart.human)

datasets.mouse = datasets[grep('Mouse', datasets$description),]
mart.mouse = useMart(biomart = marts$biomart[1], dataset = datasets.mouse$dataset[2])
filts = searchFilters(mart = mart.mouse, pattern = 'Protein stable')
filt = filts$name[1]
attributes = searchAttributes(mart = mart.mouse, pattern = 'UniProtKB Gene Name ID')
attribute = attributes$name[1]
uniprot.ids.mouse.timeshift = getBM(attributes = attribute,
                                    filters = filt,
                                    values = ensembl.prot.ids.mouse, mart = mart.mouse)

index.timeshift = unlist(lapply(unlist(uniprot.ids.mouse.timeshift), function(x) grep(pattern = x, halflives$Uniprot.IDs)))
index.timeshift = unique(index.timeshift)
halflives[,11:ncol(halflives)] = as.data.frame(sapply(halflives[,11:ncol(halflives)], function(x) gsub(pattern = ',', replacement = '.', x)))
halflives[,11:ncol(halflives)] = sapply(halflives[,11:ncol(halflives)], function(x) as.numeric(x))
timeshift.halflives = halflives[index.timeshift,]

mx = mn = 0
lsdff = c()
for (col in 11:ncol(halflives)) {
  mx.old = mx
  mn.old = mn
  # cmp = compare(dataset1 = halflives, dataset2 = timeshift.halflives, col)
  cmp = mean(timeshift.halflives[,col], na.rm = T) - mean(halflives[,col], na.rm = T)
  lsdff = c(lsdff, cmp)
  mx = max(c(mx, cmp), na.rm = T)
  if (mx.old != mx) {
    col.max = col
  }
  mn = min(c(mn, cmp), na.rm = T)
  if (mn.old != mn) {
    col.min = col
  }
  
}
setwd('plots/')
png(filename = paste(colnames(halflives)[col.max], 'all-proteins_frequency-distribution_UNTR', sep = '_'))
freq.dist(halflives[,col.max], type = 'custom', from = min(halflives[,col.max], na.rm = T), to = max(halflives[,col.max], na.rm = T), main = 'All proteins', xlab = colnames(halflives)[col.max], ylab = '# of proteins')
dev.off()
png(filename = paste(colnames(halflives)[col.max], name.timeshift.analyzed, 'frequency-distribution_UNTR', sep = '_'))
freq.dist(timeshift.halflives[,col.max], type = 'custom', from = min(halflives[,col.max], na.rm = T), to = max(halflives[,col.max], na.rm = T), main = paste(name.timeshift.analyzed, 'proteins'), xlab = colnames(halflives)[col.max], ylab = '# of proteins')
dev.off()

png(filename = paste(colnames(halflives)[col.min], 'all-proteins_frequency-distribution_UNTR', sep = '_'))
freq.dist(halflives[,col.min], type = 'custom', from = min(halflives[,col.min], na.rm = T), to = max(halflives[,col.min], na.rm = T), main = 'All proteins', xlab = colnames(halflives)[col.min], ylab = '# of proteins')
dev.off()
png(filename = paste(colnames(halflives)[col.min], name.timeshift.analyzed, 'frequency-distribution_UNTR', sep = '_'))
freq.dist(timeshift.halflives[,col.min], type = 'custom', from = min(halflives[,col.min], na.rm = T), to = max(halflives[,col.min], na.rm = T), main = paste(name.timeshift.analyzed, 'proteins'), xlab = colnames(halflives)[col.min], ylab = '# of proteins')
dev.off()

(table(timeshift.halflives[,23] < median(halflives[,23], na.rm = T)))[2] / nrow(timeshift.halflives)
mean.min.timeshift = median(timeshift.halflives[,23], na.rm = T)


