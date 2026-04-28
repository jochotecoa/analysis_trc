#' Time-shift Analysis of Correlation between TRC and Protein Expression
#' 
#' This script analyzes protein half-lives in relation to observed time-shifts
#' in transcript-protein correlation.
#'
#' Inputs:
#' - Correlation results with shift information
#' - Protein half-life data (Nature 2011)
#'
#' Outputs:
#' - Comparison plots of half-lives across different time-shifts

source("../utils.R")
forceLibrary(c('biomaRt', 'plyr', 'ggplot2', 'scales'))

# INPUT VARIABLES
comp = 'UNTR'

# Analysis
datafile.name = paste0(comp, 
                       'UNTR_correlation_results_TRC_protein_1onseveral_until072_minimumexpressedsamples-12_shift-TRUE_timeps-shifted-002.tsv')
datafile = read.table(datafile.name, header = T, stringsAsFactors = F)
shifts = as.numeric(names(table(datafile$shiftlist)))

# Open half-life data and format it
halflives = read.csv('nature10098-s5_halflife_proteins.csv', 
                     header = T, 
                     sep = ';', 
                     stringsAsFactors = F)
halflives[,10:ncol(halflives)] = as.data.frame(sapply(halflives[,10:ncol(halflives)], 
                                                      function(x) gsub(pattern = ',', 
                                                                       replacement = '.',
                                                                       x)))
halflives[,10:ncol(halflives)] = sapply(halflives[,10:ncol(halflives)], 
                                        function(x) as.numeric(x))

time = timestamp(prefix = paste('plots time-shift_analysisv6.R ', comp))
setwd('/share/analysis/hecatos/juantxo/score_protein_analysis/plots/plots_tableparams_nproteins_time-shift/')
dir.create(time)

for (s in shifts) {
  assign(x = paste0('timeshift',s), 
         datafile$protlist[grepl(s, datafile$shiftlist)])
  timeshift.analyzed = get(x = paste0('timeshift',s))
  name.timeshift.analyzed = paste('Time-shift', s)
  
  # May give errors while connecting to ENSEMBL : Error in bmRequest(request = request, ssl.verifypeer = (...)
  mart.human = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                       dataset = 'hsapiens_gene_ensembl') 
  ensembl.prot.ids.mouse = getBM(attributes = 'mmusculus_homolog_ensembl_peptide',
                                 filters = "uniprot_gn_id",
                                 values = list(timeshift.analyzed), 
                                 mart = mart.human)
  
  mart.mouse = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                       dataset = 'mmusculus_gene_ensembl')
  uniprot.ids.mouse.timeshift = getBM(attributes = 'uniprot_gn_id',
                                      filters = 'ensembl_peptide_id',
                                      values = ensembl.prot.ids.mouse, 
                                      mart = mart.mouse)
  
  uniprot.ids.mouse.timeshift = unlist(uniprot.ids.mouse.timeshift)
  index.timeshift = lapply(uniprot.ids.mouse.timeshift, function(x) grep(pattern = x, halflives$Uniprot.IDs))
  index.timeshift = unlist(index.timeshift)
  index.timeshift = unique(index.timeshift)
  timeshift.halflives = halflives[index.timeshift,]
  assign(x = paste0('timeshift', s, '.halflives'), value = timeshift.halflives)
  
  mx = mn = 0
  lsdff = c()
  for (col in 11:ncol(halflives)) {
    mx.old = mx
    mn.old = mn
    cmp = ((mean(timeshift.halflives[,col], na.rm = T) / mean(halflives[,col], na.rm = T)) - 1)*100
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
  setwd('/share/analysis/hecatos/juantxo/score_protein_analysis/plots/plots_tableparams_nproteins_time-shift/')
  setwd(time)
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
  
}
# get.double.summaries <- function(table) {
#   a = c()
#   index.names = c()
#   for (c in 1:ncol(table)) {
#     summ = summary(table[,c])
#     if (is.double(summ)) {
#       index.names = c(index.names, colnames(table)[c])
#       a = rbind(a, summ)
#     }
#   }
#   a = as.data.frame(a)
#   rownames(a) = index.names
#   return(a)
# }
# 
# 
# 
# 
# yvalues = ls()[grep(pattern = '.halflives', x = ls())][-1]
# for (col in 1:ncol(halflives)) {
#   if (is.double(x = halflives[,col])) {
#     median.diffs = c()
#     for (ts in 0:3) {
#       timesh = get(x = paste0('timeshift',ts,'.halflives'))
#       timesh.diff = (median(timesh[,col], na.rm = T) / median(halflives[,col], na.rm = T)) - 1
#       median.diffs = c(median.diffs, timesh.diff)
#     }
#     dtf = data.frame(x = c(0:3), y = median.diffs)
#     png(filename = paste0('plot_',colnames(halflives)[col], '_timeshifts_UNTR'))
#     plot = ggplot(data = dtf, aes(x, y)) +
#       geom_bar(stat = "identity", aes(fill = x)) + 
#       geom_text(aes(label = paste(round(y * 100, digits = 2), "%"),
#                     vjust = ifelse(y >= 0, 0, 1))) +
#       scale_y_continuous("% difference vs Total", labels = percent_format()) +
#       scale_x_continuous('Timeshifts') +
#       ggtitle(colnames(halflives)[col],)
#     print(plot)
#     dev.off()
#   }
# }

