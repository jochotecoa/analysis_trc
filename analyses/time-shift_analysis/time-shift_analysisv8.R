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

source("../../utils.R")
forceLibrary(c('biomaRt', 'plyr', 'ggplot2', 'scales'))

# INPUT VARIABLES
comp = if(exists("DEFAULT_COMP")) DEFAULT_COMP else 'UNTR'

# Analysis - Search for the input datafile if not explicitly provided
result_files = list.files(recursive = T, pattern = "correlation_results_TRC_protein.*\\.tsv")
if (length(result_files) > 0) {
    datafile.name = result_files[1] # Take the first one as an example or use a specific one
} else {
    datafile.name = paste0(comp, 
                       'UNTR_correlation_results_TRC_protein_1onseveral_until072_minimumexpressedsamples-12_shift-TRUE_timeps-shifted-002.tsv')
}

if (file.exists(datafile.name)) {
    datafile = read.table(datafile.name, header = T, stringsAsFactors = F)
    shifts = as.numeric(names(table(datafile$shiftlist)))
} else {
    warning("Datafile not found: ", datafile.name)
    shifts = c()
}

# Open half-life data and format it
halflife_file = 'nature10098-s5_halflife_proteins.csv'
if (file.exists(halflife_file)) {
    halflives = read.csv(halflife_file, header = T, sep = ';', stringsAsFactors = F)
} else {
    # Fallback or placeholder if file missing
    halflives = data.frame()
    warning("Half-life data file missing: ", halflife_file)
}

if (nrow(halflives) > 0) {
    halflives[,10:ncol(halflives)] = as.data.frame(sapply(halflives[,10:ncol(halflives)], 
                                                          function(x) gsub(pattern = ',', 
                                                                           replacement = '.',
                                                                           x)))
    halflives[,10:ncol(halflives)] = sapply(halflives[,10:ncol(halflives)], 
                                            function(x) as.numeric(x))
}

time_stamp = timestamp(prefix = paste('plots time-shift_analysisv6.R ', comp))
time_dir = gsub(" ", "_", gsub(":", "-", time_stamp))

# Define output base path from config if available
output_base = if(exists("BASE_ANALYSIS_PATH")) {
    file.path(BASE_ANALYSIS_PATH, 'score_protein_analysis/plots/plots_tableparams_nproteins_time-shift/', time_dir)
} else {
    file.path(getwd(), 'plots', time_dir)
}

setOrCreatewd(output_base)

# Initialize biomaRt once outside the loop
mart.human = openMart2018()
mart.mouse = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl')

for (s in shifts) {
  timeshift.analyzed = datafile$protlist[datafile$shiftlist == s]
  name.timeshift.analyzed = paste('Time-shift', s)
  
  if (length(timeshift.analyzed) == 0) next
  
  # Get homologs mapping
  ensembl.prot.ids.mouse = getBM(attributes = 'mmusculus_homolog_ensembl_peptide',
                                 filters = "uniprot_gn_id",
                                 values = timeshift.analyzed, 
                                 mart = mart.human)
  
  uniprot.ids.mouse.timeshift = getBM(attributes = 'uniprot_gn_id',
                                      filters = 'ensembl_peptide_id',
                                      values = ensembl.prot.ids.mouse$mmusculus_homolog_ensembl_peptide, 
                                      mart = mart.mouse)
  
  uniprot.ids.mouse.timeshift = unlist(uniprot.ids.mouse.timeshift)
  
  # Efficiently find matches in halflives
  matches = halflives$Uniprot.IDs %in% uniprot.ids.mouse.timeshift
  timeshift.halflives = halflives[matches, ]
  
  if (nrow(timeshift.halflives) == 0) {
      warning("No half-life data for timeshift ", s)
      next
  }

  # Analysis of difference in means
  col_indices = 11:ncol(halflives)
  means_all = colMeans(halflives[, col_indices], na.rm = T)
  means_shift = colMeans(timeshift.halflives[, col_indices], na.rm = T)
  lsdff = ((means_shift / means_all) - 1) * 100
  
  col.max = col_indices[which.max(lsdff)]
  col.min = col_indices[which.min(lsdff)]
  
  # Plots
  # Use freq.dist from utils.R
  
  create_dist_plots = function(col_idx, label) {
    col_name = colnames(halflives)[col_idx]
    
    png(filename = file.path(output_base, paste(col_name, 'all-proteins_dist', label, sep = '_')))
    freq.dist(halflives[, col_idx], type = 'custom', 
              from = min(halflives[, col_idx], na.rm = T), 
              to = max(halflives[, col_idx], na.rm = T), 
              main = 'All proteins', xlab = col_name, ylab = '# of proteins')
    dev.off()
    
    png(filename = file.path(output_base, paste(col_name, name.timeshift.analyzed, 'dist', label, sep = '_')))
    freq.dist(timeshift.halflives[, col_idx], type = 'custom', 
              from = min(halflives[, col_idx], na.rm = T), 
              to = max(halflives[, col_idx], na.rm = T), 
              main = paste(name.timeshift.analyzed, 'proteins'), 
              xlab = col_name, ylab = '# of proteins')
    dev.off()
  }
  
  create_dist_plots(col.max, "MAX")
  create_dist_plots(col.min, "MIN")
}
  
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

