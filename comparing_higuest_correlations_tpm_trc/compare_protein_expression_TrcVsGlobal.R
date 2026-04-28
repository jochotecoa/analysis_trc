read.proteindata <- function(proteomicsfilepath, 
                                        comp = NA, 
                                        samples = NA) {
  proteinvaluesfile.colnames = read.table(proteomicsfilepath, nrows = 1, 
                                          stringsAsFactors = F)
  proteinvaluesfile = read.table(proteomicsfilepath, skip = 1, 
                                 stringsAsFactors = F, 
                                 header = F, 
                                 sep = ' ', 
                                 fill = T)
  if (ncol(proteinvaluesfile) != ncol(proteinvaluesfile.colnames)) {
    proteinvaluesfile = read.table(proteomicsfilepath, 
                                   skip = 1, 
                                   stringsAsFactors = F, 
                                   header = F, 
                                   sep = '\t', 
                                   fill = T)
  }
  proteinvaluesfile = proteinvaluesfile[grep('.*\\|.*\\|.*', 
                                             proteinvaluesfile[,1]),]
  proteinvaluesfile = proteinvaluesfile[!grepl(pattern = ':', 
                                               proteinvaluesfile$V2),]
  proteinvaluesfile[,2:ncol(proteinvaluesfile)] = sapply(proteinvaluesfile[,2:ncol(proteinvaluesfile)], 
                                                         as.double)
  proteinvaluesfile = proteinvaluesfile[,c(T,!is.na(sapply(proteinvaluesfile[, -1],
                                                           mean, 
                                                           na.rm = T)))]
  proteinvaluesfile[,1] = as.character(proteinvaluesfile[,1])
  colnames(proteinvaluesfile) = proteinvaluesfile.colnames[,1:ncol(proteinvaluesfile)]
  if (!is.na(comp)) {
    protvalcomp = proteinvaluesfile[,c(1,grep(toupper(comp), 
                                              colnames(proteinvaluesfile)))] # Get all columns with our compound
  } else {
    protvalcomp = proteinvaluesfile
  }
  if (!is.na(samples)) {
    selectedtimepointcolumns = c(1,grep(pattern = paste(samples, 
                                                        collapse = '|'), 
                                        x = colnames(protvalcomp))) # Select only columns with timepoints of interest
    protvalcompsel = protvalcomp[rowSums(!is.na(protvalcomp[,selectedtimepointcolumns])) > 2, 
                                 selectedtimepointcolumns] # Filter proteins with =< 2 expressions
  } else {
    protvalcompsel = protvalcomp
  }
  return(protvalcompsel)
}
subsetIn <- function(x, x.subset = x, subset) {
  b = NULL
  for (value in subset) {
    a = grep(pattern = value, x = x.subset)
    if (!is.na(a)) {
      b = c(b, a)
    }
  }
  c = x[b, ]
}
avgRows <- function(x) {
  y = apply(X = x, MARGIN = 1, FUN = mean, na.rm = T)
}
plot.freq <- function(x, y = NULL, outliers.rm = F, ...) {
  if (outliers.rm) {
    x.old.min = min(x)
    x.old.max = max(x)
    Q1 = quantile(x, .25)
    Q3 = quantile(x, .75)
    IQ = Q3 - Q1
    lower.outer.fence = Q1 - 3*IQ
    upper.outer.fence = Q3 + 3*IQ
    x = x[x > lower.outer.fence & x < upper.outer.fence]
  }
  x.range = max(x) - min(x)
  x.breaks = seq(from = min(x), to = max(x), by = (x.range/10))
  if (outliers.rm) {
    x.breaks = c(x.old.min, x.breaks[2:10], x.old.max)
  }
  x.cut = cut(x, breaks = x.breaks)
  x.table = table(x.cut) / length(x) * 100
  if (!is.null(y)) {
    y.cut = cut(y, breaks = x.breaks)
    y.table = table(y.cut) / length(y) * 100
    xy.table = rbind(x.table, y.table)
    plt = barplot(xy.table, beside = T, xaxt = "n", ...)
    text(colMeans(plt), 
         par("usr")[3], 
         labels = names(x.table), 
         srt = 45, 
         adj = c(1.1,1.1), 
         xpd = TRUE, 
         cex = 0.6) 
  } else {
    barplot(x.table, ...)
  }
}

###########
setwd(paste0('/ngs-data/data/hecatos/Cardiac/Con_UNTR/Protein/', 
             'Proteomics_Analyses_Cardiac_UNTR_GeneData/'))
prot.expr.table = read.proteindata(paste0('Hecatos_Cardiac_Px_Untreated_',
                                          'pre-processed_renamed.txt'))
## This is a list of the proteins for which their TRC and their protein
#3 expression had the highest correlation
setwd("/share/analysis/hecatos/juantxo/score_protein_analysis/UNTR/shift-FALSE")
proteins.table = read.table('./TRC/most_correlated_proteins_trc.tsv')

sel.prot.expr.table = subsetIn(x = prot.expr.table, 
                               x.subset = prot.expr.table$Row.Names, 
                               subset = proteins.table$protlist)

subset.avg = avgRows(sel.prot.expr.table[, -1])
global.avg = avgRows(prot.expr.table[, -1])

plot.freq(x = global.avg, y = subset.avg, outliers.rm = T)

# summary(subset.avg) / summary(global.avg)
# 
# range = max(subset.avg) - min(subset.avg)
# breaks = seq(min(subset.avg), to = max(subset.avg), by = range/10)
# cut = cut(x = subset.avg, breaks = breaks)
# table = table(cut)
# table2 = rbind(table, table)
# barplot(table2, beside = T, las = 2.5)
# 
# 
# 
# rotate_x <- function(data, labels_vec, rot_angle) {
#   plt <- barplot(data, col='steelblue', xaxt="n", beside = T)
#   text(plt, par("usr")[3], labels = labels_vec, srt = rot_angle, adj = c(1.1,1.1), xpd = TRUE, cex=0.6) 
# }
# rotate_x(table2, names(table), rot_angle = 45)
# 
# x <- barplot(table(mtcars$cyl), xaxt="n")
# labs <- paste(names(table(mtcars$cyl)), "cylinders")
# text(cex=1, x=x, y=par("usr")[2], labs, xpd=TRUE, srt=45)
