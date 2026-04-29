#############

library(ggplot2)

getHalf <- function(x, order = 1) {
  y = NULL
  if (is.null(dim(x))) {
    x.length = length(x)
    x.vector = seq(from = order, to = x.length, by = 2)
    for (i in x.vector) {
      y = c(y, x[i])
    }
  } else {
    x.nrow = nrow(x)
    x.vector = seq(from = order, to = x.nrow, by = 2)
    for (i in x.vector) {
      y = rbind(y, x[i,])
    }
  }
  return(y)
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
avgRows <- function(x) {
  y = apply(X = x, MARGIN = 1, FUN = mean, na.rm = T)
}
plot.freq <- function(x, 
                      y = NULL, 
                      outliers.rm = F, 
                      x.lab = '', y.lab = '', fill = '', 
                      x.name = '', y.name = '',
                      ...) {
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
    xy.df = as.data.frame(xy.table)
    f = merge(stack(xy.df), stack(as.data.frame(t(xy.df))), by = 'values')
    f[, 3] = gsub(pattern = 'x.table', replacement = x.name, x = f[, 3])
    f[, 3] = gsub(pattern = 'y.table', replacement = y.name, x = f[, 3])
    ggplot(f, aes(x = f[, 2], y = f[, 1], fill = f[, 3])) +
      geom_bar(stat = "identity", position = "dodge") +
      theme(axis.text.x = element_text(angle = 0)) +
      labs(fill = 'Group', x = x.lab, y = y.lab)
    
    # plt = barplot(xy.table, beside = T, xaxt = "n", ...)
    # text(colMeans(plt), 
    #      par("usr")[3], 
    #      labels = names(x.table), 
    #      srt = 45, 
    #      adj = c(1.1,1.1), 
    #      xpd = TRUE, 
    #      cex = 0.6) 
  } else {
    barplot(x.table, ...)
  }
}


##### OPENING FILES #####

setwd(paste0('/share/analysis/hecatos/juantxo/score_protein_analysis/UNTR/', 
             'shift-FALSE/targetRNA_TPM/data_tables/'))
tpm_prot.table = read.table(file = 'minimumexpressedsamples3.tsv',
                       stringsAsFactors = F,
                       sep = '\t')
setwd("/share/analysis/hecatos/juantxo/score_protein_analysis/UNTR/shift-FALSE")
proteins.table = read.table('./TRC/most_correlated_proteins_trc.tsv')

tpm.table = getHalf(x = tpm_prot.table, order = 2)
rownames(tpm.table) = gsub(pattern = 'TRC', 
                           replacement = 'TPM', 
                           x = rownames(tpm.table))

subset.tpm = greppend(pattern = proteins.table$protlist, 
                      x = rownames(tpm.table), 
                      output = tpm.table)
subset.mean_tpm = avgRows(subset.tpm)
summary(subset.mean_tpm)
tpm_global.mean = avgRows(tpm.table)
summary(tpm_global.mean)

plot.freq(x = tpm_global.mean, 
          y = subset.mean_tpm, 
          outliers.rm = T, 
          x.lab = 'Mean TPM expression', 
          y.lab = '% of transcripts', 
          fill = 'Groups', 
          x.name = 'Global', 
          y.name = 'Subset')

