plot.freq <- function(x, y = NULL, outliers.rm = F, x.lab = '', y.lab = '', 
                      fill = '', x.name = '', y.name = '', angle = 0, 
                      nbreaks = 10, ...) {
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

setPaste = function(...){
  setwd(paste0(...))
}

setwd('/share/analysis/hecatos/juantxo/Score/Output_Run_mrna_SEPT2019/V3/')
setwd('Analysis/UNTR/shift_FALSE/TRC/_Tue_Sep_10_15:18:11_2019_/')
setwd('minimum_expressed_samples12/')
corrs_TRC = read.table(file = 'correlation_results_between_TRC_values_and_protein_expression.tsv', sep = '\t', header = T)
TRC = corrs_TRC[seq(2,nrow(corrs_TRC), by = 2), ]
prot_expr = corrs_TRC[seq(1,nrow(corrs_TRC), by = 2), ]
setwd('../../../targetRNA_TPM/_Tue_Sep_10_15:18:11_2019_/')
setwd('minimum_expressed_samples12/')
corrs_TPM = read.table(file = 'correlation_results_between_targetRNA_TPM_values_and_protein_expression.tsv', header = T)

plot.freq(x = corrs_TPM$corlist, y = corrs_TRC$corlist, outliers.rm = F, 
          x.name = 'TPM', y.name = 'TRC', 
          x.lab = 'Correlation values', y.lab = '% of proteins', angle = 45)
