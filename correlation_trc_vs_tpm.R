##### Functions #####
plot.freq <- function(x, y = NULL, outliers.rm = F, x.lab = '', y.lab = '', 
                      fill = '', x.name = '', y.name = '', angle = 0, 
                      nbreaks = 10, ...) {
  library(ggplot2)
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
      labs(x = x.lab, y = y.lab, ...)
  }
}

##### Analysis #####
setwd('/share/analysis/hecatos/juantxo/Score/Output_Run_mrna_SEPT2019/V3/')
setwd('TRCscore/')
trc.table = read.table('DF2_002_1_TRCscore.txt', 
                       header = T, stringsAsFactors = F)
trc.values = trc.table$TRC
tpm.values = trc.table$targetRNA_TPM

trc_tpm.cor = cor(x = trc.values, y = tpm.values)

setwd('../Analysis/UNTR/shift_FALSE/TRC/_Wed_Sep_11_14:13:37_2019_/')
setwd('minimum_expressed_samples12/')
trc_tp.table = read.table(file = 'TRC_values_and_protein_expression.tsv', 
                           header = T, stringsAsFactors = F, sep = '\t')
trc_tp.rows = grep(pattern = 'TRC', x = rownames(trc_tp.table))
trc_tp.values = trc_tp.table[trc_tp.rows, ]
summary(trc_tp.values)

setwd('../../../targetRNA_TPM/_Wed_Sep_11_14:13:37_2019_/')
setwd('minimum_expressed_samples12/')
tpm_tp.table = read.table('targetRNA_TPM_values_and_protein_expression.tsv', 
                          header = T, stringsAsFactors = F, sep = '\t')
tpm_tp.rows = grep(pattern = 'TRC', x = rownames(tpm_tp.table))
tpm_tp.values = tpm_tp.table[tpm_tp.rows, ]
rownames(tpm_tp.values) = gsub(pattern = 'TRC', replacement = 'TPM', 
                               x = rownames(tpm_tp.values))
summary(tpm_tp.values)

cors = cor(x = t(trc_tp.values), y = t(tpm_tp.values), use = "all.obs")
cors_temp = cors
cors = diag(cors)
names(cors) = rownames(cors_temp)
names(cors) = gsub(pattern = ' TRC values', replacement = '', x = names(cors))
summary(cors)
cors = cors[!is.na(cors)]

plot.freq(x = cors, y.lab = '# of proteins', x.lab = 'Correlation values', 
          title = 'Correlation between TRC and TPM values')

trc_effect = cors[cors < 0.9]
summary(trc_effect)
trc_effect.df = data.frame(trc_effect)
colnames(trc_effect.df) = 'correlation_value'

setwd('Analysis/UNTR/')
write.table(x = trc_effect.df, file = '20190916_low_correlated_TPMvsTRC.tab', 
            sep = '\t')
