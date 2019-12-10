source('/share/script/hecatos/juantxo/analysis_trc/functions.R')
forceLibrary(c('dplyr', 'tibble'))


setwd('/ngs-data/data/hecatos/Cardiac/Con_UNTR/Protein/Proteomics_Analyses_Cardiac_UNTR_GeneData/')
prot_file = read.table('Hecatos_Cardiac_Px_Untreated_pre-processed_renamed.txt', sep = '\t', header = T) %>%
  cleanProtIds() %>% 
  select(-Row.Names) %>% 
  remove_rownames() %>%
  column_to_rownames('uniprot_gn')

prot_log = prot_file %>% log2()

# Get the median of medians
expressed <- apply(prot_log, 1, function(r) all(!is.na(r)))
median_of_medians = prot_log[expressed, ] %>%
  apply(2, median, na.rm = T) %>%
  median()

# timepoints with only 1 NA will be subst by the average of the 2 others
n_corrected = 0
for (row in rownames(prot_log)) {
  tps = seq(1,ncol(prot_log), by = 3)
  for (tp in tps) {
    cols = tp:(tp + 2)
    prot_tp = prot_log[row, cols]
    n_NAs = sum(is.na(prot_tp))
    if (n_NAs == 1) {
      n_corrected = n_corrected + 1
      prot_tp = as.numeric(prot_tp)
      prot_tp[is.na(prot_tp)] = mean(prot_tp, na.rm = T)
      prot_log[row, cols] = prot_tp
    }
  }
}
print(n_corrected)
# Timepoints with only 1 or none expression will be taken away
prot_log_norm = prot_log %>% na.omit()

# Shift the dists so that all of them have the same median
for (col in colnames(prot_log_norm)) {
  fact = median_of_medians - median(prot_log_norm[, col], na.rm = T) 
  prot_log_norm[, col] = prot_log_norm[, col] + fact
}

summary(prot_log_norm)
# prot_log_norm = apply(prot_log_norm, 2, as.numeric) %>% as.data.frame()
# summary(prot_log_norm)
prot_log_norm$UNTR_The_002_1 %>% density() %>% plot()
prot_log_norm$UNTR_The_002_1 %>% qqnorm(main = sprintf('Shapiro test: %s', format(shapiro.test(.)$p.value, scientific=T, digits = 3)))
prot_log_norm$UNTR_The_002_1 %>% qqline(col = 'red')
# p.values = NULL
# prot_log_norm$p.value = NA

cols = seq(4, ncol(prot_log_norm), 3)
for (col in cols) {
  prot_log_norm[, paste0('p.value_', colnames(prot_log_norm)[col])] = NA
  for (row in rownames(prot_log_norm)) {
    X = prot_log_norm[row, 1:3] %>% as.numeric()
    Y = prot_log_norm[row, col:(col+2)] %>% as.numeric()
    t_test = t.test(x = X, y = Y, paired = T)
    prot_log_norm[row, ncol(prot_log_norm)] = t_test$p.value
  }
  
}

forceSetWd('/share/analysis/hecatos/juantxo/Score/analysis/DEPs')
write.table(x = prot_log_norm, file = 'UNTR_protein_log2_pvalues.tsv', sep = '\t')       
