source('/share/script/hecatos/juantxo/analysis_trc/functions.R')
countToTpm <- function(counts, effLen) {
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  return(exp(rate - denom + log(1e6)))
}

TrcToCount = function(TRC, effLen, counts) {
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  print(denom)
  counts = exp(log(TRC) + denom - log(10^6) + log(effLen))
}

# GEt trc values
setwd('/share/analysis/hecatos/juantxo/Score/output/UNTR/')
setwd('2019-11-11_12:04:38_UTC/TRCscore/')
trc = mergeFiles(files_patt = 'TRCscore', row_names = T, all = T)
trc_values = trc %>%
  filter(grepl('ENST', rowname)) %>%
  select(contains('TRC_'))

# Get salmon values
setwd('/share/analysis/hecatos/juantxo/Score/input/RNA_Salmon/UNTR/')
salmon = mergeFiles(all = T) %>%
  filter(grepl('ENST', salmon$Name))
count_values = salmon %>%
  select(contains('Reads'))
effLength_values = salmon %>%
  select(contains('EffectiveLength'))

# Get trc_count values
for (i in 1:ncol(trc_values)) {
  trc_values_test = trc_values[, i, drop = F] %>%
    cbind.data.frame(., trc$rowname)
  count_values_test = count_values[, i]
  effLength_values_test = effLength_values[, i] 
  
  
  salmon_test = cbind.data.frame(effLength_values_test, 
                                 count_values_test, 
                                 salmon$Name)
  salmon_test$`salmon$Name` = substr(salmon_test$`salmon$Name`, 1, 15)
  allvalues_test = merge.data.frame(x = salmon_test, 
                                    y = trc_values_test, 
                                    by.x = 'salmon$Name', 
                                    by.y = 'trc$rowname') 
  all_values_trc_test = allvalues_test %>%
    mutate(trc_counts = TrcToCount(
      TRC = allvalues_test[, 4],
      effLen = allvalues_test$effLength_values_test, 
      counts = allvalues_test$count_values_test)) 
  
  
  # Get trc_count per gene values
  all_values_genes_test = transcrToGene(table = all_values_trc_test, 
                                        aggregate = T, prot_cod = T)
  trc_counts__genes = all_values_genes_test[, c(1,5)]
  colnames(trc_counts__genes)[2] = paste0(
    'counts_', colnames(all_values_genes_test)[4])
  
  if (i == 1) {
    all__trc_counts__genes = trc_counts__genes
  } else {
    all__trc_counts__genes = merge.data.frame(
      all__trc_counts__genes, 
      trc_counts__genes, 
      by = 'ensembl_gene_id', 
      all = T)
  }
  
}

setwd('/share/analysis/hecatos/juantxo/Score/analysis/')
dir.create('TRC_counts')
setwd('TRC_counts/')
write.table(x = all__trc_counts__genes, 
            file = 'counts_TRC_UNTR_2019-11-11_12:04:38_UTC.tsv', 
            sep = '\t', row.names = F, col.names = T)
# Do DEGs