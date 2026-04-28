

setwd("/share/script/hecatos/juantxo/analysis_trc")
source('functions.R')
forceLibrary(c('dplyr', 'tibble'))
setwd('/share/analysis/hecatos/juantxo/Score/output/')
comps = list.dirs(recursive = F)
for (compound in comps) {
  setwd('/share/analysis/hecatos/juantxo/Score/output/')
  setwd(compound)
  dir_or = list.files(pattern = '03-09.*UTC') %>% sort(decreasing = T) %>% .[1] 
  if (is.na(dir_or)) {next()}
  setwd(dir_or)
  setwd('TRCscore/')
  files_trt = list.files(recursive = F)
  factor_miRNAs = seq(0.01, 1, length.out = 100)
  new_dir = paste0('../../', gsub(dir_or, pattern = 'UTC', replacement = ''), 
                   '_factor_0_1')
  if (dir.exists(new_dir) & length(list.files(new_dir)) > 0) {
    setwd(new_dir)
    
    if (length(dir(pattern = 'All_trt.rds')) == 0) {
      saveRDS(mergeFilesRds('.rds', row_names = T, all = T, progr_bar = F), 
              'All_trt.rds')
      file.remove(dir(pattern = 'TrT.rds|txt.rds'))
    }
    next()
  }
  forceSetWd(new_dir)
  p = progress_estimated(length(files_trt))
  for (file_trt in files_trt) {
    trt = read.table(paste0('../', dir_or, '/TRCscore/', file_trt))
    for (factor_miRNA in factor_miRNAs) {
      colu_name = paste0('TrT_', factor_miRNA, '_', file_trt)
      trt = trt %>% 
        rownames_to_column() %>%
        mutate(!!colu_name := targetRNA_TPM - (factor_miRNA*int_miRNA_sum)) %>% 
        column_to_rownames()
      trt[, colu_name][trt[, colu_name] < 0] = 0
      trt[, colu_name][is.na(trt[, colu_name])] = trt$targetRNA_TPM[is.na(trt[, colu_name])]
    }
    file_name = file_trt %>% gsub(pattern = 'TRCscore.txt', replacement = 'TrT') %>% 
      paste0('.rds')
    saveRDS(object = trt, file = file_name)
    p$tick()$print()
  }
  saveRDS(mergeFilesRds('.rds', row_names = T, all = T, progr_bar = F), 'All_trt.rds')
  file.remove(dir(pattern = 'TrT.rds|txt.rds'))
}
