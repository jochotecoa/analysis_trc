
# Functions ---------------------------------------------------------------


plotOmics = function(df, omics, transcript, ylim = NULL, ...) {
  if (!is.null(ylim)) {
    rownames(df) %>%  grep(pattern = transcript) %>% 
      df[., grep(omics, colnames(df))] %>% as.numeric() %>%  
      plot(main = transcript, ylim = ylim, ...)
  } else {
    rownames(df) %>%  grep(pattern = transcript) %>% 
      df[., grep(omics, colnames(df))] %>% as.numeric() %>%  
      plot(main = transcript, ylim = c(0, max(., na.rm = T)), ...)
  }
  abline(v = 12.5)
  for (l in seq(3.5, 21.5, 3)) abline(v = l, col = 'gray')
  par(new = TRUE)
}
theVsTox_fun <- function(df, omics, ...) {
  the_cols = df %>% dplyr::select(contains('The')) %>% 
    dplyr::select(starts_with(omics)) %>% colnames()
  tox_cols = df %>% dplyr::select(contains('Tox')) %>% 
    dplyr::select(starts_with(omics)) %>% colnames()
  
  res.df = df %>% select(the_cols, tox_cols) %>% zeroToNa() %>% log2() %>% 
    apply_2D(col.x = the_cols, col.y = tox_cols, ...)
}
setwd("/share/script/hecatos/juantxo/analysis_trc")
source('functions.R')
forceLibrary(c('dplyr', 'tibble'))


# Get data ----------------------------------------------------------------

comp = '5FU'

# Load proteomics data for this compound
setwd('/share/script/hecatos/juantxo/analysis_trc/')
source('differentially_expressed_proteins/DEPs_dose.R')

setwd('/share/analysis/hecatos/juantxo/Score/output/')
setwd(comp)
# ; setwd(comp)
list.files(pattern = 'factor_0_1') %>% sort(decreasing = T) %>% setwd()
# 
all_patt = 'All_trt|all_trt'

factor_miRNAs = 'TrT_' %>% paste0(seq(from = 0.01, to = 1, by = 0.01)) %>% paste0('_')


all_file = list.files(pattern = all_patt)
rm(trt_df)
# common_trt_df = mergeFiles(files_patt = 'TrT', row_names = T, progr_bar = F)
if (length(all_file) == 0) {
  if (sum(grepl(pattern = 'rds', x = list.files()))) {
    trt_df = mergeFilesRds(files_patt = 'TrT', row_names = T, progr_bar = F, 
                           all = T)
  } else {
    trt_df = mergeFiles(files_patt = 'TrT', row_names = T, progr_bar = F, 
                        all = T)
  }
  trt_df = trt_df %>% column_to_rownames() 
  colnames(trt_df) = colnames(trt_df) %>% 
    gsub(pattern = '_TrT.txt.rds', replacement = '')
  trt_df %>% select(starts_with('target')) %>% filterSamplesBySeqDepth() %>%
    ncol() %>% identical(ncol(select(trt_df, starts_with('target')))) %>% 
    stopifnot()
  saveRDS(object = trt_df, file = paste0(comp, '_all_trt.rds'))
  
} else {
  trt_df = readRDS(list.files(pattern = all_patt))
  if ('rowname' %in% colnames(trt_df)) {
    trt_df = trt_df %>% column_to_rownames() 
  }
}


comp_cas = 12
if (any(grepl(pattern = 'Dox', x = colnames(trt_df)))) {
  trt_df = trt_df %>% select(-contains('072_3'))
  comp_cas = 11
}
if (any(grepl(pattern = 'Epi', x = colnames(trt_df)))) {
  trt_df = trt_df %>% select(-contains('072'))
  comp_cas = 9
}

new_trt_df = merge.data.frame(x = rownames_to_column(trt_df), 
                              y = select(.data = proteomx_biom, 
                                         `Transcript stable ID`, `Gene stable ID`), 
                              by.x = 'rowname',
                              by.y = 'Transcript stable ID') %>% 
  column_to_rownames()

new_trt_df = new_trt_df %>% group_by(`Gene stable ID`) %>% 
  summarise_if(is.numeric, sum, na.rm = T) %>% 
  column_to_rownames('Gene stable ID')

# mutate(p.adj_trt = p.adjust(p.value_trt, method = 'BH'), 
#        p.adj_tpm = p.adjust(p.value_tpm, method = 'BH')) %>% 
#   rename(p.adj_prx = p.adj) %>% 


# t.tests -----------------------------------------------------------------
print('Calculating TPM t.test')
res.df_tpm = theVsTox_fun(df = new_trt_df, omics = 'target', FUN = t.test, 
                          complete_cases = comp_cas)
res.df_tpm$p.value_tpm = res.df_tpm$p.value %>% as.character() %>% as.numeric()
# res.df_tpm$p.value_tpm %>% na.omit() %>% hist(breaks = seq(0, 1, 0.05))




  # paste0(seq(1, 16, length.out = 100)) %>% .[-1]


results = matrix(data = -1, nrow = 1, ncol = 5) %>% as.data.frame()
print('Calculating TrT t.test')
p = progress_estimated(length(factor_miRNAs))
for (factor_miRNA in factor_miRNAs) {
  temp_res = NULL
  
  # If it takes 1h, column_to_rownames() #If result.df not found, exchange TrT/TRC
  res.df_trt = theVsTox_fun(df = new_trt_df, omics = factor_miRNA, FUN = t.test, 
                            complete_cases = comp_cas)
  res.df_trt$p.value_trt = res.df_trt$p.value %>% as.character() %>% as.numeric()
  res.df_trt$statistic.t_trt = res.df_trt$statistic.t %>% as.character() %>% as.numeric()
  
  # res.df_trt$p.value_trt %>% na.omit() %>% hist(breaks = seq(0, 1, 0.05))
  
  new_trt_df = new_trt_df %>% select(-contains('p.')) %>% 
    select(-contains('statistic'))
  
  new_trt_df = merge.data.frame(x = rownames_to_column(new_trt_df), 
                                y = rownames_to_column(select(res.df_trt, 
                                                              p.value_trt, 
                                                              statistic.t_trt)), 
                                by = 'rowname', all = T) %>% 
    merge.data.frame(y = rownames_to_column(select(res.df_tpm, p.value_tpm)), 
                     by = 'rowname', all = T) %>% 
    column_to_rownames()
  
  trt_proteomx = merge.data.frame(x = rownames_to_column(new_trt_df), 
                                  y = select(.data = proteomx_biom, 
                                             -`Transcript stable ID`), 
                                  by.x = 'rowname',
                                  by.y = 'Gene stable ID') %>% 
    rename(p.value_prx = p.value, 
           p.adj_prx = p.adj) %>% 
    mutate(p.adj_trt = p.adjust(p.value_trt, method = 'BH'), 
           p.adj_tpm = p.adjust(p.value_tpm, method = 'BH')) %>% 
    filter(!duplicated(rowname)) %>% 
    column_to_rownames()
  
  trt_proteomx$p.adj_trt %>% na.omit() %>% hist(breaks = seq(0, 1, 0.05))
  
  test = trt_proteomx %>% 
    rownames_to_column() %>% 
    filter(
      p.value_trt < 0.05,
      p.value_tpm > 0.05,
      p.value_prx < 0.05, 
      sign(statistic.t_prx) == sign(statistic.t_trt)
    ) %>% 
    column_to_rownames()
  
  temp_res = c(temp_res, nrow(test))
  
  test = trt_proteomx %>% 
    rownames_to_column() %>% 
    filter(
      p.adj_trt < 0.05,
      p.adj_tpm > 0.05,
      p.adj_prx < 0.05,
      sign(statistic.t_prx) == sign(statistic.t_trt)
    ) %>% 
    column_to_rownames()
  
  temp_res = c(temp_res, nrow(test))
  
  test = trt_proteomx %>% 
    rownames_to_column() %>% 
    filter(
      p.adj_trt < 0.05,
      p.value_tpm > 0.05,
      p.value_prx < 0.05,
      sign(statistic.t_prx) == sign(statistic.t_trt)
    ) %>% 
    column_to_rownames()
  
  temp_res = c(temp_res, nrow(test))
  
  test = trt_proteomx %>% 
    rownames_to_column() %>% 
    filter(
      p.adj_trt < 0.05,
      p.value_tpm > 0.05,
      p.adj_prx < 0.05,
      sign(statistic.t_prx) == sign(statistic.t_trt)
    ) %>% 
    column_to_rownames()
  
  temp_res = c(temp_res, nrow(test))
  
  test = trt_proteomx %>% 
    rownames_to_column() %>% 
    filter(
      p.value_trt > 0.05,
      p.adj_tpm < 0.05,
      p.value_prx > 0.05
    ) %>% 
    column_to_rownames()
  
  temp_res = c(temp_res, nrow(test))
  
  results = rbind.data.frame(results, temp_res)
  rownames(results)[nrow(results)] = factor_miRNA
  print(tail(results, 1))
  print(p$tick()$print())
}
naam_file = paste0(comp, '_diff_TrT&PrxVsTPM', '.rds', collapse = '')
saveRDS(object = results, file = naam_file)

