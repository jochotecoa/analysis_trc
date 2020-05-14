
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
# Analysis ----------------------------------------------------------------
comp = 'DOC'

setwd('/share/analysis/hecatos/juantxo/Score/output/DOC/2020-03-09_13:09:20_UTC/TRCscore/factors_DOC/')
proteomx_biom = readRDS('proteomx_biom.rds')
# ; setwd(comp)
# list.files(pattern = 'UTC') %>% sort(decreasing = T) %>% setwd()
# common_trt_df = mergeFiles(files_patt = 'TrT', row_names = T, progr_bar = F)
if (length(list.files(pattern = paste0(comp, '_all'))) == 0) {
  trt_df = mergeFilesRds(files_patt = 'TrT', row_names = T, progr_bar = F, all = T)
  trt_df = trt_df %>% column_to_rownames() 
  colnames(trt_df) = colnames(trt_df) %>% 
    gsub(pattern = '_TrT.txt.rds', replacement = '')
  trt_df %>% select(starts_with('target')) %>% filterSamplesBySeqDepth() %>%
    ncol() %>% identical(ncol(select(trt_df, starts_with('target')))) %>% 
    stopifnot()
  saveRDS(object = trt_df, file = paste0(comp, '_all_trt.rds'))
  
} else {
  trt_df = readRDS(list.files(pattern = paste0(comp, '_all')))
  if ('colname' %in% colnames(trt_df)) {
    trt_df = trt_df %>% column_to_rownames() 
  }
}


if (any(grepl(pattern = 'Dox', x = colnames(trt_df)))) {
  trt_df = trt_df %>% select(-contains('072_3'))
  comp_cas = 11
}
if (any(grepl(pattern = 'EPI', x = colnames(trt_df)))) {
  trt_df = trt_df %>% select(-contains('072'))
  comp_cas = 9
} else {
  comp_cas = 12
}
print('Calculating TPM t.test')
res.df_tpm = theVsTox_fun(df = trt_df, omics = 'target', FUN = t.test, complete_cases = comp_cas)
res.df_tpm$p.value_tpm = res.df_tpm$p.value %>% as.character() %>% as.numeric()
# res.df_tpm$p.value_tpm %>% na.omit() %>% hist(breaks = seq(0, 1, 0.05))

factor_miRNAs = 'TrT_' %>% paste0(seq(0.00078125, 0.5, length.out = 20))
# factor_miRNAs = factor_miRNAs[13:15]

results = matrix(data = -1, nrow = 1, ncol = 4) %>% as.data.frame()
print('Calculating TrT t.test')
p = progress_estimated(length(factor_miRNAs))
for (factor_miRNA in factor_miRNAs) {
  
  temp_res = NULL
  
  # If it takes 1h, column_to_rownames() #If result.df not found, exchange TrT/TRC
  res.df_trt = theVsTox_fun(df = trt_df, omics = factor_miRNA, FUN = t.test, complete_cases = comp_cas)
  res.df_trt$p.value_trt = res.df_trt$p.value %>% as.character() %>% as.numeric()
  # res.df_trt$p.value_trt %>% na.omit() %>% hist(breaks = seq(0, 1, 0.05))
  
  
  
  
  new_trt_df = merge.data.frame(x = rownames_to_column(trt_df), 
                                y = rownames_to_column(select(res.df_trt, p.value_trt)), 
                                by = 'rowname', all = T) %>% 
    merge.data.frame(y = rownames_to_column(select(res.df_tpm, p.value_tpm)), 
                     by = 'rowname', all = T) %>% 
    column_to_rownames()
  
  trt_proteomx = merge.data.frame(x = rownames_to_column(new_trt_df), 
                                  y = proteomx_biom, 
                                  by.y = 'Transcript stable ID', 
                                  by.x = 'rowname') %>% 
    rename(p.value_prx = p.value) %>% 
    mutate(p.adj_trt = p.adjust(p.value_trt, method = 'BH'), 
           p.adj_tpm = p.adjust(p.value_tpm, method = 'BH')) %>% 
    rename(p.adj_prx = p.adj) %>% 
    filter(!duplicated(rowname)) %>% 
    # filter(p.adj.tcx < 0.05, p.adj.CEL < 0.05) %>% 
    column_to_rownames()
  
  test = trt_proteomx %>% 
    rownames_to_column() %>% 
    filter(
      p.value_trt < 0.05,
      p.value_tpm > 0.05,
      p.value_prx < 0.05
    ) %>% 
    column_to_rownames()
  
  temp_res = c(temp_res, nrow(test))
  
  test = trt_proteomx %>% 
    rownames_to_column() %>% 
    filter(
      p.adj_trt < 0.05,
      p.adj_tpm > 0.05,
      p.adj_prx < 0.05
    ) %>% 
    column_to_rownames()
  
  temp_res = c(temp_res, nrow(test))
  
  test = trt_proteomx %>% 
    rownames_to_column() %>% 
    filter(
      p.adj_trt < 0.05,
      p.value_tpm > 0.05,
      p.value_prx < 0.05
    ) %>% 
    column_to_rownames()
  
  temp_res = c(temp_res, nrow(test))
  
  test = trt_proteomx %>% 
    rownames_to_column() %>% 
    filter(
      p.adj_trt < 0.05,
      p.value_tpm > 0.05,
      p.adj_prx < 0.05
    ) %>% 
    column_to_rownames()
  
  temp_res = c(temp_res, nrow(test))
  
  results = rbind.data.frame(results, temp_res)
  rownames(results)[nrow(results)] = factor_miRNA
  print(tail(results, 1))
  print(p$tick()$print())
}
naam_file = substr(x = paste0(rnorm(1), '.rds', collapse = ''), 4, 22)
saveRDS(results, file = naam_file)
