 
# Functions ---------------------------------------------------------------
pseudocount = function(x, addition = 1) {
  x = x + addition
  return(x)
}

boxplot.dose = function(df, omics, ...) {
  x = df %>% dplyr::select(matches(omics))
  the = x %>% dplyr::select(contains('The')) %>% t()
  tox = x %>% dplyr::select(contains('Tox')) %>% t()
  data = cbind.data.frame(the, tox)
  colnames(data) = c('Therapeutic', 'Toxic')
  boxplot(data, ...)
}

viewPlotProtx = function(df, row_i, p.val_col = NULL){
  r.data = df[row_i, ]
  r.name = main = rownames(r.data)
  if (!is.null(p.val_col)) {
    p.value = format.pval(r.data[, p.val_col], 3)
    main = paste0(r.name, ' ; ', 'p.value = ', p.value)
  }
  r.data %>% 
    dplyr::select(matches('Proteomics')) %>% 
    as.matrix() %>% 
    barplot(las = 2, main = main, ylab = 'Proteomics Expression')
  
  for (l in seq(3.7, 85.1, 3.6)) { 
    if (l == 14.5) {
      abline(v = 14.5, col = 'black', lwd = 3, lty = 2)
    } else {
      abline(v = l, col = 'gray')
    }
  }
}

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
  
  for (l in seq(3.5, 21.5, 3)) { 
    if (l == 12.5) {
      abline(v = 12.5, col = 'black', lwd = 6, lty = 1)
    } else {
      abline(v = l, col = 'gray')
    }
  }
  
  par(new = TRUE)
}

theVsTox_fun_limma <- function(df, omics, cond1, cond2, title = NULL, 
                               plotting = F,  ...) {
  pseudocount = function(x, addition = 1) {
    x = x + addition
    return(x)
  }
  
  cond1_cols = df %>% dplyr::select(contains(cond1)) %>% 
    dplyr::select(starts_with(omics)) %>% colnames()
  cond2_cols = df %>% dplyr::select(contains(cond2)) %>% 
    dplyr::select(starts_with(omics)) %>% colnames()
  
  stopifnot(length(cond1_cols) == length(cond2_cols))
  
  condition = c(rep(cond1, length(cond1_cols)), 
             rep(cond2, length(cond2_cols)))
  
  design = model.matrix(~ condition)
  
  res.df = df %>% 
    dplyr::select(cond1_cols, cond2_cols) %>% 
    pseudocount() %>%
    log2() %>% 
    lmFit(design) %>% 
    eBayes(trend = T) 
    # decideTests() %>% summary()
    # topTable(coef=ncol(design), p.value = 0.05) #%>% nrow()
  
  if (plotting) {
    res.df$p.value_tpm %>% .
      as.character() %>%
      as.numeric() %>%
      na.omit() %>% 
      hist(breaks = seq(0, 1, 0.05), main = 'TPM', xlab = 'p. value')
  }  
  
  res.df = res.df %>% 
    as.data.frame() %>% 
    dplyr::select(contains('condition'))
  
  res.df[, paste0('p.value_', title)] = res.df %>% 
    dplyr::select(contains('p.value')) %>% 
    .[, 1] %>% 
    as.character() %>% 
    as.numeric()
  res.df[, paste0('statistic.t_', title)] = res.df %>% 
    dplyr::select(starts_with('t')) %>% 
    .[, 1] %>% 
    as.character() %>% 
    as.numeric()
  
  res.df = res.df %>% 
    dplyr::select(matches('p.value_|statistic.t_')) %>% 
    rownames_to_column() %>% 
    mutate_all(as.character) %>% 
    mutate_at(vars(matches('p.value_|statistic.t_')), as.numeric)
  
  
  res.df[, paste0('p.adj_', title)] = 
    p.adjust(res.df[, paste0('p.value_', title)], 'BH')
  
  return(res.df)
  
}


setwd("/share/script/hecatos/juantxo/analysis_trc")
source('functions_JOA.R')
forceLibrary(c('dplyr', 'tibble', 'limma', 'edgeR'))

# Get data ----------------------------------------------------------------

comp = ''
plotting = F
miRNA_factor = 0.1
TrT_miF = paste0('TrT_', miRNA_factor, '_') 
setwd('/share/script/hecatos/juantxo/analysis_trc/')
source('differentially_expressed_proteins_JOA/DEPs_dose.R')

setwd('/share/analysis/hecatos/juantxo/Score/output/')
setwd(comp)
# ; setwd(comp)
list.files(pattern = 'factor_0_1') %>% sort(decreasing = T) %>% setwd()
# 
all_file = list.files(pattern = 'All_trt|all_trt')
rm(trt_df)
# common_trt_df = mergeFiles(files_patt = TrT_miF, row_names = T, progr_bar = F)
if (length(all_file) == 0) {
  if (sum(grepl(pattern = 'rds', x = list.files()))) {
    trt_df = mergeFilesRds(files_patt = TrT_miF, row_names = T, progr_bar = F, 
                           all = T)
  } else {
    trt_df = mergeFiles(files_patt = TrT_miF, row_names = T, progr_bar = F, 
                        all = T)
  }
  trt_df = trt_df %>% column_to_rownames() 
  colnames(trt_df) = colnames(trt_df) %>% 
    gsub(pattern = '_TrT_0.27.txt.rds', replacement = '')
  trt_df %>% dplyr::select(starts_with('target')) %>% filterSamplesBySeqDepth() %>%
    ncol() %>% identical(ncol(dplyr::select(trt_df, starts_with('target')))) %>% 
    stopifnot()
  saveRDS(object = trt_df, file = paste0(comp, '_all_trt.rds'))
  
} else {
  trt_df = readRDS(list.files(pattern = 'All|all'))
  if ('rowname' %in% colnames(trt_df)) {
    trt_df = trt_df %>% column_to_rownames() 
  }
}

trt_untr = readRDS('/share/analysis/hecatos/juantxo/Score/output/UNTR/2020-03-09_16:43:15__factor_0_1/All_trt.rds') %>% 
  column_to_rownames()

# Clean data -------------------------------------------------------------------

trt_df_ori = trt_df

colnames(trt_df) = colnames(trt_df) %>% 
  gsub(pattern = '.txt', replacement = '') %>%
  gsub(pattern = '.rds', replacement = '') %>%
  gsub(pattern = '_TrTscore', replacement = '')

colnames(trt_untr) = colnames(trt_untr) %>% 
  gsub(pattern = '.txt', replacement = '') %>%
  gsub(pattern = '.rds', replacement = '') %>%
  gsub(pattern = '_TrTscore', replacement = '')

# Remove thousands of useless columns
trt_vars = colnames(trt_df) %>% .[grep('TrT_', .)] %>% .[!grepl(TrT_miF, .)]
trt_df = trt_df[, !(colnames(trt_df) %in% trt_vars)]
trt_df = trt_df[, !grepl('hsa', colnames(trt_df))]

trt_vars = colnames(trt_untr) %>% .[grep('TrT_', .)] %>% .[!grepl(TrT_miF, .)]
trt_untr = trt_untr[, !(colnames(trt_untr) %in% trt_vars)]
trt_untr = trt_untr[, !grepl('hsa', colnames(trt_untr))]

if (plotting) {
  trt_df %>% dplyr::select(starts_with('target')) %>% 
    apply(2, median, na.rm = T) %>% 
    barplot(las = 2, main = 'Median TPM expression per Sample')
}

comp_cas = 12

# Exceptions
if (any(grepl(pattern = 'Dox', x = colnames(trt_df)))) {
  trt_df = trt_df %>% dplyr::select(-contains('072_3'))
  trt_untr = trt_untr %>% dplyr::select(-contains('072_3'))
  comp_cas = 11
}
if (any(grepl(pattern = 'Epi', x = colnames(trt_df)))) {
  trt_df = trt_df %>% dplyr::select(-contains('072'))
  trt_untr = trt_untr %>% dplyr::select(-contains('072'))
  comp_cas = 9
}

trt_comp_untr = merge.data.frame(x = rownames_to_column(trt_df), 
                                 y = rownames_to_column(trt_untr), 
                                 by = 'rowname', all = T) %>% 
  column_to_rownames()

trt_df_geneid = merge.data.frame(x = rownames_to_column(trt_comp_untr), 
                                 y = dplyr::select(.data = proteomx_biom, 
                                                   ensembl_transcript_id, ensembl_gene_id), 
                                 by.x = 'rowname',
                                 by.y = 'ensembl_transcript_id') #%>% 
# filter(!duplicated(rowname)) #%>% 
# column_to_rownames()



trt_df_geneid = trt_df_geneid %>% group_by(ensembl_gene_id) %>% 
  summarise_if(is.numeric, sum, na.rm = T) %>% 
  column_to_rownames('ensembl_gene_id')

# mutate(p.adj_theVStox_TPM = p.adjust(p.adj_theVStox_TPM, method = 'BH'), 
#        p.adj_theVStox_TPM = p.adjust(p.value_tpm, method = 'BH')) %>% 
#   rename(p.adj_prx = p.adj) %>% 



# limma -------------------------------------------------------------------


theVStox_limma_tpm = theVsTox_fun_limma(df = trt_df_geneid, omics = 'target',  
                                   cond1 = 'The', cond2 = 'Tox', 
                                   title = 'theVStox_TPM')

# If it takes 1h, column_to_rownames() #If result.df not found, exchange TrT/TRC
untVSthe_limma_tpm = theVsTox_fun_limma(df = trt_df_geneid, omics = 'target',   
                                   cond1 = 'UNTR', cond2 = 'The', 
                                   title = 'untVSthe_TPM')

untVStox_limma_tpm = theVsTox_fun_limma(df = trt_df_geneid, omics = 'target',   
                                   cond1 = 'UNTR', cond2 = 'Tox', 
                                   title = 'untVStox_TPM')

theVStox_limma_trt = theVsTox_fun_limma(df = trt_df_geneid, omics = TrT_miF,   
                                   cond1 = 'The', cond2 = 'Tox', 
                                   title = 'theVStox_TrT')

untVSthe_limma_trt = theVsTox_fun_limma(df = trt_df_geneid, omics = TrT_miF,   
                                   cond1 = 'UNTR', cond2 = 'The', 
                                   title = 'untVSthe_TrT')

untVStox_limma_trt = theVsTox_fun_limma(df = trt_df_geneid, omics = TrT_miF,   
                                   cond1 = 'UNTR', cond2 = 'Tox', 
                                   title = 'untVStox_TrT')

trt_df_limmas = merge.data.frame(x = rownames_to_column(trt_df_geneid), 
                                  y = theVStox_limma_tpm,
                                  by = 'rowname', all = T) %>% 
  merge.data.frame(y = untVSthe_limma_tpm,
                   by = 'rowname', all = T) %>% 
  merge.data.frame(y = untVStox_limma_tpm, 
                   by = 'rowname', all = T) %>% 
  merge.data.frame(y = theVStox_limma_trt, 
                   by = 'rowname', all = T) %>% 
  merge.data.frame(y = untVSthe_limma_trt, 
                   by = 'rowname', all = T) %>% 
  merge.data.frame(y = untVStox_limma_trt, 
                   by = 'rowname', all = T) %>% 
  column_to_rownames()

# Combine transcrx and protx ----------------------------------------------

proteomx_biom_filt = proteomx_biom %>% 
  mutate(id = paste0(ensembl_gene_id, '_', uniprotswissprot)) %>%
  .[!duplicated(.[, 'id']), ] %>% 
  dplyr::select(-c(id, rowname))

trt_proteomx = merge.data.frame(x = rownames_to_column(trt_df_limmas), 
                                y = dplyr::select(.data = proteomx_biom_filt, 
                                                  -ensembl_transcript_id), 
                                by.x = 'rowname',
                                by.y = 'ensembl_gene_id') %>% 
  mutate(rowname = paste0(rowname, '_', uniprotswissprot)) %>%
  column_to_rownames()

trt_proteomx$ensembl_gene_id = trt_proteomx %>% 
  rownames() %>% 
  gsub('_.*', '', .)
# rename(p.value_theVStox_prx = p.value, 
#        p.adj_prx = p.adj) %>% 
# mutate(p.adj_theVStox_TPM = p.adjust(p.adj_theVStox_TPM, method = 'BH'), 
#        p.adj_theVStox_TPM = p.adjust(p.value_tpm, method = 'BH')) %>% 
# # filter(!duplicated(rowname)) %>% 
# column_to_rownames()

if (any(grepl(pattern = 'Dox', x = colnames(trt_df)))) {
  trt_proteomx = trt_proteomx %>% dplyr::select(-contains('072_3'))
}
if (any(grepl(pattern = 'Epi', x = colnames(trt_df)))) {
  trt_proteomx = trt_proteomx %>% dplyr::select(-contains('072'))
}

# # Correlation ------------------------------------------------------------------
# 
# cor_trt_prot = apply_2D(df = dplyr::select(trt_proteomx, starts_with(TrT_miF)),
#                         FUN = cor.test, complete_cases = comp_cas*2, 
#                         col.x = 'adsf', col.y = 'asdfa',
#                         y = dplyr::select(trt_proteomx, starts_with('Proteomics'))) %>%
#   rownames_to_column() %>%
#   transmute(estimate.cor_trt = 
#               as.numeric(as.character(estimate.cor)), rowname) %>%
#   column_to_rownames()
# 
# cor_tpm_prot = apply_2D(df = dplyr::select(trt_proteomx, starts_with('target')),
#                         FUN = cor.test, complete_cases = comp_cas*2, 
#                         col.x = 'adsf', col.y = 'asdfa',
#                         y = dplyr::select(trt_proteomx, starts_with('Proteomics'))) %>%
#   rownames_to_column() %>%
#   transmute(estimate.cor_tpm = 
#               as.numeric(as.character(estimate.cor)), rowname) %>%
#   column_to_rownames()
# 
# trt_proteomx = merge.data.frame(x = rownames_to_column(trt_proteomx),
#                                 y = rownames_to_column(cor_trt_prot),
#                                 by = 'rowname', all = T) %>%
#   merge.data.frame(y = rownames_to_column(cor_tpm_prot), by = 'rowname', 
#                    all = T) %>%
#   column_to_rownames()
# 
# test = trt_proteomx %>% rownames_to_column() %>%  
#   filter(estimate.cor_trt > (0.4 + estimate.cor_tpm)) %>% column_to_rownames()
# 
# # Minimum > Maximum ------------------------------------------------------------
# 
# trt_the_range = trt_proteomx %>% dplyr::select(starts_with(TrT_miF)) %>%
#   dplyr::select(contains('The')) %>% na.omit() %>% apply(1, range, na.rm = T) %>%
#   t() %>% as.data.frame() %>% rownames_to_column() %>% 
#   filter(!is.infinite(V1)) %>% rename(min_the = V1, max_the = V2) %>% 
#   column_to_rownames()
# 
# trt_tox_range = trt_proteomx %>% dplyr::select(starts_with(TrT_miF)) %>%
#   dplyr::select(contains('Tox')) %>% na.omit() %>% apply(1, range, na.rm = T) %>%
#   t() %>% as.data.frame() %>% rownames_to_column() %>% 
#   filter(!is.infinite(V1)) %>% rename(min_tox = V1, max_tox = V2) %>% 
#   column_to_rownames()
# 
# a = merge.data.frame(rownames_to_column(trt_the_range),
#                      rownames_to_column(trt_tox_range),
#                      'rowname') %>% column_to_rownames()
# 
# b = a[a$min_the > a$max_tox,,F] %>% rownames()
# test = trt_proteomx[b,,F]
# 
# if (plotting) {
#   trt_proteomx$p.adj_theVStox_TPM %>% 
#     hist(breaks = seq(0,1,0.05), main = 'TrT', xlab = 'p. adjusted values')
#   trt_proteomx$p.adj_theVStox_TPM %>% 
#     hist(breaks = seq(0,1,0.05), main = 'TPM', xlab = 'p. adjusted values')
#   trt_proteomx$p.adj_prx %>% 
#     hist(breaks = seq(0,1,0.05), main = 'Proteomics', xlab = 'p. adjusted values')
# }
# 
# 
# Best cases --------------------------------------------------------------


# # VVV
# vvv = trt_proteomx %>% 
#   rownames_to_column() %>% 
#   filter(
#     p.adj_theVStox_TPM < 0.05,
#     p.value_tpm > 0.05,
#     p.value_theVStox_prx < 0.05, 
#     sign(statistic.t_prx) == sign(statistic.t_trt)
#   ) %>% 
#   column_to_rownames()
# # test %>% nrow() %>% print()
# 
# # AAA
# aaa = trt_proteomx %>% 
#   rownames_to_column() %>% 
#   filter(
#     p.adj_theVStox_TPM < 0.05,
#     p.adj_theVStox_TPM > 0.05,
#     p.adj_prx < 0.05,
#     sign(statistic.t_prx) == sign(statistic.t_trt)
#   ) %>% 
#   column_to_rownames()
# 

# # AVV
# avv = trt_proteomx %>% 
#   rownames_to_column() %>% 
#   filter(
#     p.adj_theVStox_TPM < 0.05,
#     p.value_tpm > 0.05,
#     p.value_theVStox_prx < 0.05,
#     sign(statistic.t_prx) == sign(statistic.t_trt)
#   ) %>% 
#   column_to_rownames()
# 
# # AVV_bad
# avv_bad = trt_proteomx %>% 
#   rownames_to_column() %>% 
#   filter(
#     p.adj_theVStox_TPM < 0.05,
#     p.value_tpm > 0.05,
#     p.value_theVStox_prx > 0.05,
#     sign(statistic.t_prx) == sign(statistic.t_trt)
#   ) %>% 
#   column_to_rownames()

# # AVA
# ava = trt_proteomx %>% 
#   rownames_to_column() %>% 
#   filter(
#     p.adj_theVStox_TPM < 0.05,
#     p.value_tpm > 0.05,
#     p.adj_prx < 0.05,
#     sign(statistic.t_prx) == sign(statistic.t_trt)
#   ) %>% 
#   column_to_rownames()

# # VAV
# vav = trt_proteomx %>% 
#   rownames_to_column() %>% 
#   filter(
#     p.adj_theVStox_TPM > 0.05,
#     p.adj_theVStox_TPM < 0.05,
#     p.value_theVStox_prx > 0.05
#   ) %>% 
#   column_to_rownames()
# 
# # VAV_bad
# vav_bad = trt_proteomx %>% 
#   rownames_to_column() %>% 
#   filter(
#     p.adj_theVStox_TPM > 0.05,
#     p.adj_theVStox_TPM < 0.05,
#     p.value_theVStox_prx < 0.05
#   ) %>% 
#   column_to_rownames()

# Sensitivity and Specificity --------------------------------------------------

#### TPM ####

conf_matrix <- function(df, pred, true, dir_pred, dir_true) {
  # Refer to column names stored as strings with the `.data` pronoun:
  tp = df %>%
    rownames_to_column() %>%
    filter(
      .data[[pred[[1]]]] < 0.05,
      .data[[true[[1]]]] < 0.05,
      sign(.data[[dir_pred[[1]]]]) == sign(.data[[dir_true[[1]]]])
    ) %>% 
    column_to_rownames() %>% 
    rownames()
  
  tn = df %>%
    rownames_to_column() %>%
    filter(
      .data[[pred[[1]]]] > 0.05,
      .data[[true[[1]]]] > 0.05
    ) %>% 
    column_to_rownames() %>% 
    rownames()
  
  fp1 = df %>%
    rownames_to_column() %>%
    filter(
      .data[[pred[[1]]]] < 0.05,
      .data[[true[[1]]]] > 0.05
    ) %>% 
    column_to_rownames() %>% 
    rownames() 
  
  fp2 =  df %>%
    rownames_to_column() %>%
    filter(
      .data[[pred[[1]]]] < 0.05,
      .data[[true[[1]]]] < 0.05,
      sign(.data[[dir_pred[[1]]]]) != sign(.data[[dir_true[[1]]]])
    ) %>% 
    column_to_rownames() %>% 
    rownames()
  
  fp = c(fp1, fp2)
  
  fn = df %>%
    rownames_to_column() %>%
    filter(
      .data[[pred[[1]]]] > 0.05,
      .data[[true[[1]]]] < 0.05
    ) %>% 
    column_to_rownames() %>% 
    rownames() 
  
  conf_matr = list(true_positive = tp,
                   true_negative = tn,
                   false_positive = fp,
                   false_negative = fn)
  
  return(conf_matr)
}
# Therapeutic versus Toxic
conf_matrix_tpm_theVStox = conf_matrix(df=trt_proteomx,
                                       pred= 'p.adj_theVStox_TPM', 
                                       true='p.value_theVStox_prx', 
                                       dir_pred='statistic.t_theVStox_TPM', 
                                       dir_true='statistic.t_theVStox_prx')

conf_matrix_trt_theVStox = conf_matrix(df=trt_proteomx,
                                       pred= 'p.adj_theVStox_TrT', 
                                       true='p.value_theVStox_prx', 
                                       dir_pred='statistic.t_theVStox_TrT', 
                                       dir_true='statistic.t_theVStox_prx')

# Untreated vs Therapeutic

conf_matrix_tpm_untVSthe = conf_matrix(df=trt_proteomx,
                                       pred= 'p.adj_untVSthe_TPM', 
                                       true='p.value_untVSthe_prx', 
                                       dir_pred='statistic.t_untVSthe_TPM', 
                                       dir_true='statistic.t_untVSthe_prx')

conf_matrix_trt_untVSthe = conf_matrix(df=trt_proteomx,
                                       pred= 'p.adj_untVSthe_TrT', 
                                       true='p.value_untVSthe_prx', 
                                       dir_pred='statistic.t_untVSthe_TrT', 
                                       dir_true='statistic.t_untVSthe_prx')

# Untreated vs Toxic

conf_matrix_tpm_untVStox = conf_matrix(df=trt_proteomx,
                                       pred= 'p.adj_untVStox_TPM', 
                                       true='p.value_untVStox_prx', 
                                       dir_pred='statistic.t_untVStox_TPM', 
                                       dir_true='statistic.t_untVStox_prx')

conf_matrix_trt_untVStox = conf_matrix(df=trt_proteomx,
                                       pred= 'p.adj_untVStox_TrT', 
                                       true='p.value_untVStox_prx', 
                                       dir_pred='statistic.t_untVStox_TrT', 
                                       dir_true='statistic.t_untVStox_prx')

# List of lists

all_conf_matrices = list(conf_matrix_trt_theVStox = conf_matrix_trt_theVStox,
                         conf_matrix_tpm_theVStox = conf_matrix_tpm_theVStox,
                         conf_matrix_trt_untVSthe = conf_matrix_trt_untVSthe,
                         conf_matrix_tpm_untVSthe = conf_matrix_tpm_untVSthe, 
                         conf_matrix_trt_untVStox = conf_matrix_trt_untVStox,
                         conf_matrix_tpm_untVStox = conf_matrix_tpm_untVStox)

i = 1
for (cf in all_conf_matrices) {
  all_conf_matrices %>% names() %>% .[i] %>% print()
  cf_summ = summary(cf)
  values = cf_summ[, 1] %>% as.numeric()
  acc = sum(values[1:2]/sum(values))
  print(acc)
  i = i + 1
}
# 
# findCorrected = function(list1, list2){
#   fn_tp = list2[['true_positive']] %in% list1[['false_negative']] %>% 
#     sum %>% 
#     print
#   fp_tn = list2[['true_negative']] %in% list1[['false_positive']] %>% 
#     sum %>% 
#     print
#   fn_tp2 = list1[['true_positive']] %in% list2[['false_negative']] %>% 
#     sum %>% 
#     print
#   fp_tn2 = list1[['true_negative']] %in% list2[['false_positive']] %>% 
#     sum %>% 
#     print
#   invisible(return(NULL))
# }
# 
# findCorrected(conf_matrix_tpm_theVStox, conf_matrix_trt_theVStox)
# findCorrected(conf_matrix_tpm_untVSthe, conf_matrix_trt_untVSthe)
# findCorrected(conf_matrix_tpm_untVStox, conf_matrix_trt_untVStox)
# 
# degs <- function(df, p.val_col, comp = comp) {
#   degs_out = df %>% 
#     rownames_to_column() %>% 
#     filter(.data[[p.val_col]] < 0.05) %>% 
#     arrange(.data[[p.val_col]]) %>% 
#     column_to_rownames() %>% 
#     rownames() %>% 
#     gsub('_.*', '', .) %>% 
#     as.data.frame()
#   
#   colnames(degs_out) = p.val_col
#   
#   return(degs_out)
# }
# 
# setwd('~/Desktop/')
# dir.create('limma')
# setwd('limma')
# dir.create(comp)
# setwd(comp)
# a = degs(df = trt_proteomx, p.val_col = 'p.adj_theVStox_TPM')
# file_naam = paste0(comp, '_', colnames(a))
# write.table(x = a, file = file_naam, quote = F, 
#             row.names = F, col.names = F)
# a = degs(df = trt_proteomx, p.val_col = 'p.adj_theVStox_TrT')
# file_naam = paste0(comp, '_', colnames(a))
# write.table(x = a, file = file_naam, quote = F, 
#             row.names = F, col.names = F)
# 
# a = degs(df = trt_proteomx, p.val_col = 'p.adj_untVSthe_TPM')
# file_naam = paste0(comp, '_', colnames(a))
# write.table(x = a, file = file_naam, quote = F, 
#             row.names = F, col.names = F)
# 
# a = degs(df = trt_proteomx, p.val_col = 'p.adj_untVSthe_TrT')
# file_naam = paste0(comp, '_', colnames(a))
# write.table(x = a, file = file_naam, quote = F, 
#             row.names = F, col.names = F)
# 
# a = degs(df = trt_proteomx, p.val_col = 'p.adj_untVStox_TPM')
# file_naam = paste0(comp, '_', colnames(a))
# write.table(x = a, file = file_naam, quote = F, 
#             row.names = F, col.names = F)
# 
# a = degs(df = trt_proteomx, p.val_col = 'p.adj_untVStox_TrT')
# file_naam = paste0(comp, '_', colnames(a))
# write.table(x = a, file = file_naam, quote = F, 
#             row.names = F, col.names = F)
# 
# #### TRT ####
# # test = avv
# 
# # colnames(test) = colnames(test) %>% gsub('_TrT', '', .)
# 
# if (plotting) {
#   for (rown in rownames(test)) {
#     layout(matrix(c(1,2), 1, 2, byrow = T))
#     row_data = test[rown, ]
#     max_tpm = test[rown, , F] %>% dplyr::select(starts_with('target')) %>% max(na.rm = T)
#     #par(mfrow = c(1, 2))
#     
#     viewPlotProtx(test, rown)
#     
#     readline(prompt = "Press [enter] to continue")
#     
#     plotOmics(df = trt_proteomx, omics = paste0('^', TrT_miF), 
#               transcript = rown, xaxt = 'n', type = 'b', 
#               ylab = 'gray: TPM | black: TrT', xlab = '',
#               ylim = c(0, max_tpm))
#     plotOmics(df = trt_proteomx, omics = '^target', transcript = rown, xaxt = 'n',
#               type = 'b', col = 'gray', yaxt = 'n', xlab = '')
#     axis(1, at = 1:24, labels = colnames(test)[grep('^target', colnames(test))], las = 2)
#     
#     readline(prompt = "Press [enter] to continue")
#     
#     layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
#     boxplot.dose(row_data, 'Proteomics', main = rown, ylab = 'Proteomics Expression')
#     boxplot.dose(row_data, '^target', main = rown, ylab = 'TPM Expression')
#     boxplot.dose(row_data, TrT_miF, main = rown, ylab = 'TrT Expression')
#     
#     
#     #test[rown, , F] %>% dplyr::select(starts_with(TrT_miF)) %>% dplyr::select(contains('The')) %>%
#     #  as.numeric() %>% mean(na.rm = T) %>% abline(h = .)
#     
#     #  test[rown, , F] %>% dplyr::select(starts_with(TrT_miF)) %>% dplyr::select(contains('Tox')) %>%
#     #   as.numeric() %>% mean(na.rm = T) %>% abline(h = .)
#     # axis(1, at = 1:24, labels = colnames(trt_proteomx)[grep('^target',
#     #colnames(trt_df))], las = 2)
#     readline(prompt = "Press [enter] to continue")
#     
#   }
# }
# 
