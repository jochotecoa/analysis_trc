
# Functions ---------------------------------------------------------------

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
theVsTox_fun <- function(df, omics, cond1, cond2, title = NULL, plot = F, ...) {
  the_cols = df %>% dplyr::select(contains(cond1)) %>% 
    dplyr::select(starts_with(omics)) %>% colnames()
  tox_cols = df %>% dplyr::select(contains(cond2)) %>% 
    dplyr::select(starts_with(omics)) %>% colnames()
  
  res.df = df %>% dplyr::select(the_cols, tox_cols) %>% zeroToNa() %>% log2() %>% 
    apply_2D(col.x = the_cols, col.y = tox_cols, ...)
  
  if (plot) {
    res.df$p.value_tpm %>% .
      as.character() %>%
      as.numeric() %>%
      na.omit() %>% 
      hist(breaks = seq(0, 1, 0.05), main = 'TPM', xlab = 'p. value')
  }  
  
  
  
  res.df[, paste0('p.value_', title)] = res.df$p.value %>% as.character() %>% as.numeric()
  res.df[, paste0('statistic.t_', title)] = res.df$statistic.t %>% as.character() %>% as.numeric()
  
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
forceLibrary(c('dplyr', 'tibble'))

# Get data ----------------------------------------------------------------
# Attention! Never run this script one after another in the same R session
# There's the danger of reusing the variable of one comp onto another
if (!exists('comp')) {
  comp = 'PTX'
}

p.value_t.test = 0.001
proteomics = F
plotting = F
miRNA_factor = 0.1
TrT_miF = paste0('TrT_', miRNA_factor, '_') 
repo_dir = '/share/script/hecatos/juantxo/analysis_trc/'

if (proteomics) {
  # This script will retrieve the DEPs of the compound of choice
  setwd(repo_dir)
  source('differentially_expressed_proteins_JOA/DEPs_dose.R')
  
}

setwd('/share/analysis/hecatos/juantxo/Score/output/')
setwd(comp)
# ; setwd(comp)
list.files(pattern = 'factor_0_1') %>% sort(decreasing = T) %>% setwd()
# 
all_file = list.files(pattern = 'All_trt|all_trt')
if (exists('trt_df')) {rm(trt_df)}

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

if (any(grepl(pattern = 'CEL|MXT|PTX', x = colnames(trt_df)))) {
  trt_df = trt_df %>% dplyr::select(-contains('002_3'))
  trt_untr = trt_untr %>% dplyr::select(-contains('002_3'))
  comp_cas = 11
}


# Add UNTR columns
trt_comp_untr = merge.data.frame(x = rownames_to_column(trt_df), 
                                 y = rownames_to_column(trt_untr), 
                                 by = 'rowname', all = T) %>% 
  column_to_rownames()

# Filter genes that do not have proteomics expression data
if (proteomics) {
  trt_df_geneid = merge.data.frame(x = rownames_to_column(trt_comp_untr), 
                                   y = dplyr::select(.data = proteomx_biom, 
                                                     ensembl_transcript_id, 
                                                     ensembl_gene_id), 
                                   by.x = 'rowname',
                                   by.y = 'ensembl_transcript_id') #%>% 
  # filter(!duplicated(rowname)) #%>% 
  # column_to_rownames()
  
} else {
  
  mart = openMart2018()
  
  biomart_table = biomaRt::getBM(attributes = c('transcript_biotype', 
                                                'uniprotswissprot', 
                                                "ensembl_gene_id", 
                                                "external_gene_name",
                                                "ensembl_transcript_id"), 
                                 filters = 'ensembl_transcript_id', 
                                 values = rownames(trt_comp_untr), mart = mart) %>% 
    dplyr::filter(transcript_biotype == 'protein_coding')
  
  
  trt_df_geneid = merge.data.frame(x = rownames_to_column(trt_comp_untr), 
                                   y = biomart_table, 
                                   by.x = 'rowname', 
                                   by.y = 'ensembl_transcript_id') %>% 
    rename(ensembl_transcript_id = rowname) %>% 
    rownames_to_column
}



trt_df_geneid = trt_df_geneid %>% group_by(ensembl_gene_id) %>% 
  summarise_if(is.numeric, sum, na.rm = T) %>% 
  column_to_rownames('ensembl_gene_id')

# mutate(p.adj_theVStox_TPM = p.adjust(p.adj_theVStox_TPM, method = 'BH'), 
#        p.adj_theVStox_TPM = p.adjust(p.value_tpm, method = 'BH')) %>% 
#   rename(p.adj_prx = p.adj) %>% 


# t.tests -----------------------------------------------------------------

theVStox_t.test_tpm = theVsTox_fun(df = trt_df_geneid, omics = 'target',  
                                   FUN = t.test, cond1 = 'The', cond2 = 'Tox', 
                                   paired = T, complete_cases = comp_cas,
                                   title = 'theVStox_TPM')

# If it takes 1h, column_to_rownames() #If result.df not found, exchange TrT/TRC
theVStox_t.test_trt = theVsTox_fun(df = trt_df_geneid, omics = TrT_miF,   
                                   FUN = t.test, cond1 = 'The', cond2 = 'Tox', 
                                   paired = T, complete_cases = comp_cas, 
                                   title = 'theVStox_TrT')

# UNTR_002_3 has low sequencing depth, so we cannot compare it to others

if (!any(grepl(pattern = 'Dox|Epi', x = colnames(trt_df)))) {
  trt_df_geneid = trt_df_geneid %>% dplyr::select(-contains('002_3'))
  comp_cas = comp_cas - 1
}





untVSthe_t.test_tpm = theVsTox_fun(df = trt_df_geneid, omics = 'target',   
                                   FUN = t.test, cond1 = 'UNTR', cond2 = 'The', 
                                   paired = T, complete_cases = comp_cas, 
                                   title = 'untVSthe_TPM')

untVStox_t.test_tpm = theVsTox_fun(df = trt_df_geneid, omics = 'target',   
                                   FUN = t.test, cond1 = 'UNTR', cond2 = 'Tox', 
                                   paired = T, complete_cases = comp_cas, 
                                   title = 'untVStox_TPM')


untVSthe_t.test_trt = theVsTox_fun(df = trt_df_geneid, omics = TrT_miF,   
                                   FUN = t.test, cond1 = 'UNTR', cond2 = 'The', 
                                   paired = T, complete_cases = comp_cas, 
                                   title = 'untVSthe_TrT')

untVStox_t.test_trt = theVsTox_fun(df = trt_df_geneid, omics = TrT_miF,   
                                   FUN = t.test, cond1 = 'UNTR', cond2 = 'Tox', 
                                   paired = T, complete_cases = comp_cas, 
                                   title = 'untVStox_TrT')

trt_df_t.tests = merge.data.frame(x = rownames_to_column(trt_df_geneid), 
                                  y = theVStox_t.test_tpm,
                                  by = 'rowname', all = T) %>% 
  merge.data.frame(y = untVSthe_t.test_tpm,
                   by = 'rowname', all = T) %>% 
  merge.data.frame(y = untVStox_t.test_tpm, 
                   by = 'rowname', all = T) %>% 
  merge.data.frame(y = theVStox_t.test_trt, 
                   by = 'rowname', all = T) %>% 
  merge.data.frame(y = untVSthe_t.test_trt, 
                   by = 'rowname', all = T) %>% 
  merge.data.frame(y = untVStox_t.test_trt, 
                   by = 'rowname', all = T) %>% 
  column_to_rownames()

# Combine transcrx and protx ----------------------------------------------

if (proteomics) {
  proteomx_biom_filt = proteomx_biom %>% 
    mutate(id = paste0(ensembl_gene_id, '_', uniprotswissprot)) %>%
    .[!duplicated(.[, 'id']), ] %>% 
    dplyr::select(-c(id, rowname))
  
  trt_proteomx = merge.data.frame(x = rownames_to_column(trt_df_t.tests), 
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
  
  
  # Old additional analyses -------------------------------------------------
  
  
  source(
    paste0(
      repo_dir, 'p.values_multiomics/correlation_minimumaximum_bestcases.R'
    )
  )
  
  # Sensitivity and Specificity --------------------------------------------------
  
  
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
    print(summary(cf))
    i = i + 1
  }
  
  i = 1
  for (cf in all_conf_matrices) {
    all_conf_matrices %>% names() %>% .[i] %>% print()
    cf_summ = summary(cf)
    values = cf_summ[, 1] %>% as.numeric()
    acc = sum(values[1:2]/sum(values))
    print(acc)
    i = i + 1
  }
  
  
  # Good cases --------------------------------------------------------------
  
  
  findCorrected = function(list1, list2){
    fn_tp = list2[['true_positive']] %in% list1[['false_negative']] %>% 
      list2[['true_positive']][.] 
    fn_tp %>% length %>% print
    fp_tn = list2[['true_negative']] %in% list1[['false_positive']] %>% 
      list2[['true_negative']][.] 
    fp_tn %>% length %>% print
    fn_tp2 = list1[['true_positive']] %in% list2[['false_negative']] %>% 
      list1[['true_positive']][.] 
    fn_tp2 %>% length %>% print
    fp_tn2 = list1[['true_negative']] %in% list2[['false_positive']] %>% 
      list1[['true_negative']][.] 
    fp_tn2 %>% length %>% print
    
    corrected = list(tp_2 = fn_tp, tn_2 = fp_tn, tp_1 = fn_tp2, tn_1 = fp_tn2)
    
    return(corrected)
  }
  
  corrctd_theVStox = findCorrected(conf_matrix_tpm_theVStox, 
                                   conf_matrix_trt_theVStox)
  corrctd_untVSthe = findCorrected(conf_matrix_tpm_untVSthe, 
                                   conf_matrix_trt_untVSthe)
  corrctd_untVStox = findCorrected(conf_matrix_tpm_untVStox,
                                   conf_matrix_trt_untVStox)
  
  names(corrctd_theVStox) = c('tp_trt', 'tn_trt', 'tp_tpm', 'tn_tpm')
  names(corrctd_untVSthe) = c('tp_trt', 'tn_trt', 'tp_tpm', 'tn_tpm')
  names(corrctd_untVStox) = c('tp_trt', 'tn_trt', 'tp_tpm', 'tn_tpm')
  
  
  
  
}
# Get DEG lists and write them as tables ----------------------------------


degs <- function(df, p.val_col, comp = comp, p.value = 0.05) {
  degs_out = df %>% 
    rownames_to_column() %>% 
    filter(.data[[p.val_col]] < p.value) %>% 
    arrange(.data[[p.val_col]]) %>% 
    column_to_rownames() %>% 
    rownames() %>% 
    gsub('_.*', '', .) %>% 
    as.data.frame()
  
  colnames(degs_out) = p.val_col
  
  return(degs_out)
}

degs_total = NULL

if (proteomics) {
  ttest_cols = 
    c('p.adj_theVStox_TPM', 'p.adj_theVStox_TrT', 'p.value_theVStox_prx',
      'p.adj_untVSthe_TPM', 'p.adj_untVSthe_TrT', 'p.value_untVSthe_prx',
      'p.adj_untVStox_TPM', 'p.adj_untVStox_TrT', 'p.value_untVStox_prx')
  
  for (ttest_col in ttest_cols) {
    setwd('/share/analysis/hecatos/juantxo/Score/analysis/t-test/')
    forceSetWd(comp)
    forceSetWd('filtered_by_protx')
    if (grepl('prx', ttest_col)) {
      forceSetWd('DEPs')
    } else {
      forceSetWd('DEGs')
    }
    degs_df = degs(df = trt_proteomx, p.val_col = ttest_col, 
                   p.value = p.value_t.test)
    file_naam = paste0(comp, '_', colnames(degs_df))
    write.table(x = degs_df, file = file_naam, quote = F, 
                row.names = F, col.names = F)
  }
  
} else {
  ttest_cols = 
    c('p.adj_theVStox_TPM', 'p.adj_theVStox_TrT',
      'p.adj_untVSthe_TPM', 'p.adj_untVSthe_TrT', 
      'p.adj_untVStox_TPM', 'p.adj_untVStox_TrT')
  
  for (ttest_col in ttest_cols) {
    setwd('/share/analysis/hecatos/juantxo/Score/analysis/t-test/')
    forceSetWd(comp)
    forceSetWd('DEGs')
    degs_df = degs(df = trt_df_t.tests, p.val_col = ttest_col, 
                   p.value = p.value_t.test)
    file_naam = paste0(comp, '_', colnames(degs_df))
    write.table(x = degs_df, file = file_naam, quote = F, 
                row.names = F, col.names = F)
    
    degs_total = c(degs_total, nrow(degs_df))
    names(degs_total)[length(degs_total)] = file_naam
  }
  
}

barplot(degs_total, las = 2)


# Write big files ---------------------------------------------------------


setwd('/share/analysis/hecatos/juantxo/Score/analysis/t-test/')
forceSetWd(comp)
if (proteomics) {
  forceSetWd('filtered_by_protx')
  saveRDS(object = trt_proteomx, file = 'whole_table_genes_filtered.rds')
  write.table(x = trt_proteomx, file = 'whole_table_genes_filtered.tsv', 
              sep = '\t')
  saveRDS(object = trt_comp_untr, file = 'whole_table_all_genes.rds')
  write.table(x = trt_comp_untr, file = 'whole_table_all_genes.tsv', 
              sep = '\t')
  
} else {
  saveRDS(object = trt_proteomx, file = 'whole_table_genes_filtered.rds')
  write.table(x = trt_proteomx, file = 'whole_table_genes_filtered.tsv', 
              sep = '\t')
  saveRDS(object = trt_comp_untr, file = 'whole_table_all_genes.rds')
  write.table(x = trt_comp_untr, file = 'whole_table_all_genes.tsv', 
              sep = '\t')
  
}


# Plotting multiomics expressions (old) -----------------------------------


source(
  paste0(
    repo_dir, 'p.values_multiomics/plot_multiomics_expressions.R'
    )
  )
