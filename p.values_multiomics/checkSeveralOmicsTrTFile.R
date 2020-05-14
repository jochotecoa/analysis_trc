
# Functions ---------------------------------------------------------------

boxplot.dose = function(df, omics, ...) {
  x = df %>% select(matches(omics))
  the = x %>% select(contains('The')) %>% t()
  tox = x %>% select(contains('Tox')) %>% t()
  data = cbind.data.frame(the, tox)
  colnames(data) = c('Therapeutic', 'Toxic')
  boxplot(data, ...)
}

viewPlotProtx = function(df, row_i, p.val_col){
  r.data = df[row_i, ]
  p.value = format.pval(r.data[, p.val_col], 3)
  r.name = rownames(r.data)
  r.data %>% 
    select(matches('Proteomics')) %>% 
    as.matrix() %>% 
    barplot(las = 2, 
            main = paste0(r.name, ' ; ', 'p.value = ', p.value),
            ylab = 'Protein Expression')
  
  abline(v = 14.5)
  for (l in seq(3.7, 85.1, 3.6)) abline(v = l, col = 'gray')
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
miRNA_factor = 0.1
TrT_miF = paste0('TrT_', miRNA_factor, '_') 
setwd('/share/script/hecatos/juantxo/analysis_trc/')
source('differentially_expressed_proteins/DEPs_dose.R')

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
  trt_df %>% select(starts_with('target')) %>% filterSamplesBySeqDepth() %>%
    ncol() %>% identical(ncol(select(trt_df, starts_with('target')))) %>% 
    stopifnot()
  saveRDS(object = trt_df, file = paste0(comp, '_all_trt.rds'))
  
} else {
  trt_df = readRDS(list.files(pattern = 'All|all'))
  if ('rowname' %in% colnames(trt_df)) {
    trt_df = trt_df %>% column_to_rownames() 
  }
}

colnames(trt_df) = colnames(trt_df) %>% 
  gsub(pattern = '.txt', replacement = '') %>%
  gsub(pattern = '.rds', replacement = '') %>%
  gsub(pattern = '_TrTscore', replacement = '')

trt_df %>% select(starts_with('target')) %>% 
  apply(2, median, na.rm = T) %>% 
  barplot(las = 2, main = 'Median TPM expression per Sample')

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
  filter(!duplicated(rowname)) %>% 
  column_to_rownames()

new_trt_df = new_trt_df %>% group_by(`Gene stable ID`) %>% 
  summarise_if(is.numeric, sum, na.rm = T) %>% 
  column_to_rownames('Gene stable ID')

# mutate(p.adj_trt = p.adjust(p.value_trt, method = 'BH'), 
#        p.adj_tpm = p.adjust(p.value_tpm, method = 'BH')) %>% 
#   rename(p.adj_prx = p.adj) %>% 


# t.tests -----------------------------------------------------------------

res.df_tpm = theVsTox_fun(df = new_trt_df, omics = 'target', FUN = t.test, 
                          complete_cases = comp_cas)
res.df_tpm$p.value_tpm = res.df_tpm$p.value %>% as.character() %>% as.numeric()
res.df_tpm$p.value_tpm %>% na.omit() %>% 
  hist(breaks = seq(0, 1, 0.05), main = 'TPM', xlab = 'p. value')

# If it takes 1h, column_to_rownames() #If result.df not found, exchange TrT/TRC
res.df_trt = theVsTox_fun(df = new_trt_df, omics = TrT_miF, FUN = t.test, 
                          complete_cases = comp_cas)
res.df_trt$p.value_trt = res.df_trt$p.value %>% as.character() %>% as.numeric()
res.df_trt$statistic.t_trt = res.df_trt$statistic.t %>% as.character() %>% as.numeric()
res.df_trt$p.value_trt %>% na.omit() %>% 
  hist(breaks = seq(0, 1, 0.05), main = 'TrT', xlab = 'p. value')

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

if (any(grepl(pattern = 'Dox', x = colnames(trt_df)))) {
  trt_proteomx = trt_proteomx %>% select(-contains('072_3'))
}
if (any(grepl(pattern = 'Epi', x = colnames(trt_df)))) {
  trt_proteomx = trt_proteomx %>% select(-contains('072'))
}

cor_trt_prot = apply_2D(df = select(trt_proteomx, starts_with(TrT_miF)),
                        FUN = cor.test, complete_cases = comp_cas*2, 
                        col.x = 'adsf', col.y = 'asdfa',
                        y = select(trt_proteomx, starts_with('Proteomics'))) %>%
  rownames_to_column() %>%
  transmute(estimate.cor_trt = 
              as.numeric(as.character(estimate.cor)), rowname) %>%
  column_to_rownames()

cor_tpm_prot = apply_2D(df = select(trt_proteomx, starts_with('target')),
                        FUN = cor.test, complete_cases = comp_cas*2, 
                        col.x = 'adsf', col.y = 'asdfa',
                        y = select(trt_proteomx, starts_with('Proteomics'))) %>%
  rownames_to_column() %>%
  transmute(estimate.cor_tpm = 
              as.numeric(as.character(estimate.cor)), rowname) %>%
  column_to_rownames()

trt_proteomx = merge.data.frame(x = rownames_to_column(trt_proteomx),
                                y = rownames_to_column(cor_trt_prot),
                                by = 'rowname', all = T) %>%
  merge.data.frame(y = rownames_to_column(cor_tpm_prot), by = 'rowname', 
                   all = T) %>%
  column_to_rownames()

test = trt_proteomx %>% rownames_to_column() %>%  
  filter(estimate.cor_trt > (0.4 + estimate.cor_tpm)) %>% column_to_rownames()

trt_the_range = trt_proteomx %>% select(starts_with(TrT_miF)) %>%
  select(contains('The')) %>% na.omit() %>% apply(1, range, na.rm = T) %>%
  t() %>% as.data.frame() %>% rownames_to_column() %>% 
  filter(!is.infinite(V1)) %>% rename(min_the = V1, max_the = V2) %>% 
  column_to_rownames()

trt_tox_range = trt_proteomx %>% select(starts_with(TrT_miF)) %>%
  select(contains('Tox')) %>% na.omit() %>% apply(1, range, na.rm = T) %>%
  t() %>% as.data.frame() %>% rownames_to_column() %>% 
  filter(!is.infinite(V1)) %>% rename(min_tox = V1, max_tox = V2) %>% 
  column_to_rownames()

a = merge.data.frame(rownames_to_column(trt_the_range),
                     rownames_to_column(trt_tox_range),
                     'rowname') %>% column_to_rownames()

b = a[a$min_the > a$max_tox,,F] %>% rownames()
test = trt_proteomx[b,,F]

trt_proteomx$p.adj_trt %>% 
  hist(breaks = seq(0,1,0.05), main = 'TrT', xlab = 'p. adjusted values')
trt_proteomx$p.adj_tpm %>% 
  hist(breaks = seq(0,1,0.05), main = 'TPM', xlab = 'p. adjusted values')
trt_proteomx$p.adj_prx %>% 
  hist(breaks = seq(0,1,0.05), main = 'Proteomics', xlab = 'p. adjusted values')


# Best cases --------------------------------------------------------------


# VVV
test = trt_proteomx %>% 
  rownames_to_column() %>% 
  filter(
    p.value_trt < 0.05,
    p.value_tpm > 0.05,
    p.value_prx < 0.05, 
    sign(statistic.t_prx) == sign(statistic.t_trt)
  ) %>% 
  column_to_rownames()
# test %>% nrow() %>% print()

# AAA
test = trt_proteomx %>% 
  rownames_to_column() %>% 
  filter(
    p.adj_trt < 0.05,
    p.adj_tpm > 0.05,
    p.adj_prx < 0.05,
    sign(statistic.t_prx) == sign(statistic.t_trt)
  ) %>% 
  column_to_rownames()
# test %>% nrow() %>% print()

# AVV
test = trt_proteomx %>% 
  rownames_to_column() %>% 
  filter(
    p.adj_trt < 0.05,
    p.value_tpm > 0.05,
    p.value_prx < 0.05,
    sign(statistic.t_prx) == sign(statistic.t_trt)
  ) %>% 
  column_to_rownames()
test %>% nrow() %>% print()

# AVA
test = trt_proteomx %>% 
  rownames_to_column() %>% 
  filter(
    p.adj_trt < 0.05,
    p.value_tpm > 0.05,
    p.adj_prx < 0.05,
    sign(statistic.t_prx) == sign(statistic.t_trt)
  ) %>% 
  column_to_rownames()
test %>% nrow() %>% print()

# VAV
test = trt_proteomx %>% 
  rownames_to_column() %>% 
  filter(
    p.value_trt > 0.05,
    p.adj_tpm < 0.05,
    p.value_prx > 0.05
  ) %>% 
  column_to_rownames()
test %>% nrow() %>% print()

for (rown in rownames(test)) {
  
  row_data = test[rown, ]
  
  par(mfrow = c(1, 2))
  viewPlotProtx(test, rown)
  max_tpm = test[rown, , F] %>% select(starts_with('target')) %>% max(na.rm = T)
  plotOmics(df = trt_proteomx, omics = paste0('^', TrT_miF), 
            transcript = rown, xaxt = 'n', type = 'b', 
            ylab = 'black: TrT // pink: TPM', xlab = 'NULL',
            ylim = c(0, max_tpm))
  plotOmics(df = trt_proteomx, omics = '^target', transcript = rown, xaxt = 'n',
            type = 'b', col = 'pink', yaxt = 'n')
  readline(prompt = "Press [enter] to continue")
  
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  boxplot.dose(row_data, 'Proteomics', main = rown, ylab = 'Proteomics Expression')
  boxplot.dose(row_data, '^target', main = rown, ylab = 'TPM Expression')
  boxplot.dose(row_data, TrT_miF, main = rown, ylab = 'TrT Expression')
  
  
  test[rown, , F] %>% select(starts_with(TrT_miF)) %>% select(contains('The')) %>%
    as.numeric() %>% mean(na.rm = T) %>% abline(h = .)
  
  test[rown, , F] %>% select(starts_with(TrT_miF)) %>% select(contains('Tox')) %>%
    as.numeric() %>% mean(na.rm = T) %>% abline(h = .)
  axis(1, at = 1:24, labels = colnames(trt_proteomx)[grep('^target',
  colnames(trt_df))], las = 2)
  readline(prompt = "Press [enter] to continue")
  
}

rown = grep('00000100345', rownames(trt_proteomx)) %>% rownames(trt_proteomx)[.]
