
# Functions ---------------------------------------------------------------

boxplot.dose = function(df, omics, ...) {
  x = df %>% select(matches(omics))
  the = x %>% select(contains('The')) %>% t()
  tox = x %>% select(contains('Tox')) %>% t()
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
    select(matches('Proteomics')) %>% 
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
  
  res.df = df %>% select(the_cols, tox_cols) %>% zeroToNa() %>% log2() %>% 
    apply_2D(col.x = the_cols, col.y = tox_cols, ...)
  
  if (plot) {
	res.df$p.value_tpm %>% 
		as.character() %>%
		as.numeric() %>%
		na.omit() %>% 
		hist(breaks = seq(0, 1, 0.05), main = 'TPM', xlab = 'p. value')
	}  
  
  res.df[, paste0('p.value_', title)] = res.df$p.value %>% as.character() %>% as.numeric()
  res.df[, paste0('statistic.t_', title)] = res.df$statistic.t %>% as.character() %>% as.numeric()
  
  return(res.df)
  
}

forceLibrary(c('dplyr', 'tibble'))
setwd("/share/script/hecatos/juantxo/analysis_trc")
source('functions.R')

# Get data ----------------------------------------------------------------

comp = 'PTX'
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



# Clean data -------------------------------------------------------------------

trt_df_ori = trt_df

colnames(trt_df) = colnames(trt_df) %>% 
  gsub(pattern = '.txt', replacement = '') %>%
  gsub(pattern = '.rds', replacement = '') %>%
  gsub(pattern = '_TrTscore', replacement = '')
  
if (plotting) {
	trt_df %>% select(starts_with('target')) %>% 
	  apply(2, median, na.rm = T) %>% 
	  barplot(las = 2, main = 'Median TPM expression per Sample')
}

comp_cas = 12

# Exceptions
if (any(grepl(pattern = 'Dox', x = colnames(trt_df)))) {
  trt_df = trt_df %>% select(-contains('072_3'))
  comp_cas = 11
}
if (any(grepl(pattern = 'Epi', x = colnames(trt_df)))) {
  trt_df = trt_df %>% select(-contains('072'))
  comp_cas = 9
}

# Remove thousands of useless columns
trt_vars = colnames(trt_df) %>% .[grep('TrT_', .)] %>% .[!grepl(TrT_miF, .)]
trt_df = trt_df[, !(colnames(trt_df) %in% trt_vars)]
trt_df = trt_df[, !grepl('hsa', colnames(trt_df))]

trt_df_geneid = merge.data.frame(x = rownames_to_column(trt_df), 
                          y = select(.data = proteomx_biom, 
                                     `Transcript stable ID`, `Gene stable ID`), 
                          by.x = 'rowname',
                          by.y = 'Transcript stable ID') %>% 
  filter(!duplicated(rowname)) %>% 
  column_to_rownames()

trt_df_geneid = trt_df_geneid %>% group_by(`Gene stable ID`) %>% 
  summarise_if(is.numeric, sum, na.rm = T) %>% 
  column_to_rownames('Gene stable ID')

# mutate(p.adj_trt = p.adjust(p.value_trt, method = 'BH'), 
#        p.adj_tpm = p.adjust(p.value_tpm, method = 'BH')) %>% 
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



trt_df_geneid = merge.data.frame(x = rownames_to_column(trt_df_geneid), 
                              y = rownames_to_column(select(theVStox_t.test_trt, 
                                                            p.value_trt, 
                                                            statistic.t_trt)), 
                              by = 'rowname', all = T) %>% 
  merge.data.frame(y = rownames_to_column(select(theVStox_t.test_tpm, p.value_tpm, statistic.t_tpm)), 
                   by = 'rowname', all = T) %>% 
  column_to_rownames()

trt_proteomx = merge.data.frame(x = rownames_to_column(trt_df_geneid), 
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

# Correlation ------------------------------------------------------------------

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

# Minimum > Maximum ------------------------------------------------------------

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

if (plotting) {
	trt_proteomx$p.adj_trt %>% 
	  hist(breaks = seq(0,1,0.05), main = 'TrT', xlab = 'p. adjusted values')
	trt_proteomx$p.adj_tpm %>% 
	  hist(breaks = seq(0,1,0.05), main = 'TPM', xlab = 'p. adjusted values')
	trt_proteomx$p.adj_prx %>% 
	  hist(breaks = seq(0,1,0.05), main = 'Proteomics', xlab = 'p. adjusted values')
}


# Best cases --------------------------------------------------------------


# VVV
vvv = trt_proteomx %>% 
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
aaa = trt_proteomx %>% 
  rownames_to_column() %>% 
  filter(
    p.adj_trt < 0.05,
    p.adj_tpm > 0.05,
    p.adj_prx < 0.05,
    sign(statistic.t_prx) == sign(statistic.t_trt)
  ) %>% 
  column_to_rownames()


# AVV
avv = trt_proteomx %>% 
  rownames_to_column() %>% 
  filter(
    p.adj_trt < 0.05,
    p.value_tpm > 0.05,
    p.value_prx < 0.05,
    sign(statistic.t_prx) == sign(statistic.t_trt)
  ) %>% 
  column_to_rownames()

# AVV_bad
avv_bad = trt_proteomx %>% 
  rownames_to_column() %>% 
  filter(
    p.adj_trt < 0.05,
    p.value_tpm > 0.05,
    p.value_prx > 0.05,
    sign(statistic.t_prx) == sign(statistic.t_trt)
  ) %>% 
  column_to_rownames()

# AVA
ava = trt_proteomx %>% 
  rownames_to_column() %>% 
  filter(
    p.adj_trt < 0.05,
    p.value_tpm > 0.05,
    p.adj_prx < 0.05,
    sign(statistic.t_prx) == sign(statistic.t_trt)
  ) %>% 
  column_to_rownames()

# VAV
vav = trt_proteomx %>% 
  rownames_to_column() %>% 
  filter(
    p.value_trt > 0.05,
    p.adj_tpm < 0.05,
    p.value_prx > 0.05
  ) %>% 
  column_to_rownames()

# VAV_bad
vav_bad = trt_proteomx %>% 
  rownames_to_column() %>% 
  filter(
    p.value_trt > 0.05,
    p.adj_tpm < 0.05,
    p.value_prx < 0.05
  ) %>% 
  column_to_rownames()
 
# Sensitivity and Specificity --------------------------------------------------
 
 #### TPM ####
# TP
tp_tpm = trt_proteomx %>% 
  rownames_to_column() %>% 
  filter(
    p.adj_tpm < 0.05,
    p.value_prx < 0.05,
    sign(statistic.t_prx) == sign(statistic.t_tpm)
  ) %>% 
  column_to_rownames()
  
# TN
tn_tpm = trt_proteomx %>% 
  rownames_to_column() %>% 
  filter(
    p.value_tpm > 0.05,
    p.value_prx > 0.05,
  ) %>% 
  column_to_rownames()

# Falses

falses_tpm = rownames(trt_proteomx) %>% 
	setdiff(rownames(tp_tpm)) %>%
	setdiff(rownames(tn_tpm)) %>% 
	trt_proteomx[., ]
	
# FP
fp = falses_tpm %>% 
  rownames_to_column() %>% 
  filter(
    p.adj_tpm < 0.05
  ) %>% 
  column_to_rownames()

# FN
fn = falses_tpm %>% 
  rownames_to_column() %>% 
  filter(
    p.value_tpm > 0.05
  ) %>% 
  column_to_rownames()

#### TRT ####
# TP
tp_trt = trt_proteomx %>% 
  rownames_to_column() %>% 
  filter(
    p.adj_trt < 0.05,
    p.value_prx < 0.05,
    sign(statistic.t_prx) == sign(statistic.t_trt)
  ) %>% 
  column_to_rownames()
  
# TN
tn_trt = trt_proteomx %>% 
  rownames_to_column() %>% 
  filter(
    p.value_trt > 0.05,
    p.value_prx > 0.05,
  ) %>% 
  column_to_rownames()

# Falses

falses_trt = setdiff(rownames(trt_proteomx), rownames(tp_trt)) %>%
	setdiff(rownames(tn_trt)) %>% trt_proteomx[., ]
	
# FP
fp_trt = falses_trt %>% 
  rownames_to_column() %>% 
  filter(
    p.adj_trt < 0.05
  ) %>% 
  column_to_rownames()

# FN
fn_trt = falses_trt %>% 
  rownames_to_column() %>% 
  filter(
    p.value_trt > 0.05
  ) %>% 
  column_to_rownames()


test = avv

colnames(test) = colnames(test) %>% gsub('_TrT', '', .)

if (plotting) {
	for (rown in rownames(test)) {
	  layout(matrix(c(1,2), 1, 2, byrow = T))
	  row_data = test[rown, ]
	  max_tpm = test[rown, , F] %>% select(starts_with('target')) %>% max(na.rm = T)
	  #par(mfrow = c(1, 2))
	  
	  viewPlotProtx(test, rown)
	  
	  readline(prompt = "Press [enter] to continue")
	  
	  plotOmics(df = trt_proteomx, omics = paste0('^', TrT_miF), 
		    transcript = rown, xaxt = 'n', type = 'b', 
		    ylab = 'gray: TPM | black: TrT', xlab = '',
		    ylim = c(0, max_tpm))
	  plotOmics(df = trt_proteomx, omics = '^target', transcript = rown, xaxt = 'n',
		    type = 'b', col = 'gray', yaxt = 'n', xlab = '')
	  axis(1, at = 1:24, labels = colnames(test)[grep('^target', colnames(test))], las = 2)
	  
	  readline(prompt = "Press [enter] to continue")

	  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
	  boxplot.dose(row_data, 'Proteomics', main = rown, ylab = 'Proteomics Expression')
	  boxplot.dose(row_data, '^target', main = rown, ylab = 'TPM Expression')
	  boxplot.dose(row_data, TrT_miF, main = rown, ylab = 'TrT Expression')
	  
	  
	  #test[rown, , F] %>% select(starts_with(TrT_miF)) %>% select(contains('The')) %>%
	  #  as.numeric() %>% mean(na.rm = T) %>% abline(h = .)
	  
	#  test[rown, , F] %>% select(starts_with(TrT_miF)) %>% select(contains('Tox')) %>%
	 #   as.numeric() %>% mean(na.rm = T) %>% abline(h = .)
	 # axis(1, at = 1:24, labels = colnames(trt_proteomx)[grep('^target',
	  #colnames(trt_df))], las = 2)
	  readline(prompt = "Press [enter] to continue")
	  
	}
}

cat(paste(comp, '\n', 'TPM_false positives:', nrow(vav), '\n', 'TPM_false negatives:', nrow(avv), '\n', 'TrT_false positives:', nrow(avv_bad), '\n', 'TrT_false negatives:', nrow(vav_bad), '\n'))
nrow(tp)
nrow(tn)
nrow(fp)
nrow(fn)
nrow(tp_trt)
nrow(tn_trt)
nrow(fp_trt)
nrow(fn_trt)

