compare_boxpl <- function(df, attrib = NULL, cols.x, cols.y, ...) {
  if (!is.null(attrib)) {
    gene_trt = df %>% 
      dplyr::select(starts_with(attrib))
    
  } else {
    gene_trt = df
  }
  
  max_y = max(gene_trt, na.rm = T)
  
  gene_trt_untr = gene_trt %>% 
    dplyr::select(matches(cols.x)) %>% 
    unlist()
  
  gene_trt_tox = gene_trt %>% 
    dplyr::select(matches(cols.y)) %>% 
    unlist()
  
  newdf = cbind.data.frame(gene_trt_untr, gene_trt_tox)
  if (is.null(attrib)) {
    colnames(newdf) = c(paste0(cols.x, attrib), paste0(cols.y, attrib))
  } else {
    colnames(newdf) = c(paste0(cols.x, '_', attrib), 
                        paste0(cols.y, '_', attrib))
    
  }
  boxplot(newdf, 
          ylim = c(0, max_y), ...)
  
}
subSet = function(x, pattern) {
  x[grep(pattern = pattern, x = x)]
}


library('VennDiagram')
setwd("/share/script/hecatos/juantxo/analysis_trc")
source('functions_JOA.R')
forceLibrary(c('dplyr', 'tibble'))



if (!exists('trt_proteomx')) {
  comp = '5FU'
  source('/share/script/hecatos/juantxo/analysis_trc/p.values_multiomics/checkSeveralOmicsTrTFile_JOA.R')
}

forceSetWd('/share/analysis/hecatos/juantxo/Score/analysis/GOrilla')

gene_details = trt_proteomx %>% 
  filter(external_gene_name == 'PAFAH1B1')



# Check attributes on TrT table -------------------------------------------


compare_boxpl(gene_details, 'targetRNA_TPM_', cols.x = 'UNTR', cols.y = 'Tox')
compare_boxpl(gene_details, 'TrT_0.1_', cols.x = 'UNTR', cols.y = 'Tox')



atts = gene_details %>% 
  dplyr::select(contains('UNTR_072_1')) %>% 
  colnames() %>% 
  gsub('_UNTR_072_1_TrT', '', .) %>% 
  .[-c(2, 17)]

for (variable in atts) {
  
  compare_boxpl(gene_details, variable, cols.x = 'UNTR', cols.y = 'Tox')
  
  readline(prompt="Press [enter] to continue")
}


# Which miRNAs are contrasted to have an effect ---------------------------


setwd('~/Downloads')
mirtarbase = read.csv('hsa_MTI.csv')

mirnas_related = mirtarbase %>% 
  filter(Target.Gene == 'PAFAH1B1',
         Support.Type != 'Non-Functional MTI',
         Support.Type != 'Non-Functional MTI (Weak)') %>% 
  dplyr::select(miRNA) %>% 
  unlist()

dim(mirnas_related)

mirna_dir = '/share/analysis/hecatos/juantxo/Score/input/miRNA_miRge2/'

setwd(mirna_dir)
setwd('5FU')
mirna_5FU = read.csv('miR.RPM.csv')

rows_exp = NULL
for (variable in mirnas_related) {
  row_exp = grep(pattern = variable, x = mirna_5FU$miRNA)
  
  if (length(row_exp) > 0) {
    rows_exp = c(rows_exp, row_exp)
  }
  
}

mirnas_related_exp_5FU = mirna_5FU[rows_exp, ] %>% 
  dplyr::select(matches('Tox|miRNA'), 
                -matches('002_3|168'))


setwd(mirna_dir)
setwd('UNTR')
mirna_UNTR = read.csv('miR.RPM.csv')

rows_exp = NULL
for (variable in mirnas_related) {
  row_exp = grep(pattern = variable, x = mirna_UNTR$miRNA)
  
  if (length(row_exp) > 0) {
    rows_exp = c(rows_exp, row_exp)
  }
  
}

mirnas_related_exp_UNTR = mirna_UNTR[rows_exp, ] %>% 
  dplyr::select(-matches('002_3|168|240|336'))

mirnas_related_exp = merge.data.frame(mirnas_related_exp_UNTR, 
                                      mirnas_related_exp_5FU, 
                                      by = 'miRNA') %>% 
  column_to_rownames('miRNA')



for (variable in rownames(mirnas_related_exp)) {
  df = mirnas_related_exp[variable,,F]
  compare_boxpl(df, cols.x = 'UNTR', cols.y = 'Tox', main = variable)
  readline(prompt="Press [enter] to continue")
  
}

