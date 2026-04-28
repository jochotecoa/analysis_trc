
# Functions ---------------------------------------------------------------


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
  
  newdf = list(gene_trt_untr, gene_trt_tox)
  if (is.null(attrib)) {
    names(newdf) = c(paste0(cols.x, attrib), paste0(cols.y, attrib))
  } else {
    names(newdf) = c(paste0(cols.x, '_', attrib), 
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



# Load data ---------------------------------------------------------------

proteomics = F

comp = '5FU'
source('/share/script/hecatos/juantxo/analysis_trc/p.values_multiomics/checkSeveralOmicsTrTFile_JOA.R')


# Find gene of interest ---------------------------------------------------

dose = 'The'

forceSetWd('/share/analysis/hecatos/juantxo/Score/analysis/GOrilla')


gene_details = trt_df_t.tests %>% 
  rownames_to_column() %>% 
  filter(grepl(118194, rowname))

dim(gene_details)


# Check expression levels -------------------------------------------------



par(mfrow = c(1, 2))

compare_boxpl(gene_details, 'targetRNA_TPM_', cols.x = 'UNTR', cols.y = dose)
compare_boxpl(gene_details, 'TrT_0.1_', cols.x = 'UNTR', cols.y = dose)

# Check attributes on TrT table -------------------------------------------


atts = gene_details %>% 
  dplyr::select(contains('UNTR_072_1')) %>% 
  colnames() %>% 
  gsub('_UNTR_072_1_TrT', '', .) %>% 
  .[-c(2, 17)]

for (variable in atts) {
  
  compare_boxpl(gene_details, variable, cols.x = 'UNTR', cols.y = dose)
  
  readline(prompt="Press [enter] to continue")
}



# miRTaRBase --------------------------------------------------------------

setwd('~/Downloads')
mirtarbase = read.csv('hsa_MTI.csv')

mirnas_related = mirtarbase %>% 
  filter(grepl('NKX2-5', Target.Gene),
         Support.Type != 'Non-Functional MTI',
         Support.Type != 'Non-Functional MTI (Weak)') %>% 
  dplyr::select(miRNA) %>% 
  unlist()

length(mirnas_related)



# Load compound miRNA expression ------------------------------------------
mirna_dir = '/share/analysis/hecatos/juantxo/Score/input/miRNA_miRge2/'

setwd(mirna_dir)
setwd(comp)
mirna_comp = read.csv('miR.RPM.csv')

rows_exp = NULL
for (variable in mirnas_related) {
  row_exp = grep(pattern = variable, x = mirna_comp$miRNA)
  
  if (length(row_exp) > 0) {
    rows_exp = c(rows_exp, row_exp)
  }
  
}

mirnas_related_exp_comp = mirna_comp[rows_exp, ] %>% 
  dplyr::select(matches('miRNA'),
                contains(dose),
                -matches('000|002_3|168|240|336'))

dim(mirnas_related_exp_comp)

# Load UNTR miRNA expression ----------------------------------------------

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

dim(mirnas_related_exp_UNTR)

mirnas_related_exp = merge.data.frame(mirnas_related_exp_UNTR, 
                                      mirnas_related_exp_comp, 
                                      by = 'miRNA') %>% 
  column_to_rownames('miRNA')


# Boxplot of miRNA expression between UNTR and compound -------------------



for (variable in rownames(mirnas_related_exp)) {
  df = mirnas_related_exp[variable,,F]
  compare_boxpl(df, cols.x = 'UNTR', cols.y = dose, main = variable)
  readline(prompt="Press [enter] to continue")
  
}

z = NULL

for (variable in rownames(mirnas_related_exp)) {
  df = mirnas_related_exp[variable,,F]
  x = df %>% dplyr::select(contains('UNTR')) %>% as.numeric() 
  y = df %>% dplyr::select(contains(dose)) %>% as.numeric() 
  
  a = NULL
  a = try(t.test(x, y, paired = T))
  if (class(a) == 'try-error') {
    a = NA
  } else {
    a = a$p.value
  }
  
  if (!is.na(a)) {
    compare_boxpl(df, cols.x = 'UNTR', cols.y = dose, main = variable)
  }
  readline(prompt="Press [enter] to continue")
  
  z = c(z, a)
}
