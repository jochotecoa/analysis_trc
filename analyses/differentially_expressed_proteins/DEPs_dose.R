# Functions & libraries ---------------------------------------------------


viewPlotsProtx = function(df, evalu, decr = F){
  if (any(grepl(pattern = 'rowname', colnames(df)))) {
    df = df %>% remove_rownames() %>% column_to_rownames()
    colnames(df)[1] %>% print()
  }
  row_data = order(as.numeric(as.character(df[, evalu])), decreasing = decr)
  for (i in row_data) {
    data = df[i, ]
    plot_data = grep('Proteomics_', colnames(df)) %>% data[.] %>% as.matrix()
    
    plot_data %>% barplot(las = 2, 
                          main = paste0(rownames(data), ' ; ',
                                        evalu, ' = ', 
                                        format.pval(data[evalu], 3)), 
                          ylab = 'Protein Expression')
    
    
    the_data = plot_data[grep('The', colnames(plot_data))]
    tox_data = plot_data[grep('Tox', colnames(plot_data))]
    mean_the_data = mean(the_data, na.rm = T)
    mean_tox_data = mean(tox_data, na.rm = T)
    sd_the_data = sd(the_data, na.rm = T)
    sd_tox_data = sd(tox_data, na.rm = T)
    thesh_up_the = mean_the_data + sd_the_data
    thesh_down_the = mean_the_data - sd_the_data
    thesh_up_tox = mean_tox_data + sd_tox_data
    thesh_down_tox = mean_tox_data - sd_tox_data
    
    abline(h = thesh_up_the, col = 'lightgreen', pch = 2)
    abline(h = mean(the_data, na.rm = T), col = 'cyan')
    abline(h = thesh_down_the, col = 'lightgreen')
    
    abline(h = thesh_up_tox, col = 'indianred')
    abline(h = mean(tox_data, na.rm = T), col = 'purple')
    abline(h = thesh_down_tox, col = 'indianred')
    
    abline(v = 14.5)
    for (l in seq(3.7, 85.1, 3.6)) abline(v = l, col = 'gray')
    
    readline(prompt = "Press [enter] to continue")
  }
  
}


setwd("/share/script/hecatos/juantxo/analysis_trc")
source('functions.R')
forceLibrary(c('tibble', 'dplyr'))

setwd('/share/analysis/hecatos/juantxo/Score/input/RNA_Salmon/')
biomart_table = readRDS('salmon_ids_biomart.rds')
colnames(biomart_table) = biomart_table[1, ] %>% unlist() %>% 
  droplevels() %>% as.character()
biomart_table = biomart_table[-1, ] 
biomart_table = biomart_table %>% filter(`Transcript type` == 'protein_coding')


### Renaming data #####
setwd('/ngs-data/data/hecatos/Cardiac/')

lett_comp = strsplit(comp, '') %>% unlist()
lett_comp[-1] = tolower(lett_comp[-1])
comps_sel = list.dirs(recursive = F)
for (letter in lett_comp) {
  comps_sel = grep(letter, comps_sel) %>% comps_sel[.]
}
if (length(comps_sel) == 0) {
  lett_comp = strsplit(comp, '') %>% unlist()
  comps_sel = list.dirs(recursive = F)
  for (letter in lett_comp) {
    comps_sel = grep(letter, comps_sel) %>% comps_sel[.]
  }
}
print(comps_sel)
setwd(comps_sel)
setwd('Protein')
prot_dir = list.files(pattern = 'Proteomics') 
stopifnot(length(prot_dir) == 1)
setwd(prot_dir)
suppressWarnings(rm(proteomx_log2, proteomx_log2_norm))
if (length(list.files(pattern = 'renamed-by-Juan.txt')) == 0) {
  proteomx_log2 = proteomx_log2_norm = 
    read.table('Results_Hecatos_Cardiac_Px_PTX_FC/Cardiac_PTX_log2_norm.txt', 
               sep = '\t', header = T) %>% 
    cleanProtIds() %>% 
    remove_rownames() %>% 
    column_to_rownames('uniprot_gn') %>% 
    select(-Row.Names)
  # select(-X, -contains('.1')) %>% 
  # select(contains('PTX')) 
  
  proteomx_log2 %>% naToZero() %>% as.matrix.data.frame() %>% heatmap()
  
  # proteomx_log2 %>% dist() %>% naToZero() %>% hclust() %>% plot()
  
  colnames(proteomx_log2) %>% 
    gsub(pattern = '_rep', replacement = '') %>% strsplit('_') %>% 
    as.data.frame() %>% t() %>% as.data.frame() %>% .[, 'V4'] %>% sort() %>% head()
  
  sample_order = colnames(proteomx_log2) %>% 
    gsub(pattern = '_rep', replacement = '') %>% strsplit('_') %>% 
    as.data.frame() %>% t() %>% as.data.frame() %>% .[, 'V4'] %>% order()
  
  
  proteomx_log2 = proteomx_log2[, sample_order]
  
  the_samples = data.frame(dose = rep('PTX_The', 24))
  the_samples$timepoint = 
    lapply(c('000', '002', '008', '024', '072', '168', '240', '336'), rep, 3) %>% 
    unlist()
  the_samples$replicate = rep(c('1', '2', '3'), 8)
  
  tox_samples = data.frame(dose = rep('PTX_Tox', 18))
  tox_samples$timepoint = 
    lapply(c('002', '008', '024', '072', '168', '240'), rep, 3) %>% 
    unlist()
  tox_samples$replicate = rep(c('1', '2', '3'), 6)
  
  sample_names = rbind.data.frame(the_samples, tox_samples) %>% 
    apply(1, paste, collapse = '_')
  
  colnames(proteomx_log2) = sample_names
  
  proteomx_log2 %>% naToZero() %>% as.matrix.data.frame() %>% heatmap()
  
  write.table(x = proteomx_log2, 
              file = 
                'Hecatos_Cardiac_Px_PTX_Ther_Tox_log2_norm-renamed-by-Juan.txt')
  
} else {
  proteomx_log2 = 
    read.table(list.files(pattern = 'renamed-by-Juan.txt'))
}

##### Analysis #####
colnames(proteomx_log2) = colnames(proteomx_log2) %>% gsub(pattern = 'XPTX', replacement = 'PTX')

prot_log_norm = filterSamplesBySeqDepth(proteomx_log2)
proteomx_log2 %>% apply(2, median, na.rm = T) %>% 
  barplot(las = 2, main = 'Median Proteomics Expression per Sample')
### Get differentially exp %>% sed proteins
the_cols = prot_log_norm %>% select(contains('The')) %>% 
  dplyr::select(matches('002|008|024|072')) %>% colnames()
tox_cols = prot_log_norm %>% select(contains('Tox')) %>% 
  dplyr::select(matches('002|008|024|072')) %>% colnames()
rm(res.df)
res.df = prot_log_norm %>% apply_2D(FUN = t.test, col.x = the_cols, 
                  col.y = tox_cols, paired = T, complete_cases = 12)
res.df$p.value = res.df$p.value %>% as.character() %>% as.numeric()
res.df$statistic.t_prx = res.df$statistic.t %>% as.character() %>% as.numeric()
res.df$p.value %>% na.omit() %>% hist(breaks = seq(0, 1, 0.05))
### Bind the p-values to the proteomics data 
proteomx = 2^prot_log_norm %>% 
  dplyr::select(matches('002|008|024|072'))
colnames(proteomx) = paste0('Proteomics_', colnames(proteomx))
proteomx_filt = filterSamplesBySeqDepth(proteomx) 
proteomx_pval = merge.data.frame(x = rownames_to_column(proteomx_filt), 
                                y = rownames_to_column(select(res.df, p.value, statistic.t_prx)), 
                                by = 'rowname', all = T) %>% 
  mutate(p.adj = p.adjust(p.value, "BH")) %>% 
  column_to_rownames() 

proteomx_pval$p.adj %>% as.numeric() %>% na.omit() %>% hist(breaks = seq(0, 1, 0.05))

proteomx_biom = merge.data.frame(x = rownames_to_column(proteomx_pval), 
                                 y = biomart_table, by.x = 'rowname', by.y = 'UniProtKB/Swiss-Prot ID') %>% 
  filter(`Transcript type` == 'protein_coding') %>% 
  rename(`UniProtKB/Swiss-Prot ID` = rowname)

proteomx_sign = proteomx_biom %>% filter(p.value < 0.05) 


# setwd('/share/analysis/hecatos/juantxo/Score/output/DOC/2020-03-09_13:09:20_UTC/TRCscore/factors_DOC/')
#saveRDS(proteomx_biom, file = 'proteomx_biom.rds')
# viewPlotsProtx(filter(proteomx_pval, p.adj < 0.05), evalu = 'p.adj', T)
