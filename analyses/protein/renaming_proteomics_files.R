source('/share/script/hecatos/juantxo/analysis_trc/functions.R')
forceLibrary(c('dplyr', 'tibble'))


setwd('/ngs-data/data/hecatos/Cardiac/Paclitaxel/Protein/Proteomics_Analysis_Cardiac_FC_TxvsT0_revised_workflow_Genedata/Results_Hecatos_Cardiac_Px_PTX_FC/')
prot_file = read.table('Cardiac_PTX_log2.txt', sep = '\t', header = T) %>%
  cleanProtIds() %>% 
  select(-Row.Names) %>% 
  remove_rownames() %>%
  column_to_rownames('uniprot_gn')

prot_file = colnames(prot_file) %>% strsplit('_') %>% as.data.frame() %>% 
  .[4, ] %>% order() %>% prot_file[, .]
colnames(prot_file) = colnames(prot_file) %>% strsplit('_') %>% as.data.frame() %>% t() %>% 
  as.data.frame() %>% mutate(V4 = rep(1:3, nrow(.)/3)) %>% 
  apply(1, paste, collapse = '_') %>% 
  gsub(pattern = 'Ther', replacement = 'The')
colnames(prot_file) = colnames(prot_file) %>% strsplit('_') %>% 
  as.data.frame() %>% t() %>% as.data.frame() %>%
  mutate(V3 = gsub('T', '00', V3)) %>% 
  mutate(V3 = gsub('00240', '240', V3)) %>% 
  mutate(V3 = gsub('0024', '024', V3)) %>% 
  mutate(V3 = gsub('0072', '072', V3)) %>% 
  mutate(V3 = gsub('00168', '168', V3)) %>% 
  mutate(V3 = gsub('00336', '336', V3)) %>% 
  apply(1, paste, collapse = '_')
setwd('../../../Protein_JOA/')
write.table(x = prot_file, file = 'Cardiac_PTX_log2_renamedbyJuan.txt', sep = '\t')
