source('/share/script/hecatos/juantxo/analysis_trc/functions.R')

setwd('/share/analysis/hecatos/juantxo/Score/analysis/DESeq2/')


counts_DEGs = read.table('Salmon_counts/list_DEGs.txt') 
TPM_DEGs = read.table('Salmon_TPM/list_DEGs.txt')
voom_TRC_DEGs = read.table('TRC_voom/list_DEGs.txt')
salmon_TRC_DEGs = read.table('TRC_no_voom/list_DEGs.txt') 


counts_DEGs[is.na.data.frame(counts_DEGs)] = F
TPM_DEGs[is.na.data.frame(TPM_DEGs)] = F
voom_TRC_DEGs[is.na.data.frame(voom_TRC_DEGs)] = F
salmon_TRC_DEGs[is.na.data.frame(salmon_TRC_DEGs)] = F

# mart = openMart2018()
# counts_ENSG = getBM(attributes = 'ensembl_gene_id', 
#                     filters = 'ensembl_gene_id_version', 
#                     values = counts, mart = mart)

proteins_DEPs = read.table('protein/list_DEGs.txt') 
proteins_DEPs_uniprot = cleanProtIds(proteins_DEPs)
mart = openMart2018()
proteins_DEGs = getBM(attributes = c('uniprot_gn', 'ensembl_gene_id'), filters = 'uniprot_gn', 
                  values = proteins_DEPs_uniprot$uniprot_gn, mart = mart)
proteins_DEGs = merge.data.frame(proteins_DEPs_uniprot, proteins_DEGs, by = 'uniprot_gn')
proteins_DEGs = proteins_DEGs %>% .[, -(1:2)] %>% select(ensembl_gene_id.y, everything())
proteins_DEGs = proteins_DEGs %>% mutate(ensembl_gene_id = ensembl_gene_id.y)
proteins_DEGs[is.na.data.frame(proteins_DEGs)] = F
lis = list(proteins_DEGs$ensembl_gene_id, counts_DEGs$V1, TPM_DEGs$V1, voom_TRC_DEGs$V1, salmon_TRC_DEGs$V1)

sapply(lis, length)

venn.diagram(x = list(counts_DEGs$ensembl_gene_id[counts_DEGs$timepoints_UNTR_008_vs_UNTR_002], 
                      TPM_DEGs$ensembl_gene_id[TPM_DEGs$timepoints_UNTR_008_vs_UNTR_002], 
                      voom_TRC_DEGs$ensembl_gene_id[voom_TRC_DEGs$timepoints_UNTR_008_vs_UNTR_002],
                      salmon_TRC_DEGs$ensembl_gene_id[salmon_TRC_DEGs$timepoints_UNTR_008_vs_UNTR_002], 
                      proteins_DEGs$ensembl_gene_id[proteins_DEGs$timepoints_UNTR_008_vs_UNTR_002]), 
             filename = 'vennon52.png', imagetype = 'png', 
             category.names = c('counts_DEGs', 'TPM_DEGs', 'voom_TRC_DEGs', 'salmon_TRC_DEGs', 'proteins_DEGs'), 
             width = 3000)



counts_n_DGEs = t(counts_DEGs) %>% .[-1, ] %>% apply(., 1, as.logical) %>%t() %>% apply(., 1, sum, na.rm = T) %>% data.frame(DGEs_n = .) %>%
  mutate(quant = 'Salmon_counts', comparisons = c('002vs008', '002vs024', '002vs072')) 
counts_n_DGEs_tmp = t(TPM_DEGs) %>% .[-1, ] %>% apply(., 1, as.logical) %>%t() %>% apply(., 1, sum, na.rm = T) %>% data.frame(DGEs_n = .) %>%
  mutate(quant = 'Salmon_TPM', comparisons = c('002vs008', '002vs024', '002vs072')) 
counts_n_DGEs = rbind.data.frame(counts_n_DGEs, counts_n_DGEs_tmp)
counts_n_DGEs_tmp = t(voom_TRC_DEGs) %>% .[-1, ] %>% apply(., 1, as.logical) %>%t() %>% apply(., 1, sum, na.rm = T) %>% data.frame(DGEs_n = .) %>%
  mutate(quant = 'voom_TRC', comparisons = c('002vs008', '002vs024', '002vs072')) 
counts_n_DGEs = rbind.data.frame(counts_n_DGEs, counts_n_DGEs_tmp)
counts_n_DGEs_tmp = t(salmon_TRC_DEGs) %>% .[-1, ] %>% apply(., 1, as.logical) %>%t() %>% apply(., 1, sum, na.rm = T) %>% data.frame(DGEs_n = .) %>%
  mutate(quant = 'no_voom_TRC', comparisons = c('002vs008', '002vs024', '002vs072')) 
counts_n_DGEs = rbind.data.frame(counts_n_DGEs, counts_n_DGEs_tmp)
counts_n_DGEs_tmp = t(proteins_DEGs) %>% .[-1, ] %>% apply(., 1, as.logical) %>%t() %>% apply(., 1, sum, na.rm = T) %>% data.frame(DGEs_n = .) %>%
  mutate(quant = 'protein', comparisons = c('002vs008', '002vs024', '002vs072')) 
counts_n_DGEs = rbind.data.frame(counts_n_DGEs, counts_n_DGEs_tmp)
ggpubr::ggpaired(counts_n_DGEs, x = 'quant', y = 'DGEs_n', id = 'comparisons')
