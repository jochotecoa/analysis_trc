source('/share/script/hecatos/juantxo/analysis_trc/functions.R')

setwd('/share/analysis/hecatos/juantxo/Score/analysis/DESeq2/')


counts_DEGs = read.table('Salmon_counts/list_DEGs.txt') 
TPM_DEGs = read.table('Salmon_TPM/list_DEGs.txt')
voom_TRC_DEGs = read.table('TRC_voom/list_DEGs.txt')
salmon_TRC_DEGs = read.table('TRC_no_voom/list_DEGs.txt') 
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
lis = list(proteins_DEGs$ensembl_gene_id, counts_DEGs$V1, TPM_DEGs$V1, voom_TRC_DEGs$V1, salmon_TRC_DEGs$V1)

sapply(lis, length)

venn.diagram(x = list(proteins_DEGs$ensembl_gene_id, counts_DEGs$V1, TPM_DEGs$V1, voom_TRC_DEGs$V1, salmon_TRC_DEGs$V1), 
             filename = 'vennon5.png', imagetype = 'png', 
             category.names = c('Protein', 'Counts', 'TPM', 'Voom_TRC', 'Salmon_TRC'))

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
