source('/share/script/hecatos/juantxo/analysis_trc/functions.R')

setwd('/share/analysis/hecatos/juantxo/Score/analysis/DESeq2/')

forceLibrary(c('VennDiagram', 'dplyr'))

counts_DEGs_deseq2 = read.table(
  '/share/analysis/hecatos/juantxo/Score/analysis/DESeq2/Salmon_counts/list_DEGs.txt') 
TPM_DEGs_deseq2 = read.table(
  '/share/analysis/hecatos/juantxo/Score/analysis/DESeq2/Salmon_TPM_nofilt/list_DEGs.txt')
TRC_DEGs_deseq2 = read.table(
  '/share/analysis/hecatos/juantxo/Score/analysis/DESeq2/TRC_nofilt/list_DEGs.txt') 
TRC_counts_DEGs_deseq2 = read.table(
  '/share/analysis/hecatos/juantxo/Score/analysis/DESeq2/TRC_counts/list_DEGs.txt')

TRC_DEGs_edgeR = read.table('/share/analysis/hecatos/juantxo/Score/analysis/edgeR/TRC_nofilt/DEG_table_BH_timepointUNTR_008_glm_.tsv')
TPM_DEGs_edgeR = read.table('/share/analysis/hecatos/juantxo/Score/analysis/edgeR/Salmon_TPM/DEG_table_BH_timepointUNTR_008_glm_.tsv')


counts_DEGs_deseq2[is.na.data.frame(counts_DEGs_deseq2)] = F
TPM_DEGs_deseq2[is.na.data.frame(TPM_DEGs_deseq2)] = F
TRC_DEGs_deseq2[is.na.data.frame(TRC_DEGs_deseq2)] = F
TRC_counts_DEGs_deseq2[is.na.data.frame(TRC_counts_DEGs_deseq2)] = F

counts_DEGs_deseq2_2_8 = counts_DEGs_deseq2 %>% 
  dplyr::select(ensembl_gene_id, timepoints_UNTR_008_vs_UNTR_002) %>% 
  filter(timepoints_UNTR_008_vs_UNTR_002) %>%
   dplyr::select(-timepoints_UNTR_008_vs_UNTR_002)

TRC_counts_DEGs_deseq2_2_8 = TRC_counts_DEGs_deseq2 %>%
  dplyr::select(ensembl_gene_id, timepoints_UNTR_008_vs_UNTR_002) %>% 
  filter(timepoints_UNTR_008_vs_UNTR_002) %>%
  dplyr::select(-timepoints_UNTR_008_vs_UNTR_002)

TPM_DEGs_deseq2_2_8 = TPM_DEGs_deseq2 %>% 
   dplyr::select(ensembl_gene_id, timepoints_UNTR_008_vs_UNTR_002) %>% 
  filter(timepoints_UNTR_008_vs_UNTR_002) %>%
   dplyr::select(-timepoints_UNTR_008_vs_UNTR_002)

TRC_DEGs_deseq2_2_8 = TRC_DEGs_deseq2 %>%
   dplyr::select(ensembl_gene_id, timepoints_UNTR_008_vs_UNTR_002) %>% 
  filter(timepoints_UNTR_008_vs_UNTR_002) %>%
   dplyr::select(-timepoints_UNTR_008_vs_UNTR_002)


TRC_DEGs_edgeR_2_8 = TRC_DEGs_edgeR %>%
  rownames_to_column() %>%
  filter(PValue < 0.05) %>%
   dplyr::select(rowname)

TPM_DEGs_edgeR_2_8 = TPM_DEGs_edgeR %>%
  rownames_to_column() %>%
  filter(PValue < 0.05) %>%
  dplyr::select(rowname)


# mart = openMart2018()
# counts_ENSG = getBM(attributes = 'ensembl_gene_id', 
#                     filters = 'ensembl_gene_id_version', 
#                     values = counts, mart = mart)

proteins_DEPs = read.table('/share/analysis/hecatos/juantxo/Score/analysis/DEPs/UNTR_protein_log2_pvalues.tsv') %>%
  tibble::rownames_to_column() %>% mutate(uniprot_gn = rowname)


# proteins_DEPs_uniprot = cleanProtIds(proteins_DEPs)
proteins_DEGs = openMart2018() %>%
  getBM(attributes = c('uniprot_gn', 'ensembl_gene_id'), 
        filters = 'uniprot_gn', 
        values = proteins_DEPs$uniprot_gn, mart = .) %>%
  merge.data.frame(proteins_DEPs, by = 'uniprot_gn')

proteins_DEGs_2_8 = proteins_DEGs %>% 
  filter(p.value_UNTR_The_008_1 < 0.05) %>%
  dplyr::select(ensembl_gene_id)
  
proteins_DEGs[is.na.data.frame(proteins_DEGs)] = F
# lis = list(proteins_DEGs$ensembl_gene_id, counts_DEGs_deseq2$V1, TPM_DEGs_deseq2$V1, TRC_DEGs_edgeR$V1, TRC_DEGs_deseq2$V1)
# 
# sapply(lis, length)

venn.diagram(x = list(TPM_DEGs_deseq2_2_8$ensembl_gene_id, 
                      TRC_DEGs_deseq2_2_8$ensembl_gene_id,
                      TPM_DEGs_edgeR_2_8$rowname,
                      TRC_DEGs_edgeR_2_8$rowname),
             filename = 'TPMvsTRC&DESeq2vsedgeR.png', imagetype = 'png', 
             category.names = c('TPM_DESeq2', 'TRC_DESeq2', 'TPM_edgeR', 'TRC_edgeR'), 
             width = 3000)



counts_n_DGEs = t(counts_DEGs_deseq2) %>% 
  .[-1, ] %>% 
  apply(., 1, as.logical) %>%
  t() %>% apply(., 1, sum, na.rm = T) %>% 
  data.frame(DGEs_n = .) %>%
  mutate(quant = 'Salmon_counts', comparisons = c('002vs008', '002vs024', '002vs072')) 
counts_n_DGEs_tmp = t(TPM_DEGs_deseq2) %>% .[-1, ] %>% apply(., 1, as.logical) %>%t() %>% apply(., 1, sum, na.rm = T) %>% data.frame(DGEs_n = .) %>%
  mutate(quant = 'Salmon_TPM', comparisons = c('002vs008', '002vs024', '002vs072')) 
counts_n_DGEs = rbind.data.frame(counts_n_DGEs, counts_n_DGEs_tmp)
# counts_n_DGEs_tmp = t(TRC_DEGs_edgeR) %>% .[-1, ] %>% apply(., 1, as.logical) %>%t() %>% apply(., 1, sum, na.rm = T) %>% data.frame(DGEs_n = .) %>%
#   mutate(quant = 'voom_TRC', comparisons = c('002vs008', '002vs024', '002vs072')) 
# counts_n_DGEs = rbind.data.frame(counts_n_DGEs, counts_n_DGEs_tmp)
counts_n_DGEs_tmp = t(TRC_DEGs_deseq2) %>% .[-1, ] %>% apply(., 1, as.logical) %>%t() %>% apply(., 1, sum, na.rm = T) %>% data.frame(DGEs_n = .) %>%
  mutate(quant = 'no_voom_TRC', comparisons = c('002vs008', '002vs024', '002vs072'))
counts_n_DGEs = rbind.data.frame(counts_n_DGEs, counts_n_DGEs_tmp)
# counts_n_DGEs_tmp = t(proteins_DEGs) %>% .[-1, ] %>% apply(., 1, as.logical) %>%t() %>% apply(., 1, sum, na.rm = T) %>% data.frame(DGEs_n = .) %>%
#   mutate(quant = 'protein', comparisons = c('002vs008', '002vs024', '002vs072')) 
# counts_n_DGEs = rbind.data.frame(counts_n_DGEs, counts_n_DGEs_tmp)
ggpubr::ggpaired(counts_n_DGEs, x = 'quant', y = 'DGEs_n', id = 'comparisons')

DEGs = c(nrow(counts_DEGs_deseq2_2_8), 
         nrow(TPM_DEGs_deseq2_2_8), 
         nrow(TRC_DEGs_deseq2_2_8), 
         nrow(TPM_DEGs_edgeR_2_8), 
         nrow(TRC_DEGs_edgeR_2_8), 
         nrow(proteins_DEGs_2_8)) %>%
  data.frame(n_DEGs = .) %>%
  mutate(names = c('Counts_DESeq2', 
                   'TPM_DESeq2', 
                   'TRC_DESeq2', 
                   'TPM_edgeR', 
                   'TRC_edgeR', 
                   'Protein_DEGs'))

ggplot(DEGs, aes(names, n_DEGs)) + 
  geom_bar(stat = "identity") + 
  labs(x = 'Methods', y = 'Number of DEGs')
