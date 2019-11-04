source('/share/script/hecatos/juantxo/analysis_trc/functions.R')
library(dplyr)
library(biomaRt)

setwd('/share/analysis/hecatos/juantxo/Score/input/miRNA_miRge2/UNTR/')
mirna = read.csv('miR.Counts_voom.csv', row.names = 1)
mir451_expr = mirna %>% 
  filter(grepl(pattern = '451a', x = rownames(mirna)))

setwd('/share/analysis/hecatos/juantxo/Score/input/RNA_Salmon/UNTR/')
salmon = getSalmonCols(cols = 'TPM', salfiles = '/quant_voom.sf')
salmon$enst = read.table('UNTR_002_1_quant/quant_voom.sf', header = T) %>%
  pull(Name)
salmon_crna = salmon[grepl(x = salmon$enst, pattern = 'ENST'), ]
rownames(salmon_crna) = salmon_crna$enst
salmon_crna = salmon_crna[, -ncol(salmon_crna)]
salmon_gene = transcrToGene(table = salmon_crna, aggregate = T)
# salmon_as = aggregate.data.frame(x = salmon_gene[, 2:13], by = list(salmon_gene$ensembl_gene_id), FUN = sum)
# head(salmon_gene)

salmon_MIF = salmon_gene %>% filter(grepl(pattern = 'ENSG00000276701', x = ensembl_gene_id))

mart = openMart2018()

MIF_uniprot = getBM(attributes = 'uniprot_gn', filters = 'ensembl_gene_id', 
                    values = 'ENSG00000276701', mart = mart)

prot = read.delim('/ngs-data/data/hecatos/Cardiac/Con_UNTR/Protein/Proteomics_Analyses_Cardiac_UNTR_GeneData/Hecatos_Cardiac_Px_Untreated_pre-processed_renamed.txt') %>%
  cleanProtIds()
head(prot$uniprot_gn)

MIF_prot = prot[grep(MIF_uniprot$uniprot_gn[1], prot$uniprot_gn), ]

salmon_MIF_num = salmon_MIF[, -1] %>% .[, seq(1, 12, by = 3)] %>% as.numeric()
plot(salmon_MIF_num)
ggplot2::ggsave('/share/script/hecatos/juantxo/analysis_trc/integrated_timeline/plots/MIF_Salmon_TPM.png')
mir451_expr_num = mir451_expr[, seq(1, 12, by = 3)] %>% as.numeric()
plot(mir451_expr_num)
cor(tpm_MIF_2, mir451_expr_num)

# running integrated_timeline...

trc = read.table('/share/analysis/hecatos/juantxo/Score/output/Output_Run_mrna_SEPT2019/V3/output/UNTR/2019-11-04_09:29:30_UTC/TRCscore/UNTR_002_2_TRCscore.txt')
trc = trc %>% dplyr::select(contains('targetRNA_'))

trc_MIF = global_gene.table[grep(global_gene.table$ensembl_gene_id, pattern = 'ENSG00000276701'), ]
trc_MIF_2 =  trc_MIF[, grepl(pattern = 'TRC\\.', x = colnames(trc_MIF))]
trc_MIF_3 = trc_MIF_2[, seq(1, 12, by = 3)] %>% as.numeric()

tpm_MIF_2 = trc_MIF[, grepl(pattern = 'targetRNA_TPM', x = colnames(trc_MIF))] %>%
  .[, seq(1, 12, by = 3)] %>% as.numeric()
plot(trc_MIF_3, ylim = c(max(tpm_MIF_2), min(trc_MIF_3)))
plot(tpm_MIF_2, ylim = c(max(tpm_MIF_2), min(trc_MIF_3)))

# View(listAttributes(mart))
# 
# salmon_gene_names = getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
#                           filters = 'ensembl_gene_id', 
#                           values = salmon_gene$ensembl_gene_id, 
#                           mart = mart)
# 
# salmon_gene_names[grep('MIF', salmon_gene_names$external_gene_name), ]
