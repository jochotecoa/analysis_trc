source('/share/script/hecatos/juantxo/analysis_trc/functions.R')
library(dplyr)
library(biomaRt)

# miRNA expression of 451 -------------------------------------------------

setwd('/share/analysis/hecatos/juantxo/Score/input/miRNA_miRge2/UNTR/')
mirna = read.csv('miR.Counts_voom.csv', row.names = 1)
mir451_expr = mirna %>% 
  filter(grepl(pattern = '451a', x = rownames(mirna)))

# Gene expression of MIF --------------------------------------------------

setwd('/share/analysis/hecatos/juantxo/Score/input/RNA_Salmon/UNTR/')
salmon = read.table(file = 'total_quant_voom.sf', header = T)
salmon_crna = salmon[grepl(x = salmon$Name, pattern = 'ENST'), ]
rownames(salmon_crna) = salmon_crna$Name
salmon_crna = salmon_crna[, -1]
salmon_gene = transcrToGene(table = salmon_crna, aggregate = T)
salmon_MIF = salmon_gene %>% filter(grepl(pattern = 'ENSG00000276701', 
                                          x = ensembl_gene_id))

# Protein expression of MIF -----------------------------------------------

mart = openMart2018()
MIF_uniprot = getBM(attributes = 'uniprot_gn', filters = 'ensembl_gene_id', 
                    values = 'ENSG00000276701', mart = mart)
setwd('/ngs-data/data/hecatos/Cardiac/Con_UNTR/Protein/')
setwd('Proteomics_Analyses_Cardiac_UNTR_GeneData/')
prot = read.delim('Hecatos_Cardiac_Px_Untreated_pre-processed_renamed.txt') %>%
  cleanProtIds()
MIF_prot = prot[grep(MIF_uniprot$uniprot_gn[1], prot$uniprot_gn), ]

# TRC values --------------------------------------------------------------

setwd('/share/analysis/hecatos/juantxo/Score/output/Output_Run_mrna_SEPT2019/')
setwd('V3/output/UNTR/2019-11-04_11:56:17_UTC/TRCscore/')
trc = read.table(file = 'UNTR_002_1_TRCscore.txt')


# Plotting ----------------------------------------------------------------

salmon_MIF_num = salmon_MIF[, -1] %>% .[, seq(1, 12, by = 3)] %>% as.numeric()
plot(salmon_MIF_num)
ggplot2::ggsave('/share/script/hecatos/juantxo/analysis_trc/integrated_timeline/plots/MIF_Salmon_TPM.png')

mir451_expr_num = mir451_expr[, seq(1, 12, by = 3)] %>% as.numeric()
plot(mir451_expr_num)

cor(tpm_MIF_2, mir451_expr_num)



# running integrated_timeline...

read.table('/share/analysis/hecatos/juantxo/Score/output/Output_Run_mrna_SEPT2019/V3/output/UNTR/2019-11-04_11:56:17_UTC/TRCscore/UNTR_002_2_TRCscore.txt') %>% dplyr::select(contains('targetRNA_')) %>% head(, n = 1)

read.table('UNTR_002_1_quant/quant_voom.sf', header = T) %>% .[grep('ENST00000419783', .[, 1]), ]

trc_MIF = global_gene.table[grep(global_gene.table$ensembl_gene_id, pattern = 'ENSG00000276701'), ]
trc_MIF_2 =  trc_MIF[, grepl(pattern = 'TRC\\.', x = colnames(trc_MIF))]
trc_MIF_3 = trc_MIF_2[, seq(1, 12, by = 3)] %>% as.numeric()

tpm_MIF_2 = trc_MIF[, grepl(pattern = 'targetRNA_TPM', x = colnames(trc_MIF))] %>%
  .[, seq(1, 12, by = 3)] %>% as.numeric()
plot(trc_MIF_3, ylim = c(max(tpm_MIF_2), min(trc_MIF_3)))
plot(tpm_MIF_2, ylim = c(max(tpm_MIF_2), min(trc_MIF_3)))


