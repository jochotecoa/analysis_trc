source('/share/script/hecatos/juantxo/analysis_trc/functions.R')
library(dplyr)
library(biomaRt)

# Functions ---------------------------------------------------------------

getMirnaExpr <- function(file, mirna) {
  mirna_voom = read.csv(file, row.names = 1)
  mir451_expr_voom = mirna_voom %>% 
    filter(grepl(pattern = mirna, x = rownames(mirna_voom)))
}

getGeneExpr = function(file, gene) {
  if (typeof(file) == 'character') {
    salmon = read.table(file = file, header = T)
  } else {
    salmon = file
  }
  salmon_crna = salmon[grepl(x = salmon$Name, pattern = 'ENST'), ]
  rownames(salmon_crna) = salmon_crna$Name
  salmon_crna = salmon_crna[, -1]
  salmon_gene = transcrToGene(table = salmon_crna, aggregate = T)
  salmon_MIF = salmon_gene %>% filter(grepl(pattern = gene, 
                                            x = ensembl_gene_id))
}
cleanCircRna = function(x) {
  y = x[grepl(x = x[, 'Name'], pattern = 'ENST'), ]
}
onlyTPM = function(x, names = F) {
  if (any(grepl('targetRNA_TPM', colnames(x)))) {
    TPM_col = 'targetRNA_TPM'
  } else {
    TPM_col = 'TPM'
  }
  if (!names) {
    y = x[, grep(TPM_col, colnames(x))]
  } else {
    y = x[, c(1, grep(TPM_col, colnames(x)))]
  }
  
}
is.expressed = function(x) {
  rows = apply(x[, -1] >= 1, 1, any, na.rm = T)
  y = x[rows, ]
}

# miRNA expression of 451 -------------------------------------------------

setwd('/share/analysis/hecatos/juantxo/Score/input/miRNA_miRge2/UNTR/')
voom_451 = getMirnaExpr(file = 'miR.Counts_voom.csv', mirna = '451a')
orig_451 = getMirnaExpr(file = 'miR.Counts.csv', mirna = '451a') 

voom_451_triplicate_1 = voom_451 %>% 
  .[seq(1, ncol(.), 3)] %>% .[1:4] %>% as.numeric()
orig_451_triplicate_1 = orig_451 %>% 
  .[seq(1, ncol(.), 3)] %>% .[1:4] %>% as.numeric()
# Gene expression of MIF --------------------------------------------------

setwd('/share/analysis/hecatos/juantxo/Score/input/RNA_Salmon/UNTR/')
all_voom_2 = mergeFiles(files_patt = 'quant_voom_2.sf', by_col = "Name", all = T)
all_salmon = mergeFiles(files_patt = 'quant.sf', by_col = 'Name', all = T)

# Give me the gene IDs of these transcripts, but do not aggregate them
crna_voom_2_genes = transcrToGene(table = all_voom_2, aggregate = F)

# What is the expression in the fist sample of the MIF transcripts?
MIF_voom_row = grep('ENSG00000276701', crna_voom_2_genes$ensembl_gene_id)
MIF_rna_voom_2_1 = crna_voom_2_genes$`TPM_UNTR_002_1_quant/quant_voom_2.sf` %>% 
  .[MIF_row]
print(MIF_rna_voom_2_1)

# crna_voom_2_genes = crna_voom_2_genes %>%
#   mutate(Name = ensembl_gene_id) %>% rename(ensembl_gene_id = Name)
# 
# crna_voom_2 = onlyTPM(crna_voom_2_genes, names = T) #gene names lost
# grep('ENSG00000276701', crna_voom_2$ensembl_gene_id)
# 
# crna_voom_2 = is.expressed(crna_voom_2)
# print(grep('ENSG00000276701', crna_voom_2$ensembl_gene_id))
# 
# searchgene(crna_voom_2)
  
# Give me the gene IDs of these transcripts, but do not aggregate them
crna_salmon_genes = transcrToGene(table = all_salmon, aggregate = F)

# What is the expression in the fist sample of the MIF transcripts?
MIF_salmon_row = grep('ENSG00000276701', crna_salmon_genes$ensembl_gene_id)
MIF_salmon_2_1 = crna_salmon_genes$`TPM_UNTR_002_1_quant/quant.sf`[MIF_salmon_row]
print(MIF_salmon_2_1)

# crna_salmon_genes$Name = crna_salmon_genes$ensembl_gene_id
# colnames(crna_salmon_genes)[1] = 'ensembl_gene_id'
# crna_salmon = onlyTPM(crna_salmon_genes, names = T) #gene names lost
# print(grep('ENSG00000276701', crna_salmon$ensembl_gene_id))
# crna_salmon = is.expressed(crna_salmon)
# print(grep('ENSG00000276701', crna_salmon$ensembl_gene_id))
# searchgene(crna_salmon)

# crna_salmon = cleanCircRna(all_salmon) %>% transcrToGene(F) %>% onlyTPM(names = T) %>% is.expressed() %>% searchgene()
# tpm_voom_2 = onlyTPM(crna_voom_2, T) %>% transcrToGene()
# expr_voom_2 = is.expressed(tpm_voom_2) %>% transcrToGene()
# searchgene = function(x) {a = x[grep('ENSG00000276701', x[, 'ensembl_gene_id']), ];  return(a)}
# 
# b = onlyTPM(trc_all, T)
# 
# 
# summary(rowSums(crna_voom_2[, -1], na.rm = T))
# 
# MIF_voom_2 = getGeneExpr(file = 'total_quant_voom_2.sf', gene = 'ENSG00000276701')
# MIF_orig = getGeneExpr(file = all_salmon, gene = 'ENSG00000276701')
# 
# MIF_voom_2_expr_1 = MIF_voom_2[, grepl('TPM', colnames(MIF_voom_2))] %>% 
#   .[seq(1, ncol(.), 3)] %>% as.numeric()
# MIF_orig_expr_1 = MIF_orig[, grepl('TPM', colnames(MIF_orig))] %>% 
#   .[seq(1, ncol(.), 3)] %>% as.numeric()
# 
# cor(MIF_orig_expr_1, orig_451_expr_1)
# cor(MIF_voom_2_expr_1, voom_451_expr_1)
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
setwd('V3/output/UNTR/2019-11-06_17:29:23_UTC/')
setwd('TRCscore/')
trc_all = mergeFiles('TRCscore', row_names = T, all = T)
write.table(trc_all, 'trc_all.tsv', sep = '\t')
trc_all_genes = transcrToGene(trc_all, F)
genenames = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), 
                  filters = 'ensembl_gene_id', 
                  values = salmon_gene$ensembl_gene_id, 
                  mart = mart)
trc_all_genes2 = merge(trc_all_genes, genenames, by = 'ensembl_gene_id')
trc_all_genes2$external_gene_name[grep('MIF', trc_all_genes2$external_gene_name)]
trc_MIF = trc_all_genes %>% 
  filter(grepl(pattern = 'ENSG00000276701', x = ensembl_gene_id))

# Plotting ----------------------------------------------------------------

salmon_MIF_num = salmon_MIF[, -1] %>% .[, seq(1, 12, by = 3)] %>% as.numeric()
plot(salmon_MIF_num)
ggplot2::ggsave('/share/script/hecatos/juantxo/analysis_trc/integrated_timeline/plots/MIF_Salmon_TPM.png')

mir451_expr_num = mir451_expr[, seq(1, 12, by = 3)] %>% as.numeric()
plot(mir451_expr_num)

cor(salmon_MIF_num, mir451_expr_num)



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


