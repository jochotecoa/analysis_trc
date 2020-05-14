#### Transcriptomics ####
source('/share/script/hecatos/juantxo/analysis_trc/functions.R')
forceLibrary(c('corrplot', 'RColorBrewer'))

# Load transcriptomics data
setwd('/share/analysis/hecatos/juantxo/Score/input/RNA_Salmon/UNTR/')
tpm_all = read.table('total_quant.sf', header = T) %>% 
  dplyr::select(ensembl_transcript_id_version = Name, contains('TPM'))

# Only protein-coding transcripts
mart = openMart2018()
enst_pc = getBM(attributes = c('ensembl_transcript_id_version', 'transcript_biotype'), 
                filters = 'ensembl_transcript_id_version', 
                values = tpm_all$ensembl_transcript_id_version, 
                mart = mart)
enst_pc = enst_pc[enst_pc$transcript_biotype == 'protein_coding', 
                  'ensembl_transcript_id_version', 
                  drop = F]
tpm_pc = merge.data.frame(x = tpm_all, y = enst_pc, by = 'ensembl_transcript_id_version')

# Aggregate transcripts by protein
enst_unip = getBM(attributes = c('ensembl_transcript_id_version', 'uniprot_gn'), 
                  filters = 'ensembl_transcript_id_version', 
                  values = tpm_all$ensembl_transcript_id_version, 
                  mart = mart)
tpm_uniprot = merge.data.frame(x = tpm_pc, y = enst_unip, 
                               by = 'ensembl_transcript_id_version') %>%
  dplyr::select(-ensembl_transcript_id_version) 
tpm_uniprot = aggregate.data.frame(x = tpm_uniprot[, -13], by = list(tpm_uniprot[, 13]), FUN = sum) %>%
  mutate(uniprot_gn = Group.1) %>%
  dplyr::select(-Group.1)

#### Proteomics ####
setwd('/ngs-data/data/hecatos/Cardiac/Con_UNTR/Protein/')
setwd('Proteomics_Analyses_Cardiac_UNTR_GeneData/')
protein_table = read.table(
  file = 'Hecatos_Cardiac_Px_Untreated_pre-processed_renamed.txt', 
  sep = '\t', header = T) %>%
  cleanProtIds() %>%
  dplyr::select(-Row.Names)
 
#protein_table[, -ncol(protein_table)] = protein_table[, -ncol(protein_table)] %>% filterSamplesBySeqDepth
 
 median_of_medians = protein_table %>% 
 	remove_rownames %>% 
 	column_to_rownames('uniprot_gn') %>% 
 	apply(2, median, na.rm = T) %>% 
 	median(na.rm = T)
 
 for (col in colnames(protein_table)[-ncol(protein_table)]) {
 	median_col = protein_table[, col] %>% median(na.rm = T)
 	correction = median_of_medians / median_col
 #	protein_table[, col] = protein_table[, col] * correction
 }

#### Correlation ####
tpm_prot = merge.data.frame(x = tpm_uniprot, y = protein_table, 
                            by = 'uniprot_gn') %>% 
  column_to_rownames('uniprot_gn')

colnames(tpm_prot) = gsub('_quant.*', '', colnames(tpm_prot))		

estimates = pvalues = NULL
for (col.x in seq(1, 12, 3)) {
  for (row.x in 1:nrow(tpm_prot)) {
    x = tpm_prot[row.x, col.x:(col.x+2)] %>%
      as.numeric() %>%
      na.omit() 
    
    y = tpm_prot[row.x, 22:24] %>%
      as.numeric() %>%
      na.omit() 
    if (length(x) == 0 | length(y) == 0) {
      estimates = c(estimates, NA)
      pvalues = c(pvalues, NA)
    } else {
      corr = cor.test(x, y)
      estimates = c(estimates, corr$estimate)
      pvalues = c(pvalues, corr$p.value)
    }
    estimates_total = 
  }
  
}
estimates072 = estimates

tpm_prot_2 = data.frame(TPM_002 = unlist(tpm_prot[, 1:3]))
tpm_prot_2$TPM_008 = unlist(tpm_prot[, 4:6])
tpm_prot_2$TPM_024 = unlist(tpm_prot[, 7:9])
tpm_prot_2$TPM_072 = unlist(tpm_prot[, 10:12])
tpm_prot_2$Protein_002 = unlist(tpm_prot[, 13:15])
tpm_prot_2$Protein_008 = unlist(tpm_prot[, 16:18])
tpm_prot_2$Protein_024 = unlist(tpm_prot[, 19:21])
tpm_prot_2$Protein_072 = unlist(tpm_prot[, 22:24])
tpm_prot_2$Protein_168 = unlist(tpm_prot[, 25:27])
tpm_prot_2$Protein_240 = unlist(tpm_prot[, 28:30])
tpm_prot_2$Protein_336 = unlist(tpm_prot[, 31:33])

M = tpm_prot_2 %>% na.omit() %>% cor()
corrplot(M, type="upper", order="original",
         col=brewer.pal(n=8, name="RdYlBu"), addCoef.col = "black")

corrplot(cor(na.omit(tpm_prot)), type="upper", order="original",
         col=brewer.pal(n=8, name="RdYlBu"))


cor.test(c(tpm_prot[, 10], tpm_prot[, 11], tpm_prot[, 12]), c(tpm_prot[, 22], tpm_prot[, 23], tpm_prot[, 24]))

