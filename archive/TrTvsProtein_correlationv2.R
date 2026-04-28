#### Transcriptomics ####
source('/share/script/hecatos/juantxo/analysis_trc/functions.R')
forceLibrary(c('corrplot', 'RColorBrewer'))

# Load transcriptomics data
setwd('/share/analysis/hecatos/juantxo/Score/output/UNTR/2020-03-09_16:43:15__factor_0_1/')
TPM_all = readRDS('All_trt.rds') %>% 
	dplyr::select(ensembl_transcript_id = rowname, 
			contains('targetRNA_TPM')) #targetRNA_TPM TPM_0.1_UNTR_0

# Only protein-coding transcripts
mart = openMart2018()
enst_pc = getBM(attributes = c('ensembl_transcript_id', 'transcript_biotype'), 
                filters = 'ensembl_transcript_id', 
                values = TPM_all$ensembl_transcript_id, 
                mart = mart) %>%
                filter(transcript_biotype == 'protein_coding') %>%
                dplyr::select(ensembl_transcript_id)

TPM_pc = merge.data.frame(x = TPM_all, y = enst_pc, by = 'ensembl_transcript_id')

# Aggregate transcripts by protein
enst_unip = getBM(attributes = c('ensembl_transcript_id', 'uniprot_gn'), 
                  filters = 'ensembl_transcript_id', 
                  values = TPM_all$ensembl_transcript_id, 
                  mart = mart)
TPM_uniprot = merge.data.frame(x = TPM_pc, y = enst_unip, 
                               by = 'ensembl_transcript_id') %>%
  dplyr::select(-ensembl_transcript_id) 
TPM_uniprot = aggregate.data.frame(x = TPM_uniprot[, -13], by = list(TPM_uniprot[, 13]), FUN = sum, na.rm = T) %>%
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
# Proteomics normalization
 median_of_medians = protein_table %>% 
 	remove_rownames %>% 
 	column_to_rownames('uniprot_gn') %>% 
 	apply(2, median, na.rm = T) %>% 
 	median(na.rm = T)
 
 for (col in colnames(protein_table)[-ncol(protein_table)]) {
 	median_col = protein_table[, col] %>% median(na.rm = T)
 	correction = median_of_medians / median_col
 	protein_table[, col] = protein_table[, col] * correction
 }

#### Correlation ####
TPM_prot = merge.data.frame(x = TPM_uniprot, y = protein_table, 
                            by = 'uniprot_gn') %>% 
  column_to_rownames('uniprot_gn')

colnames(TPM_prot) = gsub('_quant.*', '', colnames(TPM_prot))		

estimates = pvalues = NULL
for (col.x in seq(1, 12, 3)) {
  for (row.x in 1:nrow(TPM_prot)) {
    x = TPM_prot[row.x, col.x:(col.x+2)] %>%
      as.numeric() %>%
      na.omit() 
    
    y = TPM_prot[row.x, 22:24] %>%
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

TPM_prot_2 = data.frame(TPM_002 = unlist(TPM_prot[, 1:3]))
TPM_prot_2$TPM_008 = unlist(TPM_prot[, 4:6])
TPM_prot_2$TPM_024 = unlist(TPM_prot[, 7:9])
TPM_prot_2$TPM_072 = unlist(TPM_prot[, 10:12])
TPM_prot_2$Proteomics_002 = unlist(TPM_prot[, 13:15])
TPM_prot_2$Proteomics_008 = unlist(TPM_prot[, 16:18])
TPM_prot_2$Proteomics_024 = unlist(TPM_prot[, 19:21])
TPM_prot_2$Proteomics_072 = unlist(TPM_prot[, 22:24])
TPM_prot_2$Proteomics_168 = unlist(TPM_prot[, 25:27])
TPM_prot_2$Proteomics_240 = unlist(TPM_prot[, 28:30])
TPM_prot_2$Proteomics_336 = unlist(TPM_prot[, 31:33])

M = TPM_prot_2 %>% zeroToNa %>% na.omit() %>% cor()
M_168 = TPM_prot_2 %>% .[nrow(TPM_prot)+1:nrow(TPM_prot_2), ] %>% zeroToNa %>% na.omit() %>% cor()
M[, 'Proteomics_168'] = M_168[, 'Proteomics_168']
M['Proteomics_168', ] = M_168['Proteomics_168', ]
png(
filename=
'/share/script/hecatos/juantxo/analysis_trc/correlation_transcr_prot/plots/correlation_plot_TPM_norm_corrected.png',
width = 1440, height = 1440)
corrplot(M, type="upper", order="original",
         col=brewer.pal(n=8, name="RdYlBu"), tl.cex = 3, cl.cex = 2.4)
dev.off()
#corrplot(cor(na.omit(TPM_prot)), type="upper", order="original",
#         col=brewer.pal(n=8, name="RdYlBu"))


#cor.test(c(TPM_prot[, 10], TPM_prot[, 11], TPM_prot[, 12]), c(TPM_prot[, 22], TPM_prot[, 23], TPM_prot[, 24]))
TPM_prot_2 %>% zeroToNa %>% na.omit() %>% nrow()

