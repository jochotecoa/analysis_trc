library('biomaRt')

#### Parameters ####

level = 'ensembl_gene_id' # ensembl_transcript_id ensembl_gene_id
n_top = 5
compound = 'UNTR'
comp_id = 'Con_UNTR'

#### Functions ####
source('/share/script/hecatos/juantxo/analysis_trc/functions.R')

#### Analysis ####
setwd('/ngs-data/data/hecatos/Cardiac/')
setwd(comp_id)
setwd('Protein/')
proteomics_dir = list.dirs()[grep(pattern = 'Proteomics', list.dirs())]
setwd(proteomics_dir)
proteomics_file = list.files()[grep(pattern = 'renamed', list.files())]
protein_table = read.table(proteomics_file, header = T, sep = '\t')
protein_table = cleanProtIds(protein_table)
protein_table[is.na(protein_table)] = 0

protein_table_num = as.numeric(protein_table[, c(-1, -(ncol(protein_table)))])
min_max_diff = NULL
for (row in rownames(protein_table)) {
  vector = as.numeric(protein_table[row, c(-1, -(ncol(protein_table)))])
  min_max_vect = max(vector, T) - min(vector, T)
  min_max_diff = rbind(min_max_diff, min_max_vect)
}

protein_table$min_max_diff = min_max_diff
top_diff = max(protein_table$min_max_diff, T) 

top_prot = protein_table[grep(top_diff, protein_table$min_max_diff), ]


# #### Get the most drastic miRNA change ####
# 
# setwd('/share/analysis/hecatos/juantxo/Score/input/miRNA_miRge2/')
# setwd(compound)
# 
# mirna_expr.table = read.csv(file = 'miR.RPM.csv', stringsAsFactors = F)
# mirna_expr.table$diff_2_8 = abs(mirna_expr.table[,2] - mirna_expr.table[,5])
# top_diff = sort(mirna_expr.table$diff_2_8, decreasing = T)[n_top]
# top_mirna = mirna_expr.table[grep(top_diff, mirna_expr.table$diff_2_8), ]
# # top_diff2 = sort(mirna_expr.table$diff_2_8, decreasing = T)[3]
# # top_mirna2 = mirna_expr.table[grep(top_diff2, mirna_expr.table$diff_2_8), ]
# 
# #### Get which transcripts are targeted by this miRNA ####
# 
# setwd('../../Interactions/')
# if (length(ls(pattern = 'miranda.table')) == 0) {
#   miranda.table = read.table('ALLmRNAtranscripts_vs_ALLmiRNA3.parsed', 
#                              stringsAsFactors = F)
#   colnames(miranda.table)[1:2] = c('miRNA_id', 'ensembl_transcript_id')
# }
# miranda.subset = miranda.table[grep(top_mirna$miRNA, miranda.table$miRNA_id), ]
# miranda.subset = transcrToGene(miranda.subset)
# miranda.ids = unique(miranda.subset$ensembl_gene_id)


# miranda.subset2 = miranda.table[grep(top_mirna2$miRNA, 
#                                      miranda.table$miRNA_id), ]
#### How do these transcripts behave? ####
setwd('/share/analysis/hecatos/juantxo/Score/output/Output_Run_mrna_SEPT2019/')
setwd('V3/output/TRCscore/')
# setwd(compound)
if (length(ls(pattern = 'untr_.\\.table')) == 2) {
  untr_2.table = read.table('UNTR_002_1_TRCscore.txt', stringsAsFactors = F)
  untr_8.table = read.table('UNTR_008_1_TRCscore.txt', stringsAsFactors = F)
  untr_2.table[, 'ensembl_transcript_id'] = rownames(untr_2.table)
  untr_8.table[, 'ensembl_transcript_id'] = rownames(untr_8.table)
  untr_2.agg = untr_2.table[, !grepl('hsa', colnames(untr_2.table))]
  untr_8.agg = untr_8.table[, !grepl('hsa', colnames(untr_8.table))]
}
if (level == 'ensembl_gene_id') {
  untr_2.agg = transcrToGene(untr_2.agg, aggregate = T)
  untr_8.agg = transcrToGene(untr_8.agg, aggregate = T)
}

#### Subset by proteomics ids ####



utr_2.subset = merge.data.frame(x = untr_2.agg, 
                                y = gene_prot, 
                                by = level)
utr_8.subset = merge.data.frame(x = untr_8.agg, 
                                y = gene_prot, 
                                by = level)




utr_2.subset = merge(x = utr_2.subset, y = protein_table, by = 'uniprot_gn')
utr_8.subset = merge(x = utr_8.subset, y = protein_table, by = 'uniprot_gn')

utr_2_8.subset = merge(utr_2.subset, utr_8.subset, by = level)
colnames(utr_2_8.subset) = gsub('.x', '.2', colnames(utr_2_8.subset))
colnames(utr_2_8.subset) = gsub('.y', '.8', colnames(utr_2_8.subset))

utr_2_8.subset$diff_prot = utr_2_8.subset$UNTR_The_002_1.2 - utr_2_8.subset$UNTR_The_008_1.2
utr_2_8.subset$diff_prot = abs(utr_2_8.subset$diff_prot)
utr_2_8.subset = utr_2_8.subset[utr_2_8.subset$diff_prot > 0, ]
diff_prot = as.data.frame(sort(unique(utr_2_8.subset$diff_prot), T))
diff_prot$order_prot = as.numeric(rownames(diff_prot))
colnames(diff_prot)[1] = 'diff_prot'
utr_2_8.subset = merge(x = utr_2_8.subset, y = diff_prot, by = 'diff_prot')

utr_2_8.subset$diff_tpm = utr_2_8.subset$targetRNA_TPM.2 - utr_2_8.subset$targetRNA_TPM.8
utr_2_8.subset$diff_tpm = abs(utr_2_8.subset$diff_tpm)
utr_2_8.subset = utr_2_8.subset[utr_2_8.subset$diff_tpm > 0, ]
diff_tpm = as.data.frame(sort(unique(utr_2_8.subset$diff_tpm), T))
diff_tpm$order_tpm = as.numeric(rownames(diff_tpm))
colnames(diff_tpm)[1] = 'diff_tpm'
utr_2_8.subset = merge(x = utr_2_8.subset, y = diff_tpm, by = 'diff_tpm')

utr_2_8.subset$diff_prot_tpm = abs(utr_2_8.subset$order_prot - utr_2_8.subset$order_tpm)

maxmin = utr_2_8.subset[grep(max(utr_2_8.subset$diff_prot_tpm), utr_2_8.subset$diff_prot_tpm), ]

#### Subset by miRanda targets ####

utr_2.subset = utr_2.subset[, level] %in% miranda.ids
utr_8.subset = utr_8.subset[, level] %in% miranda.ids
utr_2.subset = utr_2.subset[utr_2.subset, ]
utr_8.subset = utr_8.subset[utr_8.subset, ]

#### Subset by common in both timepoints ####

comm_transcr = utr_2.subset[, level] %in% utr_8.subset[, level]
comm_transcr = utr_2.subset[, level][comm_transcr]
utr_2.subset = utr_2.subset[utr_2.subset[, level] %in% comm_transcr, ]
utr_8.subset = utr_8.subset[utr_8.subset[, level] %in% comm_transcr, ]

#### Get the least affected gene in TPM terms ####

diffs.tpm = abs(utr_2.subset$targetRNA_TPM - utr_8.subset$targetRNA_TPM)
utr_8.subset$diff_tpm_2_8 = diffs.tpm
utr_2.subset$diff_tpm_2_8 = diffs.tpm

min_tpm = diff.min.tpm = min(utr_8.subset$diff_tpm_2_8)
min_tpm_2 = utr_2.subset[grep(min_tpm, utr_2.subset$diff_tpm_2_8), ]
min_tpm_8 = utr_8.subset[grep(min_tpm, utr_8.subset$diff_tpm_2_8), ]
min_tpm_2.row = rownames(min_tpm_2)
min_tpm_8.row = rownames(min_tpm_8)

#### Let's check the difference in TRC ####
diff.min.trc = abs(utr_2.subset[min_tpm_2.row, 
                                'TRC'] - utr_8.subset[min_tpm_8.row, 'TRC'])
diff.min.prot = abs(sum(utr_2.subset$UNTR_The_002_1[min_tpm_2.row], 
                        -(utr_2.subset$UNTR_The_008_1[min_tpm_2.row]), 
                        na.rm = T))

print(c(top_diff, diff.min.tpm, diff.min.trc, diff.min.prot))
