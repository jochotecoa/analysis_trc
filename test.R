library('biomaRt')

#### Parameters ####

level = 'ensembl_gene_id' # ensembl_transcript_id
n_top = 5
compound = 'Dox'

#### Functions ####
cleanProtIds = function(protein_table) {
  protein_table = protein_table[!grepl(protein_table[, 1], pattern = ':'), ]
  names = strsplit(as.character(protein_table[, 1]), '\\|')
  names = as.character(lapply(names, '[', 2))
  protein_table$uniprot_gn = names
  return(protein_table)
}

transcrToGene = function(table, aggregate = F) {
  table[, 'rownames'] = rownames(table)
  sampl = table[nrow(table), ]
  enst_col = grep(pattern = 'ENST', x = sampl)[1]
  
  version = grepl('\\.', sampl[, enst_col])
  if (length(version) == 0) {version = F}
  if (version) {
    transcript_id = 'ensembl_transcript_id_version'
  } else {
    transcript_id = 'ensembl_transcript_id'
  }
  
  values = table[, enst_col]
  mart.human = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                       dataset = 'hsapiens_gene_ensembl',
                       host = 'http://apr2018.archive.ensembl.org') 
  
  new_cols = getBM(attributes = c(transcript_id, 'ensembl_gene_id'), 
                   filters = transcript_id, values = values, mart = mart.human)
  
  table = merge.data.frame(x = table, y = new_cols, 
                           by.x = colnames(table)[enst_col], 
                           by.y = transcript_id)
  if (aggregate) {
    int_cols = grepl('integer', sapply(X = table[1, ], FUN = typeof))
    int_cols = int_cols + grepl('double', sapply(X = table[1, ], 
                                                 FUN = typeof))
    int_cols = as.logical(int_cols)
    table = aggregate(x = table[, int_cols], by = list(table$ensembl_gene_id), 
                      FUN = sum)
    colnames(table)[1] = 'ensembl_gene_id'
  }
  return(table)
}

#### Analysis ####

setwd('/share/analysis/hecatos/juantxo/Score/output/Output_Run_mrna_SEPT2019/')
setwd('V3/output/TRCscore/')
if (length(ls(pattern = 'untr_.\\.table')) != 2) {
  untr_2.table = read.table('UNTR_002_1_TRCscore.txt', stringsAsFactors = F)
  untr_8.table = read.table('UNTR_008_1_TRCscore.txt', stringsAsFactors = F)
  untr_2.table[, 'ensembl_transcript_id'] = rownames(untr_2.table)
  untr_8.table[, 'ensembl_transcript_id'] = rownames(untr_8.table)
  untr_2.agg = untr_2.table
  untr_8.agg = untr_8.table
}
if (level == 'ensembl_gene_id') {
  untr_2.agg = transcrToGene(untr_2.table, aggregate = T)
  untr_8.agg = transcrToGene(untr_8.table, aggregate = T)
}


utr_2.subset = untr_2.agg[, level] %in% miranda.ids
utr_8.subset = untr_8.agg[, level] %in% miranda.ids
utr_2.subset = untr_2.agg[utr_2.subset, ]
utr_8.subset = untr_8.agg[utr_8.subset, ]

comm_transcr = utr_2.subset[, level] %in% utr_8.subset[, level]
comm_transcr = utr_2.subset[, level][comm_transcr]
utr_2.subset = utr_2.subset[utr_2.subset[, level] %in% comm_transcr, ]
utr_8.subset = utr_8.subset[utr_8.subset[, level] %in% comm_transcr, ]

# Incorporate protein column

mart.human = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                     dataset = 'hsapiens_gene_ensembl',
                     host = 'http://apr2018.archive.ensembl.org') 
gene_prot = getBM(attributes = c('ensembl_gene_id', 'uniprot_gn'), 
                  filters = 'ensembl_gene_id', 
                  values = utr_2.subset$ensembl_gene_id, 
                  mart = mart.human)

utr_2.subset = merge.data.frame(x = utr_2.subset, 
                                y = gene_prot, 
                                by = 'ensembl_gene_id')
utr_8.subset = merge.data.frame(x = utr_8.subset, 
                                y = gene_prot, 
                                by = 'ensembl_gene_id')

comp_id = 'Con_UNTR'
setwd('/ngs-data/data/hecatos/Cardiac/')
setwd(comp_id)
setwd('Protein/')
proteomics_dir = list.dirs()[grep(pattern = 'Proteomics', list.dirs())]
setwd(proteomics_dir)
proteomics_file = list.files()[grep(pattern = 'renamed', list.files())]
protein_table = read.table(proteomics_file, header = T, sep = '\t')
protein_table = cleanProtIds(protein_table)

utr_2.subset = merge(x = utr_2.subset, y = protein_table, by = 'uniprot_gn')
utr_8.subset = merge(x = utr_8.subset, y = protein_table, by = 'uniprot_gn')

diffs.tpm = abs(utr_2.subset$targetRNA_TPM - utr_8.subset$targetRNA_TPM)
utr_8.subset$diff_tpm_2_8 = diffs.tpm
utr_2.subset$diff_tpm_2_8 = diffs.tpm

min_tpm = diff.min.tpm = min(utr_8.subset$diff_tpm_2_8)
min_tpm_2 = utr_2.subset[grep(min_tpm, utr_2.subset$diff_tpm_2_8), ]
min_tpm_8 = utr_8.subset[grep(min_tpm, utr_8.subset$diff_tpm_2_8), ]
min_tpm_2.row = rownames(min_tpm_2)
min_tpm_8.row = rownames(min_tpm_8)

# #### Let's check the difference in TRC ####
diff.min.trc = abs(utr_2.subset[min_tpm_2.row, 
                                'TRC'] - utr_8.subset[min_tpm_8.row, 'TRC'])
diff.min.prot = abs(sum(utr_2.subset$DOX_The_002_1[min_tpm_2.row], 
                        -(utr_2.subset$DOX_The_008_1[min_tpm_2.row]), 
                        na.rm = T))



