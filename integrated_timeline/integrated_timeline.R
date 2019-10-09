library('biomaRt')

#### Parameters ####

level = 'ensembl_gene_id' # ensembl_transcript_id ensembl_gene_id
n_top = 5
compound = 'UNTR'
comp_id = 'Con_UNTR'

#### Functions ####
source('/share/script/hecatos/juantxo/analysis_trc/functions.R')

#### Get the top protein ####
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

#### Get all the 
setwd("/share/analysis/hecatos/juantxo/Score/output/Output_Run_mrna_SEPT2019/")
setwd('V3/output/UNTR/TRCscore')

trc.files = list.files(pattern = 'UNTR')

first_trc.table = read.table(trc.files[1], stringsAsFactors = F)
first_trc.table = rmMirnas(first_trc.table)
colnames(first_trc.table) = paste(colnames(first_trc.table), trc.files[1], 
                                  sep = '.')
first_trc.table[, 'ensembl_transcript_id'] = rownames(first_trc.table)
global_trc.table = first_trc.table

for (f in trc.files[-1]) {
  ind_trc.table = read.table(f, stringsAsFactors = F)
  ind_trc.table = rmMirnas(ind_trc.table)
  colnames(ind_trc.table) = paste(colnames(ind_trc.table), f, sep = '.')
  ind_trc.table[, 'ensembl_transcript_id'] = rownames(ind_trc.table)
  
  global_trc.table = merge.data.frame(x = global_trc.table, y = ind_trc.table, 
                                      by = 'ensembl_transcript_id')
}

global_gene.table = transcrToGene(global_trc.table, T)

mart.human = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                     dataset = 'hsapiens_gene_ensembl',
                     host = 'http://apr2018.archive.ensembl.org') 
gene_prot = getBM(attributes = c('ensembl_gene_id', 'uniprot_gn'), 
                  filters = 'uniprot_gn', values = protein_table$uniprot_gn, 
                  mart = mart.human)

protein_table = merge(protein_table, gene_prot, by = 'uniprot_gn')
prot_expr.cols = grep('UNTR', colnames(protein_table))
new.cols = paste0(colnames(protein_table)[prot_expr.cols], '.protein')
colnames(protein_table)[prot_expr.cols] = new.cols
global.table = merge(global_gene.table, protein_table, by = 'ensembl_gene_id')
# Protein expression
top_prot_naam = top_prot$uniprot_gn
top.all = global.table[global.table$uniprot_gn == top_prot_naam, ]
top_all.prot = top.all[, grep(x = colnames(top.all), pattern = '.protein')]
trpl_1.prot = top_all.prot[, grep(x = colnames(top_all.prot), pattern = '_1\\.')]
trpl_1.prot = as.data.frame(t(trpl_1.prot))
trpl_1.prot$timepoints = as.numeric(substr(rownames(trpl_1.prot), 10, 12))
# TPM expression
top_all.tpm = top.all[, grep(x = colnames(top.all), pattern = 'targetRNA_TPM')]
trpl_1.tpm = top_all.tpm[, grep(x = colnames(top_all.tpm), pattern = '_1_')]
trpl_1.tpm = as.data.frame(t(trpl_1.tpm))
ind = gregexpr(pattern = '002', text = rownames(trpl_1.tpm)[1])[[1]][1]
trpl_1.tpm$timepoints = as.numeric(substr(rownames(trpl_1.tpm), ind, ind + 2))
# TRC expression
top_all.trc = top.all[, grep(x = colnames(top.all), pattern = 'TRC\\.')]
trpl_1.trc = top_all.trc[, grep(x = colnames(top_all.trc), pattern = '_1_')]
trpl_1.trc = as.data.frame(t(trpl_1.trc))
ind = gregexpr(pattern = '002', text = rownames(trpl_1.trc)[1])[[1]][1]
trpl_1.trc$timepoints = as.numeric(substr(rownames(trpl_1.trc), ind, ind + 2))
# miRNA effect
top_all.sp = top.all[, grep(x = colnames(top.all), pattern = 'miRNA_sum')]
trpl_1.sp = top_all.sp[, grep(x = colnames(top_all.sp), pattern = '_1_')]
trpl_1.sp = as.data.frame(t(trpl_1.sp))
ind = gregexpr(pattern = '002', text = rownames(trpl_1.sp)[1])[[1]][1]
trpl_1.sp$timepoints = as.numeric(substr(rownames(trpl_1.sp), ind, ind + 2))

# circRNA effect
top_all.circ = top.all[, grep(x = colnames(top.all), pattern = 'E_circ_sum')]
trpl_1.circ = top_all.circ[, grep(x = colnames(top_all.circ), pattern = '_1_')]
trpl_1.circ = as.data.frame(t(trpl_1.circ))
ind = gregexpr(pattern = '002', text = rownames(trpl_1.circ)[1])[[1]][1]
trpl_1.circ$timepoints = as.numeric(substr(rownames(trpl_1.circ), ind, ind + 2))

#### Plotting ####
## add extra space to right margin of plot within frame
par(mar=c(5, 4, 4, 6) + 0.1)
## Plot first set of data and draw its axis
plot(x = trpl_1.prot$timepoints, y = trpl_1.prot$`160`, pch=16, axes=FALSE, xlab="", ylab="", 
     type="b",col="black", xlim = c(0, 336))
# axis(2, col="black",las=1)  ## las=1 makes horizontal labels
mtext("Expression",side=2,line=2.5)
mtext("Timepoint",side=1,line=2.5)
axis(1,pretty(range(trpl_1.prot$timepoints),10))
box()

## Allow a second plot on the same graph
par(new=TRUE)

## Plot the second plot and put axis scale on right
plot(x = trpl_1.tpm$timepoints, y = trpl_1.tpm$`160`, pch=15,  xlab="", ylab="", 
     axes=FALSE, type="b", col="red", xlim = c(0, 336))
axis(4, col="red",col.axis="red",las=1)

## Allow a second plot on the same graph
par(new=TRUE)

## Plot the second plot and put axis scale on right
plot(x = trpl_1.trc$timepoints, y = trpl_1.trc$`160`, pch=15,  xlab="", ylab="", 
     axes=FALSE, type="b", col="blue", xlim = c(0, 336))
# axis(4, col="blue",col.axis="red",las=1)

## Allow a second plot on the same graph
par(new=TRUE)

## Plot the second plot and put axis scale on right
plot(x = trpl_1.sp$timepoints, y = trpl_1.sp$`160`, pch=15,  xlab="", ylab="", 
     axes=FALSE, type="b", col="green", xlim = c(0, 336))
## Allow a second plot on the same graph
par(new=TRUE)

## Plot the second plot and put axis scale on right
plot(x = trpl_1.circ$timepoints, y = trpl_1.circ$`160`, pch=15,  xlab="", ylab="", 
     axes=FALSE, type="b", col="orange", xlim = c(0, 336))

# axis(4, col="blue",col.axis="red",las=1)
legend('bottomright', 
       legend = c('Protein', 'TPM', 'TRC', 'miRNA effect', 'circRNA'), 
       col = c('black', 'red', 'blue', 'yellow', 'orange'), 
       text.col = c('black', 'red', 'blue', 'green', 'orange'))


#### Get the most drastic miRNA change ####

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
