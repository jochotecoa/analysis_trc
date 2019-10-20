#### Functions ####
source('/share/script/hecatos/juantxo/analysis_trc/functions.R')

#### Libraries ####
forceLibrary('biomaRt')

#### Parameters ####

level = 'ensembl_gene_id' # ensembl_transcript_id ensembl_gene_id
n_top_prot = 5
compound = 'UNTR'
comp_id = 'Con_UNTR'


#### Get the top protein ####
# Open protein file
setwd('/ngs-data/data/hecatos/Cardiac/')
setwd(comp_id)
setwd('Protein/')
proteomics_dir = list.dirs()[grep(pattern = 'Proteomics', list.dirs())]
setwd(proteomics_dir)
proteomics_file = list.files()[grep(pattern = 'renamed', list.files())]
protein_table = read.table(proteomics_file, header = T, sep = '\t')
# Clean data
protein_table = cleanProtIds(protein_table)
protein_table[is.na(protein_table)] = 0

# Get the range of protein expressions
protein_table_num = lapply(protein_table[, c(-1, -(ncol(protein_table)))], 
                           as.numeric)
protein_table_num = as.data.frame(protein_table_num)
min_max_diff = NULL
for (row in rownames(protein_table)) {
  vector = as.numeric(protein_table[row, c(-1, -(ncol(protein_table)))])
  min_max_vect = max(vector, T) - min(vector, T)
  min_max_diff = rbind(min_max_diff, min_max_vect)
}

protein_table$min_max_diff = min_max_diff

# Get the protein with the biggest/maximum range

top_prot = protein_table[order(protein_table$min_max_diff, 
                                 decreasing = T)[n_top_prot], ]

#### Get TPM, TRC & miRNA parameters
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
# If that protein is not related to any OMICs, comparison cannot be done
stopifnot(any(global.table$uniprot_gn == top_prot_naam))
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
mirna_effect = 'int_miRNA_sum'
top_all.sp = top.all[, grep(x = colnames(top.all), pattern = mirna_effect)]
trpl_1.sp = top_all.sp[, grep(x = colnames(top_all.sp), pattern = '_1_')]
trpl_1.sp = as.data.frame(t(trpl_1.sp))
ind = gregexpr(pattern = '002', text = rownames(trpl_1.sp)[1])[[1]][1]
trpl_1.sp$timepoints = as.numeric(substr(rownames(trpl_1.sp), ind, ind + 2))
trpl_1.sp[, 1] = naToZero(trpl_1.sp[, 1])

# circRNA effect
circRNA_effect = 'E_circ_sum'
top_all.circ = top.all[, grep(x = colnames(top.all), pattern = circRNA_effect)]
trpl_1.circ = top_all.circ[, grep(x = colnames(top_all.circ), pattern = '_1_')]
trpl_1.circ = as.data.frame(t(trpl_1.circ))
ind = gregexpr(pattern = '002', text = rownames(trpl_1.circ)[1])[[1]][1]
trpl_1.circ$timepoints = as.numeric(substr(rownames(trpl_1.circ), ind, ind + 2))
trpl_1.circ[, 1] = naToZero(trpl_1.circ[, 1])

#### Plotting ####
setwd("/share/script/hecatos/juantxo/analysis_trc/integrated_timeline/plots")
png(paste('protein', 'TPM', 'TRC', mirna_effect, circRNA_effect, 
          paste0('top_prot_', n_top_prot, '.png'), sep = '&'))
## add extra space to right margin of plot within frame
par(mar=c(5, 4, 4, 6) + 0.1)
## Plot first set of data and draw its axis
plot(x = trpl_1.prot$timepoints, y = trpl_1.prot[, 1], pch=16, axes=FALSE, xlab="", ylab="", 
     type="b",col="black", xlim = c(0, 336))
# axis(2, col="black",las=1)  ## las=1 makes horizontal labels
mtext("Expression",side=2,line=2.5)
mtext("Timepoint",side=1,line=2.5)
axis(1,pretty(range(trpl_1.prot$timepoints),10))
box()

## Allow a second plot on the same graph
par(new=TRUE)

## Plot the second plot and put axis scale on right
plot(x = trpl_1.tpm$timepoints, y = trpl_1.tpm[, 1], pch=15,  xlab="", ylab="", 
     axes=FALSE, type="b", col="red", xlim = c(0, 336))
axis(4, col="red",col.axis="red",las=1)

## Allow a second plot on the same graph
par(new=TRUE)

## Plot the second plot and put axis scale on right
plot(x = trpl_1.trc$timepoints, y = trpl_1.trc[, 1], pch=15,  xlab="", ylab="", 
     axes=FALSE, type="b", col="blue", xlim = c(0, 336))
# axis(4, col="blue",col.axis="red",las=1)

if(sum(trpl_1.sp[, 1]) != 0) {
  ## Allow a second plot on the same graph
  par(new=TRUE)
  
  ## Plot the second plot and put axis scale on right
  plot(x = trpl_1.sp$timepoints, y = trpl_1.sp[, 1], pch=15,  xlab="", ylab="", 
       axes=FALSE, type="b", col="green", xlim = c(0, 336))
}



if (sum(trpl_1.circ[, 1]) != 0) {
  ## Allow a second plot on the same graph
  par(new=TRUE)
  ## Plot the second plot and put axis scale on right
  plot(x = trpl_1.circ$timepoints, y = trpl_1.circ[, 1], pch=15,  xlab="", ylab="", 
       axes=FALSE, type="b", col="orange", xlim = c(0, 336))
  
}

# axis(4, col="blue",col.axis="red",las=1)
legend('bottomright', 
       legend = c('Protein', 'TPM', 'TRC', 'miRNA effect', 'circRNA'), 
       col = c('black', 'red', 'blue', 'yellow', 'orange'), 
       text.col = c('black', 'red', 'blue', 'green', 'orange'))

dev.off()
