read.salmon.quant <- function(sample) {
  quantfile.path = paste(sample, 'quant.sf', sep = '/')
  sample.file = read.table(file = quantfile.path, header = T)
}

get.expressed.rows <- function(dataset) {
  dataset.means = rowMeans(x = dataset)
  names(dataset.means) = rownames(dataset)
  datasetexpressed.means = dataset.means[dataset.means>0]
  datasetexpressed.rows = names(datasetexpressed.means)
  dataset.expressed = dataset[datasetexpressed.rows,]
  return(dataset.expressed)
}

setwd('/share/analysis/hecatos/juantxo/Score/RNA_Salmon/UNTR/')
listsamples.dirs = list.dirs(full.names = F, recursive = F)
rna = matrix()
for (sample in listsamples.dirs) {
  sample.file = read.salmon.quant(sample = sample)
  rna = cbind.data.frame(rna, sample.file$TPM)
  colnames(rna)[length(colnames(rna))] = sample
  rownames(rna) = sample.file[,1]
}
rna = rna[,-1]
rna.expressed = get.expressed.rows(dataset = rna)
circRNAexpressed.rows = grep(pattern = ':', x = rownames(rna.expressed))
circRNA.expressed = rna.expressed[circRNAexpressed.rows,]
noncircrna.expressed = rna.expressed[-circRNAexpressed.rows,]
noncircrna.expressed.means = rowMeans(noncircrna.expressed)
names(noncircrna.expressed.means) = rownames(noncircrna.expressed)
noncircrna.expressed.means.top = sort(noncircrna.expressed.means)[1:35000]
noncircrna.expressed.means.top.rows = names(noncircrna.expressed.means.top)
noncircrna.expressed.top = noncircrna.expressed[noncircrna.expressed.means.top.rows,]


all.expressed.list = c()
for (r in 1:nrow(noncircrna.expressed)) {
  value = all(noncircrna.expressed[r,]>0, na.rm = T)
  all.expressed.list = c(all.expressed.list, value)
}

setwd('/share/analysis/hecatos/juantxo/score_protein_analysis/plots/heatmaps/')
png(filename = 'heatmap_non-circRNAs_expressed_UNTR.png')
heatmap(as.matrix(noncircrna.expressed[all.expressed.list,]))
dev.off()
