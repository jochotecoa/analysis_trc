# load required packages
require("DESeq2")
require("ggplot2")
require("pheatmap")
require("grid")
library("gridExtra")

options(scipen=999) # turn off scientific notation like 1e+06

# Setup your working directory
predir <- 'C:/D/florian/andrea/'
WORK.DIR <- paste(predir,"Batch1/", sep="")
setwd(WORK.DIR)

# read the raw data table
dataFile <- "Andrea_raw_count.txt"
mRNAtot <- read.delim(dataFile, header=TRUE, row.names=1, sep="\t")
colnames(mRNAtot) <- gsub("X","",colnames(mRNAtot))

rnaseq_id <- read.table(file = "all_genes.txt", header = TRUE, sep = "\t")

## Read the metadata table
samplekey <- read.csv("samplekey_r_nr_beforehormonaltreatment.csv", header=TRUE)
rownames(samplekey) <- sub(pattern="X", replacement="", x=rownames(samplekey))

mRNAtot <- mRNAtot[,colnames(mRNAtot) %in% rownames(samplekey)]

## Look for the number of reads per samples: 
raw_read_count <- colSums(mRNAtot)

barplot(raw_read_count, las=2, cex.names=0.6)
barplot(log2(raw_read_count), las=2, cex.names=0.6)

# filter samples with less than 1000000 reads

samplekey <- samplekey[rownames(samplekey) %in% colnames(mRNAtot),]



## Start the statistical analysis made on the Compounds variable from the metadata file. You could change this by modifying the Compound by another field from Metadata
dds <- DESeqDataSetFromMatrix(countData = round(mRNAtot),
                              colData = samplekey,
                              design = ~ status)

# If you want DeSEQ NOT to normalize you data, do : 
# sizeFactors(dds) <- rep(1, ncol(dds))


dds<-DESeq(dds)

# perform the normalization of the read count 
norm_data <- counts(dds,normalized=TRUE)               

norm_data_name <- merge(norm_data,rnaseq_id["gene_name"],by="row.names",all.x=TRUE)
rownames(norm_data_name) <- norm_data_name$Row.names
norm_data_name <- norm_data_name[,-1]

norm_read_count <- colSums(norm_data)
barplot(norm_read_count, las=2, cex.names=0.6)
barplot(log2(norm_read_count), las=2, cex.names=0.6)

## Set the pvalue to 0.05
res <- results(dds, alpha=0.05, contrast=c("status", "Response" , "No_response"))
res <- res[order(res$padj),]

# perform a PCA plot. You can change the condition assessed by the pca plot by changing the value between [] from 1 to 3 (corresponding to the column number of the metadata file)
rld <- vst(dds)
plotPCA(rld,intgroup=c(colnames(samplekey)[5]))
savePlot(filename="localisation.png",type="png")

## extract the differentially expressed genes passing FDR multiple testing
DEmRNA_fdr <-  subset(res,res$padj < 0.05)
DEmRNA_fdr <- DEmRNA_fdr[order(DEmRNA_fdr$padj),]

print(paste("there is",nrow(DEmRNA_fdr), "gene(s) passing FDR correction"))

DEmRNA_all <-  subset(res,res$pvalue < 0.05)

print(paste("there is",nrow(DEmRNA_all), "gene(s) have a pvalue < 0.05"))

# save the normalized count table: 
write.table(norm_data,file="batch2_norm_counts.txt", sep="\t", quote=FALSE)

# save the pvalue and fdr corrected value for each gene: 
write.table(res,file="batch2_DE_Analysis.txt", sep="\t", quote=FALSE)

# save the pvalue and fdr corrected value for each gene: 
write(rownames(DEmRNA_fdr),file="list_DEGs.txt")

# perform a heat Map:
#norm_data_heatmap <- subset(norm_data_name, rownames(norm_data_name) %in% rownames(DEmRNA_fdr)[1:100])
#pheatmap(log2(norm_data_heatmap+0.1), show_rownames=FALSE)            

norm_data_heatmap <- subset(norm_data_name, rownames(norm_data_name) %in% rownames(DEmRNA_fdr)[1:100])
rownames(norm_data_heatmap) <- norm_data_heatmap$gene_name)
pheatmap(log2(norm_data_heatmap[,1:ncol(norm_data_heatmap)-1]+1), show_rownames=TRUE,labels_row=norm_data_heatmap$gene_name, cellwidth=30 , cellheight= 6.5 , fontsize= 7)


# Look for the normalized read count of a gene. You can freely change the ensembl ID to plot your gene of interest
gene_to_plot <- "ENSG00000091831"
gene_count <- subset(norm_data, rownames(norm_data) == gene_to_plot)
barplot(gene_count, las=2, cex.names = 0.6, main="RUNX1T1")
savePlot(filename="DMD.png",type="png")

corm = cor(norm_data)

for (i in 1:(ncol(norm_data)-1)) {
  for (j in (i+1):ncol(norm_data)) {
    
    DF <- data.frame(norm_data[,i], norm_data[,j])
    colnames(DF) <- c(colnames(norm_data)[i],colnames(norm_data)[j])
    DF_filt <- subset(DF,rowSums(DF) > 0)
    DF_filt$diff <- DF_filt[,2]-DF_filt[,1]
    
    
    my_grob = grid.text(paste("R2 =",round(corm[i,j],digit=2)), x=0.1,  y=0.95, gp=gpar(col="red3", fontsize=12, fontface="bold"))
    g1 <- ggplot(DF_filt, aes(x=log2(DF_filt[,1]), y=log2(DF_filt[,2]))) + geom_point() + geom_smooth(method="lm",col="red3") + annotation_custom(my_grob) + labs(y=colnames(DF_filt[2]), x=colnames(DF_filt[1]))
    g2 <- ggplot(DF_filt, aes(x=log2(DF_filt[,1]), y=log2(DF_filt[,2]))) + geom_point() + geom_smooth(method="lm",col="red3") + annotation_custom(my_grob) + labs(y=colnames(DF_filt[2]), x=colnames(DF_filt[1]))
    
    plot(g)
    
    
  }}
