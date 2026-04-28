# load required packages
source("../../utils.R")
forceLibrary(c("DESeq2", "ggplot2", "pheatmap", "grid", "gridExtra"))

options(scipen=999) # turn off scientific notation like 1e+06

# Use paths from config or relative paths
WORK.DIR <- getwd()
setwd(WORK.DIR)

# read the raw data table
dataFile <- "Andrea_raw_count.txt"
if (file.exists(dataFile)) {
    mRNAtot <- read.delim(dataFile, header=TRUE, row.names=1, sep="\t")
    colnames(mRNAtot) <- gsub("X","",colnames(mRNAtot))
} else {
    warning("Data file not found: ", dataFile)
    mRNAtot <- data.frame()
}

rnaseq_id_file = "all_genes.txt"
if (file.exists(rnaseq_id_file)) {
    rnaseq_id <- read.table(file = rnaseq_id_file, header = TRUE, sep = "\t")
} else {
    rnaseq_id <- data.frame()
}

## Read the metadata table
metadata_file = "samplekey_r_nr_beforehormonaltreatment.csv"
if (file.exists(metadata_file)) {
    samplekey <- read.csv(metadata_file, header=TRUE)
    rownames(samplekey) <- sub(pattern="X", replacement="", x=rownames(samplekey))
} else {
    samplekey <- data.frame()
}

if (nrow(mRNAtot) > 0 && nrow(samplekey) > 0) {
    mRNAtot <- mRNAtot[,colnames(mRNAtot) %in% rownames(samplekey)]
    
    ## Look for the number of reads per samples
    raw_read_count <- colSums(mRNAtot)
    
    # filter samples with less than 1000000 reads
    samplekey <- samplekey[rownames(samplekey) %in% colnames(mRNAtot),]
    
    ## DESeq2 Analysis
    dds <- DESeqDataSetFromMatrix(countData = round(mRNAtot),
                                  colData = samplekey,
                                  design = ~ status)
    
    dds <- DESeq(dds)
    
    # perform the normalization of the read count 
    norm_data <- counts(dds, normalized=TRUE)               
    
    if ("gene_name" %in% colnames(rnaseq_id)) {
        norm_data_name <- merge(norm_data, rnaseq_id["gene_name"], by="row.names", all.x=TRUE)
        rownames(norm_data_name) <- norm_data_name$Row.names
        norm_data_name <- norm_data_name[,-1]
    } else {
        norm_data_name <- as.data.frame(norm_data)
    }
    
    ## Set the pvalue to 0.05
    res <- results(dds, alpha=0.05, contrast=c("status", "Response" , "No_response"))
    res <- res[order(res$padj),]
    
    # perform a PCA plot
    rld <- vst(dds)
    plotPCA(rld, intgroup=c(colnames(samplekey)[min(5, ncol(samplekey))]))
    
    ## extract the differentially expressed genes passing FDR multiple testing
    DEmRNA_fdr <- subset(res, res$padj < 0.05)
    DEmRNA_fdr <- DEmRNA_fdr[order(DEmRNA_fdr$padj),]
    
    print(paste("there is", nrow(DEmRNA_fdr), "gene(s) passing FDR correction"))
    
    # save the results
    write.table(norm_data, file="batch2_norm_counts.txt", sep="\t", quote=FALSE)
    write.table(res, file="batch2_DE_Analysis.txt", sep="\t", quote=FALSE)
    write(rownames(DEmRNA_fdr), file="list_DEGs.txt")
    
    # Heatmap of top 100 DEGs
    if (nrow(DEmRNA_fdr) > 0) {
        top_genes = rownames(DEmRNA_fdr)[1:min(100, nrow(DEmRNA_fdr))]
        norm_data_heatmap <- subset(norm_data_name, rownames(norm_data_name) %in% top_genes)
        
        if ("gene_name" %in% colnames(norm_data_heatmap)) {
            # Fix syntax error and handle mapping
            plot_data = norm_data_heatmap[, setdiff(colnames(norm_data_heatmap), "gene_name")]
            pheatmap(log2(plot_data + 1), show_rownames=TRUE, 
                     labels_row=norm_data_heatmap$gene_name, 
                     cellwidth=30, cellheight=6.5, fontsize=7)
        } else {
            pheatmap(log2(norm_data_heatmap + 1), show_rownames=TRUE, 
                     cellwidth=30, cellheight=6.5, fontsize=7)
        }
    }
}


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
