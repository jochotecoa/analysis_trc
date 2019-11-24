# Libraries ---------------------------------------------------------------

source('/share/script/hecatos/juantxo/analysis_trc/functions.R')
forceLibrary(c('biomaRt', "tximport", "dplyr", "DESeq2", "grid", "ggplot2", 
               "pheatmap", "BiocParallel"))
register(MulticoreParam(20))
getCts <- function(project_name) {
  if (grepl('TRC', project_name)) {
    setwd('/share/analysis/hecatos/juantxo/Score/output/UNTR')
    transcript_id = 'ensembl_transcript_id'
    col_quant = 'TRC_'
    if (grepl('voom', project_name)) {
      setwd('2019-11-11_09:41:21_UTC/')
    } else {
      setwd('2019-11-11_12:04:38_UTC/')
    }
    quant_file = mergeFiles(files_patt = 'TRCscore', row_names = T, all = T)
  }
  
  if (grepl('Salmon', project_name)) {
    transcript_id = 'ensembl_transcript_id_version'
    setwd('/share/analysis/hecatos/juantxo/Score/input/RNA_Salmon/UNTR/')
    if (grepl('counts', project_name)) {
      col_quant = 'NumReads'
    } else {
      col_quant = 'TPM'
    }
    quant_file = mergeFiles(files_patt = 'quant.sf', 
                            by_col = 'Name', all = T)
  }
  
  if (grepl('protein', project_name)) {
    source('/share/script/hecatos/juantxo/analysis_trc/functions.R')
    setwd('/ngs-data/data/hecatos/Cardiac/Con_UNTR/Protein/Proteomics_Analyses_Cardiac_UNTR_GeneData/')
    prot = read.table('Hecatos_Cardiac_Px_Untreated_pre-processed_renamed.txt', header = T, sep = '\t', row.names = 1)
    prot_num = apply(prot, MARGIN = 2, as.numeric) %>% .[, 1:12]
    rownames(prot_num) = rownames(prot)
    quant_file_counts = na.omit(prot_num)
    cts = apply(quant_file_counts, 2, round) %>% as.data.frame()
    
  } else {
    quant_file = quant_file %>% 
      filter(grepl(quant_file[, 1], pattern = 'ENST'))
    mart = openMart2018()
    tx2gene = getBM(attributes = c(transcript_id, 'ensembl_gene_id'), 
                    filters = transcript_id, values = quant_file[, 1],
                    mart = mart)
    nrow(tx2gene) %>% all(. > 0) %>% stopifnot()
    quant_file_genes = merge.data.frame(quant_file, tx2gene, 
                                        by.x = colnames(quant_file)[1], 
                                        by.y = transcript_id)
    quant_file_counts = quant_file_genes  %>% 
      dplyr::select(contains(col_quant))
    cts = apply(quant_file_counts, 2, round) %>% as.data.frame()
    # cts$ensembl_gene_id = quant_file_genes$ensembl_gene_id
    cts = aggregate.data.frame(x = cts, 
                               by = list(quant_file_genes$ensembl_gene_id), 
                               FUN = sum, na.rm = T) %>%
      tibble::column_to_rownames(var = 'Group.1')
  }
}


# Input data --------------------------------------------------------------

# Parameters

project_name = 'Salmon_counts' # Salmon_counts TRC_no_voom
setSizeFactorToOne = F # Default: FALSE
filtering = T # Default: TRUE

# Analysis 
cts = getCts(project_name)

# Metadata ----------------------------------------------------------------


timepoints = NULL
for (tp in c('UNTR_002', 'UNTR_008','UNTR_024', 'UNTR_072')) {
    repls = rep(tp, 3)
    timepoints = c(timepoints, repls)
}

coldata = data.frame(row.names = colnames(cts), timepoints = timepoints)

rownames(coldata) %in% colnames(cts) %>% all() %>% stopifnot()
all(rownames(coldata) == colnames(cts)) %>% stopifnot()


# DESeq2 pipeline ---------------------------------------------------------

# txi = tximport(files = list.files(pattern = 'quant.sf', recursive = T)[-1], type = "salmon", tx2gene = tx2gene)
dds = DESeqDataSetFromMatrix(countData = cts, 
                             colData = coldata, 
                             design = ~ timepoints)

if (setSizeFactorToOne) {
  sizeFactors(dds) <- rep(1, ncol(dds)) # this is the important line
}


dds <- DESeq(dds)

resultsNames(dds) # lists the coefficients

# Independent filtering can be turned off by setting independentFiltering to FALSE.
# cooksCutoff: count outliers


# perform the normalization of the read count 
norm_data <- counts(dds,normalized=TRUE)               

# norm_data_name <- merge(norm_data,rnaseq_id["gene_name"],by="row.names",all.x=TRUE)
# rownames(norm_data_name) <- norm_data_name$Row.names
# norm_data_name <- norm_data_name[,-1]

norm_read_count <- colSums(norm_data)
# barplot(norm_read_count, las=2, cex.names=0.6)
# barplot(log2(norm_read_count), las=2, cex.names=0.6)

## Set the pvalue to 0.05
# res <- results(dds, alpha=0.05, contrast=c("status", "Response" , "No_response"))
# res <- res[order(res$padj),]

# perform a PCA plot. You can change the condition assessed by the pca plot by changing the value between [] from 1 to 3 (corresponding to the column number of the metadata file)
rld <- vst(dds)
pca = plotPCA(rld, intgroup = names(colData(dds)))
setwd("/share/analysis/hecatos/juantxo/Score/analysis/DESeq2/")
forceSetWd(project_name)
ggsave(filename="localisation.png", device = "png", plot = pca)


DEmRNA_fdr = data.frame(ensembl_gene_id = character())
for (name in resultsNames(dds)[-1]) {
  if (filtering) {
    res <- results(dds, name = name)
  } else {
    res <- results(dds, 
                   independentFiltering=F, cooksCutoff = F, 
                   name = name) 
    }
  res <- results(dds, name = name)
  DEmRNA_fdr_tmp = data.frame(ensembl_gene_id = rownames(res), res$padj < 0.05)
  colnames(DEmRNA_fdr_tmp)[2] = name
  print(paste("There are", sum(DEmRNA_fdr_tmp[, 2], na.rm = T), 
              "gene(s) that have a pvalue < 0.05"))
  DEmRNA_fdr = merge.data.frame(DEmRNA_fdr, DEmRNA_fdr_tmp, 
                                by = 'ensembl_gene_id', all = T)
}
  
# save the normalized count table: 
write.table(norm_data,file="batch2_norm_counts.txt", sep="\t", quote=FALSE)

# save the pvalue and fdr corrected value for each gene: 
write.table(res,file="batch2_DE_Analysis.txt", sep="\t", quote=FALSE)

# save the pvalue and fdr corrected value for each gene: 
write.table(DEmRNA_fdr, file="list_DEGs.txt", )

# perform a heat Map:
#norm_data_heatmap <- subset(norm_data_name, rownames(norm_data_name) %in% rownames(DEmRNA_fdr)[1:100])
#pheatmap(log2(norm_data_heatmap+0.1), show_rownames=FALSE)            

# norm_data_heatmap <- subset(norm_data, rownames(norm_data) %in% rownames(DEmRNA_fdr)[1:100]) %>% as.data.frame()
# rownames(norm_data_heatmap) <- subset(rownames(norm_data), rownames(norm_data) %in% rownames(DEmRNA_fdr)[1:100])
# phm = pheatmap(log2(norm_data_heatmap[,1:ncol(norm_data_heatmap)-1]+1), show_rownames=TRUE,labels_row=norm_data_heatmap$gene_name, cellwidth=30 , cellheight= 6.5 , fontsize= 7)
# ggsave('heatmap_DEGs.png', plot = phm)

# Look for the normalized read count of a gene. You can freely change the ensembl ID to plot your gene of interest
# gene_to_plot <- "ENSG00000091831"
# gene_count <- subset(norm_data, rownames(norm_data) == gene_to_plot)
# barplot(gene_count, las=2, cex.names = 0.6, main="RUNX1T1")
# savePlot(filename="DMD.png",type="png")

corm = cor(norm_data)

forceSetWd('plots')
forceSetWd('log2_correlation')
is.log2 = c(T,F)

for (l in is.log2) {
  for (i in 1:(ncol(norm_data)-1)) {
    for (j in (i+1):ncol(norm_data)) {
      
      DF <- data.frame(norm_data[,i], norm_data[,j])
      colnames(DF) <- c(colnames(norm_data)[i],colnames(norm_data)[j])
      DF_filt <- subset(DF,rowSums(DF) > 0)
      DF_filt$diff <- DF_filt[,2]-DF_filt[,1]
      
      
      my_grob = grid.text(paste("R2 =", round(cor(DF_filt[,1], 
                                                  DF_filt[,2]), 
                                              digit=2)),
                          x=0.1,  y=0.95, 
                          gp=gpar(col="red3", fontsize=12, fontface="bold"))
      g1 <- ggplot(DF_filt, aes(x=log2(DF_filt[,1]), y=log2(DF_filt[,2]))) + 
        geom_point() + geom_smooth(method="lm",col="red3") + 
        annotation_custom(my_grob) + 
        labs(y=colnames(DF_filt[2]), x=colnames(DF_filt[1]))
      g2 <- ggplot(DF_filt, aes(x=DF_filt[,1], y=DF_filt[,2])) + 
        geom_point() + 
        geom_smooth(method="lm",col="red3") + 
        annotation_custom(my_grob) + 
        labs(y=colnames(DF_filt[2]), x=colnames(DF_filt[1]))
      
      if (l) {
        forceSetWd('../log2_correlation')
        ggsave(paste0('correlation_log2_between_', i, '_and_', j, '.png'), 
               plot = g1, width = 7, height = 7)
      } else {
        forceSetWd('../correlation')
        ggsave(paste0('correlation_between_', i, '_and_', j, '.png'), 
               plot = g2, width = 7, height = 7)
      }
    }
  }
}
