source('/share/script/hecatos/juantxo/analysis_trc/functions.R')
forceLibrary(c('limma', 'dplyr'))
getCts <- function(project_name) {
  if (grepl('TRC', project_name)) {
    setwd('/share/analysis/hecatos/juantxo/Score/output/UNTR')
    if (grepl('counts', project_name)) {
      setwd('/share/analysis/hecatos/juantxo/Score/analysis/TRC_counts/')
      quant_file = read.table(
        'counts_TRC_UNTR_2019-11-11_12:04:38_UTC.tsv', header = T)
    }
    if (grepl('voom', project_name)) {
      setwd('2019-11-11_09:41:21_UTC/')
      quant_file = mergeFiles(files_patt = 'TRCscore', row_names = T, all = T)
    } else {
      setwd('2019-11-11_12:04:38_UTC/')
      quant_file = read.table('TRCscore_total.txt', header = T)
    }
  }
  
  if (grepl('Salmon', project_name)) {
    setwd('/share/analysis/hecatos/juantxo/Score/input/RNA_Salmon/UNTR/')
    quant_file = read.table('total_quant.sf', header = T)
  }
  return(quant_file)
}
transformLogTPM = function(x) {
  x[is.na.data.frame(x)] = 0.25
  x[x == 0] = 0.25
  y = log2(x)
}

#########################
salmon = getCts('Salmon_TPM')
# salmon_tpm = salmon %>% 
#   dplyr::select(contains('TPM')) %>%
#   mutate(Name = salmon$Name) %>%
#   filter(grepl('ENST', Name)) %>%
#   transcrToGene(aggregate = T) %>%
#   tibble::column_to_rownames('ensembl_gene_id')



trc_output = getCts('TRC')

trc = trc_output %>%
  dplyr::select(contains('TRC_')) %>%
  mutate(Name = x$rowname) %>%
  filter(grepl('ENST', Name)) %>%
  tibble::column_to_rownames('Name') 

processLima = function(x, prot_cod = F, DETs = F) {
  if (prot_cod) {
    if (DETs) {
      trc = trc %>%
        filterProtCod() %>% 
        dplyr::select(-rownames)
    } else {
      trc =  trc %>% 
        transcrToGene(aggregate = T, prot_cod = T) %>% 
        tibble::column_to_rownames('ensembl_gene_id')
    }
    
  } else {
    if (!DETs) {
      trc =  trc %>% 
        transcrToGene(aggregate = T, prot_cod = F) %>% 
        column_to_rownames('ensembl_gene_id')
    }
  }
  return(trc)
}

  log_trc = transformLogTPM(trc)

  timepoints = NULL
  for (tp in c('UNTR_002', 'UNTR_008','UNTR_024', 'UNTR_072')) {
    repls = rep(tp, 3)
    timepoints = c(timepoints, repls)
  }
  
  coldata = data.frame(row.names = colnames(log_tpm), 
                       timepoints = timepoints)
  all(rownames(coldata) == colnames(log_tpm)) %>% stopifnot()
  
  tps_2 = seq(4, 10, 3)
  
  i = 0; rm(DEGs, DEGs_final); x = log_trc
  
  for (tp_2 in tps_2) {
    tp_2 = tp_2:(tp_2 + 2)
    design = model.matrix(~ timepoints[c(1:3, tp_2)])
    
    fit = lmFit(x[, c(1:3, tp_2)], design = design[, 2])
    fit = eBayes(fit, trend = T)
    DEGs = topTable(fit, p.value = 0.05, number = Inf)
    if (length(DEGs) > 0) {
      # if (DETs) {
      #   DEGs = transcrToGene(DEGs, aggregate = T) %>% 
      #     column_to_rownames('ensembl_gene_id')
      # }
      colnames(DEGs) = paste(colnames(DEGs), timepoints[tp_2[1]], sep = '_')
      DEGs = rownames_to_column(DEGs)
      
      first_DEG = i == 0
      if (first_DEG) {
        DEGs_final = DEGs
      } else {
        DEGs_final = merge.data.frame(x = DEGs_final, 
                                      y = DEGs, 
                                      by = 'rowname', 
                                      all = T)
      }
      i = i + 1
      print(length(rownames(DEGs)))
    }
  }

  



setwd("/share/script/hecatos/juantxo/analysis_trc/DEGs/limma/output")
write.table(x = DEGs, file = 'DEGs_trc_withoutNCRNA.tsv', sep = '\t')
DETs = DEGs
DEGs = transcrToGene(DEGs, aggregate = T) %>% 
  column_to_rownames('ensembl_gene_id')
paste(rownames(DEGs), collapse = ' ')
