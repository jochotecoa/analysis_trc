source('/share/script/hecatos/juantxo/analysis_trc/DEGs/limma/limma_JOA.R')
doLimmaSingleTp = function(x, timepoints) {
  i = 0; rm(DEGs, DEGs_final)
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
  return(DEGs_final)
}

trc = processLima(x = trc, prot_cod = F, DETs = T)
tpm = processLima(x = salmon_tpm, prot_cod = F, DETs = T)

log_trc = transformLogTPM(trc)
log_tpm = transformLogTPM(tpm)

timepoints = NULL
for (tp in c('UNTR_002', 'UNTR_008','UNTR_024', 'UNTR_072')) {
  repls = rep(tp, 3)
  timepoints = c(timepoints, repls)
}
coldata = data.frame(row.names = colnames(log_tpm), 
                     timepoints = timepoints)
all(rownames(coldata) == colnames(log_tpm)) %>% stopifnot()
treatment = seq(4, 10, 3)




trc_degs = doLimmaSingleTp(x = log_trc, timepoints = timepoints)

names = rownames(log_tpm) %>% substr(1, 15) 
log_tpm = log_tpm[names %in% rownames(log_trc), ]

tpm_degs = doLimmaSingleTp(x = log_tpm, timepoints = timepoints)

names = tpm_degs$rowname %>% substr(1, 15) %>% 
  setdiff(trc_degs$rowname) 

getBM(attributes = c('ensembl_transcript_id', 'external_gene_name'), 
      filters = 'ensembl_transcript_id', 
      values = names, mart = mart.human)

  





# setwd("/share/script/hecatos/juantxo/analysis_trc/DEGs/limma/output")
# write.table(x = DEGs, file = 'DEGs_trc_withoutNCRNA.tsv', sep = '\t')
# DETs = DEGs
# DEGs = transcrToGene(DEGs, aggregate = T) %>% 
#   column_to_rownames('ensembl_gene_id')
# paste(rownames(DEGs), collapse = ' ')

results = data.frame()

for (thresh in seq(0,0.00001,0.000001)) {
  salmon_tpm = salmon %>%
    dplyr::select(contains('TPM')) %>%
    mutate(Name = salmon$Name) %>%
    filter(grepl('ENST', Name)) %>%
    tibble::column_to_rownames('Name')
  salmon_tpm[salmon_tpm < thresh] = 0
  log_tpm_2 = transformLogTPM(salmon_tpm)
  names = rownames(log_tpm_2) %>% substr(1, 15) 
  log_tpm_2 = log_tpm_2[names %in% rownames(log_trc), ] 
  fit = lmFit(log_tpm_2, design = model.matrix(~ timepoints))
  fit = eBayes(fit, trend = T)
  DEGs = topTable(fit, p.value = 0.05, number = Inf)
  results = rbind(results, c(thresh, nrow(DEGs)))
}
decideTests(fit) %>% summary()
plot(results[, 1], results[, 2], type = 'b', 
     xlab = 'Threshold', ylab = 'Number of DETs', 
     main = 'Number of DETs at different cutoff values', 
     sub = 'TPM_UNTR_002 vs TPM_UNTR_072')


tpm_degs$rowname = tpm_degs$rowname %>% substr(1,15) 

degs_trc = topBy(x = trc_degs, by = 'adj.P.Val_samplesUNTR_008', n = 10) %>% 
  merge(topBy(x = trc_degs, by = 'adj.P.Val_samplesUNTR_024', n = 10), by = 'rowname', all = T) %>%
  merge(topBy(x = trc_degs, by = 'adj.P.Val_samplesUNTR_072', n = 10), by = 'rowname', all = T) 

degs_tpm = topBy(x = tpm_degs, by = 'adj.P.Val_samplesUNTR_008', n = 10) %>% 
  merge(topBy(x = tpm_degs, by = 'adj.P.Val_samplesUNTR_024', n = 10), by = 'rowname', all = T) %>%
  merge(topBy(x = tpm_degs, by = 'adj.P.Val_samplesUNTR_072', n = 10), by = 'rowname', all = T) 


degs = merge.data.frame(x = degs_tpm, y = degs_trc, by = 'rowname', all = T) %>%
  select(contains('adj.P')) %>% makeBinaryMatrix()
topBy = function(x,by,n) {
  x = x[order(x[, by]), ]
  x = x[1:n, c(1, grep(by, colnames(x)))]
  return(x)
}

degs = merge.data.frame(x = trc_degs, y = tpm_degs, by = 'rowname', all = T) %>% 
  select(contains('Adj.P')) 
colnames(degs) = gsub('.x', '.TRC', colnames(degs))
colnames(degs) = gsub('.y', '.TPM', colnames(degs))
makeBinaryMatrix = function(x) {
  x[!is.na.data.frame(x)] = 1
  x[is.na.data.frame(x)] = 0
  return(x)
}

upset(data = degs, nsets = 6, matrix.color = 'blue', sets.bar.color = 'red')
