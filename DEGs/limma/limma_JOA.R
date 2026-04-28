#### Functions ####
source('/share/script/hecatos/juantxo/analysis_trc/functions.R')
forceLibrary(c('limma', 'dplyr', 'tibble', 'VennDiagram', 'UpSetR'))
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
processLima = function(x, prot_cod = F, DETs = F) {
  if (prot_cod) {
    if (DETs) {
      x = x %>%
        filterProtCod() %>% 
        dplyr::select(-rownames)
    } else {
      x =  x %>% 
        transcrToGene(aggregate = T, prot_cod = T) %>% 
        tibble::column_to_rownames('ensembl_gene_id')
    }
    
  } else {
    if (!DETs) {
      x =  x %>% 
        transcrToGene(aggregate = T, prot_cod = F) %>% 
        column_to_rownames('ensembl_gene_id')
    }
  }
  return(x)
}

doLimma = function(x, control = NULL, samples, all_genes = F) {
  samples_factor = as.factor(samples)
  if (is.null(control)) {control = levels(samples_factor)[1]}
  treatment = levels(samples_factor) %>% .[. != control]
  i = 0
  design = model.matrix(~ samples)
  treatment = paste0('samples', treatment)
  fit = lmFit(x, design = design)
  fit = eBayes(fit, trend = T)
  for (treatm in treatment) {
    DEGs = topTable(fit, p.value = 0.05, number = Inf, coef = treatm)
    if (length(DEGs) > 0 & all_genes == F) {
      colnames(DEGs) = paste(colnames(DEGs), treatm, sep = '_')
      DEGs = rownames_to_column(DEGs)
      
      first_DEG = i == 0
      if (!exists('DEGs_final')) {
        DEGs_final = DEGs
        i = i + 1
      } else {
        DEGs_final = merge.data.frame(x = DEGs_final, 
                                      y = DEGs, 
                                      by = 'rowname', 
                                      all = T)
      }
      
      print(length(rownames(DEGs)))
    } else {
      non_DEGs = topTable(fit, number = Inf, coef = treatm)
      first_DEG = i == 0
      
      colnames(non_DEGs) = paste(colnames(non_DEGs), treatm, sep = '_')
      non_DEGs = rownames_to_column(non_DEGs)
      if (!exists('non_DEGs_final')) {
        non_DEGs_final = non_DEGs
        i = i + 1
      } else {
        non_DEGs_final = merge.data.frame(x = non_DEGs_final, 
                                          y = non_DEGs, 
                                          by = 'rowname', 
                                          all = T)
      }
    }
  }
  if (exists('DEGs_final') & exists('non_DEGs_final')) {
    return(list(DEGs_final, non_DEGs_final))
  }
  if (exists('DEGs_final')) {return(DEGs_final)} else {
    print('No DEGs found')
    return(non_DEGs_final)
  }
}
doLimmaSingleTp = function(x, timepoints, all_genes = F) {
  i = 0
  if (exists('DEGs_final')) {rm('DEGs_final')}
  tps_2 = seq(4, 10, 3)
  for (tp_2 in tps_2) {
    tp_2 = tp_2:(tp_2 + 2)
    design = model.matrix(~ timepoints[c(1:3, tp_2)])
    
    fit = lmFit(x[, c(1:3, tp_2)], design = design[, 2])
    fit = eBayes(fit, trend = T)
    if (all_genes) {
      DEGs = topTable(fit, number = Inf)
    } else {
      DEGs = topTable(fit, p.value = 0.05, number = Inf)
    }
    if (length(DEGs) > 0) {
      # if (DETs) {
      #   DEGs = transcrToGene(DEGs, aggregate = T) %>% 
      #     column_to_rownames('ensembl_gene_id')
      # }
      colnames(DEGs) = paste(colnames(DEGs), timepoints[tp_2[1]], sep = '_')
      DEGs = rownames_to_column(DEGs)
      
      if (!exists('DEGs_final')) {
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
  if (i != 0) {return(DEGs_final)} else {print('No DEGs')}
}

#########################
setwd('/share/analysis/hecatos/juantxo/Score/input/RNA_Salmon/UNTR/')
salmon = read.table('total_quant.sf', header = T)
salmon_tpm = salmon %>%
  dplyr::select(contains('TPM')) %>%
  mutate(Name = salmon$Name) %>%
  filter(grepl('ENST', Name)) %>%
  tibble::column_to_rownames('Name')

# salmon_tpm[salmon_tpm < 1] = 0

setwd('/share/analysis/hecatos/juantxo/Score/output/UNTR/2019-11-11_12:04:38_UTC/')
trc_output = read.table('TRCscore_total.txt', header = T)
trc = trc_output %>%
  dplyr::select(contains('TRC_')) %>%
  mutate(Name = trc_output$rowname) %>%
  filter(grepl('ENST', Name)) %>%
  tibble::column_to_rownames('Name') 

trc = processLima(x = trc, prot_cod = T, DETs = F)
tpm = processLima(x = salmon_tpm, prot_cod = T, DETs = F)
stopifnot(!identical(tpm, trc))

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



names = rownames(log_tpm) %>% substr(1, 15) 
log_tpm = log_tpm[names %in% rownames(log_trc), ]
log_trc = log_trc[rownames(log_trc) %in% names, ]
trc_degs = doLimmaSingleTp(x = log_trc, timepoints = timepoints, all_genes = T)
tpm_degs = doLimmaSingleTp(x = log_tpm, timepoints = timepoints, all_genes = T)
# lmFit(log_tpm, model.matrix(~ timepoints)) %>% eBayes(trend = T) %>% decideTests() %>% summary()

forceSetWd('/share/analysis/hecatos/juantxo/Score/analysis/DEPs')
prot_log_norm = read.table(file = 'UNTR_protein_log2_pvalues.tsv', sep = '\t')       
prot_DEGs_072 = rownames(prot_log_norm) %>% 
  .[prot_log_norm$p.value_UNTR_The_072_1 < 0.05] %>%
  getBM(attributes = c('ensembl_gene_id', 'uniprot_gn'), 
        filters = 'uniprot_gn', values = ., mart = mart.human)
degs_names = getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                   filters = 'ensembl_gene_id', values = tpm_degs$rowname, 
                   mart = mart.human)
trc_degs = trc_degs %>%
  merge.data.frame(y = degs_names, by.x = 'rowname', by.y = 'ensembl_gene_id')
tpm_degs = tpm_degs %>%
  merge.data.frame(y = degs_names, by.x = 'rowname', by.y = 'ensembl_gene_id')

prot_DEGs_072 = trc_degs %>% dplyr::select(rowname, contains('adj.P')) %>% 
  merge.data.frame(y = ., x = prot_DEGs_072, by.x = 'ensembl_gene_id', 
                   by.y = 'rowname')
colnames(prot_DEGs_072) = colnames(prot_DEGs_072) %>% 
  gsub(pattern = 'adj.P', replacement = 'TRC_adj.P')
prot_DEGs_072 = tpm_degs %>% dplyr::select(rowname, contains('adj.P')) %>% 
  merge.data.frame(y = ., x = prot_DEGs_072, by.x = 'ensembl_gene_id', 
                   by.y = 'rowname')
colnames(prot_DEGs_072) = colnames(prot_DEGs_072) %>% 
  gsub(pattern = '^adj.P', replacement = 'TPM_adj.P')

ggprot_DEGs_072 = melt(prot_DEGs_072, id = c('ensembl_gene_id', 'uniprot_gn'))
UnitandTimepoint = gsub('_adj.P.Val', '', ggprot_DEGs_072$variable)
ggplot(data = ggprot_DEGs_072) + 
  geom_density(mapping = aes(x = value, colour = UnitandTimepoint)) + 
  label_value(labels = 'Adj.P.value')
'Transcriptomics Adj.P.Values for DEPs at 72h'