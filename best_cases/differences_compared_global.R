forceLibrary('VennDiagram')
comp = '5FU'
source('/share/script/hecatos/juantxo/analysis_trc/p.values_multiomics/checkSeveralOmicsTrTFile_JOA.R')

names(corrctd_theVStox) = c('tp_trt', 'tn_trt', 'tp_tpm', 'tn_tpm')
names(corrctd_untVSthe) = c('tp_trt', 'tn_trt', 'tp_tpm', 'tn_tpm')
names(corrctd_untVStox) = c('tp_trt', 'tn_trt', 'tp_tpm', 'tn_tpm')

a = c(corrctd_untVSthe$tp_trt, corrctd_untVSthe$tn_trt,
      corrctd_untVStox$tp_trt, corrctd_untVStox$tn_trt) %>% 
  gsub('_.*', '', .) %>% 
  unique()

setwd('/share/analysis/hecatos/juantxo/Score/analysis/')


# write.table(x = a, file = 'trt_corrected', quote = F, 
#             row.names = F, col.names = F)

trt_proteomx_best = c(corrctd_untVSthe$tp_trt, corrctd_untVSthe$tn_trt,
      corrctd_untVStox$tp_trt, corrctd_untVStox$tn_trt) %>% 
  trt_proteomx[., ]


compare_trt_all_means = function(x, y, comp_col){
  x = x %>%
    dplyr::select(starts_with(comp_col)) %>% 
    rowMeans(na.rm = T) %>% 
    mean(na.rm = T)
  
  y = y %>% 
    dplyr::select(starts_with(comp_col)) %>% 
    rowMeans(na.rm = T) %>% 
    mean(na.rm = T) 
  
  z = print(x/y)
  
}

compare_trt_all_medians = function(x, y, comp_col){
  x = x %>%
    dplyr::select(starts_with(comp_col)) %>% 
    apply(1, median, na.rm = T) %>% 
    median(na.rm = T)
  
  y = y %>% 
    dplyr::select(starts_with(comp_col)) %>% 
    apply(1, median, na.rm = T) %>% 
    median(na.rm = T)
  
  z = print(x/y)
  
}


compare_trt_all_medians(trt_proteomx_best, trt_proteomx, "sum_Seed_proportion_")
compare_trt_all_medians(trt_proteomx_best, trt_proteomx, "int_miRNA_sum_")
compare_trt_all_medians(trt_proteomx_best, trt_proteomx, "int_miRNA_min_")
compare_trt_all_medians(trt_proteomx_best, trt_proteomx, "int_miRNA_max_")
compare_trt_all_medians(trt_proteomx_best, trt_proteomx, "int_miRNA_mean_")
compare_trt_all_medians(trt_proteomx_best, trt_proteomx, "int_miRNA_median_")
compare_trt_all_medians(trt_proteomx_best, trt_proteomx, "int_miRNA_sd_")
compare_trt_all_medians(trt_proteomx_best, trt_proteomx, "targetRNA_TPM_")
compare_trt_all_medians(trt_proteomx_best, trt_proteomx, "SUMseeds.targetRNA_")
compare_trt_all_medians(trt_proteomx_best, trt_proteomx, "X.miRNAs_")
compare_trt_all_medians(trt_proteomx_best, trt_proteomx, "sum_miRNA_TPM_")
compare_trt_all_medians(trt_proteomx_best, trt_proteomx, "E_circ_sum_")
compare_trt_all_medians(trt_proteomx_best, trt_proteomx, "E_circ_N_")
compare_trt_all_medians(trt_proteomx_best, trt_proteomx, "E_decoy_sum_")
compare_trt_all_medians(trt_proteomx_best, trt_proteomx, "E_decoy_N_")
compare_trt_all_medians(trt_proteomx_best, trt_proteomx, "TrT_0.1_")


