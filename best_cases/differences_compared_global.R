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


write.table(x = a, file = 'trt_corrected', quote = F, 
            row.names = F, col.names = F)

trt_proteomx_best = c(corrctd_untVSthe$tp_trt, corrctd_untVSthe$tn_trt,
      corrctd_untVStox$tp_trt, corrctd_untVStox$tn_trt) %>% 
  trt_proteomx[., ]


compare_trt_all = function(x, y, comp_col){
  x %>%
    dplyr::select(starts_with(comp_col)) %>% 
    rowMeans(na.rm = T) %>% 
    mean(na.rm = T) %>% 
    print
  
  y %>% 
    dplyr::select(starts_with(comp_col)) %>% 
    rowMeans(na.rm = T) %>% 
    mean(na.rm = T) %>% 
    print
  
}


compare_trt_all(trt_proteomx_best, trt_proteomx, 'E_decoy_N')

X.miRNAs_