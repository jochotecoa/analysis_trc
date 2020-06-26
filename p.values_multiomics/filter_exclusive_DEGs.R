library('VennDiagram')
setwd("/share/script/hecatos/juantxo/analysis_trc")
source('functions_JOA.R')
forceLibrary(c('dplyr', 'tibble'))



if (!exists('trt_df')) {
  comp = '5FU'
  # source('/share/script/hecatos/juantxo/analysis_trc/p.values_multiomics/checkSeveralOmicsTrTFile_JOA.R')
}

setwd('/share/analysis/hecatos/juantxo/Score/analysis/t-test/')

comps = 
  list.dirs(recursive=F,
            path = "/share/analysis/hecatos/juantxo/Score/analysis/t-test")

for (comp in comps) {
  forceSetWd(comp)
  forceSetWd('DEGs')
  
  output.files = list.files(pattern = 'p.adj.*_TPM', recursive = T)
  items_excl = output.files %>% grepl('exclusive', .) 
  excl.files = output.files[items_excl]
  for (excl.file in excl.files) {
    file.remove(excl.file)
  }
  output.files = list.files(pattern = 'p.adj.*_TPM', recursive = T)
  
  
  for (output.file in output.files) {
    tpm_naam = output.file
    trt_naam = gsub('TPM', 'TrT', tpm_naam)
    
    if (file.exists(trt_naam)) {
      print(output.file)
      tpm_file = read.table(tpm_naam)
      trt_file = read.table(trt_naam)
      
      excl_tpm = tpm_file[!(tpm_file[, 1] %in% trt_file[, 1]), ]
      excl_trt = trt_file[!(trt_file[, 1] %in% tpm_file[, 1]), ]
      
      nieuw_naam_tpm = tpm_naam %>% 
        strsplit('/') %>% 
        unlist 
      
      dir.create(paste(nieuw_naam_tpm[1], 'exclusive_DEGs', sep = '/'))
      
      nieuw_naam_tpm = 
        paste(nieuw_naam_tpm[1], 'exclusive_DEGs', nieuw_naam_tpm[2], sep = '/') 
      
      
      nieuw_naam_trt = trt_naam %>% 
        strsplit('/') %>% 
        unlist  
      nieuw_naam_trt = 
        paste(nieuw_naam_trt[1], 'exclusive_DEGs', nieuw_naam_trt[2], sep = '/') 
      
      forceSetWd(paste(nieuw_naam_tpm[1], 'exclusive_DEGs', sep = '/'))
      
      write.table(x = excl_tpm, file = nieuw_naam_tpm, quote = F, 
                  row.names = F, col.names = F)
      write.table(x = excl_trt, file = nieuw_naam_trt, quote = F, 
                  row.names = F, col.names = F)
      
    }
    
  }
  
}

