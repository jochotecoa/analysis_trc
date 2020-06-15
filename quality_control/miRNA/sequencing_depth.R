library(dplyr)
source("/share/script/hecatos/juantxo/analysis_trc/functions_JOA.R")

setwd('/share/analysis/hecatos/juantxo/Score/input/miRNA_miRge2/')
count_files = list.files(pattern='Counts', recursive=T)
count_files = count_files[!count_files %in% "UNTR/miR.Counts_voom.csv"]

for (count_file in count_files) {
  count_file %>% 
    read.csv(row.names = 1) %>% 
    .["miRNAtotal", , F] %>% 
    filterSamplesBySeqDepth()
}
