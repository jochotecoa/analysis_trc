library(dplyr)
source("/share/script/hecatos/juantxo/analysis_trc/functions_JOA.R")

setwd('/share/analysis/hecatos/juantxo/Score/input/RNA_Salmon/')
count_files = list.files(pattern='total_quant_2.sf', recursive=T)

for (count_file in count_files) {
  read.table(count_file) %>% 
    select(contains('NumReads')) %>% 
    filterSamplesBySeqDepth()
}


