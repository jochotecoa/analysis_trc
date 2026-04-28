list.of.packages <- c("biomaRt")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(new.packages)
# setRepositories()

library(biomaRt) # setRepositories() # install.packages('biomaRt')


setwd('/share/analysis/hecatos/juantxo/score_protein_analysis/data_extra')
mirtarbase.table = read.csv(file = 'hsa_MTI.csv', 
                            header = T, 
                            stringsAsFactors = F)
mart.human = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                     dataset = 'hsapiens_gene_ensembl') 
head(searchFilters(mart.human, 'Gene Name'))
ensembl.prot.ids.mouse = getBM(attributes = c('ensembl_gene_id', 'uniprot_gn_id'),
                               filters = 'external_gene_name',
                               values = list(mirtarbase.table$Target.Gene), 
                               mart = mart.human)
