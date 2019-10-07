list.of.packages <- c("biomaRt")
inst.packages = installed.packages()[,"Package"]
new.packages <- list.of.packages[!(list.of.packages %in% inst.packages)]
if (length(new.packages)) install.packages(new.packages)
# setRepositories()

library(biomaRt) # setRepositories() # install.packages('biomaRt')

tablePct <- function(x) {
  x.table = table(x)
  x.proportion = x.table / length(x)
  x.percent = x.proportion * 100
}

setwd('/share/analysis/hecatos/juantxo/score_protein_analysis/data_extra')
mirtarbase.table = read.csv(file = 'hsa_MTI.csv', 
                            header = T, 
                            stringsAsFactors = F)
setwd("/share/analysis/hecatos/juantxo/Score/output/Output_Run_mrna_SEPT2019/")
setwd('V3/analysis/UNTR/')
proteins.table = read.table('20190916_low_correlated_TPMvsTRC.tab')
proteins.naam = rownames(proteins.table)
mart.human = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                     dataset = 'hsapiens_gene_ensembl',
                     host = 'http://apr2018.archive.ensembl.org') 

# View(listAttributes(mart = mart.human))
mirtarbase.Gn.UpId = getBM(attributes = c('external_gene_name', 
                                                'uniprot_gn'),
                           filters = 'uniprot_gn',
                           values = proteins.naam, 
                           mart = mart.human)

rows = mirtarbase.table$Target.Gene %in% mirtarbase.Gn.UpId[, 1]
mirtarbase.table.trc = mirtarbase.table[rows, ]

trc.supptype.pct = tablePct(mirtarbase.table.trc$Support.Type) 
mirtar.supptype.pct = tablePct(mirtarbase.table$Support.Type) 
mirtar.TrcVsAll.diff = trc.supptype.pct - mirtar.trc.supptype.pct

