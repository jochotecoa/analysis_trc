##### LIBRARIES #####
list.of.packages <- c("pbmcapply", "foreach", "doParallel", "biomaRt")
new.packages.log = !(list.of.packages %in% installed.packages()[,"Package"])
new.packages <- list.of.packages[new.packages.log]
if (length(new.packages)) install.packages(new.packages)
new.packages.log = !(list.of.packages %in% installed.packages()[,"Package"])
new.packages <- list.of.packages[new.packages.log]
if (length(new.packages)) {
  setRepositories(graphics = F, ind = c(1,2))
  install.packages(new.packages)
} 
# Progress bar library
library(pbmcapply)
# Parallel computing libraries
library(foreach)
library(doParallel)
library(biomaRt)

##### FUNCTIONS #####

greppend <- function(pattern, 
                     x, 
                     x.col.input = colnames(x),
                     x.col.output = colnames(x)) {
  pb = progressBar(min = 0, max = length(pattern))
  new.table = data.frame()
  for (variable in 1:length(pattern)) {
    pos = grep(pattern = pattern[variable], x = x[, x.col.input])
    new.table = rbind.data.frame(new.table, x[pos, x.col.output])
    setTxtProgressBar(pb, value = variable)
  }
  close(pb)
  return(new.table)
}


##### OPENING FILES #####

setwd('/share/analysis/hecatos/juantxo/Score/Interactions/')
miRanda.table = read.table(file = 'ALLmRNAtranscripts_vs_ALLmiRNA3.parsed', 
                           stringsAsFactors = F)
setwd("/share/analysis/hecatos/juantxo/score_protein_analysis/UNTR/shift-FALSE")
proteins.table = read.table('./TRC/most_correlated_proteins_trc.tsv')
colnames(miRanda.table)[1:2] = c('miRNA.ID', 'mRNA.ID')
# setwd('/share/analysis/hecatos/juantxo/tableomics/')
# protrans = read.csv('./uniprot-yourlist_M20190501.csv2', header = T, 
#                     stringsAsFactors = F)


##### GET PROTEIN NAMES TO UNIQUE MRNA NAMES #####


mart.human = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                     dataset = 'hsapiens_gene_ensembl') 

mart.attributes = searchAttributes(mart = mart.human)
mart.filters = searchFilters(mart.human)

trc.mrnas.table = getBM(attributes = c('ensembl_transcript_id', 
                                       'transcript_biotype',
                                       'uniprot_gn_id'),
                        filters = 'uniprot_gn_id',
                        values = list(proteins.table$protlist), 
                        mart = mart.human)

# We only want the interactions of the transcripts that translate to the protein
prot_cod.log = trc.mrnas.table$transcript_biotype == 'protein_coding'
trc.mrnas.table = trc.mrnas.table[prot_cod.log, ]

inter.trc = greppend(pattern = trc.mrnas.table[, 1], 
                     x = miRanda.table, 
                     x.col.input = 'mRNA.ID')
setwd('/share/analysis/hecatos/juantxo/score_protein_analysis/UNTR/')
write.table(x = inter.trc, 
            file = 'topTrcInteractionsMirna.tsv', 
            sep = '\t',
            row.names = F)



# mrnas.unique = as.data.frame(unique(miRanda.table$mRNA.ID))
# colnames(mrnas.unique) = 'transcripts'
# mrnas.unique$protein = NA

# pb = progressBar(min = 0, max = nrow(mrnas.unique))
# 
# mrnas.notfound = 0
# 
# for (i in 1:nrow(mrnas.unique)) {
#   mrna.protrans.pos = grep(pattern = mrnas.unique[i, 1], 
#                           x = protrans[, 1])
#   if (length(mrna.protrans.pos) > 0) {
#     mrnas.unique$protein[i] = protrans$Entry[mrna.protrans.pos]
#   } else {
#     mrnas.notfound = mrnas.notfound + 1
#   }
#   if (i == nrow(mrnas.unique)) {
#     print(paste('Number of mRNAs not found:', mrnas.notfound))
#   }
#   setTxtProgressBar(pb, value = i)
# }
# close(pb)

# a = greppend(pattern = proteins.table$protlist,
#              x = protrans, 
#              x.col.input = 'Entry',
#              x.col.output = 1)

# trc.mrnas.table = greppend(pattern = proteins.table$protlist, 
#                            x = mrnas.unique, x.col.input = 'protein')


# miRanda.df = as.data.frame(miRanda.table)
# miRanda.df$protein = NA
# 
# 
# ################# PARALLEL COMPUTING #################
# no_cores <- detectCores() / 2
# cl <- makeCluster(no_cores)
# registerDoParallel(cl)
# 
# # pb = progressBar(min = 0, max = nrow(mrnas.unique))
# 
# foreach(i = 1:nrow(mrnas.unique), .combine = rbind, .verbose = T) %dopar% {
#   pos = grep(pattern = mrnas.unique[i, 1], x = miRanda.df$mRNA.ID)
#   miRanda.df$protein[pos] = mrnas.unique$protein[i]
#   miRanda.df
#   # setTxtProgressBar(pb, value = i)
# }
# # close(pb)
# stopCluster(cl)
# 
# ##### ASSIGN PROTEINS TO MIRANDA DF ####
# pb = progressBar(min = 0, max = nrow(mrnas.unique))
# 
# for (i in 1:nrow(mrnas.unique)) {
#   pos = grep(pattern = mrnas.unique[i, 1], x = miRanda.df$mRNA.ID)
#   miRanda.df$protein[pos] = mrnas.unique$protein[i]
#   setTxtProgressBar(pb, value = i)
# }
# close(pb)
# 
# # ##### COUNT SEEDS #####
# # 
# # proteins.table$seeds = NA
# # 
# # for (i in 1:nrow(proteins.table)) {
# #   mrna.miranda.pos = grep(pattern = proteins.table$protlist[i], 
# #                           x = miRanda.df$protein)
# #   mrna.miranda.seeds = length(mrna.miranda.pos)
# #   proteins.table$seeds[i] = mrna.miranda.seeds
# #   print(mrna.miranda.seeds)
# # }
# # 
# # summary(miRanda.table$mRNA.ID)
