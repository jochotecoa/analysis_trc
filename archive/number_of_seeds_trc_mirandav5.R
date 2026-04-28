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

greppend <- function(pattern, x, x.col.input = colnames(x),
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

setwd('/share/analysis/hecatos/juantxo/Score/input/Interactions/')
miRanda.table = read.table(file = 'ALLmRNAtranscripts_vs_ALLmiRNA3.parsed', 
                           stringsAsFactors = F)
colnames(miRanda.table)[1:2] = c('miRNA.ID', 'mRNA.ID')

setwd("/share/analysis/hecatos/juantxo/Score/output/Output_Run_mrna_SEPT2019/")
setwd('V3/analysis/UNTR/other/')
proteins.table = read.table('20190916_low_correlated_TPMvsTRC.tab')

# setwd('/share/analysis/hecatos/juantxo/tableomics/')
# protrans = read.csv('./uniprot-yourlist_M20190501.csv2', header = T, 
#                     stringsAsFactors = F)


##### GET PROTEIN NAMES TO UNIQUE MRNA NAMES #####


mart.human = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                     dataset = 'hsapiens_gene_ensembl',
                     host = 'http://apr2018.archive.ensembl.org') 

# mart.attributes = searchAttributes(mart = mart.human)
# mart.filters = searchFilters(mart.human)

trc.mrnas.table = getBM(attributes = c('ensembl_transcript_id', 
                                       'transcript_biotype',
                                       'uniprot_gn'),
                        filters = 'uniprot_gn',
                        values = list(rownames(proteins.table)), 
                        mart = mart.human)

# We only want the interactions of the transcripts that translate to the protein
prot_cod.log = trc.mrnas.table$transcript_biotype == 'protein_coding'
trc.mrnas.table = trc.mrnas.table[prot_cod.log, ]

inter.trc = greppend(pattern = trc.mrnas.table[, 1], 
                     x = miRanda.table, 
                     x.col.input = 'mRNA.ID')

mRNA.IDs.subset = unique(inter.trc$mRNA.ID)
mRNA.IDs.global = unique(miRanda.table$mRNA.ID)
seeds_per_mRNA.subset = nrow(inter.trc) / length(mRNA.IDs.subset)
seeds_per_mRNA.global = nrow(miRanda.table) / length(mRNA.IDs.global)
seeds.diff = seeds_per_mRNA.subset / seeds_per_mRNA.global
print(seeds.diff)

ex = inter.trc$mRNA.ID[1]
n_seeds_ex = length(grep(pattern = ex, x = miRanda.table$mRNA.ID))
setwd('../../../output/TRCscore/')
trc.table = read.table(file = 'UNTR_002_1_TRCscore.txt')
trc_ex.row = grep(pattern = ex, x = rownames(trc.table))
n_seeds_ex.trc_table = trc.table$SUMseeds.targetRNA[trc_ex.row]
print('Are the number of seeds in miranda the same as in the output table?')
print(n_seeds_ex.trc_table == n_seeds_ex)
print('Are the number of seeds in the output table equal across samples?')
trc.table = read.table(file = 'UNTR_002_2_TRCscore.txt')
trc_ex.row = grep(pattern = ex, x = rownames(trc.table))
n_seeds_ex.trc_table2 = trc.table$SUMseeds.targetRNA[trc_ex.row]
print(n_seeds_ex.trc_table == n_seeds_ex.trc_table2)




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
