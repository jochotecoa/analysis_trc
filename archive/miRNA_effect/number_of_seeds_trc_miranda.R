##### LIBRARIES #####
list.of.packages <- c("pbmcapply", "foreach", "doParallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(new.packages)

library(pbmcapply)
library(foreach)
library(doParallel)

##### OPENING FILES

setwd('/share/analysis/hecatos/juantxo/Score/Interactions/')
miRanda.table = read.table(file = 'ALLmRNAtranscripts_vs_ALLmiRNA3.parsed', 
                           stringsAsFactors = F)
setwd("/share/analysis/hecatos/juantxo/score_protein_analysis/UNTR/shift-FALSE")
proteins.table = read.table('./TRC/most_correlated_proteins_trc.tsv')
colnames(miRanda.table)[1:2] = c('miRNA.ID', 'mRNA.ID')
setwd('/share/analysis/hecatos/juantxo/tableomics/')
protrans = read.csv('./uniprot-yourlist_M20190501.csv2', header = T, 
                    stringsAsFactors = F)


##### GET PROTEIN NAMES TO UNIQUE MRNA NAMES #####

mrnas.unique = as.data.frame(unique(miRanda.table$mRNA.ID))
mrnas.unique$protein = NA

pb = progressBar(min = 0, max = nrow(mrnas.unique))

for (i in 1:nrow(mrnas.unique)) {
  mrnas.notfound = 0
  mrna.protrans.pos = grep(pattern = mrnas.unique[i, 1], 
                          x = protrans[, 1])
  if (length(mrna.protrans.pos) > 0) {
    mrnas.unique$protein[i] = protrans$Entry[mrna.protrans.pos]
  } else {
    mrnas.notfound = mrnas.notfound + 1
  }
  if (i == nrow(mrnas.unique)) {
    print(paste('Number of mRNAs not found:', mrnas.notfound))
  }
  setTxtProgressBar(pb, value = i)
}
close(pb)

miRanda.df = as.data.frame(miRanda.table)
miRanda.df$protein = NA


# ################# PARALLEL COMPUTING #################
# no_cores <- detectCores() / 2
# cl<-makeCluster(no_cores)
# registerDoParallel(cl)
# 
# # pb = progressBar(min = 0, max = nrow(mrnas.unique))
# 
# foreach(i = 1:nrow(mrnas.unique), .combine = rbind) %dopar% {
#   pos = grep(pattern = mrnas.unique[i, 1], x = miRanda.df$mRNA.ID)
#   miRanda.df$protein[pos] = mrnas.unique$protein[i]
#   miRanda.df
#   # setTxtProgressBar(pb, value = i)
# }
# # close(pb)
# stopCluster(cl)
# 
##### ASSIGN PROTEINS TO MIRANDA DF ####
pb = progressBar(min = 0, max = nrow(mrnas.unique))

for (i in 1:nrow(mrnas.unique)) {
  pos = grep(pattern = mrnas.unique[i, 1], x = miRanda.df$mRNA.ID)
  miRanda.df$protein[pos] = mrnas.unique$protein[i]
  setTxtProgressBar(pb, value = i)
}
close(pb)

# ##### COUNT SEEDS #####
# 
# proteins.table$seeds = NA
# 
# for (i in 1:nrow(proteins.table)) {
#   mrna.miranda.pos = grep(pattern = proteins.table$protlist[i], 
#                           x = miRanda.df$protein)
#   mrna.miranda.seeds = length(mrna.miranda.pos)
#   proteins.table$seeds[i] = mrna.miranda.seeds
#   print(mrna.miranda.seeds)
# }
# 
# summary(miRanda.table$mRNA.ID)
