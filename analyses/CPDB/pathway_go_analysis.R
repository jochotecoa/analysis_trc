
library(dplyr)
source('/share/script/hecatos/juantxo/analysis_trc/functions_JOA.R')
forceLibrary('UpSetR')

Subset <- function(x, pattern) {
  
  if (is.null(dim(x))) {
    ind <- grepl(pattern, x)
    return(x[ind])
  } else {
    ind <- grepl(pattern, names(x))
    return(x[, ind])
  }
}

setwd('~/Downloads/CPDB/')

compounds = list.dirs(recursive=F)

for (compound in compounds) {
  setwd(compound)
  comparisons = list.dirs(recursive=F)
  for (comparison in comparisons) {
    setwd(comparison)
    quantifiers = list.dirs(recursive=F)
    for (quantifier in quantifiers) {
      setwd(quantifier)
      analyses = list.dirs(recursive=F)
      for (analysis in analyses) {
        setwd(analysis)
        condit = list.files() %>% length() > 0
        if (condit) {
          temp_file = read.table('ORA_results.tab', sep='\t', header=T)
          name = paste('df', compound, comparison, quantifier, analysis, 
                       sep = '_') %>% 
            gsub('./', '', .) %>% 
            gsub('-', '', .)
          assign(name,temp_file)
        } else {
          paste(getwd(), 'has no result file') %>% warning()
        }
        setwd('../')
      }
      setwd('../')
    }
    setwd('../')
  }
  setwd('../')
}
 
path_dataframes = ls(pattern = 'df_') %>% 
  Subset(pattern = 'pathway')

go_dataframes = ls(pattern = 'df_') %>% 
  Subset(pattern = 'geneonthology')

# trt_dfs_names = dataframes %>% 
#   Subset(pattern = 'TRT') %>% 
#   Subset(pattern = 'pathway')
path_dfs = get(path_dataframes[1])
for (variable in path_dataframes[-1]) {
  temp_df = get(variable)
  path_dfs = rbind.data.frame(path_dfs, temp_df)
}

all_pathways = path_dfs$pathway %>% unique()

path_binary = data.frame(row.names = all_pathways)

for (variable in path_dataframes) {
  temp_paths = get(variable) %>% .[, 'pathway'] %>% as.character()
  temp_bi = rownames(path_binary) %in% temp_paths
  path_binary = cbind.data.frame(path_binary, temp_bi)
  colnames(path_binary)[ncol(path_binary)] = variable
}

path_binary = path_binary %>% lapply(as.integer) %>% as.data.frame
rownames(path_binary) = all_pathways
colnames(path_binary) = colnames(path_binary) %>% 
  gsub(pattern = '_pathway_analysis', 
       replacement = '')



upset(path_binary, 
      order.by = "freq", 
      nsets = 6, 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 3, 
      point.size = 2.8, 
      line.size = 1
)



go_dataframes = ls(pattern = 'df_') %>% 
  Subset(pattern = 'geneonthology')

go_dfs = get(go_dataframes[1])
for (variable in go_dataframes[-1]) {
  temp_df = get(variable)
  go_dfs = rbind.data.frame(go_dfs, temp_df)
}

all_term_names = go_dfs$term_name %>% unique()

go_binary = data.frame(row.names = all_term_names)

for (variable in go_dataframes) {
  temp_gos = get(variable) %>% .[, 'term_name'] %>% as.character()
  temp_bi = rownames(go_binary) %in% temp_gos
  go_binary = cbind.data.frame(go_binary, temp_bi)
  colnames(go_binary)[ncol(go_binary)] = variable
}

go_binary = go_binary %>% lapply(as.integer) %>% as.data.frame
rownames(go_binary) = all_term_names
colnames(go_binary) = colnames(go_binary) %>% 
  gsub(pattern = '_geneonthology_analysis', 
       replacement = '')



upset(go_binary, 
      order.by = "freq", 
      nsets = 6, 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 3, 
      point.size = 2.8, 
      line.size = 1
)



# tpm_dfs_names = dataframes %>% 
#   Subset(pattern = 'TPM') %>% 
#   Subset(pattern = 'pathway')
# tpm_dfs = get(tpm_dfs_names[1])
# for (variable in tpm_dfs_names[-1]) {
#   temp_df = get(variable)
#   tpm_dfs = merge.data.frame(tpm_dfs, temp_df, 'pathway', all = T)
# }
# 
# setdiff(unique(trt_dfs$pathway), unique(tpm_dfs$pathway))
# 
