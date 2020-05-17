# Functions & libraries ---------------------------------------------------


viewPlotsProtx = function(df, evalu, decr = F){
  if (any(grepl(pattern = 'rowname', colnames(df)))) {
    df = df %>% remove_rownames() %>% column_to_rownames()
    colnames(df)[1] %>% print()
  }
  row_data = order(as.numeric(as.character(df[, evalu])), decreasing = decr)
  for (i in row_data) {
    data = df[i, ]
    plot_data = grep('Proteomics_', colnames(df)) %>% data[.] %>% as.matrix()
    
    plot_data %>% barplot(las = 2, 
                          main = paste0(rownames(data), ' ; ',
                                        evalu, ' = ', 
                                        format.pval(data[evalu], 3)), 
                          ylab = 'Protein Expression')
    
    
    the_data = plot_data[grep('The', colnames(plot_data))]
    tox_data = plot_data[grep('Tox', colnames(plot_data))]
    mean_the_data = mean(the_data, na.rm = T)
    mean_tox_data = mean(tox_data, na.rm = T)
    sd_the_data = sd(the_data, na.rm = T)
    sd_tox_data = sd(tox_data, na.rm = T)
    thesh_up_the = mean_the_data + sd_the_data
    thesh_down_the = mean_the_data - sd_the_data
    thesh_up_tox = mean_tox_data + sd_tox_data
    thesh_down_tox = mean_tox_data - sd_tox_data
    
    abline(h = thesh_up_the, col = 'lightgreen', pch = 2)
    abline(h = mean(the_data, na.rm = T), col = 'cyan')
    abline(h = thesh_down_the, col = 'lightgreen')
    
    abline(h = thesh_up_tox, col = 'indianred')
    abline(h = mean(tox_data, na.rm = T), col = 'purple')
    abline(h = thesh_down_tox, col = 'indianred')
    
    abline(v = 14.5)
    for (l in seq(3.7, 85.1, 3.6)) abline(v = l, col = 'gray')
    
    readline(prompt = "Press [enter] to continue")
  }
  
}

shift_median = function(df, median_of_medians) {
	for (i in colnames(df)) {
		sample = df[, i]
		corr_factor = median_of_medians - median(sample, na.rm = T)
		df[, i] = df[, i] + corr_factor
	}
	return(df)
}

normalizeProteomics = function(df) {
	medians = apply(df, 2, median, na.rm = T)
	median_of_medians = median(medians, na.rm = T)
	df = shift_median(df, median_of_medians)
	return(df)
}

setwd("/share/script/hecatos/juantxo/analysis_trc")
source('functions_JOA.R')
forceLibrary(c('tibble', 'dplyr'))

setwd('/share/analysis/hecatos/juantxo/Score/input/RNA_Salmon/')
biomart_table = readRDS('salmon_ids_biomart.rds')
colnames(biomart_table) = biomart_table[1, ] %>% unlist() %>% 
  droplevels() %>% as.character()
biomart_table = biomart_table[-1, ] 
biomart_table = biomart_table %>% filter(`Transcript type` == 'protein_coding')


### Renaming data #####
setwd('/ngs-data/data/hecatos/Cardiac/')

lett_comp = strsplit(comp, '') %>% unlist()
lett_comp[-1] = tolower(lett_comp[-1])
comps_sel = list.dirs(recursive = F)
for (letter in lett_comp) {
  comps_sel = grep(letter, comps_sel) %>% comps_sel[.]
}
if (length(comps_sel) == 0) {
  lett_comp = strsplit(comp, '') %>% unlist()
  comps_sel = list.dirs(recursive = F)
  for (letter in lett_comp) {
    comps_sel = grep(letter, comps_sel) %>% comps_sel[.]
  }
}
if (length(comps_sel) > 1) {
	for (comp_sel in comps_sel) {
		comp_sel_ori = comp_sel
		comp_sel = comp_sel %>% strsplit('') %>% unlist() %>% .[3:5]
		iscomp = NULL
		for (letter in lett_comp) {
			iscomp = comp_sel %>% grepl(letter, .) %>% sum %>% as.logical
			if (!iscomp) {break}
		}
		if(iscomp) {
			comps_sel = comp_sel_ori
			break
		}
	}
}
print(comps_sel)
setwd(comps_sel)
setwd('Protein')
prot_dir = list.files(pattern = 'Proteomics') 
stopifnot(length(prot_dir) == 1)
setwd(prot_dir)
suppressWarnings(rm(proteomx_log2, proteomx_log2_norm))
if (length(list.files(pattern = 'renamed-by-Juan.txt')) == 0) {
  proteomx_log2 = proteomx_log2_norm = 
    read.table('Results_Hecatos_Cardiac_Px_PTX_FC/Cardiac_PTX_log2_norm.txt', 
               sep = '\t', header = T) %>% 
    cleanProtIds() %>% 
    remove_rownames() %>% 
    column_to_rownames('uniprot_gn') %>% 
    dplyr::select(-Row.Names)
  # dplyr::select(-X, -contains('.1')) %>% 
  # dplyr::select(contains('PTX')) 
  if (plotting) {
    proteomx_log2 %>% naToZero() %>% as.matrix.data.frame() %>% heatmap()
  }
  
  # proteomx_log2 %>% dist() %>% naToZero() %>% hclust() %>% plot()
  
  colnames(proteomx_log2) %>% 
    gsub(pattern = '_rep', replacement = '') %>% strsplit('_') %>% 
    as.data.frame() %>% t() %>% as.data.frame() %>% .[, 'V4'] %>% sort() %>% head()
  
  sample_order = colnames(proteomx_log2) %>% 
    gsub(pattern = '_rep', replacement = '') %>% strsplit('_') %>% 
    as.data.frame() %>% t() %>% as.data.frame() %>% .[, 'V4'] %>% order()
  
  
  proteomx_log2 = proteomx_log2[, sample_order]
  
  the_samples = data.frame(dose = rep('PTX_The', 24))
  the_samples$timepoint = 
    lapply(c('000', '002', '008', '024', '072', '168', '240', '336'), rep, 3) %>% 
    unlist()
  the_samples$replicate = rep(c('1', '2', '3'), 8)
  
  tox_samples = data.frame(dose = rep('PTX_Tox', 18))
  tox_samples$timepoint = 
    lapply(c('002', '008', '024', '072', '168', '240'), rep, 3) %>% 
    unlist()
  tox_samples$replicate = rep(c('1', '2', '3'), 6)
  
  sample_names = rbind.data.frame(the_samples, tox_samples) %>% 
    apply(1, paste, collapse = '_')
  
  colnames(proteomx_log2) = sample_names
  
  if (plotting) {
    proteomx_log2 %>% naToZero() %>% as.matrix.data.frame() %>% heatmap()
  }
  
  try(write.table(x = proteomx_log2, 
              file = 
                'Hecatos_Cardiac_Px_PTX_Ther_Tox_log2_norm-renamed-by-Juan.txt'))
  
} else {
  proteomx_log2 = 
    read.table(list.files(pattern = 'renamed-by-Juan.txt'))
}

##### Loading Control data #####

setwd('/ngs-data/data/hecatos/Cardiac/Con_UNTR/Protein/Proteomics_Analyses_Cardiac_UNTR_GeneData/')
proteomx_untr = read.table('Hecatos_Cardiac_Px_Untreated_pre-processed_renamed.txt', header = T, sep = '\t') %>%
	cleanProtIds %>%
	remove_rownames() %>% 
	column_to_rownames('uniprot_gn') %>% 
	dplyr::select(-Row.Names) %>%
	filterSamplesBySeqDepth %>%
  log2() %>%
	normalizeProteomics %>%
	dplyr::select(matches('002|008|024|072'))

##### Analysis #####
colnames(proteomx_log2) = colnames(proteomx_log2) %>% gsub(pattern = 'XPTX', replacement = 'PTX')

prot_log_norm = filterSamplesBySeqDepth(proteomx_log2)
if (plotting) {
  proteomx_log2 %>% apply(2, median, na.rm = T) %>% 
    barplot(las = 2, main = 'Median Proteomics Expression per Sample')
}
### Get differentially exp %>% sed proteins

the_cols = prot_log_norm %>% dplyr::select(contains('The')) %>% 
  dplyr::select(matches('002|008|024|072')) %>% colnames()
tox_cols = prot_log_norm %>% dplyr::select(contains('Tox')) %>% 
  dplyr::select(matches('002|008|024|072')) %>% colnames()

# ii) for every treatment/dose combination and for each time-point 
# the common protein set between the controls and the treatment samples was determined

proteomx_untVSthe = proteomx_log2[, the_cols] %>% 
	rownames_to_column %>% 
	merge.data.frame(rownames_to_column(proteomx_untr), by = 'rowname') %>% 
	column_to_rownames

proteomx_untVStox = proteomx_log2[, tox_cols] %>% 
	rownames_to_column %>% 
	merge.data.frame(rownames_to_column(proteomx_untr), by = 'rowname') %>% 
	column_to_rownames

# iii)
# the median of the medians of the (in general 3) normalized control samples was determined
# using this common protein set between the controls and the treatment samples

proteomx_untVSthe[, 13:24] = normalizeProteomics(proteomx_untVSthe[, 13:24])
proteomx_untVStox[, 13:24] = normalizeProteomics(proteomx_untVStox[, 13:24])

mom_unthe = median(proteomx_untVSthe[, 13], na.rm = T)
mom_untox = median(proteomx_untVStox[, 13], na.rm = T)

# iv) the
# data from the samples of the treatments were shifted to these medians.

proteomx_untVSthe[, 1:12] = shift_median(proteomx_untVSthe[, 1:12], mom_unthe)
proteomx_untVStox[, 1:12] = shift_median(proteomx_untVStox[, 1:12], mom_untox)


rm(theVStox_t.test)
theVStox_t.test = prot_log_norm %>% apply_2D(FUN = t.test, col.x = the_cols, 
                  col.y = tox_cols, paired = T, complete_cases = 12) 
theVStox_t.test$p.value_theVStox_prx = theVStox_t.test$p.value %>% as.character() %>% as.numeric()
theVStox_t.test$statistic.t_theVStox_prx = theVStox_t.test$statistic.t %>% as.character() %>% as.numeric()
if (plotting) {theVStox_t.test$p.value %>% na.omit() %>% hist(breaks = seq(0, 1, 0.05))}

rm(untVSthe_t.test)
untVSthe_t.test = prot_log_norm %>% apply_2D(FUN = t.test, col.x = 13:24, 
                  col.y = the_cols, paired = T, complete_cases = 12)
untVSthe_t.test$p.value_untVSthe_prx = untVSthe_t.test$p.value %>% as.character() %>% as.numeric()
untVSthe_t.test$statistic.t_untVSthe_prx = untVSthe_t.test$statistic.t %>% as.character() %>% as.numeric()
if (plotting) {untVSthe_t.test$p.value %>% na.omit() %>% hist(breaks = seq(0, 1, 0.05))}

rm(untVStox_t.test)
untVStox_t.test = prot_log_norm %>% apply_2D(FUN = t.test, col.x = 13:24, 
                  col.y = tox_cols, paired = T, complete_cases = 12)
untVStox_t.test$p.value_untVStox_prx = untVStox_t.test$p.value %>% as.character() %>% as.numeric()
untVStox_t.test$statistic.t_untVStox_prx = untVStox_t.test$statistic.t %>% as.character() %>% as.numeric()
if (plotting) {untVStox_t.test$p.value %>% na.omit() %>% hist(breaks = seq(0, 1, 0.05))}


### Bind the p-values to the proteomics data 
proteomx = 2^prot_log_norm %>% 
  dplyr::select(matches('002|008|024|072'))
colnames(proteomx) = paste0('Proteomics_', colnames(proteomx))

proteomx_pval = merge.data.frame(x = rownames_to_column(proteomx), 
                                y = rownames_to_column(dplyr::select(theVStox_t.test, matches('p.value_|statistic.t_'))), 
                                by = 'rowname', all = T) %>% 
                	merge.data.frame(y = rownames_to_column(dplyr::select(untVSthe_t.test, 
                				matches('p.value_|statistic.t_'))), 
                                by = 'rowname', all = T) %>% 
                        merge.data.frame(y = rownames_to_column(dplyr::select(untVStox_t.test, 
                				matches('p.value_|statistic.t_'))), 
                                by = 'rowname', all = T) %>% 
  			column_to_rownames() 

# proteomx_pval$p.adj %>% as.numeric() %>% na.omit() %>% hist(breaks = seq(0, 1, 0.05))
p.val_vars = colnames(proteomx_pval)%>% .[grep('p.value', .)]

mart = openMart2018()

biomart_table = biomaRt::getBM(attributes = c('transcript_biotype', 
                                  'uniprotswissprot', 
                                  "ensembl_gene_id", 
                                  "external_gene_name",
                                  "ensembl_transcript_id"), 
               filters = 'uniprotswissprot', 
               values = rownames(proteomx_pval), mart = mart) %>% 
  dplyr::filter(transcript_biotype == 'protein_coding')

proteomx_biom = merge.data.frame(x = rownames_to_column(proteomx_pval), 
                                 y = biomart_table, by.x = 'rowname', 
                                 by.y = 'uniprotswissprot') %>% 
  rename(uniprotswissprot = rowname) %>% 
  rownames_to_column

proteomx_sign = proteomx_biom %>% try(filter(p.value < 0.05))


# setwd('/share/analysis/hecatos/juantxo/Score/output/DOC/2020-03-09_13:09:20_UTC/TRCscore/factors_DOC/')
#saveRDS(proteomx_biom, file = 'proteomx_biom.rds')
# viewPlotsProtx(filter(proteomx_pval, p.adj < 0.05), evalu = 'p.adj', T)
