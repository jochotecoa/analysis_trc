library('VennDiagram')

setwd('/share/analysis/hecatos/juantxo/Score/analysis/CPDB/gene_sets/')

geneontologyterm = 'GO:0044295    axonal growth cone    (CC level 4)'
axon_entrez_genes = read.table(
  file = geneontologyterm,
  header = F, sep = '\t')
colnames(axon_entrez_genes) = c('entrezgene', 'genename')

if (!exists('trt_proteomx')) {
  comp = '5FU'
  source('/share/script/hecatos/juantxo/analysis_trc/p.values_multiomics/checkSeveralOmicsTrTFile_JOA.R')
}

if (!any(grepl('entrez', colnames(trt_proteomx)))) {
  mart = openMart2018()
  
  biomaRt::listAttributes(mart) %>% .[, 'name'] %>% subset(., grepl('entrez', .))
  
  ens_entr = biomaRt::getBM(attributes = c("entrezgene", "ensembl_gene_id"), 
                            filters = "ensembl_gene_id", 
                            values = trt_proteomx$ensembl_gene_id, 
                            mart = mart)
  
  trt_proteomx = merge.data.frame(x = trt_proteomx, y = ens_entr, 
                                  by = 'ensembl_gene_id')
}


# axon_entrez_genes$user.provided.identifier = axon_degs$user.provided.identifier %>% 
#   as.character %>% 
#   gsub('     ', '', .)

axon_rows = trt_proteomx$entrezgene %in% axon_entrez_genes$entrezgene

axon_trt_prot = trt_proteomx[axon_rows, ]

axon_tpm = axon_trt_prot %>% 
  filter(p.adj_untVStox_TPM < 0.05) %>% 
  dplyr::select(ensembl_gene_id) %>% 
  unlist() %>% 
  as.character()

axon_trt = axon_trt_prot %>% 
  filter(p.adj_untVStox_TrT < 0.05) %>% 
  dplyr::select(ensembl_gene_id) %>% 
  unlist() %>% 
  as.character()

axon_prot = axon_trt_prot %>% 
  filter(p.value_untVStox_prx < 0.05) %>% 
  dplyr::select(ensembl_gene_id) %>% 
  unlist() %>% 
  as.character()




# # Chart
# venn.diagram(
#   x = list(set1, set2, set3),
#   category.names = c("Set 1" , "Set 2 " , "Set 3"),
#   filename = '#14_venn_diagramm.png',
#   output=TRUE
# )

setwd('vennDiagrams')

venn.diagram(
  x = list(axon_tpm, axon_trt, axon_prot),
  category.names = c("TPM" , "TRT" , "Proteomics"),
  filename = paste0(geneontologyterm, '.png', collapse = ''),
  output=TRUE
)