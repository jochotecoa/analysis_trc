axon_degs = read.table('~/Documents/5FU_p.adj_untVStox_TPM/GO:0048675    axon extension    (BP level 5)', header = T, sep = '\t')

axon_degs$user.provided.identifier = axon_degs$user.provided.identifier %>% 
  as.character %>% 
  gsub('     ', '', .)

axon_rows = trt_proteomx$ensembl_gene_id %in% axon_degs$user.provided.identifier

axon_trt_prot = trt_proteomx[axon_rows, ]

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


forceLibrary('VennDiagram')

# # Chart
# venn.diagram(
#   x = list(set1, set2, set3),
#   category.names = c("Set 1" , "Set 2 " , "Set 3"),
#   filename = '#14_venn_diagramm.png',
#   output=TRUE
# )

venn.diagram(
  x = list(axon_trt_prot$ensembl_gene_id, axon_trt, axon_prot),
  category.names = c("TPM" , "TRT" , "Proteomics"),
  filename = 'axonextension_venn_diagramm.png',
  output=TRUE
)
