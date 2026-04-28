#' Comparison of Transcript (TRC) and Protein Expression for Several Transcripts
#' 
#' This script analyzes the correlation between transcriptomics (TRC/TPM) and proteomics data.
#' It handles time-shifted correlations and minimum sample filtering.
#'
#' Inputs:
#' - Proteomics pre-processed files
#' - TRC score files
#' - BiomaRt for transcript-protein mapping
#'
#' Outputs:
#' - Correlation tables (TSV)
#' - Distribution plots (PNG)

source("../utils.R")
forceLibrary(c('pbmcapply', 'biomaRt'))

##### Input Variables #####

compound = 'Con_UNTR'
comp = 'UNTR'
control = T
doses = c('The','Tox')
timepoints = c('002', '008', '024', '072')
triplicates = 1:3
MITdir = '/share/analysis/hecatos/juantxo/Score/output/Output_Run_mrna_SEPT2019/V3/output/UNTR/2019-10-29_16:38:03_UTC/TRCscore/'
organpath = '/ngs-data/data/hecatos/Cardiac'
predictive.values.list = c('targetRNA_TPM', 'TRC') # TRC / targetRNA_TPM
minimum_samples.list = c(3, 6, 9, 12)
# F: no shift, T: shift, c(F,T): both
is.shift.list = c(F, T)
shift = c(0,1,2,3)
transcr_sample_shifted = 1:3
timepoints.shifted = length(transcr_sample_shifted) / 3

##### Analysis #####
# Save metadata
script.name <- basename(strsplit(commandArgs(trailingOnly = FALSE)[4],                                 "=")[[1]][2])
script.path = paste(getwd(), script.name, sep = '/')
metadata = ls()

# Generate sample names
samples = getSamples(timepoints = timepoints, 
                     triplicates = triplicates, 
                     control = control, 
                     doses = doses)


# Open protein expression file
proteinpath = paste(organpath, compound, 'Protein', sep = '/')
proteindirs = list.dirs(path = proteinpath)
proteomicspath = subset(x = proteindirs, 
                        grepl(pattern = 'Proteomics', x = proteindirs))
proteomicsfiles = list.files(path = proteomicspath)
proteomicsfile = subset(x = proteomicsfiles, 
                        grepl(pattern = 'pre-processed_renamed', 
                              x = proteomicsfiles))
proteomicsfilepath = paste(proteomicspath, proteomicsfile, sep = '/')
protvalcompsel = read.and.format.proteindata(proteomicsfilepath = proteomicsfilepath, 
                                             comp = comp, 
                                             samples = samples)
# Format protein IDs to only UniProtIDs 
protvalcompsel = protvalcompsel[!grepl(pattern = ':', protvalcompsel$Row.Names), ]
protnaams = strsplit(x = protvalcompsel$Row.Names, split = '\\|') # Split IDs
protnaams = as.character(lapply(X = protnaams, FUN = '[', 2)) # Only UniProt ID
protvalcompsel$Row.Names = protnaams

# Open Protein-Transcript table
mart.human = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                     dataset = 'hsapiens_gene_ensembl', 
                     host = 'http://apr2018.archive.ensembl.org') 

transcr_all.list = getBM(attributes = c('uniprot_gn', 
                                        'ensembl_transcript_id', 
                                        'transcript_biotype'),
                         filters = 'uniprot_gn',
                         values = protnaams, 
                         mart = mart.human)
prot_cod.rows = transcr_all.list$transcript_biotype == 'protein_coding'
protrans = transcr_all.list[prot_cod.rows, -3]


setwd('/share/analysis/hecatos/juantxo/tableomics/enst_uniprot_tables/')
protrans2 = read.csv(file = 'uniprot-yourlist_M20190501.csv2',
                    header = T,
                    stringsAsFactors = F)

pb1 = progressBar(min = 0, max = length(minimum_samples.list), style = "ETA")
# prefix = paste(comp, script.name, '', 
               # sep = '_')
time = timestamp()
time = gsub(pattern = '#', replacement = '', x = time)
time = gsub(pattern = '-', replacement = '', x = time)
time = gsub(pattern = ' ', replacement = '_', x = time)


for (is.shift in is.shift.list) {
  for (predictive.values in predictive.values.list) {
    for (minimum_samples in minimum_samples.list) {
      # Generate empty variables needed
      protwouttrans = valuesnotfound = protnotminexpr = naprotcorlist = NULL
      transcr_not_found_trcscorefile = non_expressed_shifted_tp = NULL
      allprotTRCvalsel = data.frame(numeric())
      corlist = protlist = translist = shiftlist = NULL
      
      pb2 = progressBar(min = 0, max = nrow(protvalcompsel), style = "ETA")
      for (protrow in 1:nrow(protvalcompsel)) {
        # Get the Protein ID
        protvalsel = protvalcompsel[protrow,]
        protnaam = protvalsel[1,1]
        
        # Look for transcripts related to our protein
        transcr.rows = grep(protrans$uniprot_gn, pattern = protnaam)
        transel = protrans$ensembl_transcript_id[transcr.rows] # Which transcripts are related to this UniProtID
        
        # View(searchAttributes(mart = mart.human))
                              
        if (length(transel) == 0) { 
          protwouttrans[length(protwouttrans) + 1] = protnaam
        } else {
          transel = unlist(strsplit(transel, ';'))
          
          version = 0
          # Now that we have the transcript(s), we can get the data out of the score
          protTRCvalsel = rbind(protvalsel, rep(x = 0, ncol(protvalsel))) # 1 row for prot expression, and another row to be filled with TRC values
          
          for (t in transel) {
            # Get the TRC value for each sample
            for (s in samples) {
              setwd(MITdir)
              files = list.files()
              samplefilepath = files[grep(pattern = s, x = files)]
              if (length(samplefilepath) == 0) {
                samplefilepath = paste0(comp,'_',s, '_TRCscore.txt')
              }
              if (file.exists(samplefilepath)) {
                if (length(grep(pattern = s, x = ls())) > 0) {
                  TRCscorefile = get(s)
                } else {
                  TRCscorefile = read.table(samplefilepath)
                  assign(x = s, value = TRCscorefile)
                }
                row.transcr.TRC = grep(x = rownames(TRCscorefile), pattern = t)
                if (length(row.transcr.TRC) > 0) {
                  TRCsel = TRCscorefile[row.transcr.TRC,predictive.values] 
                  col.sample = grep(pattern = s, x = colnames(protTRCvalsel))
                  if (length(transel) > 1) {
                    protTRCvalsel[2, col.sample] = protTRCvalsel[2, col.sample] + TRCsel
                    
                  } else {
                    protTRCvalsel[2, col.sample] = TRCsel
                  }
                } else {
                  transcr_not_found_trcscorefile = c(transcr_not_found_trcscorefile, t)
                }
              } else {
                print(paste0('Warning: ', samplefilepath, ' not found'))
              }
            }
          }
          protTRCvalsel[2,1] = paste(protTRCvalsel[1,1],'TRC values')
          rownames(protTRCvalsel) = protTRCvalsel[,1]
          protTRCvalsel = protTRCvalsel[,-1]
          colsbothTRCprot = !is.na(colSums(protTRCvalsel))
          numcolsavail = length(grep(T,colsbothTRCprot))
          numcolstransexpr = sum(protTRCvalsel[2,] != 0, na.rm = T)
          if (is.na(numcolsavail)) {
            numcolsavail = 0
          }
          if (numcolsavail >= minimum_samples) {
            protTRCvalsel2 = protTRCvalsel[,colsbothTRCprot]
            if (is.shift) {
              if (all(as.numeric(protTRCvalsel[2,transcr_sample_shifted]) != 0)) {
                tempcorlist = NULL
                for (s in shift) {
                  singcor = cor(x = as.numeric(protTRCvalsel[1,(transcr_sample_shifted + 3*s)]), 
                                y = as.numeric(protTRCvalsel[2,transcr_sample_shifted]))
                  tempcorlist = rbind(tempcorlist, c(singcor, s))
                }
                max.cor = max(tempcorlist[,1], na.rm = T)
                if (max.cor >= -1) {
                  corlist[length(corlist) + 1] = max.cor
                  protlist[length(protlist) + 1] = protnaam
                  translist[length(translist) + 1] = t[1]
                  shiftlist[length(shiftlist) + 1] = tempcorlist[grep(pattern = max.cor, 
                                                                      x = tempcorlist[,1]), 2][1]
                } else {
                  naprotcorlist = c(naprotcorlist, max.cor)
                }
                
              } else {
                non_expressed_shifted_tp = c(non_expressed_shifted_tp, protnaam)
              }
            } else {
              singcor = cor(x = as.numeric(protTRCvalsel2[1,]), 
                            y = as.numeric(protTRCvalsel2[2,]))
              corlist[length(corlist) + 1] = singcor
              protlist[length(protlist) + 1] = protnaam
              translist[length(translist) + 1] = t[1]
              
            }
          } else {
            protnotminexpr = c(protnotminexpr, protnaam)
          }
          allprotTRCvalsel = rbind(allprotTRCvalsel, protTRCvalsel)
        }
        
        # print(paste0('Progress: ', protrow/nrow(protvalcompsel)*100,'%'))
        setTxtProgressBar(pb2, protrow)
        # print(rep)
      }
      
      close(pb2)
      # Warnings
      if (length(protwouttrans)) {
        print(paste0('Warning: there are ', 
                     length(protwouttrans), 
                     ' proteins for which a transcript has not been found'))
      }
      if (length(protnotminexpr)) {
        print(paste0('Warning: ', 
                     length(protnotminexpr) ,
                     ' proteins had < ', 
                     minimum_samples, 
                     ' samples expressed'))
      }
      if (nrow(protvalcompsel) != (length(corlist) + 
                                   length(naprotcorlist) + 
                                   length(protwouttrans) + 
                                   length(protnotminexpr) + 
                                   length(non_expressed_shifted_tp))) {
        print(paste0('Error: there are ', 
                     nrow(protvalcompsel) - 
                       (length(corlist) + 
                          length(naprotcorlist) + 
                          length(protwouttrans) + 
                          length(protnotminexpr) + 
                          length(non_expressed_shifted_tp)), 
                     ' proteins that have not been analyzed'))
      }
      if (length(transcr_not_found_trcscorefile)) {
        print(paste0('Warning: ', length(transcr_not_found_trcscorefile), 
                     ' transcripts not found'))
      }
      
      corlistnames = cbind(translist, protlist, corlist, shiftlist)
      setwd(MITdir)
      setOrCreatewd('Analysis')
      setOrCreatewd(comp)
      setOrCreatewd(paste0('shift_', is.shift))
      setOrCreatewd(predictive.values)
      setOrCreatewd(time)
      setOrCreatewd(paste0('minimum_expressed_samples', minimum_samples))
      corlistnames.filenaam = paste0('correlation_results_between_',
                                     predictive.values,
                                     '_values_and_protein_expression.tsv')
      saveMetadata(x = metadata)
      write.table(x = corlistnames, 
                  file = corlistnames.filenaam, 
                  sep = '\t', 
                  quote = F, 
                  row.names = F)
      write.table(x = allprotTRCvalsel,
                  file = paste0(predictive.values, 
                                '_values_and_protein_expression.tsv'), 
                  sep = '\t', quote = F, row.names = T)
      
      setOrCreatewd('plots/')
      if (is.shift) {
        setOrCreatewd('plots_correlations_nproteins_stringencies_optimalshift/')
      } else {
        setOrCreatewd('plots_correlations_nproteins_stringencies/')
      }
      png(paste0('frequency_distribution_correlation_results_between',
                 predictive.values,
                 '_and_protein.png'))
      if (is.shift) {
        sub = paste0('Shift: max(', paste0(shift, collapse = ','), ')')
      } else {
        sub = ''
      }
      x = table(cut(corlist, 
                    breaks = seq(from = -1, 
                                 to = 1, 
                                 by = 0.1)))
      x = (x / length(corlist)) * 100
      barplot(x, 
              las = 2, 
              main = paste0('Correlation between ', 
                            predictive.values, 
                            ' and protein (Min.samples = ', 
                            minimum_samples, ')'), 
              sub = sub, 
              ylim = c(0,length(corlist)), 
              xlab = 'Correlation values', ylab = '% of proteins')
      dev.off()
      setTxtProgressBar(pb1, minimum_samples)
    }
    
  }
  
}

close(pb1)
