#' Correlation of TRC and Protein Expression
#' 
#' This script analyzes the correlation between transcriptomics (TRC) 
#' and proteomics data.
#'
#' Inputs:
#' - Proteomics pre-processed files (defined in config.R)
#' - TRC score files (defined in config.R)
#' - BiomaRt for transcript-protein mapping
#'
#' Outputs:
#' - Correlation tables (TSV)
#' - Distribution plots (PNG)

if (file.exists("utils.R")) { source("utils.R") } else if (file.exists("../utils.R")) { source("../utils.R") }
forceLibrary(c('pbmcapply', 'biomaRt', 'dplyr', 'tibble'))

##### Input Variables #####

compound = if(exists("DEFAULT_COMPOUND")) DEFAULT_COMPOUND else 'Con_UNTR'
comp = if(exists("DEFAULT_COMP")) DEFAULT_COMP else 'UNTR'
control = T
doses = c('The','Tox')
timepoints = if(exists("DEFAULT_TIMEPOINTS")) DEFAULT_TIMEPOINTS else c('002', '008', '024', '072')
triplicates = if(exists("DEFAULT_TRIPLICATES")) DEFAULT_TRIPLICATES else 1:3

# Construct paths using config variables if available
if (exists("SCORE_OUTPUT_PATH")) {
  MITdir = getLatestTimestampDir(file.path(SCORE_OUTPUT_PATH, "V3/output/UNTR"), pattern = "TRCscore")
  if (is.null(MITdir)) {
    # Fallback to specific known timestamp if latest not found or structure differs
    MITdir = file.path(SCORE_OUTPUT_PATH, 'V3/output/UNTR/2019-10-29_16:38:03_UTC/TRCscore/')
  }
} else {
  MITdir = '/share/analysis/hecatos/juantxo/Score/output/Output_Run_mrna_SEPT2019/V3/output/UNTR/2019-10-29_16:38:03_UTC/TRCscore/'
}

organpath = if(exists("PROTEOMICS_BASE_PATH")) PROTEOMICS_BASE_PATH else '/ngs-data/data/hecatos/Cardiac'

predictive.values = 'TRC' 
minimum_samples = 12
is.shift = FALSE
shift = 0
transcr_sample_shifted = 1:3

##### Analysis #####

# Generate sample names
samples = getSamples(timepoints = timepoints, 
                     triplicates = triplicates, 
                     control = control, 
                     doses = doses)


# Open protein expression file
proteinpath = file.path(organpath, compound, 'Protein')
proteindirs = list.dirs(path = proteinpath)
proteomicspath = subset(x = proteindirs, 
                        grepl(pattern = 'Proteomics', x = proteindirs))
proteomicsfiles = list.files(path = proteomicspath)
proteomicsfile = subset(x = proteomicsfiles, 
                        grepl(pattern = 'pre-processed_renamed', 
                              x = proteomicsfiles))
proteomicsfilepath = file.path(proteomicspath, proteomicsfile)

if (file.exists(proteomicsfilepath)) {
  protvalcompsel = read.and.format.proteindata(proteomicsfilepath = proteomicsfilepath, 
                                               comp = comp, 
                                               samples = samples)

  # Format protein IDs to only UniProtIDs using utility function
  protvalcompsel = cleanProtIds(protvalcompsel)
  protnaams = protvalcompsel$uniprot_gn
  rownames(protvalcompsel) = protnaams

  # Open Protein-Transcript table
  mart.human = openMart2018()

  transcr_all.list = getBM(attributes = c('uniprot_gn', 
                                          'ensembl_transcript_id', 
                                          'transcript_biotype'),
                           filters = 'uniprot_gn',
                           values = protnaams, 
                           mart = mart.human)
  prot_cod.rows = transcr_all.list$transcript_biotype == 'protein_coding'
  protrans = transcr_all.list[prot_cod.rows, -3]

  # Cache TRC files to avoid repeated setwd and listing
  if (dir.exists(MITdir)) {
    trc_files = list.files(MITdir)
  } else {
    trc_files = character()
  }

  # Generate empty variables needed
  protwouttrans = valuesnotfound = protnotminexpr = naprotcorlist = NULL
  transcr_not_found_trcscorefile = non_expressed_shifted_tp = NULL
  allprotTRCvalsel = data.frame()
  corlist = protlist = translist = NULL
  
  pb2 = txtProgressBar(min = 0, max = nrow(protvalcompsel), style = 3)
  for (protrow in 1:nrow(protvalcompsel)) {
    # Get the Protein ID
    protvalsel = protvalcompsel[protrow, -ncol(protvalcompsel)] # exclude uniprot_gn column
    protnaam = rownames(protvalcompsel)[protrow]
    
    # Look for transcripts related to our protein
    transcr.rows = grep(protrans$uniprot_gn, pattern = protnaam)
    transel = protrans$ensembl_transcript_id[transcr.rows] # Which transcripts are related to this UniProtID
    
    if (length(transel) == 0) { 
      protwouttrans[length(protwouttrans) + 1] = protnaam
    } else {
      transel = unlist(strsplit(transel, ';'))
      
      # Now that we have the transcript(s), we can get the data out of the score
      protTRCvalsel = rbind(protvalsel, rep(x = 0, ncol(protvalsel))) 
      
      for (t in transel) {
        # Get the TRC value for each sample
        for (s in samples) {
          samplefilepath = trc_files[grep(pattern = s, x = trc_files)]
          if (length(samplefilepath) == 0) {
            samplefilepath = paste0(comp,'_',s, '_TRCscore.txt')
          }
          
          full_sample_path = file.path(MITdir, samplefilepath)
          
          if (file.exists(full_sample_path)) {
            if (exists(s, envir = .GlobalEnv)) {
              TRCscorefile = get(s, envir = .GlobalEnv)
            } else {
              TRCscorefile = read.table(full_sample_path)
              assign(x = s, value = TRCscorefile, envir = .GlobalEnv)
            }
            row.transcr.TRC = grep(x = rownames(TRCscorefile), pattern = t)
            if (length(row.transcr.TRC) > 0) {
              TRCsel = TRCscorefile[row.transcr.TRC, predictive.values] 
              col.sample = grep(pattern = s, x = colnames(protTRCvalsel))
              if (length(transel) > 1) {
                protTRCvalsel[2, col.sample] = protTRCvalsel[2, col.sample] + TRCsel
              } else {
                protTRCvalsel[2, col.sample] = TRCsel
              }
            } else {
              transcr_not_found_trcscorefile = c(transcr_not_found_trcscorefile, t)
            }
          }
        }
      }
      
      protTRCvalsel[2,1] = paste(protTRCvalsel[1,1],'TRC values')
      rownames(protTRCvalsel) = protTRCvalsel[,1]
      protTRCvalsel = protTRCvalsel[,-1]
      colsbothTRCprot = !is.na(colSums(protTRCvalsel))
      numcolsavail = sum(colsbothTRCprot)
      
      if (numcolsavail >= minimum_samples) {
        protTRCvalsel2 = protTRCvalsel[,colsbothTRCprot]
        singcor = cor(x = as.numeric(protTRCvalsel2[1,]), 
                      y = as.numeric(protTRCvalsel2[2,]))
        corlist[length(corlist) + 1] = singcor
        protlist[length(protlist) + 1] = protnaam
        translist[length(translist) + 1] = transel[1]
      } else {
        protnotminexpr = c(protnotminexpr, protnaam)
      }
      allprotTRCvalsel = rbind(allprotTRCvalsel, protTRCvalsel)
    }
    setTxtProgressBar(pb2, protrow)
  }
  close(pb2)
  
  # Results saving
  corlistnames = data.frame(translist, protlist, corlist)
  output_base = file.path('results', comp, predictive.values)
  setOrCreatewd(output_base)
  
  write.table(x = corlistnames, 
              file = paste0('correlation_results_between_', predictive.values, '_and_protein.tsv'), 
              sep = '\t', quote = F, row.names = F)
  
  png(paste0('frequency_distribution_', predictive.values, '.png'))
  hist(corlist, main = paste0('Correlation ', predictive.values, ' vs Protein'), 
       xlab = 'Correlation', col = 'steelblue')
  dev.off()
} else {
  warning("Proteomics file not found: ", proteomicsfilepath)
}
