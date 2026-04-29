#' Correlation of TRC and Protein Expression with Time-Shifts
#' 
#' This script analyzes the correlation between TRC and proteomics data,
#' considering potential time-shifts (e.g., protein expression lagging behind mRNA).
#'
#' Inputs:
#' - Proteomics pre-processed files (defined in config.R)
#' - TRC score files (defined in config.R)
#' - BiomaRt for transcript-protein mapping
#'
#' Outputs:
#' - Correlation tables with optimal shifts (TSV)
#' - Distribution plots (PNG)

source("utils.R")
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
is.shift = TRUE
shift_range = 0:3
transcr_sample_shifted = 1:3 # First 3 samples of TRC are used for all shifts

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
  protwouttrans = NULL
  corlist = protlist = translist = shiftlist = NULL
  
  pb2 = txtProgressBar(min = 0, max = nrow(protvalcompsel), style = 3)
  for (protrow in 1:nrow(protvalcompsel)) {
    # Get the Protein ID
    protvalsel = protvalcompsel[protrow, -ncol(protvalcompsel)] 
    protnaam = rownames(protvalcompsel)[protrow]
    
    # Look for transcripts related to our protein
    transcr.rows = grep(protrans$uniprot_gn, pattern = protnaam)
    transel = protrans$ensembl_transcript_id[transcr.rows] 
    
    if (length(transel) == 0) { 
      protwouttrans = c(protwouttrans, protnaam)
    } else {
      transel = unlist(strsplit(transel, ';'))
      protTRCvalsel = rbind(protvalsel, rep(x = 0, ncol(protvalsel))) 
      
      for (t in transel) {
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
              protTRCvalsel[2, col.sample] = protTRCvalsel[2, col.sample] + TRCsel
            }
          }
        }
      }
      
      protTRCvalsel = protTRCvalsel[,-1]
      
      if (all(as.numeric(protTRCvalsel[2, transcr_sample_shifted]) != 0)) {
        tempcorlist = NULL
        for (sh in shift_range) {
          # Compare shifted protein samples with first mRNA samples
          # (e.g. Protein at 2h, 8h, 24h vs mRNA at 2h)
          # Note: The logic in v33 was slightly different, let's stick to v33's logic:
          # singcor = cor(x = as.numeric(protTRCvalsel[1,(transcr_sample_shifted + 3*sh)]), 
          #               y = as.numeric(protTRCvalsel[2,transcr_sample_shifted]))
          
          idx_prot = transcr_sample_shifted + 3*sh
          if (max(idx_prot) <= ncol(protTRCvalsel)) {
             singcor = cor(x = as.numeric(protTRCvalsel[1, idx_prot]), 
                           y = as.numeric(protTRCvalsel[2, transcr_sample_shifted]))
             tempcorlist = rbind(tempcorlist, c(singcor, sh))
          }
        }
        
        if (!is.null(tempcorlist)) {
          max_idx = which.max(tempcorlist[,1])
          corlist = c(corlist, tempcorlist[max_idx, 1])
          protlist = c(protlist, protnaam)
          translist = c(translist, transel[1])
          shiftlist = c(shiftlist, tempcorlist[max_idx, 2])
        }
      }
    }
    setTxtProgressBar(pb2, protrow)
  }
  close(pb2)
  
  # Results saving
  corlistnames = data.frame(translist, protlist, corlist, shiftlist)
  output_base = file.path('results', comp, paste0(predictive.values, '_shifted'))
  setOrCreatewd(output_base)
  
  write.table(x = corlistnames, 
              file = paste0('correlation_results_TRC_protein_shifted.tsv'), 
              sep = '\t', quote = F, row.names = F)
  
  png('frequency_distribution_shifts.png')
  barplot(table(shiftlist), main = 'Optimal Time-Shifts Distribution', 
          xlab = 'Shift', ylab = '# of proteins', col = 'orange')
  dev.off()
}
