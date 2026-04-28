##### Personal functions ####

getSamples <- function(timepoints, triplicates, doses = '', control = F) {
  samples = c()
  if (control) {
    for (tp in timepoints) {
      for (tr in triplicates) {
        samples = c(samples, paste(tp, tr, sep = '_'))
      }
    }
  } else {
    for (d in doses) {
      for (tp in timepoints) {
        for (tr in triplicates) {
          samples = c(samples, paste(d, tp, tr, sep = '_'))
        }
      }
    }
  }
  return(samples)
}
get.or.assign <- function(var.name, filepath) {
  if (length(ls(pattern = var.name)) > 0) {
    file = get(var.name)
  } else {
    file = read.table(filepath)
    assign(x = var.name, value = file)
  }
  return(file)
}
setOrCreatewd <- function(dir.name) {
  if (dir.exists(dir.name)) {
    setwd(dir.name)
  } else {
    dir.create(dir.name)
    if (dir.exists(dir.name)) {
      setwd(dir.name)
    } else {
      stop(paste('Permission denied to create', dir.name, 'in', getwd()))
    }
    
  }
}
read.and.format.proteindata <- function(proteomicsfilepath, comp, samples) {
  proteinvaluesfile.colnames = read.table(proteomicsfilepath, nrows = 1, stringsAsFactors = F)
  proteinvaluesfile = read.table(proteomicsfilepath, skip = 1, stringsAsFactors = F, header = F, sep = ' ', fill = T)
  if (ncol(proteinvaluesfile) != ncol(proteinvaluesfile.colnames)) {
    proteinvaluesfile = read.table(proteomicsfilepath, skip = 1, stringsAsFactors = F, header = F, sep = '\t', fill = T)
  }
  proteinvaluesfile = proteinvaluesfile[grep('.*\\|.*\\|.*', proteinvaluesfile[,1]),]
  proteinvaluesfile = proteinvaluesfile[!grepl(pattern = ':', proteinvaluesfile$V2),]
  proteinvaluesfile[,2:ncol(proteinvaluesfile)] = sapply(proteinvaluesfile[,2:ncol(proteinvaluesfile)], as.double)
  proteinvaluesfile = cbind(proteinvaluesfile[,1], proteinvaluesfile[,!is.na(sapply(proteinvaluesfile, mean, na.rm = T))])
  proteinvaluesfile[,1] = as.character(proteinvaluesfile[,1])
  colnames(proteinvaluesfile) = proteinvaluesfile.colnames[,1:ncol(proteinvaluesfile)]
  protvalcomp = proteinvaluesfile[,c(1,grep(toupper(comp), colnames(proteinvaluesfile)))] # Get all columns with our compound
  selectedtimepointcolumns = c(1,grep(pattern=paste(samples, collapse='|'), x=colnames(protvalcomp))) # Select only columns with timepoints of interest
  protvalcompsel = protvalcomp[rowSums(!is.na(protvalcomp[,selectedtimepointcolumns])) > 2, selectedtimepointcolumns] # Filter proteins with =< 2 expressions
}
##### Input Variables #####
compound = 'Doxorubicin'
comp = 'Dox'
control = F
doses = c('The','Tox')
timepoints = c('002', '008', '024', '072')
triplicates = 1:3
MITdir = '/share/analysis/hecatos/juantxo/Score/Output_RunJUN2019_MV/'
organpath = '/ngs-data/data/hecatos/Cardiac'
predictive.values.list = c('TRC', 'targetRNA_TPM') # TRC / targetRNA_TPM
minimum_samples.list = c(3, 6, 9, 12)
is.shift.list = c(F,T)
shift = c(0,1,2,3)
transcr_sample_shifted = c(1,2,3)
timepoints.shifted = length(transcr_sample_shifted) / 3


# Generate sample names
samples = getSamples(timepoints = timepoints, triplicates = triplicates, control = control, doses = doses)


# Open protein expression file
proteinpath = paste(organpath,compound,'Protein', sep = '/')
proteindirs = list.dirs(path = proteinpath)
proteomicspath = subset(x = proteindirs, grepl(pattern = 'Proteomics', x = proteindirs))
proteomicsfiles = list.files(path = proteomicspath)
proteomicsfile = subset(x = proteomicsfiles, grepl(pattern = 'pre-processed_renamed', x = proteomicsfiles))
proteomicsfilepath = paste(proteomicspath, proteomicsfile, sep = '/')
protvalcompsel = read.and.format.proteindata(proteomicsfilepath = proteomicsfilepath, comp = comp, samples = samples)

# Open Protein-Transcript table
protrans = read.csv('/share/analysis/hecatos/juantxo/tableomics/uniprot-yourlist_M20190501.csv2', header = T, stringsAsFactors = F)

pb1 = txtProgressBar(min = 0, max = length(minimum_samples.list), style = 3)
time = timestamp(prefix = paste(comp, 'trc_protein_comparison_severaltranscriptsv27.R', sep = '_'))

for (is.shift in is.shift.list) {
  for (predictive.values in predictive.values.list) {
    for (minimum_samples in minimum_samples.list) {
      # Generate empty variables needed
      protwouttrans = NULL
      corlist = NULL
      protlist = NULL
      translist = NULL
      allprotTRCvalsel = data.frame(numeric())
      valuesnotfound = NULL
      protnotminexpr = NULL
      shiftlist = NULL
      naprotcorlist = NULL
      non_expressed_shifted_tp = NULL
      transcr_not_found_trcscorefile = NULL
      pb2 = txtProgressBar(min = 0, max = nrow(protvalcompsel), style = 3)
      for (protrow in 1:nrow(protvalcompsel)) {
        # Get the Protein ID
        protvalsel = protvalcompsel[protrow,]
        protnaam = strsplit(x = protvalsel[1,1], split = '\\|') # Only 1 ID
        protnaam = as.character(lapply(X = protnaam, FUN = '[', 2)) # Only UniProt ID
        
        
        # Look for transcripts related to our protein
        transel = protrans[grep(protrans$Entry, pattern = protnaam), 1] # Which transcripts are related to this UniProtID
        if (length(transel)==0) { 
          protwouttrans[length(protwouttrans)+1] = protnaam
        } else {
          transel = unlist(strsplit(transel, ';'))
          
          version = 0
          # Now that we have the transcript(s), we can get the data out of the score
          protTRCvalsel = rbind(protvalsel, rep(x = 0, ncol(protvalsel))) # 1 row for prot expression, and another row to be filled with TRC values
          
          for (t in transel) {
            # Get the TRC value for each sample
            for (s in samples) {
              setwd(MITdir)
              samplefilepath = paste0(comp,'_',s, '_TRCscore.txt') 
              if (file.exists(samplefilepath)) {
                if (length(grep(pattern = s, x = ls())) > 0) {
                  TRCscorefile = get(s)
                } else {
                  TRCscorefile = read.table(samplefilepath)
                  assign(x = s, value = TRCscorefile)
                }
                row.transcr.TRC = grep(x = rownames(TRCscorefile), pattern = t)
                if(length(row.transcr.TRC) > 0){
                  TRCsel = TRCscorefile[row.transcr.TRC,predictive.values] 
                  col.sample = grep(pattern = s, x = colnames(protTRCvalsel))
                  if (length(transel)>1) {
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
          numcolstransexpr = sum(protTRCvalsel[2,]!=0, na.rm=T)
          if (is.na(numcolsavail)) {
            numcolsavail = 0
          }
          if (numcolsavail >= minimum_samples) {
            protTRCvalsel2 = protTRCvalsel[,colsbothTRCprot]
            if (is.shift) {
              if (all(as.numeric(protTRCvalsel[2,transcr_sample_shifted])!=0)) {
                tempcorlist = NULL
                for (s in shift) {
                  singcor = cor(x = as.numeric(protTRCvalsel[1,(transcr_sample_shifted+3*s)]), y = as.numeric(protTRCvalsel[2,transcr_sample_shifted]))
                  tempcorlist = rbind(tempcorlist, c(singcor, s))
                }
                max.cor = max(tempcorlist[,1], na.rm = T)
                if (max.cor>=-1) {
                  corlist[length(corlist)+1] = max.cor
                  protlist[length(protlist)+1] = protnaam
                  translist[length(translist)+1] = t[1]
                  shiftlist[length(shiftlist)+1] = tempcorlist[grep(pattern = max.cor, x = tempcorlist[,1]), 2][1]
                } else {
                  naprotcorlist = c(naprotcorlist, max.cor)
                }
                
              } else {
                non_expressed_shifted_tp = c(non_expressed_shifted_tp, protnaam)
              }
            } else {
              singcor = cor(x = as.numeric(protTRCvalsel2[1,]), y = as.numeric(protTRCvalsel2[2,]))
              corlist[length(corlist)+1] = singcor
              protlist[length(protlist)+1] = protnaam
              translist[length(translist)+1] = t[1]
              
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
      if (length(protwouttrans)>0) {
        print(paste0('Warning: transcript not found for these proteins: ', paste(protwouttrans, collapse = ',')))
      }
      if (length(protnotminexpr)>0) {
        print(paste0('Warning: ', length(protnotminexpr) ,' proteins had < ', minimum_samples, ' samples expressed'))
      }
      if (nrow(protvalcompsel) != (length(corlist) + length(naprotcorlist) + length(protwouttrans) + length(protnotminexpr) + length(non_expressed_shifted_tp))) {
        print(paste0('Error: there are ', nrow(protvalcompsel) - (length(corlist) + length(naprotcorlist) + length(protwouttrans) + length(protnotminexpr) + length(non_expressed_shifted_tp)), ' proteins that have not been analyzed'))
      }
      
      corlistnames = cbind(translist, protlist, corlist, shiftlist)
      setwd('/share/analysis/hecatos/juantxo/score_protein_analysis/')
      setOrCreatewd(comp)
      setOrCreatewd(paste0('shift_', is.shift))
      setOrCreatewd(predictive.values)
      setOrCreatewd(time)
      corlistnames.filenaam = paste0(comp,'_correlation_results_',predictive.values,'_protein_1onseveral_until072_minimumexpressedsamples-', minimum_samples, '_shift-', is.shift, '_timeps-shifted-', timepoints[timepoints.shifted], '.tsv')
      write.table(x=corlistnames,file=corlistnames.filenaam, sep='\t', quote=F, row.names=F)
      write.table(x=allprotTRCvalsel,file=paste0(comp, '_dataframe_containing_', predictive.values, '_protein_1onseveral_until072_minimumexpressedsamples', minimum_samples,'_shift-', is.shift, '_timeps-shifted-', timepoints[timepoints.shifted], '.tsv'), sep='\t', quote=F, row.names=T)
      
      setwd('/share/analysis/hecatos/juantxo/score_protein_analysis/plots/')
      if (is.shift) {
        setwd('plots_correlations_nproteins_stringencies_optimalshift/')
      } else {
        setwd('plots_correlations_nproteins_stringencies/')
      }
      
      setOrCreatewd(comp)
      setOrCreatewd(predictive.values)
      setOrCreatewd(time)
      png(paste0('frequency_distribution_',comp,'_correlation_results_',predictive.values,'_protein_1onseveral_until072_minimumexpressedsamples-', minimum_samples,'_shift-', is.shift, '_timeps-shifted-', timepoints[timepoints.shifted], '.png'))
      barplot(table(cut(corlist, breaks = seq(from = -1, to = 1, by = 0.1))), las = 2, main = paste0('Correlation between ', predictive.values, ' and protein (Min.samples = ', minimum_samples, ')'), sub = paste0('Shift: max(', paste0(shift, collapse=','), ')'), ylim = c(0,120))
      dev.off()
      setTxtProgressBar(pb1, minimum_samples)
    }
    
  }
  
}

close(pb1)