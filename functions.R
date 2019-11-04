cleanProtIds = function(protein_table) {
  protein_table = protein_table[!grepl(protein_table[, 1], pattern = ':'), ]
  names = strsplit(as.character(protein_table[, 1]), '\\|')
  names = as.character(lapply(names, '[', 2))
  protein_table$uniprot_gn = names
  return(protein_table)
}

naToZero = function(x) {
  x[is.na(x)] = 0
  return(x)
}

transcrToGene = function(table, aggregate = F) {
  library('biomaRt')
  table[, 'rownames'] = rownames(table)
  sampl = table[nrow(table), ]
  enst_col = grep(pattern = 'ENST', x = sampl)[1]
  
  version = grepl('\\.', sampl[, enst_col])
  if (length(version) == 0) {version = F}
  if (version) {
    transcript_id = 'ensembl_transcript_id_version'
  } else {
    transcript_id = 'ensembl_transcript_id'
  }
  
  values = table[, enst_col]
  mart.human = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                       dataset = 'hsapiens_gene_ensembl',
                       host = 'http://apr2018.archive.ensembl.org') 
  
  new_cols = getBM(attributes = c(transcript_id, 'ensembl_gene_id'), 
                   filters = transcript_id, values = values, mart = mart.human)
  
  table = merge.data.frame(x = table, y = new_cols, 
                           by.x = colnames(table)[enst_col], 
                           by.y = transcript_id)
  if (aggregate) {
    int_cols = grepl('integer', sapply(X = table[1, ], FUN = typeof))
    int_cols = int_cols + grepl('double', sapply(X = table[1, ], 
                                                 FUN = typeof))
    int_cols = as.logical(int_cols)
    table = aggregate(x = table[, int_cols], by = list(table$ensembl_gene_id), 
                      FUN = sum, na.rm = T)
    colnames(table)[1] = 'ensembl_gene_id'
  }
  return(table)
}

rmMirnas = function(x) {
  mirna.cols = grep(pattern = 'hsa', x = colnames(x))
  y = x[, -mirna.cols]
  return(y)
}
forceLibrary <- function(list.of.packages) {
  checkNewPackages <- function(list.of.packages) {
    new.packages.log = !(list.of.packages %in% installed.packages()[,"Package"])
    new.packages <- list.of.packages[new.packages.log]
    return(new.packages)
  }
  new.packages = checkNewPackages(list.of.packages)
  if (length(new.packages)) {
    print(paste('Trying to install the following packages:', new.packages))
    install.packages(new.packages)
    new.packages = checkNewPackages(list.of.packages)
    if (length(new.packages)) {
      print(paste(new.packages, 'were not installed through the easy way'))
      print("Let's try the hard way then")
      setRepositories(graphics = F, ind = 1:8)
      install.packages(new.packages)
      new.packages = checkNewPackages(list.of.packages)
      if (length(new.packages)) {
        stop('forceLibrary was not able to install the following packages: ', 
             new.packages)
      }
    }
  } 
  
  lapply(list.of.packages, library, character.only = T)
  
  print(NULL)
}

forceSetWd = function(x) {
  if (dir.exists(x)) {
    setwd(x)
  } else {
    dir.create(x)
    if (dir.exists(x)) {
      setwd(x)
    } else {
      warning(c('Warning: ', x, 
                ' could not be created as a dir due to permission issues'))
    }
  }
}

getSalmonCols = function(cols =  NULL, salfiles = '/quant.sf') {
  samples = list.dirs(recursive = F)
  filename = paste0(samples[1], salfiles)
  rna = read.table(filename, header = T)
  if (is.null(cols)) {
    cols = 1:ncol(rna)
  } else {
    cols = grep(cols, colnames(rna))
  }
  rna_num = data.frame(rna[, cols])
  rownames(rna_num) = rownames(rna)
  colnames(rna_num)[ncol(rna_num)] = paste0(colnames(rna_num)[ncol(rna_num)], samples[1])
  for (sample in samples[-1]) {
    filename = paste0(sample, '/quant.sf')
    rna = read.table(filename, header = T)
    rna_num = cbind(rna_num, rna[, cols])
    colnames(rna_num)[ncol(rna_num)] = paste0(colnames(rna_num)[ncol(rna_num)], sample)
  }
  return(rna_num)
} 

openMart2018 <- function(variables) {
  library(biomaRt)
  mart.human = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                       dataset = 'hsapiens_gene_ensembl',
                       host = 'http://apr2018.archive.ensembl.org') 
}