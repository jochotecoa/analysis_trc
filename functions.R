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
  
  sampl = table[nrow(table), ]
  enst_col = grep(pattern = 'ENST', x = sampl)[1]
  if (length(enst_col) == 0) {
    table[, 'rownames'] = rownames(table)
    enst_col = grep(pattern = 'ENST', x = sampl)[1]
  }
  
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
  
  invisible()
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

mergeFiles = function(files_patt =  'quant.sf', by_col = 'Name', row_names = F, ...) {
  
  forceLibrary(c('pbmcapply', 'dplyr'))
  files = list.files(pattern = files_patt, recursive = T)
  files = files[!grepl('total', files)]
  # files = files[-1]
  print(paste('Number of files found:', length(files)))
  file = files[1]
  stopifnot(file.exists(file))
  voom_file = read.table(file, header = T, stringsAsFactors = F)
  if (row_names) {
    voom_file = voom_file %>% tibble::rownames_to_column() %>% 
      dplyr::select(rowname, everything())
    by_col = 'rowname'
  }
  colnames(voom_file)[-1] = paste(colnames(voom_file)[-1], file, sep = '_')
  big_quant_voom = voom_file
  pb = progressBar(max = length(files[-1]))
  for (file in files[-1]) {
    voom_file = read.table(file, header = T, stringsAsFactors = F)
    if (row_names) {
      voom_file = voom_file %>% tibble::rownames_to_column() %>% 
        dplyr::select(rowname, everything())
    }
    colnames(voom_file)[-1] = paste(colnames(voom_file)[-1], file, sep = '_')
    big_quant_voom = merge.data.frame(big_quant_voom, voom_file, by = by_col, ...)
    setTxtProgressBar(pb, grep(file, files[-1]))
  }
  close(pb)
  return(big_quant_voom)
} 

openMart2018 <- function(...) {
  library(biomaRt)
  mart.human = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                       dataset = 'hsapiens_gene_ensembl',
                       host = 'http://apr2018.archive.ensembl.org', ...) 
}