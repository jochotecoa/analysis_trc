#' Utility functions for TRC and Proteomics analysis
#' Consolidated from functions.R, functions_JOA.R and analysis scripts.

if (file.exists("config.R")) {
  source("config.R")
} else if (file.exists("../config.R")) {
  source("../config.R")
}

#' Install and load a list of packages
#' @param list.of.packages Vector of package names
forceLibrary <- function(list.of.packages) {
  checkNewPackages <- function(list.of.packages) {
    new.packages.log = !(list.of.packages %in% installed.packages()[,"Package"])
    new.packages <- list.of.packages[new.packages.log]
    return(new.packages)
  }
  new.packages = checkNewPackages(list.of.packages)
  if (length(new.packages)) {
    print(paste('Trying to install the following packages:', paste(new.packages)))
    install.packages(new.packages)
    new.packages = checkNewPackages(list.of.packages)
    if (length(new.packages)) {
      print(paste(paste(new.packages), 'were not installed through the easy way'))
      print("Let's try the hard way then")
      setRepositories(graphics = F, ind = 1:8)
      install.packages(new.packages)
      new.packages = checkNewPackages(list.of.packages)
      if (length(new.packages)) {
        stop('forceLibrary was not able to install the following packages: ', 
             paste(new.packages))
      }
    }
  } 
  
  lapply(list.of.packages, library, character.only = T)
  invisible()
}

#' Clean Protein IDs from UniProt format
#' @param protein_table Data frame with IDs in the first column
cleanProtIds = function(protein_table) {
  protein_table = protein_table[!grepl(protein_table[, 1], pattern = ':'), ]
  if (class(protein_table) != "data.frame") {
    protein_table = as.data.frame(protein_table)
  }
  names = strsplit(as.character(protein_table[, 1]), '\|')
  names = as.character(lapply(names, '[', 2))
  protein_table$uniprot_gn = names
  return(protein_table)
}

#' Convert NA values to Zero
naToZero = function(x) {
  x[is.na(x)] = 0
  return(x)
}

#' Convert Zero values to NA
zeroToNa = function(x) {
  x[x == 0] = NA
  return(x)
}

#' Add a pseudocount to values
pseudocount = function(x, addition = 1) {
  x = x + addition
  return(x)
}

#' Convert Transcript IDs to Gene IDs using biomaRt
transcrToGene = function(table, aggregate = F, prot_cod = F) {
  forceLibrary(c('biomaRt', 'dplyr'))
  
  enst_col = apply(table, 2, grepl, pattern = 'ENST') %>% 
    apply(2, any) %>%
    as.logical()
  if (sum(enst_col) == 0) {
    table[, 'rownames'] = rownames(table)
    enst.rown = grepl(pattern = 'ENST', x = table[, 'rownames'])
    if (!sum(enst.rown)) {stop(print(table[1, ]))}
    enst_col = grepl(pattern = 'rownames', x = colnames(table))
  }
  
  table[, enst_col] = as.character(table[, enst_col])
  version = grepl('\.', table[, enst_col]) %>% sum()
  if (version > 0) {
    transcript_id = 'ensembl_transcript_id_version'
  } else {
    transcript_id = 'ensembl_transcript_id'
  }
  
  values = table[, enst_col]
  mart.human = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                       dataset = 'hsapiens_gene_ensembl',
                       host = 'http://apr2018.archive.ensembl.org') 
  if (prot_cod) {
    transcr_biotypes = getBM(attributes = c(transcript_id, 'transcript_biotype'), 
                             filters = transcript_id, values = values, 
                             mart = mart.human)
    isProtCod = transcr_biotypes$transcript_biotype == 'protein_coding'
    values = values[isProtCod]
    print(paste0(sum(!isProtCod), ' transcripts were not protein_coding'))
  }
  
  new_cols = getBM(attributes = c(transcript_id, 'ensembl_gene_id'), 
                   filters = transcript_id, values = values, mart = mart.human)
  
  table = merge.data.frame(x = table, y = new_cols, 
                           by.x = colnames(table)[enst_col], 
                           by.y = transcript_id)
  if (nrow(table) < length(values)) {
    print(paste0(length(values) - nrow(table), 
                 ' transcripts had no gene matched'))
  }
  if (prot_cod & !aggregate) {
    table = tibble::rownames_to_column(table)
  }
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

#' Remove miRNA columns (those matching 'hsa')
rmMirnas = function(x) {
  mirna.cols = grep(pattern = 'hsa', x = colnames(x))
  if (length(mirna.cols) > 0) {
    x = x[, -mirna.cols]
  }
  return(x)
}

#' Set working directory and create if it doesn't exist
forceSetWd = function(x) {
  if (dir.exists(x)) {
    setwd(x)
  } else {
    dir.create(x, recursive = T)
    if (dir.exists(x)) {
      setwd(x)
    } else {
      warning(c('Warning: ', x, 
                ' could not be created as a dir due to permission issues'))
    }
  }
}

#' Alias for forceSetWd
setOrCreatewd <- forceSetWd

#' Merge multiple files into a single data frame
mergeFiles = function(files_patt =  'quant.sf', by_col = 'Name', 
                      row_names = F, progr_bar = T, ...) {
  if (progr_bar) {
    forceLibrary('pbmcapply')
  }
  forceLibrary('dplyr')
  files = list.files(pattern = files_patt, recursive = T)
  files = files[!grepl('total', files)]
  print(paste('Number of files found:', length(files)))
  if (length(files) == 0) return(NULL)
  
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
  
  if (progr_bar) {
    pb = txtProgressBar(max = length(files))
  }
  
  for (i in 2:length(files)) {
    file = files[i]
    voom_file = read.table(file, header = T, stringsAsFactors = F)
    if (row_names) {
      voom_file = voom_file %>% tibble::rownames_to_column() %>% 
        dplyr::select(rowname, everything())
    }
    colnames(voom_file)[-1] = paste(colnames(voom_file)[-1], file, sep = '_')
    big_quant_voom = merge.data.frame(big_quant_voom, voom_file, by = by_col, ...)
    if (progr_bar) {
      setTxtProgressBar(pb, i)
    }
  }
  if (progr_bar) {
    close(pb)
  }
  return(big_quant_voom)
} 

#' Open Ensembl Mart (2018 archive)
openMart2018 <- function(...) {
  forceLibrary('biomaRt')
  mart.human = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                       dataset = 'hsapiens_gene_ensembl',
                       host = 'http://apr2018.archive.ensembl.org', ...) 
  return(mart.human)
}

#' Filter for protein coding transcripts
filterProtCod = function(table) {
  forceLibrary('biomaRt')
  
  sampl = table[nrow(table), ]
  enst_col = grep(pattern = 'ENST', x = sampl)
  if (length(enst_col) > 1) {enst_col = enst_col[1]}
  if (length(enst_col) == 0) {
    table[, 'rownames'] = rownames(table)
    sampl = table[nrow(table), ]
    enst_col = grep(pattern = 'ENST', x = sampl)[1]
  }
  
  version = grepl('\.', sampl[, enst_col])
  if (length(version) == 0) {version = F}
  if (version) {
    transcript_id = 'ensembl_transcript_id_version'
  } else {
    transcript_id = 'ensembl_transcript_id'
  }
  
  values = table[, enst_col]
  mart.human = openMart2018()
  
  transcr_biotypes = getBM(attributes = c(transcript_id, 'transcript_biotype'), 
                           filters = transcript_id, values = values, 
                           mart = mart.human)
  isProtCod = transcr_biotypes$transcript_biotype == 'protein_coding'
  values = values[isProtCod]
  print(paste0(sum(!isProtCod), ' transcripts were not protein_coding'))
  new_table = table[table[, enst_col] %in% values, ]
  return(new_table)
}

#' Filter samples by sequencing depth
filterSamplesBySeqDepth = function(df) {
  seq_depth_ratio <- df %>% 
    colSums(na.rm = T) %>% 
    `/` (mean(.)) %>% 
    log2() 
  
  if (!all(seq_depth_ratio > -2 & seq_depth_ratio < 2)) {
    bad_samples = names(seq_depth_ratio)[seq_depth_ratio <= -2 | seq_depth_ratio >= 2]
    warning(length(bad_samples), 
            ' sample(s) filtered out due to sequencing depth: ', 
            paste(bad_samples, collapse = ", "), immediate. = T)
  }
  
  df = df[, seq_depth_ratio > -2 & seq_depth_ratio < 2]
  return(df)
}

#' Apply a function to rows of a data frame comparing two groups of columns
apply_2D = function(df, FUN, col.x = NULL, col.y = NULL, complete_cases = NULL, 
                    y = NULL, ...) {
  forceLibrary('dplyr')
  if (is.null(complete_cases)) {
    complete_cases = length(col.x)
  }
  
  result.list = list()
  
  for (row in rownames(df)) {
    if (!is.null(y)) {
      X = df[row, ] %>% as.numeric()
      row.y = grep(row, rownames(df))[1]
      Y = y[row.y, ] %>% as.numeric()
    } else {
      X = df[row, col.x] %>% as.numeric()
      Y = df[row, col.y] %>% as.numeric()
    }
    
    if (sum(complete.cases(X, Y)) >= complete_cases) {
      res = FUN(x = X, y = Y, ...)
      result.list[[row]] = unlist(res)
    } 
  }
  
  if (length(result.list) > 0) {
    result.df = do.call(rbind, lapply(result.list, as.data.frame))
    colnames(result.df) = names(result.list[[1]])
    return(result.df)
  } else {
    return(NULL)
  }
}

#' Interquartile Range
iq = function(x, na.rm = F) {
	if (na.rm) {x = na.omit(x)}
	first_q = as.numeric(quantile(x, 0.25))
	third_q = as.numeric(quantile(x, 0.75))
	y = third_q - first_q
	return(y)
}

#' Identify outliers
outlier = function(x, y=NULL, na.rm = F, onlyExtreme = F) {
	if (na.rm) {x = na.omit(x)}
	q1 = quantile(x, 0.25)
	q3 = quantile(x, 0.75)
	multiplier = if(onlyExtreme) 3 else 1.5
	
	if (!is.null(y)) {
		isOutlier = y < q1 - multiplier*iq(x) | y > q3 + multiplier*iq(x)
		return(as.logical(isOutlier))
	}
	
	outliers = x[x < q1 - multiplier*iq(x) | x > q3 + multiplier*iq(x)]
	return(outliers)
}

#' Shift values by median
shift_median = function(df, median_of_medians) {
	for (i in colnames(df)) {
		sample = df[, i]
		corr_factor = median_of_medians - median(sample, na.rm = T)
		df[, i] = df[, i] + corr_factor
	}
	return(df)
}

#' Normalize proteomics data
normalizeProteomics = function(df) {
	medians = apply(df, 2, median, na.rm = T)
	median_of_medians = median(medians, na.rm = T)
	df = shift_median(df, median_of_medians)
	return(df)
}

#' Compare two datasets by column summary
compare <- function(dataset1, dataset2, column) {
  difference = summary(dataset2[,column]) - summary(dataset1[,column])
  return(abs(difference))
}

#' Generate frequency distribution barplot
freq.dist <- function(column, type = 'sd', from = '', to = '', rot_angle = 45, main ='', xlab = '', ylab = '') {
  sd_val = sd(column, na.rm = T)
  mn = mean(column, na.rm = T)
  mxm = max(column, na.rm = T)
  mnm = min(column, na.rm = T)
  if (type == 'sd') {
    bye = (2*sd_val + mn - (mn - 2*sd_val))/10
    brks = seq(from = mnm, to = mxm, by = bye)
  }
  if (type == 'minmax') {
    bye = (mxm - mnm)/10
    brks = seq(from = mnm, to = mxm, by = bye)
  }
  if (type == 'custom') {
    bye = (to - from)/10
    brks = seq(from = from, to = to, by = bye)
  }
  
  if (type == 'sd') {
    if (mnm < (mn - 2*sd_val)) {
      brks = c(mnm, brks)
    }
    if (mxm > (mn + 2*sd_val)) {
      brks = c(brks, mxm)
    }
  }
  ct = cut(column, breaks = brks)
  tbl = table(ct)
  plt <- barplot(tbl, col='steelblue', xaxt="n", main = main, xlab = xlab, ylab = ylab)
  text(plt, par("usr")[3], labels = names(tbl), srt = rot_angle, adj = c(1.1,1.1), xpd = TRUE, cex=0.6) 
  return(tbl)
}

#' Calculate table percentages
tablePct <- function(x) {
  x.table = table(x)
  x.proportion = x.table / length(x)
  x.percent = x.proportion * 100
  return(x.percent)
}

#' Plot frequency distribution using ggplot2
plot.freq <- function(x, y = NULL, outliers.rm = F, x.lab = '', y.lab = '', 
                      fill = '', x.name = '', y.name = '', angle = 0, 
                      nbreaks = 10, ...) {
  forceLibrary('ggplot2')
  
  x = naToZero(x)
  if (outliers.rm) {
    x.old.min = min(x)
    x.old.max = max(x)
    Q1 = quantile(x, .25, na.rm = T)
    Q3 = quantile(x, .75, na.rm = T)
    IQ = Q3 - Q1
    lower.outer.fence = Q1 - 3*IQ
    upper.outer.fence = Q3 + 3*IQ
    x = x[x > lower.outer.fence & x < upper.outer.fence]
  }
  x.range = max(x) - min(x)
  x.breaks = seq(from = min(x), to = max(x), by = (x.range/nbreaks))
  if (outliers.rm) {
    x.breaks = c(x.old.min, x.breaks[2:nbreaks], x.old.max)
  }
  x.cut = cut(x, breaks = x.breaks, include.lowest = T)
  x.table = table(x.cut) / length(x) * 100
  if (!is.null(y)) {
    y = naToZero(y)
    y.cut = cut(y, breaks = x.breaks, include.lowest = T)
    y.table = table(y.cut) / length(y) * 100
    xy.table = rbind(x.table, y.table)
    xy.df = as.data.frame(xy.table)
    f = merge(stack(xy.df), stack(as.data.frame(t(xy.df))), by = 'values')
    if (x.name != '') {
      f[, 3] = gsub(pattern = 'x.table', replacement = x.name, x = f[, 3])
    }
    if (y.name != '') {
      f[, 3] = gsub(pattern = 'y.table', replacement = y.name, x = f[, 3])
    }
    ggplot(f, aes(x = f[, 2], y = f[, 1], fill = f[, 3])) +
      geom_bar(stat = "identity", position = "dodge") +
      theme(axis.text.x = element_text(angle = angle)) +
      labs(fill = fill, x = x.lab, y = y.lab)
    
  } else {
    x.df = as.data.frame(x.table)
    x.df$names = rownames(x.df)
    ggplot(x.df, aes(x = x.cut, y = x.df[, 2])) +
      geom_bar(stat = 'identity') +
      theme(axis.text.x = element_text(angle = angle)) +
      labs(x = x.lab, y = y.lab, ...)
  }
}

#' Generate sample names based on experimental design
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

#' Read a table and assign it to a variable name if not already in memory
get.or.assign <- function(var.name, filepath) {
  if (exists(var.name, envir = .GlobalEnv)) {
    file = get(var.name, envir = .GlobalEnv)
  } else {
    file = read.table(filepath)
    assign(x = var.name, value = file, envir = .GlobalEnv)
  }
  return(file)
}

#' Read and format proteomics data
read.and.format.proteindata <- function(proteomicsfilepath, comp, samples) {
  proteinvaluesfile.colnames = read.table(proteomicsfilepath, nrows = 1, 
                                          stringsAsFactors = F)
  # Attempt reading with space separator first
  proteinvaluesfile = tryCatch({
    read.table(proteomicsfilepath, skip = 1, stringsAsFactors = F, 
               header = F, sep = ' ', fill = T)
  }, error = function(e) NULL)
  
  if (is.null(proteinvaluesfile) || ncol(proteinvaluesfile) != ncol(proteinvaluesfile.colnames)) {
    proteinvaluesfile = read.table(proteomicsfilepath, skip = 1, stringsAsFactors = F, 
                                   header = F, sep = '	', fill = T)
  }
  
  proteinvaluesfile = proteinvaluesfile[grep('.*\|.*\|.*', proteinvaluesfile[,1]),]
  proteinvaluesfile = proteinvaluesfile[!grepl(pattern = ':', proteinvaluesfile$V2),]
  
  data_cols = 2:ncol(proteinvaluesfile)
  proteinvaluesfile[,data_cols] = sapply(proteinvaluesfile[,data_cols], as.double)
  
  means = sapply(proteinvaluesfile, mean, na.rm = T)
  proteinvaluesfile = cbind(proteinvaluesfile[,1], proteinvaluesfile[,!is.na(means)])
  
  colnames(proteinvaluesfile) = proteinvaluesfile.colnames[1, 1:ncol(proteinvaluesfile)]
  
  protvalcomp = proteinvaluesfile[,c(1, grep(toupper(comp), colnames(proteinvaluesfile)))]
  selectedtimepointcolumns = c(1, grep(pattern = paste(samples, collapse = '|'), 
                                      x = colnames(protvalcomp)))
  
  protvalcompsel = protvalcomp[rowSums(!is.na(protvalcomp[,selectedtimepointcolumns])) > 2, 
                               selectedtimepointcolumns]
  return(protvalcompsel)
}

#' Save metadata for a list of objects
saveMetadata <- function(x, dir = 'metadata') {
  original_wd = getwd()
  forceSetWd(dir)
  for (obj.name in x) {
    if (exists(obj.name)) {
      obj.data = get(obj.name)
      # Skip functions and closures
      if (!is.function(obj.data) && !is.primitive(obj.data)) {
        tryCatch({
          write.table(x = obj.data, file = paste0(obj.name, '.tsv'), sep = "\t", quote = F)
        }, error = function(e) {
          warning(paste("Could not save metadata for", obj.name, ":", e$message))
        })
      }
    }
  }
  setwd(original_wd)
}
