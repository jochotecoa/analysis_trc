#' Time-shift Analysis of Correlation between TRC and Protein Expression
#' 
#' This script analyzes protein half-lives in relation to observed time-shifts
#' in transcript-protein correlation.
#'
#' Inputs:
#' - Correlation results with shift information (from 06_Correlation_With_Shifts.R)
#' - Protein half-life data (Nature 2011)

if (file.exists("utils.R")) { source("utils.R") } else if (file.exists("../utils.R")) { source("../utils.R") }
forceLibrary(c('biomaRt', 'plyr', 'ggplot2', 'scales'))

# INPUT VARIABLES
comp = if(exists("DEFAULT_COMP")) DEFAULT_COMP else 'UNTR'
results_dir = 'results'
datafile.path = file.path(results_dir, comp, 'TRC_shifted', 'correlation_results_TRC_protein_shifted.tsv')

if (file.exists(datafile.path)) {
    datafile = read.table(datafile.path, header = T, stringsAsFactors = F)
    shifts = as.numeric(names(table(datafile$shiftlist)))
} else {
    stop("Input datafile not found: ", datafile.path)
}

# Open half-life data and format it
# Assuming the file is in a data directory or current dir
halflife_file = 'nature10098-s5_halflife_proteins.csv'
if (!file.exists(halflife_file)) {
    # Try looking in archive or other known locations
    halflife_file = 'archive/nature10098-s5_halflife_proteins.csv' # Placeholder if it exists there
}

if (file.exists(halflife_file)) {
    halflives = read.csv(halflife_file, header = T, sep = ';', stringsAsFactors = F)
} else {
    warning("Half-life data file missing: ", halflife_file)
    halflives = data.frame()
}

if (nrow(halflives) > 0) {
    # Clean half-life data
    col_start = 10
    if (ncol(halflives) >= col_start) {
        halflives[, col_start:ncol(halflives)] = as.data.frame(sapply(halflives[, col_start:ncol(halflives)], 
                                                              function(x) as.numeric(gsub(',', '.', x))))
    }
}

output_dir = file.path(results_dir, comp, 'time_shift_analysis')
setOrCreatewd(output_dir)

# Initialize biomaRt
mart.human = openMart2018()
mart.mouse = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                     dataset = 'mmusculus_gene_ensembl',
                     host = 'http://apr2018.archive.ensembl.org')

for (s in shifts) {
  timeshift.analyzed = datafile$protlist[datafile$shiftlist == s]
  if (length(timeshift.analyzed) == 0) next
  
  # Get homologs mapping (UniProt Human to Mouse)
  ensembl.prot.ids.mouse = getBM(attributes = 'mmusculus_homolog_ensembl_peptide',
                                 filters = "uniprot_gn_id",
                                 values = timeshift.analyzed, 
                                 mart = mart.human)
  
  uniprot.ids.mouse.timeshift = getBM(attributes = 'uniprot_gn_id',
                                      filters = 'ensembl_peptide_id',
                                      values = ensembl.prot.ids.mouse$mmusculus_homolog_ensembl_peptide, 
                                      mart = mart.mouse)
  
  uniprot.ids.mouse.timeshift = unlist(uniprot.ids.mouse.timeshift)
  
  # Match with halflives
  matches = halflives$Uniprot.IDs %in% uniprot.ids.mouse.timeshift
  timeshift.halflives = halflives[matches, ]
  
  if (nrow(timeshift.halflives) > 0) {
      # Plot half-life distribution for this shift
      # Focus on a representative column, e.g., the 11th column as in v8
      if (ncol(halflives) >= 11) {
          col_idx = 11
          col_name = colnames(halflives)[col_idx]
          
          png(file.path(output_dir, paste0('halflife_dist_shift_', s, '.png')))
          hist(timeshift.halflives[, col_idx], 
               main = paste('Half-life Dist - Shift', s),
               xlab = col_name, col = 'lightgreen')
          dev.off()
      }
  }
}
