##### Library and functions #####
list.of.packages <- c("biomaRt")
inst.packages = installed.packages()[,"Package"]
new.packages <- list.of.packages[!(list.of.packages %in% inst.packages)]
if (length(new.packages)) install.packages(new.packages)
# setRepositories()

library(biomaRt) # setRepositories() # install.packages('biomaRt')

tablePct <- function(x) {
  x.table = table(x)
  x.proportion = x.table / length(x)
  x.percent = x.proportion * 100
}

plot.freq <- function(x, y = NULL, outliers.rm = F, x.lab = '', y.lab = '', 
                      fill = '', x.name = '', y.name = '', angle = 0, 
                      nbreaks = 10, ...) {
  library(ggplot2)
  naToZero <- function(x) {
    x[is.na(x)] = 0
    return(x)
  }
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

##### Analysis #####

setwd('/share/analysis/hecatos/juantxo/score_protein_analysis/data_extra')
mirtarbase.table = read.csv(file = 'hsa_MTI.csv', 
                            header = T, 
                            stringsAsFactors = F)
setwd("/share/analysis/hecatos/juantxo/Score/output/Output_Run_mrna_SEPT2019/")
setwd('V3/analysis/UNTR/other/')
proteins.table = read.table('20190916_low_correlated_TPMvsTRC.tab')
proteins.naam = rownames(proteins.table)
mart.human = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                     dataset = 'hsapiens_gene_ensembl',
                     host = 'http://apr2018.archive.ensembl.org') 

# View(listAttributes(mart = mart.human))
mirtarbase.Gn.UpId = getBM(attributes = c('external_gene_name', 
                                                'uniprot_gn'),
                           filters = 'uniprot_gn',
                           values = proteins.naam, 
                           mart = mart.human)

rows = mirtarbase.table$Target.Gene %in% mirtarbase.Gn.UpId[, 1]
mirtarbase.table.trc = mirtarbase.table[rows, ]

# How different are the support types between both groups?
mirtar.trc.supptype.pct = tablePct(mirtarbase.table.trc$Support.Type) 
mirtar.supptype.pct = tablePct(mirtarbase.table$Support.Type) 
mirtar.TrcVsAll.diff = mirtar.trc.supptype.pct - mirtar.supptype.pct
print(mirtar.TrcVsAll.diff)

# Are validated relationships having a bigger impact on TRC?
# Let's first select only the best validated interactions
mtb.validated.rows = mirtarbase.table$Support.Type == 'Functional MTI'
mtb.validated = mirtarbase.table[mtb.validated.rows, ]

# Let's check the expression of the miRNAs involved, and subset only the most
# expressed ones
setwd('/share/analysis/hecatos/juantxo/Score/input/miRNA_miRge2/UNTR/')
mirna_expr = read.csv(file = 'miR.RPM.csv', stringsAsFactors = F)
mirna_expr$RPM.avg = apply(X = mirna_expr[, -1], MARGIN = 1, mean, na.rm = T)
quant_90 = quantile(x = mirna_expr$RPM.avg, probs = 0.9, na.rm = T)
mirna_expr = mirna_expr[mirna_expr$RPM.avg >= quant_90, ]

mtb_most_expr = merge(x = mtb.validated, y = mirna_expr, by = 'miRNA')

# To see whether they have an effect on TRC, we need to know which transcripts
# are involved in these inhibitory interactions

gn_eti = getBM(attributes = c('external_gene_name',
                              'uniprot_gn',
                              'ensembl_transcript_id', 
                              'transcript_biotype'), 
               filters = 'external_gene_name', 
               values =  unique(mtb_most_expr$Target.Gene), 
               mart = mart.human)
gn_eti.pc = grep(pattern = 'protein_coding', x = gn_eti$transcript_biotype)
gn_eti = gn_eti[gn_eti.pc, ]

mtb_gn_uniprot = merge.data.frame(x = mtb_most_expr, by.x = 'Target.Gene',
                               y = gn_eti, by.y = 'external_gene_name')

# First let's check on the TRC values of these interactions
setwd('/share/analysis/hecatos/juantxo/Score/output/Output_Run_mrna_SEPT2019/')
setwd('V3/analysis/UNTR/shift_FALSE/TRC/_Wed_Sep_11_14:13:37_2019_/')
setwd('minimum_expressed_samples12/')
trc.table = read.table(file = 'TRC_values_and_protein_expression.tsv', 
                       header = T, sep = '\t')
trc.values = trc.table[seq(from = 2, to = nrow(trc.table), by = 2), ]
rownames(trc.values) = gsub(pattern = ' TRC values', replacement = '', 
                             x = rownames(trc.values))
colnames(trc.values) = Map(f = paste0, colnames(trc.values), '.TRC')
trc.values$TRC.avg = apply(X = trc.values, MARGIN = 1, FUN = mean, na.rm = T)
trc.values$uniprot_gn = rownames( trc.values)

mtb_trc = merge.data.frame(x = mtb_gn_uniprot, y = trc.values, 
                           by = 'uniprot_gn')
colnames(mtb_trc) = gsub(pattern = '.fastq', replacement = '.miRNA.RPM', 
                         x = colnames(mtb_trc))

summary(mtb_trc$TRC.avg[!duplicated(mtb_trc$uniprot_gn)])
mtb.subset = mtb_trc[!duplicated(mtb_trc$uniprot_gn), ]
summary(trc.values$TRC.avg)

plot.freq(x = trc.values$TRC.avg, y = mtb.subset$TRC.avg, outliers.rm = T, 
          x.name = 'Global', y.name = 'Validated interactions', 
          x.lab = 'TRC values', y.lab = '% proteins')
