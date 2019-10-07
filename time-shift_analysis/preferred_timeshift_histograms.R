setwd('/share/analysis/hecatos/juantxo/score_protein_analysis/')
all.files = list.files()
shiftlist.files.index = grep(pattern = 'UNTR_correlation_results_.*_protein_1onseveral_until072_minimumexpressedsamples-.*_shift-TRUE_timeps-shifted-002.tsv', x = all.files)
shiftlist.files = all.files[shiftlist.files.index]
for (shiftlist in shiftlist.files) {
  setwd('/share/analysis/hecatos/juantxo/score_protein_analysis/')
  data = read.table(shiftlist, header = T)
  data.freq = table(data$shiftlist)
  splitted.strings = unlist(strsplit(shiftlist, split = '_'))
  strings.naam = splitted.strings[c(1,4,length(splitted.strings)-2)]
  strings.naam.pasted = paste(strings.naam[1], strings.naam[2], strings.naam[3], 'amount of proteins per preferred time shift')
  setwd('plots/')
  png(paste(strings.naam[1], strings.naam[2], strings.naam[3], 'frequency-preferred-timeshift_histogram.png', sep = '_'))
  barplot(data.freq, main = strings.naam.pasted)
  dev.off()
}

