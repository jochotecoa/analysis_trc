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
  png(paste(strings.naam[1], strings.naam[2], strings.naam[3], 'frequency-preferred-timeshift_histogram.png', sep = '_'))
  strings.naam[3] = unlist(strsplit(strings.naam[3], split = '-'))[2]
  strings.naam.pasted = paste('Amount of proteins per preferred time shift\n', strings.naam[1], strings.naam[2], '(Min. samples =', strings.naam[3], ')')
  setwd('plots/plots_timeshifts_nproteins_stringencies_optimalshift/')
  barplot(data.freq, main = strings.naam.pasted, ylim = c(0,400))
  dev.off()
}

