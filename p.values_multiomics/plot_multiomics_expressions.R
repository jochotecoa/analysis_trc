test = avv

colnames(test) = colnames(test) %>% gsub('_TrT', '', .)

if (plotting) {
  for (rown in rownames(test)) {
    layout(matrix(c(1,2), 1, 2, byrow = T))
    row_data = test[rown, ]
    max_tpm = test[rown, , F] %>% dplyr::select(starts_with('target')) %>% max(na.rm = T)
    par(mfrow = c(1, 2))
    
    viewPlotProtx(test, rown)
    
    readline(prompt = "Press [enter] to continue")
    
    plotOmics(df = trt_proteomx, omics = paste0('^', TrT_miF),
              transcript = rown, xaxt = 'n', type = 'b',
              ylab = 'gray: TPM | black: TrT', xlab = '',
              ylim = c(0, max_tpm))
    plotOmics(df = trt_proteomx, omics = '^target', transcript = rown, xaxt = 'n',
              type = 'b', col = 'gray', yaxt = 'n', xlab = '')
    axis(1, at = 1:24, labels = colnames(test)[grep('^target', colnames(test))], las = 2)
    
    readline(prompt = "Press [enter] to continue")
    
    layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
    boxplot.dose(row_data, 'Proteomics', main = rown, ylab = 'Proteomics Expression')
    boxplot.dose(row_data, '^target', main = rown, ylab = 'TPM Expression')
    boxplot.dose(row_data, TrT_miF, main = rown, ylab = 'TrT Expression')
    
    
    test[rown, , F] %>% dplyr::select(starts_with(TrT_miF)) %>% dplyr::select(contains('The')) %>%
      as.numeric() %>% mean(na.rm = T) %>% abline(h = .)
    
    test[rown, , F] %>% dplyr::select(starts_with(TrT_miF)) %>% dplyr::select(contains('Tox')) %>%
      as.numeric() %>% mean(na.rm = T) %>% abline(h = .)
    axis(1, at = 1:24, labels = colnames(trt_proteomx)[grep('^target',
                                                            colnames(trt_df))], las = 2)
    readline(prompt = "Press [enter] to continue")
    
  }
}

