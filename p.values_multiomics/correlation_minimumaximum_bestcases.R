# # Correlation ------------------------------------------------------------------
# 
# cor_trt_prot = apply_2D(df = dplyr::select(trt_proteomx, starts_with(TrT_miF)),
#                         FUN = cor.test, complete_cases = comp_cas*2, 
#                         col.x = 'adsf', col.y = 'asdfa',
#                         y = dplyr::select(trt_proteomx, starts_with('Proteomics'))) %>%
#   rownames_to_column() %>%
#   transmute(estimate.cor_trt = 
#               as.numeric(as.character(estimate.cor)), rowname) %>%
#   column_to_rownames()
# 
# cor_tpm_prot = apply_2D(df = dplyr::select(trt_proteomx, starts_with('target')),
#                         FUN = cor.test, complete_cases = comp_cas*2, 
#                         col.x = 'adsf', col.y = 'asdfa',
#                         y = dplyr::select(trt_proteomx, starts_with('Proteomics'))) %>%
#   rownames_to_column() %>%
#   transmute(estimate.cor_tpm = 
#               as.numeric(as.character(estimate.cor)), rowname) %>%
#   column_to_rownames()
# 
# trt_proteomx = merge.data.frame(x = rownames_to_column(trt_proteomx),
#                                 y = rownames_to_column(cor_trt_prot),
#                                 by = 'rowname', all = T) %>%
#   merge.data.frame(y = rownames_to_column(cor_tpm_prot), by = 'rowname', 
#                    all = T) %>%
#   column_to_rownames()
# 
# test = trt_proteomx %>% rownames_to_column() %>%  
#   filter(estimate.cor_trt > (0.4 + estimate.cor_tpm)) %>% column_to_rownames()
# 
# # Minimum > Maximum ------------------------------------------------------------
# 
# trt_the_range = trt_proteomx %>% dplyr::select(starts_with(TrT_miF)) %>%
#   dplyr::select(contains('The')) %>% na.omit() %>% apply(1, range, na.rm = T) %>%
#   t() %>% as.data.frame() %>% rownames_to_column() %>% 
#   filter(!is.infinite(V1)) %>% rename(min_the = V1, max_the = V2) %>% 
#   column_to_rownames()
# 
# trt_tox_range = trt_proteomx %>% dplyr::select(starts_with(TrT_miF)) %>%
#   dplyr::select(contains('Tox')) %>% na.omit() %>% apply(1, range, na.rm = T) %>%
#   t() %>% as.data.frame() %>% rownames_to_column() %>% 
#   filter(!is.infinite(V1)) %>% rename(min_tox = V1, max_tox = V2) %>% 
#   column_to_rownames()
# 
# a = merge.data.frame(rownames_to_column(trt_the_range),
#                      rownames_to_column(trt_tox_range),
#                      'rowname') %>% column_to_rownames()
# 
# b = a[a$min_the > a$max_tox,,F] %>% rownames()
# test = trt_proteomx[b,,F]
# 
# if (plotting) {
#   trt_proteomx$p.adj_theVStox_TPM %>% 
#     hist(breaks = seq(0,1,0.05), main = 'TrT', xlab = 'p. adjusted values')
#   trt_proteomx$p.adj_theVStox_TPM %>% 
#     hist(breaks = seq(0,1,0.05), main = 'TPM', xlab = 'p. adjusted values')
#   trt_proteomx$p.adj_prx %>% 
#     hist(breaks = seq(0,1,0.05), main = 'Proteomics', xlab = 'p. adjusted values')
# }
# 
# 
# Best cases --------------------------------------------------------------


# # VVV
# vvv = trt_proteomx %>%
#   rownames_to_column() %>%
#   filter(
#     p.adj_theVStox_TPM < 0.05,
#     p.value_tpm > 0.05,
#     p.value_theVStox_prx < 0.05,
#     sign(statistic.t_prx) == sign(statistic.t_trt)
#   ) %>%
#   column_to_rownames()
# # test %>% nrow() %>% print()
#
# # AAA
# aaa = trt_proteomx %>%
#   rownames_to_column() %>%
#   filter(
#     p.adj_theVStox_TPM < 0.05,
#     p.adj_theVStox_TPM > 0.05,
#     p.adj_prx < 0.05,
#     sign(statistic.t_prx) == sign(statistic.t_trt)
#   ) %>%
#   column_to_rownames()
#

# # AVV
# avv = trt_proteomx %>%
#   rownames_to_column() %>%
#   filter(
#     p.adj_theVStox_TPM < 0.05,
#     p.value_tpm > 0.05,
#     p.value_theVStox_prx < 0.05,
#     sign(statistic.t_prx) == sign(statistic.t_trt)
#   ) %>%
#   column_to_rownames()
#
# # AVV_bad
# avv_bad = trt_proteomx %>%
#   rownames_to_column() %>%
#   filter(
#     p.adj_theVStox_TPM < 0.05,
#     p.value_tpm > 0.05,
#     p.value_theVStox_prx > 0.05,
#     sign(statistic.t_prx) == sign(statistic.t_trt)
#   ) %>%
#   column_to_rownames()

# # AVA
# ava = trt_proteomx %>%
#   rownames_to_column() %>%
#   filter(
#     p.adj_theVStox_TPM < 0.05,
#     p.value_tpm > 0.05,
#     p.adj_prx < 0.05,
#     sign(statistic.t_prx) == sign(statistic.t_trt)
#   ) %>%
#   column_to_rownames()

# # VAV
# vav = trt_proteomx %>%
#   rownames_to_column() %>%
#   filter(
#     p.adj_theVStox_TPM > 0.05,
#     p.adj_theVStox_TPM < 0.05,
#     p.value_theVStox_prx > 0.05
#   ) %>%
#   column_to_rownames()
#
# # VAV_bad
# vav_bad = trt_proteomx %>%
#   rownames_to_column() %>%
#   filter(
#     p.adj_theVStox_TPM > 0.05,
#     p.adj_theVStox_TPM < 0.05,
#     p.value_theVStox_prx < 0.05
#   ) %>%
#   column_to_rownames()
