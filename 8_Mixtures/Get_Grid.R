compounds  <- c("pest_3pba_ilog2", "pest_tcpy_ilog2", "pest_dmp_ms_ilog2", "pest_dep_ms_ilog2")
main_cov   <- c("mo_age", "mo_edu", "mo_smoke", "gest_age", "ch_eth", "ch_sex", "Meanlog2oddsContamination")
cell_types <- c("Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC")
probs_grid <- seq(0.01, 0.99, by = 0.005)   # 1st-99th percentile (no extremes)

sep_df <- readRDS("db/subsets/df_main.rds")
sep_grid_long <- NULL

for (tp in c("t2", "t3", "avg")) {
  exp_names <- paste0(compounds, "_", tp)
  sep_cc    <- sep_df[complete.cases(sep_df[, c(exp_names, main_cov, cell_types)]), ]
  cat(tp, "| SEPAGES N (complete cases):", nrow(sep_cc), "\n")   # check: should be 312
  for (k in seq_along(compounds)) {
    sep_vals <- quantile(sep_cc[[exp_names[k]]], probs = probs_grid, names = FALSE)
    sep_grid_long <- rbind(sep_grid_long,
                           data.frame(tp = tp, compound = compounds[k], prob = probs_grid,
                                      value = round(sep_vals, 4), n = nrow(sep_cc)))
  }
}

write.csv(sep_grid_long, "results/2_Gcomp_Results/quantile_cutoffs/sepages_pctgrid.csv", row.names = FALSE)