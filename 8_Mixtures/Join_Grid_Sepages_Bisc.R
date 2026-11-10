compounds  <- c("pest_3pba_ilog2", "pest_tcpy_ilog2", "pest_dmp_ms_ilog2", "pest_dep_ms_ilog2")
main_cov   <- c("mo_age", "mo_edu", "mo_smoke", "gest_age", "ch_eth", "ch_sex", "Meanlog2oddsContamination", "hospital")
cell_types <- c("Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC")
q_mix      <- 3                                              # 3 = tertiles, 4 = quartiles
probs_cut  <- seq(0, 1, by = 1 / q_mix)[-c(1, q_mix + 1)]    # interior cut points

sep_grid_long  <- read.csv("results/sepages/sepages_pctgrid.csv")
shared_cutoffs <- list()

for (tp in c("t2", "t3", "avg")) {
  bisc_df   <- readRDS(paste0("db/subsets/df_main_", tp, ".rds"))   # BiSC 'tp' sample
  exp_names <- paste0(compounds, "_", tp)
  bisc_cc   <- bisc_df[complete.cases(bisc_df[, c(exp_names, main_cov, cell_types)]), ]
  cat("\n===", tp, "| BiSC N:", nrow(bisc_cc), "\n")               # check: 378 / 437 / 503
  
  shared_cutoffs[[tp]] <- list()
  for (k in seq_along(compounds)) {
    sep_rows  <- sep_grid_long$tp == tp & sep_grid_long$compound == compounds[k]
    sep_vals  <- sep_grid_long$value[sep_rows]
    sep_probs <- sep_grid_long$prob[sep_rows]
    n_sep     <- sep_grid_long$n[sep_rows][1]
    bisc_vals <- bisc_cc[[exp_names[k]]]
    n_bisc    <- length(bisc_vals)
    
    # Pooled CDF = N-weighted mix of the two cohort CDFs, then invert it
    x_grid   <- sort(unique(c(sep_vals, bisc_vals)))
    F_sep    <- approx(sep_vals, sep_probs, xout = x_grid, rule = 2, ties = "ordered")$y
    F_bisc   <- ecdf(bisc_vals)(x_grid)
    F_pooled <- (n_sep * F_sep + n_bisc * F_bisc) / (n_sep + n_bisc)
    cut_shared <- approx(F_pooled, x_grid, xout = probs_cut, ties = "ordered")$y
    
    # Check: SEPAGES-only cut-offs from the grid (should match your SEPAGES output)
    cut_sep_check <- approx(sep_probs, sep_vals, xout = probs_cut, ties = "ordered")$y
    cat(compounds[k], "| shared:", round(cut_shared, 2),
        "| SEPAGES (grid):", round(cut_sep_check, 2),
        "| BiSC:", round(quantile(bisc_vals, probs_cut, names = FALSE), 2), "\n")
    
    shared_cutoffs[[tp]][[exp_names[k]]] <- cut_shared
  }
}

saveRDS(shared_cutoffs, paste0("results/2_Gcomp_Results/quantile_cutoffs/shared_bisc_sepages_q", q_mix, ".rds"))
