#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the
# BiSC and SEPAGES cohorts.
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org
# Date: November 2026
#==============================================================================#
# [DESCRIPTION]: Shared BiSC + SEPAGES quantile cut-offs of the 4 mixture
# compounds, per time point: pooled CDF = N-weighted mix of the SEPAGES CDF
# (from its percentile grid) and the BiSC empirical CDF, then inverted at the
# interior quantile probabilities. Used by 4_Run_Gcomp.R with cut-offs = shared.
# Run on the BiSC PC:
#   Rscript script/8_Mixtures/3_Join_Grid_Sepages_Bisc.R
#
# [INPUT]:  results/sepages/8_Mixtures/2_Get_Grid_Sepages/sepages_pctgrid.csv
#           db/subsets/df_main_<tp>.rds (BiSC tp sample)
# [OUTPUT]: results/8_Mixtures/3_Shared_Cutoffs/shared_bisc_sepages_q<q>.rds
#           (copy it to the SEPAGES PC to run SEPAGES with shared cut-offs)
#==============================================================================#
source("script/00_paths.R")

compounds  <- c("pest_3pba_ilog2", "pest_tcpy_ilog2", "pest_dmp_ms_ilog2", "pest_dep_ms_ilog2")
main_cov   <- c("mo_age", "mo_edu", "mo_smoke", "gest_age", "ch_eth", "ch_sex", "Meanlog2oddsContamination", "hospital")
cell_types <- c("Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC")
q_mix      <- 3                                              # 3 = tertiles, 4 = quartiles
probs_cut  <- seq(0, 1, by = 1 / q_mix)[-c(1, q_mix + 1)]    # interior cut points

grid_file <- results_dir("sepages", "8_Mixtures", "2_Get_Grid_Sepages", "sepages_pctgrid.csv")
stopifnot(file.exists(grid_file))   # copied with results/sepages/ from the SEPAGES PC
sep_grid_long  <- read.csv(grid_file)
shared_cutoffs <- list()

for (tp in c("t2", "t3", "avg")) {
  bisc_df   <- readRDS(db_subset("bisc", "main", "tp", tp))   # BiSC 'tp' sample
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

outdir <- shared_cutoffs_dir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
saveRDS(shared_cutoffs, file.path(outdir, paste0("shared_bisc_sepages_q", q_mix, ".rds")))
