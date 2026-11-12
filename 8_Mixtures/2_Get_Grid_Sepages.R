#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the
# BiSC and SEPAGES cohorts.
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: SEPAGES percentile grid (P1-P99, step 0.5) of the 4 mixture
# compounds per time point, on the main analytic sample (complete cases).
# Only summary percentiles leave the SEPAGES PC: 3_Join_Grid_Sepages_Bisc.R
# uses them on the BiSC PC to build shared BiSC + SEPAGES quantile cut-offs.
# Run on the SEPAGES PC:
#   Rscript script/8_Mixtures/2_Get_Grid_Sepages.R
#
# [INPUT]:  db/subsets/df_main.rds (SEPAGES)
# [OUTPUT]: results/sepages/8_Mixtures/2_Get_Grid_Sepages/sepages_pctgrid.csv
#==============================================================================#
source("script/00_paths.R")

compounds  <- c("pest_3pba_ilog2", "pest_tcpy_ilog2", "pest_dmp_ms_ilog2", "pest_dep_ms_ilog2")
main_cov   <- c("mo_age", "mo_edu", "mo_smoke", "gest_age", "ch_eth", "ch_sex", "Meanlog2oddsContamination")
cell_types <- c("Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC")
probs_grid <- seq(0.01, 0.99, by = 0.005)   # 1st-99th percentile (no extremes)

sep_df <- readRDS(db_subset("sepages", "main", "2tp"))
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

outdir <- results_dir("sepages", "8_Mixtures", "2_Get_Grid_Sepages", create = TRUE)
# Plain csv on purpose: read back by later scripts (not a results table)
write.csv(sep_grid_long, file.path(outdir, "sepages_pctgrid.csv"), row.names = FALSE)

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----