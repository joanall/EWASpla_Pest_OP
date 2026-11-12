#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the
# BiSC and SEPAGES cohorts.
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: Tertile and quartile cut-offs of the 4 mixture compounds
# (3-PBA, TCPy, ΣDMP, ΣDEP), as used in 4_Run_Gcomp.R: same analytic sample
# (complete cases for the compounds + covariates of the analysis), same
# variables (LOD-imputed, log2) and same function (qgcomp::quantize).
# Plot 1: per time point, histogram of each compound with tertile (solid) and
#         quartile (dashed) cut-offs; one panel per compound (2 x 2)
# Plot 2: the same, one PNG per compound and time point (solo/)
# Table 1: cut-offs on the log2 scale and back-transformed (2^x, natural units)
# Table 2: N women per tertile and per quartile. If more than 1/3 (or 1/4) of
#          the women share the lowest (LOD-imputed) value, quantize merges the
#          tied groups: n_groups_tert < 3 / n_groups_quart < 4, and the last
#          group column(s) are 0
# Arguments (optional):
#   1) cohort:   sepages / bisc   (default sepages)
#   2) analysis: main / s1 ... s5 (default main)
#   3) sample:   2tp / tp         (default 2tp; BiSC: tp)
#   Rscript script/8_Mixtures/1_Plot_Mixture_Quantiles.R sepages main 2tp

# [INPUT]:  db/subsets/df_<analysis>.rds (2tp) or df_<analysis>_<tp>.rds (tp)
# [OUTPUT]: results/<cohort>/8_Mixtures/1_Quantile_Cutoffs/<cohort>_<sample>_<analysis>/
#           hist_cutoffs_<tp>.png, solo/hist_cutoffs_<compound>_<tp>.png,
#           table_cutoffs.csv/.png, table_n_per_group.csv/.png
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_paths.R")
source("script/00_template_plot.R")
library(qgcomp)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
cohort     <- if (length(args) >= 1) args[1] else "sepages"
analysis   <- if (length(args) >= 2) args[2] else "main"
sample_set <- if (length(args) >= 3) args[3] else "2tp"

check_cohort(cohort, c("sepages", "bisc"))
stopifnot(analysis %in% c("main", paste0("s", 1:5)))
stopifnot(sample_set %in% c("2tp", "tp"))

time_points <- c("t2", "t3", "avg")
q_values    <- c(3, 4)    # 3 = tertiles, 4 = quartiles

# Same subsets as 4_Run_Gcomp.R
subset_dir <- db_subset_dir(cohort, sample_set)
extra_cov  <- if (cohort == "bisc") "hospital" else character(0)

# Same covariates as 4_Run_Gcomp.R (complete cases are taken on these too)
main_cov   <- c("mo_age", "mo_edu", "mo_smoke", "gest_age",
                "ch_eth", "ch_sex", "Meanlog2oddsContamination", extra_cov)
cell_types <- c("Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC")

comp_code   <- c("3pba",  "tcpy", "dmp",  "dep")
comp_labels <- c("3-PBA", "TCPy", "ΣDMP", "ΣDEP")
comp_units  <- c("ng/mL", "ng/mL", "µmol/L", "µmol/L")

outdir <- results_dir(cohort, "8_Mixtures", "1_Quantile_Cutoffs",
                      paste(cohort, sample_set, analysis, sep = "_"))
dir.create(file.path(outdir, "solo"), recursive = TRUE, showWarnings = FALSE)

#==============================================================================#----
# CUT-OFFS PER TIME POINT                                                      #----
#==============================================================================#----
cutoff_tab <- NULL   # Table 1
n_tab      <- NULL   # Table 2

for (tp in time_points) {
  
  cat("\n=== Time point:", tp, "\n")
  
  mix_exposures <- paste0(c("pest_3pba_ilog2", "pest_tcpy_ilog2",
                            "pest_dmp_ms_ilog2", "pest_dep_ms_ilog2"), "_", tp)
  covars <- list(main = c(main_cov, cell_types),
                 s1   = c(main_cov, cell_types),
                 s2   = c(setdiff(main_cov, "ch_eth"), cell_types),
                 s3   = c(main_cov, cell_types, paste0("pool_sg_", tp)),
                 s4   = c(main_cov, cell_types, "mo_fv_g1000kcal"),
                 s5   = setdiff(main_cov, "gest_age"))[[analysis]]
  
  file_suffix <- if (sample_set == "2tp") "" else paste0("_", tp)
  subset_path <- paste0(subset_dir, "df_", analysis, file_suffix, ".rds")
  stopifnot(file.exists(subset_path))
  metadata <- readRDS(subset_path)
  
  setdiff(c(mix_exposures, covars), colnames(metadata))   # should be character(0)
  stopifnot(all(c(mix_exposures, covars) %in% colnames(metadata)))
  
  # Analytic sample as in 4_Run_Gcomp.R: complete cases for compounds + covariates
  metadata <- metadata[complete.cases(metadata[, c(mix_exposures, covars)]), ]
  n_women  <- nrow(metadata)
  cat("Women (complete cases):", n_women, "\n")
  
  plot_list <- list()
  
  for (k in seq_along(comp_code)) {
    
    var_k <- mix_exposures[k]
    x_k   <- as.numeric(metadata[[var_k]])
    
    # Cut-offs exactly as qgcomp::quantize: inner breaks only (the outer ones are
    # -Inf/Inf or -1e64/1e64 depending on the qgcomp version)
    breaks_q3 <- quantize(metadata, expnms = var_k, q = 3)$breaks[[1]]
    breaks_q4 <- quantize(metadata, expnms = var_k, q = 4)$breaks[[1]]
    cut_q3 <- breaks_q3[-c(1, length(breaks_q3))]
    cut_q4 <- breaks_q4[-c(1, length(breaks_q4))]
    cat(comp_labels[k], "| tertile cut-offs:", round(cut_q3, 2),
        "| quartile cut-offs:", round(cut_q4, 2), "\n")
    
    # Percentile label of each cut-off that survives (tied percentiles are merged by quantize)
    q3_full <- quantile(x_k, probs = seq(0, 1, by = 1 / 3))
    q4_full <- quantile(x_k, probs = seq(0, 1, by = 1 / 4))
    q3_keep <- q3_full[!duplicated(q3_full)]
    q4_keep <- q4_full[!duplicated(q4_full)]
    lab_q3  <- paste0("P", round(as.numeric(sub("%", "", names(q3_keep)))))[-c(1, length(q3_keep))]
    lab_q4  <- paste0("P", round(as.numeric(sub("%", "", names(q4_keep)))))[-c(1, length(q4_keep))]
    stopifnot(all.equal(unname(q3_keep[-c(1, length(q3_keep))]), cut_q3),
              all.equal(unname(q4_keep[-c(1, length(q4_keep))]), cut_q4))   # same as quantize
    if (length(cut_q3) < 2) warning(var_k, ": fewer than 3 tertile groups (ties)")
    if (length(cut_q4) < 3) warning(var_k, ": fewer than 4 quartile groups (ties)")
    
    # Group sizes (same coding as quantize: 0, 1, ..., q-1)
    grp_q3 <- quantize(metadata, expnms = var_k, q = 3)$data[[var_k]]
    grp_q4 <- quantize(metadata, expnms = var_k, q = 4)$data[[var_k]]
    n_q3 <- table(factor(grp_q3, levels = 0:2))
    n_q4 <- table(factor(grp_q4, levels = 0:3))
    stopifnot(sum(n_q3) == n_women, sum(n_q4) == n_women)
    
    # Share of women at the lowest value (LOD-imputed values are tied)
    pct_min <- round(100 * mean(x_k == min(x_k)), 1)
    
    # Table 1: one row per cut-off
    cutoff_tab <- rbind(cutoff_tab, data.frame(
      time_point = tp, compound = comp_labels[k], n = n_women,
      quantiles  = c(rep("tertiles", length(cut_q3)), rep("quartiles", length(cut_q4))),
      cutoff     = c(lab_q3, lab_q4),
      log2_value = round(c(cut_q3, cut_q4), 3),
      natural_value = signif(2^c(cut_q3, cut_q4), 3),
      unit       = comp_units[k],
      stringsAsFactors = FALSE))
    
    # Table 2: one row per compound and time point
    n_tab <- rbind(n_tab, data.frame(
      time_point = tp, compound = comp_labels[k], n = n_women,
      pct_at_min = pct_min,
      n_groups_tert = length(cut_q3) + 1, n_groups_quart = length(cut_q4) + 1,   # < 3 / < 4: ties merged groups
      T1 = n_q3[[1]], T2 = n_q3[[2]], T3 = n_q3[[3]],
      Q1 = n_q4[[1]], Q2 = n_q4[[2]], Q3 = n_q4[[3]], Q4 = n_q4[[4]],
      stringsAsFactors = FALSE))
    
    # Plot: histogram + tertile (solid) and quartile (dashed) cut-offs
    df_k <- data.frame(x = x_k)
    df_lines <- data.frame(value = c(cut_q3, cut_q4),
                           type  = c(rep("Tertile cut-off", length(cut_q3)),
                                     rep("Quartile cut-off", length(cut_q4))))
    df_lines$type <- factor(df_lines$type, levels = c("Tertile cut-off", "Quartile cut-off"))
    
    p_k <- ggplot(df_k, aes(x = x)) +
      geom_histogram(bins = 30, fill = "grey80", colour = "grey50", linewidth = 0.2) +
      geom_vline(data = df_lines, aes(xintercept = value, linetype = type, colour = type),
                 linewidth = 0.7) +
      scale_linetype_manual(values = c("Tertile cut-off" = "solid", "Quartile cut-off" = "dashed"),
                            name = NULL) +
      scale_colour_manual(values = pal_cutoff,
                          name = NULL) +
      labs(title = comp_labels[k],
           subtitle = paste0("N = ", n_women, " | ", pct_min, "% at the lowest value | groups: ",
                             length(cut_q3) + 1, " tertile, ", length(cut_q4) + 1, " quartile"),
           x = paste0("log2 concentration (", comp_units[k], ", LOD-imputed)"),
           y = "Number of women") +
      theme(legend.position = "bottom")
    
    plot_list[[k]] <- p_k
    save_plot(p_k, file.path(outdir, "solo", paste0("hist_cutoffs_", comp_code[k], "_", tp, ".png")),
              width = 14, height = 12)
  }
  
  # Plot 1: 4 compounds of one time point (2 x 2)
  tp_title <- c(t2 = "second trimester", t3 = "third trimester",
                avg = "average of second and third trimester")[[tp]]
  p_panel <- wrap_plots(plot_list, ncol = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(title = paste0(toupper(cohort), " - mixture compounds, ", tp_title),
                    subtitle = paste0("Tertile and quartile cut-offs (qgcomp::quantize), analysis: ",
                                      analysis)) &
    theme(legend.position = "bottom")
  save_plot(p_panel, file.path(outdir, paste0("hist_cutoffs_", tp, ".png")),
            width = 24, height = 20)
}

#==============================================================================#----
# SAVE TABLES                                                                  #----
#==============================================================================#----
nrow(cutoff_tab)   # check: 3 time points x 4 compounds x 5 cut-offs = 60 (fewer if ties)
nrow(n_tab)        # check: 3 x 4 = 12
n_tab

save_table(cutoff_tab, file.path(outdir, "table_cutoffs"), base_size = 8)
save_table(n_tab,      file.path(outdir, "table_n_per_group"), base_size = 9)
cat("\nSaved in:", outdir, "\n")

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----