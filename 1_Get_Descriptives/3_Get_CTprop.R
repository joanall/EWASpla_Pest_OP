#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the
# BiSC and SEPAGES cohorts.
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org
# Date: November 2026
#==============================================================================#
# [DESCRIPTION]: Placental cell type proportions (estimated from DNAm) in the
# analytic sample(s) of one cohort:
# 1) Stacked bar plot: one bar per sample, the 6 cell types stacked (sum = 1),
#    samples ordered by Syncytiotrophoblast proportion (largest cell type)
# 2) Table: mean, SD, median, P5, P95, min and max of each cell type
# Analytic samples:
# - SEPAGES: main (women with both time points)
# - BiSC:    main t2 / main t3 / both time points / average (>= 1 time point)
#
# Argument: cohort (sepages / bisc); run on the PC of that cohort
#   Rscript script/1_Get_Descriptives/3_Get_CT_prop.R sepages
#   Rscript script/1_Get_Descriptives/3_Get_CT_prop.R bisc
#
# [INPUT SEPAGES]: db/subsets/df_main.rds
# [INPUT BiSC]:    db/subsets/df_main_t2.rds, df_main_t3.rds, df_main_avg.rds;
#                  db/subsets_2tp/df_main.rds
# [OUTPUT]: results/<cohort>/1_Get_Descriptives/cell_types/
#           ct_prop_<cohort>_<sample>.png   (one plot per analytic sample)
#           ct_prop_summary_<cohort>.xlsx   (one sheet per analytic sample)
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_paths.R")
source("script/00_template_plot.R")   # ggplot2, theme_ewas(), palettes, save_plot(), save_table()

args   <- commandArgs(trailingOnly = TRUE)
cohort <- if (length(args) >= 1) args[1] else "bisc"
check_cohort(cohort, c("sepages", "bisc"))
cohort_label <- cohort_labels[[cohort]]

outdir <- results_dir(cohort, "1_Get_Descriptives", "cell_types", create = TRUE)

# Cell types, in stacking order (bottom to top: Syncytiotrophoblast on top)
cell_types <- c("Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC",
                "Syncytiotrophoblast")

# Colours: pal_cell_type from 00_template_plot.R (same colours in every figure)
stopifnot(exists("pal_cell_type"), all(cell_types %in% names(pal_cell_type)))

#==============================================================================#----
# DATA                                                                         #----
#==============================================================================#----
# Analytic samples of this cohort (paths from 00_paths.R)
if (cohort == "sepages") {
  samples_list <- list(main = readRDS(db_subset("sepages", "main", "2tp")))
  sample_labels <- c(main = "main sample (both time points)")
}
if (cohort == "bisc") {
  samples_list <- list(t2  = readRDS(db_subset("bisc", "main", "tp", "t2")),
                       t3  = readRDS(db_subset("bisc", "main", "tp", "t3")),
                       "2tp" = readRDS(db_subset("bisc", "main", "2tp")),
                       avg = readRDS(db_subset("bisc", "main", "tp", "avg")))
  sample_labels <- c(t2    = "main sample, second trimester",
                     t3    = "main sample, third trimester",
                     "2tp" = "main sample, both time points",
                     avg   = "main sample, average (>= 1 time point)")
}

# Check: cell type columns present and N per sample
for (s in names(samples_list)) {
  stopifnot(all(c("Basename", cell_types) %in% colnames(samples_list[[s]])))
  cat(s, "| N =", nrow(samples_list[[s]]), "\n")
}

#==============================================================================#----
# PLOT AND SUMMARY PER SAMPLE                                                  #----
#==============================================================================#----
summary_list <- list()

for (s in names(samples_list)) {
  ct_wide <- samples_list[[s]][, c("Basename", cell_types)]
  ct_wide <- ct_wide[complete.cases(ct_wide), ]
  n_samples <- nrow(ct_wide)
  
  # Check: the 6 proportions sum to ~1 in every sample
  ct_sum <- rowSums(ct_wide[, cell_types])
  cat("\n", s, "| N =", n_samples, "| sum of proportions: min", round(min(ct_sum), 3),
      "max", round(max(ct_sum), 3), "\n")
  
  # Order samples by Syncytiotrophoblast proportion (lowest first)
  sample_order <- ct_wide$Basename[order(ct_wide$Syncytiotrophoblast)]
  
  # Long format (base R): one row per sample x cell type
  ct_long <- data.frame(
    Basename   = rep(ct_wide$Basename, times = length(cell_types)),
    cell_type  = rep(cell_types, each = n_samples),
    proportion = unlist(ct_wide[, cell_types], use.names = FALSE),
    stringsAsFactors = FALSE
  )
  ct_long$Basename  <- factor(ct_long$Basename, levels = sample_order)
  ct_long$cell_type <- factor(ct_long$cell_type, levels = rev(cell_types))   # Syncytiotrophoblast on top
  stopifnot(nrow(ct_long) == n_samples * length(cell_types))
  
  p_ct <- ggplot(ct_long, aes(x = Basename, y = proportion, fill = cell_type)) +
    geom_col(width = 1) +
    scale_fill_manual(values = pal_cell_type, breaks = cell_types) +
    scale_y_continuous(limits = c(0, 1.0001), expand = c(0, 0)) +
    labs(title = paste0(cohort_label, " - cell type proportions, ", sample_labels[[s]],
                        " (N = ", n_samples, ")"),
         x = "Samples (ordered by Syncytiotrophoblast proportion)",
         y = "Proportion", fill = "Cell type") +
    theme_ewas() +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = "bottom")
  
  save_plot(p_ct, file.path(outdir, paste0("ct_prop_", cohort, "_", s, ".png")),
            width = 18, height = 12)
  
  # Summary table: one row per cell type
  ct_values <- ct_wide[, cell_types]
  summary_list[[s]] <- data.frame(
    cell_type = cell_types,
    N         = colSums(!is.na(ct_values)),
    mean      = round(colMeans(ct_values), 3),
    sd        = round(apply(ct_values, 2, sd), 3),
    median    = round(apply(ct_values, 2, median), 3),
    p5        = round(apply(ct_values, 2, quantile, probs = 0.05), 3),
    p95       = round(apply(ct_values, 2, quantile, probs = 0.95), 3),
    min       = round(apply(ct_values, 2, min), 3),
    max       = round(apply(ct_values, 2, max), 3),
    row.names = NULL
  )
  print(summary_list[[s]])
}

# One .xlsx, one sheet per analytic sample
save_table(summary_list, file.path(outdir, paste0("ct_prop_summary_", cohort)))
cat("\nSaved in", outdir, "\n")

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----