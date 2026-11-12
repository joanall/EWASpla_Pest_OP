#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the
# BiSC and SEPAGES cohorts.
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: OBJECTIVE 1 - Overview figure of the number of suggestive CpGs
# (p < 1e-5) per exposure in the main models, one panel per time point.
# - Plot:  horizontal stacked bars per exposure, split by direction (higher /
#          lower methylation with higher exposure); panel header with total
#          suggestive CpGs, n passing Bonferroni and median |effect| of the
#          suggestive CpGs at that time point.
# - Table: counts per exposure x time point (+ totals per time point).
# - ΣDAP (composite of ΣDMP + ΣDEP) is reported in the table only, not plotted.
# Self-contained: counts are taken from the annotated main results.
# Arguments: 1) cohort (bisc / sepages / meta / bisc_mix / sepages_mix / meta_mix), 2) pattern for the runs (one cohort, sample and run date), e.g.
#   Rscript script/6_Analyze_DMPs/3d_Get_Suggestive_Counts_Barplot.R sepages "sepages_2tp_t[23]_20261021"
#   Rscript script/6_Analyze_DMPs/3d_Get_Suggestive_Counts_Barplot.R sepages "sepages_2tp_(t2|t3|avg)_20261021"
#
# (<cohort> = bisc / sepages; meta = results/7_Meta/)
# [INPUT]:  results/<cohort>/5_Annotation/files_annotated_index.csv
#           -> results/<cohort>/5_Annotation/<run>/main/<exposure>_<tp>_annotated.rds
# [OUTPUT]: results/<cohort>/6_Analyze_DMPs/summary_tables/<cohort>_<sample>_<date>/
#           barplot_suggestive_counts_<tps>.png, table_suggestive_counts_<tps>.xlsx
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_paths.R")

args <- commandArgs(trailingOnly = TRUE)
cohort      <- if (length(args) >= 1) args[1] else "bisc"
run_pattern <- if (length(args) >= 2 && args[2] != "") args[2] else "."
check_cohort(cohort)

source("script/00_template_plot.R")

index_file  <- results_dir(cohort, "5_Annotation", "files_annotated_index.csv")
cutoff_sugg <- 1e-5    # suggestive threshold
cutoff_bn   <- 0.05    # Bonferroni-adjusted p threshold

# run_pattern <- "sepages_2tp_t[23]_20261021"   # TEST

files_index <- read.csv(index_file, stringsAsFactors = FALSE)
files_index <- files_index[grepl(run_pattern, files_index$run) & files_index$analysis == "main", ]
stopifnot(nrow(files_index) > 0)
stopifnot(all(file.exists(files_index$path)))

# Check: one cohort, sample and run date (counts must not mix runs)
stopifnot(length(unique(files_index$cohort)) == 1,
          length(unique(files_index$sample)) == 1,
          length(unique(files_index$run_date)) == 1)

# Check: one file per exposure x time point
table(files_index$exposure, files_index$time_point)
stopifnot(all(table(files_index$exposure, files_index$time_point) == 1))

tag    <- paste(files_index$cohort[1], files_index$sample[1], files_index$run_date[1], sep = "_")
outdir <- file.path(results_dir(cohort, "6_Analyze_DMPs", "summary_tables"), tag)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Exposure labels and order (top to bottom in the plot)
exp_code   <- c("24d",   "3pba",  "tcpy", "dmp",  "dep",  "dap",  "mix")
exp_labels <- c("2,4-D", "3-PBA", "TCPy", "ΣDMP", "ΣDEP", "ΣDAP", "Mixture")
stopifnot(all(files_index$exposure %in% exp_code))

# Exposures counted in the per-time-point totals.
exp_in_totals <- c("24d", "3pba", "tcpy", "dmp", "dep")
stopifnot(all(exp_in_totals %in% exp_code))

# Time points present, in order, and panel titles
tp_order    <- c("t2", "t3", "avg")
tp_titles   <- c(t2 = "Second trimester", t3 = "Third trimester", avg = "Pregnancy average")
time_points <- tp_order[tp_order %in% files_index$time_point]
time_points

# Direction colours (same as volcano / Manhattan plots)
pal_dir <- c("Higher methylation" = pal_direction[["Hypermethylated"]],
             "Lower methylation"  = pal_direction[["Hypomethylated"]])

#==============================================================================#----
# 1. COUNTS OF SUGGESTIVE CpGs PER EXPOSURE x TIME POINT                       #----
#==============================================================================#----

counts   <- NULL
eff_sugg <- list()   # effects of the suggestive CpGs, pooled per time point (excl. ΣDAP)
cpg_sugg <- list()   # CpG IDs of the suggestive CpGs, pooled per time point (excl. ΣDAP)

for (i in seq_len(nrow(files_index))) {
  f   <- files_index[i, ]
  res <- readRDS(f$path)
  sugg <- res[!is.na(res$pvalue) & res$pvalue < cutoff_sugg, ]
  
  counts <- rbind(counts, data.frame(
    exposure          = f$exposure,
    time_point        = f$time_point,
    n_cpgs            = nrow(res),
    n_sugg            = nrow(sugg),
    n_higher          = sum(sugg$Effect > 0, na.rm = TRUE),
    n_lower           = sum(sugg$Effect < 0, na.rm = TRUE),
    n_bonferroni      = sum(!is.na(sugg$pval_bn) & sugg$pval_bn < cutoff_bn),
    median_abs_effect = if (nrow(sugg) > 0) median(abs(sugg$Effect), na.rm = TRUE) else NA,
    stringsAsFactors  = FALSE))
  if (f$exposure %in% exp_in_totals) {
    eff_sugg[[f$time_point]] <- c(eff_sugg[[f$time_point]], sugg$Effect)
    cpg_sugg[[f$time_point]] <- c(cpg_sugg[[f$time_point]], sugg$CpG)
  }
  cat(f$exposure, f$time_point, "| suggestive:", nrow(sugg),
      "| Bonferroni:", counts$n_bonferroni[nrow(counts)], "\n")
}

# Check: higher + lower adds up to the suggestive CpGs (no Effect == 0 / NA)
stopifnot(all(counts$n_higher + counts$n_lower == counts$n_sugg))

counts <- counts[order(match(counts$time_point, time_points), match(counts$exposure, exp_code)), ]
counts$exposure_label <- exp_labels[match(counts$exposure, exp_code)]

# Totals per time point (panel header)
tp_summary <- NULL
for (tp in time_points) {
  c_tp <- counts[counts$time_point == tp & counts$exposure %in% exp_in_totals, ]
  c_tp_dap <- counts[counts$time_point == tp & counts$exposure == "dap", ]
  tp_summary <- rbind(tp_summary, data.frame(
    time_point        = tp,
    n_sugg            = sum(c_tp$n_sugg),
    n_unique_cpgs     = length(unique(cpg_sugg[[tp]])),
    n_higher          = sum(c_tp$n_higher),
    n_lower           = sum(c_tp$n_lower),
    n_bonferroni      = sum(c_tp$n_bonferroni),
    n_bonferroni_dap  = sum(c_tp_dap$n_bonferroni),
    median_abs_effect = if (length(eff_sugg[[tp]]) > 0) median(abs(eff_sugg[[tp]]), na.rm = TRUE) else NA,
    stringsAsFactors  = FALSE))
  # Check: pooled vectors match the summed counts
  stopifnot(length(cpg_sugg[[tp]]) == sum(c_tp$n_sugg),
            length(eff_sugg[[tp]]) == sum(c_tp$n_sugg))
}
tp_summary

#==============================================================================#----
# TABLE: COUNTS                                                                #----
#==============================================================================#----

tps_tag <- paste(time_points, collapse = "_")
save_table(list(per_exposure   = counts,
                per_time_point = tp_summary),
           file.path(outdir, paste0("table_suggestive_counts_", tps_tag)))

#==============================================================================#----
# PLOT: BARPLOT OF SUGGESTIVE CpGs                                             #----
#==============================================================================#----

# Panel header: time point + totals
tp_summary$panel <- paste0(
  tp_titles[tp_summary$time_point], "\n",
  tp_summary$n_sugg, " suggestive CpGs (", tp_summary$n_unique_cpgs, " unique) | ",
  ifelse(tp_summary$n_bonferroni == 0, "none", tp_summary$n_bonferroni), " passed Bonferroni | ",
  "median |effect| ", signif(tp_summary$median_abs_effect, 2))
counts$panel <- tp_summary$panel[match(counts$time_point, tp_summary$time_point)]

# Plot only the exposures counted in the totals (ΣDAP kept in the table, not plotted)
counts_plot <- counts[counts$exposure %in% exp_in_totals, ]
stopifnot(!any(counts_plot$exposure == "dap"))

# Long format: one row per exposure x time point x direction
df_bar <- rbind(data.frame(panel = counts_plot$panel, exposure = counts_plot$exposure_label,
                           direction = "Higher methylation", n = counts_plot$n_higher),
                data.frame(panel = counts_plot$panel, exposure = counts_plot$exposure_label,
                           direction = "Lower methylation",  n = counts_plot$n_lower))
df_bar$panel     <- factor(df_bar$panel, levels = tp_summary$panel)
df_bar$exposure  <- factor(df_bar$exposure, levels = rev(exp_labels[exp_labels %in% df_bar$exposure]))
df_bar$direction <- factor(df_bar$direction, levels = names(pal_dir))

# Totals at the end of each bar
df_tot <- data.frame(panel    = factor(counts_plot$panel, levels = tp_summary$panel),
                     exposure = factor(counts_plot$exposure_label, levels = levels(df_bar$exposure)),
                     n_sugg   = counts_plot$n_sugg)

# Check: plotted counts add up to the header totals per time point
n_plotted <- tapply(df_bar$n, df_bar$panel, sum)
n_plotted
stopifnot(all(n_plotted[tp_summary$panel] == tp_summary$n_sugg))

# Plot 1. Suggestive CpGs per exposure (y) and number of CpGs (x), stacked by
# direction (colour: higher / lower methylation); one panel per time point,
# header with totals; number at the end of each bar = total suggestive CpGs.
p_bar <- ggplot(df_bar, aes(x = n, y = exposure, fill = direction)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.65) +
  geom_text(data = df_tot, aes(x = n_sugg, y = exposure, label = n_sugg),
            inherit.aes = FALSE, hjust = -0.4, size = 3.2, fontface = "bold") +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_manual(values = pal_dir) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title    = "Prenatal pesticide exposure and placental DNA methylation",
       subtitle = paste0(tag, " | main models, suggestive CpGs per biomarker"),
       x = "Number of suggestive CpGs", y = NULL, fill = NULL) +
  theme_ewas() +
  theme(axis.text.x     = element_text(angle = 0, hjust = 0.5),
        strip.text      = element_text(hjust = 0, lineheight = 1.1),
        legend.position = "bottom")

save_plot(p_bar, file.path(outdir, paste0("barplot_suggestive_counts_", tps_tag, ".png")),
          width = 10 * length(time_points), height = 10)
cat("Saved barplot_suggestive_counts:", tps_tag, "\n")

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----