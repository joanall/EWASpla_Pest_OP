#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the
# BiSC and SEPAGES cohorts.
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: Panel of Manhattan plots (main models): the 6 exposures of one
# time point (t2, t3, average) stacked in one figure (6 x 1), with a SHARED
# y-axis (0 to max -log10(p) across the 6 plots) so the exposures are directly
# comparable. The x-axis (genomic position) is already the same in all plots.
# Same as 2b_Get_Manhattan_Panel.R except for the shared y-axis.
# Uses the ggplot objects saved by the Manhattan plot script (run it first).
# Argument (optional): pattern for the runs (default: all), e.g.
#   Rscript script/6_Analyze_Results/5c_Get_Manhattan_Panel_Shared_Axes.R "sepages_2tp_t[23]_20261021"
#
# [INPUT]:  results/5_Annotation/files_annotated_index.csv
#           results/6_Analyze_Results/<run>/main/<exposure>_<tp>_manhplot.rds
# [OUTPUT]: results/6_Analyze_Results/<run>/main/manhattan_panel_main_<tp>_shared_axes.png
#==============================================================================#


#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_template_plot.R")
source("script/00_functions.R")
library(patchwork)

index_file    <- "results/5_Annotation/files_annotated_index.csv"
output_folder <- "results/6_Analyze_Results/"

# Select runs (argument: pattern, default all)
args <- commandArgs(trailingOnly = TRUE)
run_pattern <- if (length(args) >= 1 && args[1] != "") args[1] else "."
# run_pattern <- "sepages_2tp_t2_20261021"   # TEST

files_index <- read.csv(index_file, stringsAsFactors = FALSE)
files_index <- files_index[grepl(run_pattern, files_index$run) & files_index$analysis == "main", ]
stopifnot(nrow(files_index) > 0)

# Path of the single Manhattan plot (ggplot object) saved by the Manhattan script
files_index$plot_rds <- file.path(output_folder, files_index$run, files_index$analysis,
                                  paste0(files_index$exposure, "_", files_index$time_point,
                                         "_manhplot.rds"))
files_index$plot_rds[!file.exists(files_index$plot_rds)]   # should be character(0)
stopifnot(all(file.exists(files_index$plot_rds)))           # run the Manhattan script first

# Exposure labels and panel order
exp_code   <- c("24d",   "3pba",  "tcpy", "dmp",  "dep",  "dap")
exp_labels <- c("2,4-D", "3-PBA", "TCPy", "ΣDMP", "ΣDEP", "ΣDAP")

# One panel per run x time point
files_index$panel <- paste(files_index$run, files_index$time_point, sep = "|")
table(files_index$panel)   # 6 per panel


#==============================================================================#----
# MANHATTAN PLOTS - PANEL WITH SHARED Y-AXIS                                   #----
#==============================================================================#----

for (pn in unique(files_index$panel)) {
  
  f_panel <- files_index[files_index$panel == pn, ]
  f_panel <- f_panel[order(match(f_panel$exposure, exp_code)), ]
  run_i   <- f_panel$run[1]
  tp_i    <- f_panel$time_point[1]
  cat("\nPanel:", run_i, "|", tp_i, "| n plots:", nrow(f_panel), "\n")
  if (nrow(f_panel) != 6) warning("Panel ", pn, " has ", nrow(f_panel), " plots (expected 6)")
  
  # Read the single plots and the y range each one shows
  # (range in plot units: -log10(p), incl. the threshold lines)
  plot_list <- list()
  y_range   <- matrix(NA, nrow = nrow(f_panel), ncol = 2)
  for (j in seq_len(nrow(f_panel))) {
    p_single  <- readRDS(f_panel$plot_rds[j])
    exp_label <- exp_labels[match(f_panel$exposure[j], exp_code)]
    plot_list[[j]] <- p_single + labs(title = exp_label, subtitle = NULL)
    y_range[j, ]   <- layer_scales(p_single)$y$range$range
  }
  
  # Check: y range of each plot before sharing
  print(data.frame(exposure = f_panel$exposure, y_max = y_range[, 2]))
  
  # Shared limit: y from 0 to the highest -log10(p) of the 6 plots
  y_lim_shared <- max(y_range[, 2], na.rm = TRUE)
  cat("  Shared y max:", signif(y_lim_shared, 3), "\n")
  
  for (j in seq_along(plot_list)) {
    plot_list[[j]] <- plot_list[[j]] + coord_cartesian(ylim = c(0, y_lim_shared))
  }
  
  # Plot 1. Manhattan plots of the 6 exposures for one time point (main models),
  # shared y-axis: genomic position (x) vs -log10(p) (y); suggestive CpGs
  # coloured by direction and labelled; dotted = suggestive, dashed = Bonferroni.
  p_panel <- wrap_plots(plot_list, ncol = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(title = paste0("Main models - ", tp_i, " (shared y-axis)"),
                    subtitle = run_i)
  
  plot_path <- file.path(output_folder, run_i, "main")
  save_plot(p_panel, file.path(plot_path, paste0("manhattan_panel_main_", tp_i, "_shared_axes.png")),
            width = 24, height = 50)
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----