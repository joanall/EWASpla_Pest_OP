#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the
# BiSC and SEPAGES cohorts.
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: Panel of Manhattan plots (main models): the 6 exposures of one
# time point (t2, t3, average) stacked in one figure (6 x 1).
# Uses the ggplot objects saved by the Manhattan plot script (run it first).
# Argument (optional): pattern for the runs (default: all), e.g.
#   Rscript script/6_Analyze_Results/5b_Get_Manhattan_Panel.R "sepages_2tp_t[23]_20261021"
#
# [INPUT]:  results/5_Annotation/files_annotated_index.csv
#           results/6_Analyze_Results/<run>/main/<exposure>_<tp>_manhplot.rds
# [OUTPUT]: results/6_Analyze_Results/<run>/main/manhattan_panel_main_<tp>.png
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
# MANHATTAN PLOTS - PANEL (6 EXPOSURES PER TIME POINT)                         #----
#==============================================================================#----

for (pn in unique(files_index$panel)) {
  
  f_panel <- files_index[files_index$panel == pn, ]
  f_panel <- f_panel[order(match(f_panel$exposure, exp_code)), ]
  run_i   <- f_panel$run[1]
  tp_i    <- f_panel$time_point[1]
  cat("\nPanel:", run_i, "|", tp_i, "| n plots:", nrow(f_panel), "\n")
  if (nrow(f_panel) != 6) warning("Panel ", pn, " has ", nrow(f_panel), " plots (expected 6)")
  
  # Single plots: title = exposure only (time point and analysis go in the panel title)
  plot_list <- list()
  for (j in seq_len(nrow(f_panel))) {
    p_single  <- readRDS(f_panel$plot_rds[j])
    exp_label <- exp_labels[match(f_panel$exposure[j], exp_code)]
    plot_list[[j]] <- p_single + labs(title = exp_label, subtitle = NULL)
  }
  
  # Plot 1. Manhattan plots of the 6 exposures for one time point (main models):
  # genomic position (x) vs -log10(p) (y); suggestive CpGs coloured by direction
  # and labelled; dotted line = suggestive, dashed line = Bonferroni.
  p_panel <- wrap_plots(plot_list, ncol = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(title = paste0("Main models - ", tp_i),
                    subtitle = run_i)
  
  plot_path <- file.path(output_folder, run_i, "main")
  save_plot(p_panel, file.path(plot_path, paste0("manhattan_panel_main_", tp_i, ".png")),
            width = 24, height = 50)
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----