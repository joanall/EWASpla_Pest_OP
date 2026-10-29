#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: Volcano plot for each annotated EWAS result (main + sensitivity).
# Bonferroni-significant CpGs coloured by direction: hyper- vs hypomethylated.
# Argument (optional): pattern for the runs to plot (default: all), e.g.
#   Rscript script/6_Analyze_Results/1_Get_Volcano_Plot.R "sepages_2tp_t[23]_20261021"
#
# [INPUT]:  results/5_Annotation/files_annotated_index.csv
#           -> results/5_Annotation/<run>/<analysis>/<exposure>_<tp>_annotated.rds
# [OUTPUT]: results/6_Analyze_Results/<run>/<analysis>/<exposure>_<tp>_volcanoplot.png
#           results/6_Analyze_Results/<run>/<analysis>/<exposure>_<tp>_volcanoplot.rds (ggplot object)
#==============================================================================#


#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_template_plot.R")
source("script/00_functions.R")

list_file  <- "results/5_Annotation/files_annotated.txt"
index_file <- "results/5_Annotation/files_annotated_index.csv"
output_folder  <- "results/6_Analyze_Results/"
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# Select runs (argument: pattern, default all)
args <- commandArgs(trailingOnly = TRUE)
run_pattern <- if (length(args) >= 1 && args[1] != "") args[1] else "."
# files_index <- files_index[1, ]   # TEST - remove for the real run
run_pattern <- "sepages_2tp_t2_20261021"   # all 36 files of that run (6 analyses × 6 exposures)
files_index <- read.csv(index_file, stringsAsFactors = FALSE)
files_index <- files_index[grepl(run_pattern, files_index$run), ]
cat("Files to plot:", nrow(files_index), "\n")
stopifnot(nrow(files_index) > 0)

files_in <- files_index$path
files_in[!file.exists(files_in)]   # should be character(0)
stopifnot(all(file.exists(files_in)))

# Exposure labels for titles
exp_code   <- c("24d",   "3pba",  "tcpy", "dmp",  "dep",  "dap")
exp_labels <- c("2,4-D", "3-PBA", "TCPy", "ΣDMP", "ΣDEP", "ΣDAP")
cutoff_sugg <- 1E-05

#==============================================================================#----
# VOLCANO PLOTS - SINGLE                                                       #----
#==============================================================================#----

for (i in seq_len(nrow(files_index))) {
  
  f <- files_index[i, ]
  cat("\n[", i, "/", nrow(files_index), "]", f$run, "|", f$analysis, "|", f$exposure, "\n")
  
  ewas_results <- readRDS(f$path)
  exp_label    <- exp_labels[match(f$exposure, exp_code)]
  
  # Check: Bonferroni-significant CpGs
  cat("  Bonferroni < 0.05:", sum(ewas_results$pval_bn < 0.05, na.rm = TRUE), "\n")
  
  # Plot 1. Volcano plot for one exposure, time point and analysis: effect (x)
  # vs -log10(p) (y); Bonferroni CpGs coloured by direction.
  p_volcano <- volcano_plot(ewas_results,
                            title = paste0(exp_label, " - ", f$time_point, " (", f$analysis, ")"),
                            p_cutoff = cutoff_sugg)
  
  # Save: same <run>/<analysis>/ structure as the results
  plot_path <- file.path(output_folder, f$run, f$analysis)
  dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
  base_name <- paste0(f$exposure, "_", f$time_point, "_volcanoplot")
  
  save_plot(p_volcano, file.path(plot_path, paste0(base_name, ".png")), width = 14, height = 12)
  saveRDS(p_volcano, file.path(plot_path, paste0(base_name, ".rds")))
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----