#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: October - November 2026
#==============================================================================#

# [DESCRIPTION]: Manhattan plot for each annotated EWAS result (main +
# sensitivity). Suggestive CpGs (p < 1e-5) coloured by direction and labelled
# as "CpG (gene)"; dashed line = Bonferroni, dotted = suggestive.
# Argument (optional): pattern for the runs to plot (default: all), e.g.
#   Rscript script/6_Analyze_Results/5_Get_Manhattan_Plot.R "sepages_2tp_t[23]_20261021"
#
# [INPUT]:  results/5_Annotation/files_annotated_index.csv
#           -> results/5_Annotation/<run>/<analysis>/<exposure>_<tp>_annotated.rds
#           (chr and pos already in the annotated files)
# [OUTPUT]: results/6_Analyze_Results/<run>/<analysis>/<exposure>_<tp>_manhplot.png
#           results/6_Analyze_Results/<run>/<analysis>/<exposure>_<tp>_manhplot.rds (ggplot object)

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_template_plot.R")
source("script/00_functions.R")

index_file    <- "results/5_Annotation/files_annotated_index.csv"
output_folder <- "results/6_Analyze_Results/"
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

cutoff_sugg <- 1e-5    # suggestive threshold: dotted line, coloured and labelled CpGs

# Select runs (argument: pattern, default all)
args <- commandArgs(trailingOnly = TRUE)
run_pattern <- if (length(args) >= 1 && args[1] != "") args[1] else "."
# run_pattern <- "sepages_2tp_t2_20261021"   # TEST - one run

files_index <- read.csv(index_file, stringsAsFactors = FALSE)
files_index <- files_index[grepl(run_pattern, files_index$run), ]
# files_index <- files_index[1, ]            # TEST - one file
cat("Files to plot:", nrow(files_index), "\n")
stopifnot(nrow(files_index) > 0)

files_index$path[!file.exists(files_index$path)]   # should be character(0)
stopifnot(all(file.exists(files_index$path)))

#==============================================================================#----
# MANHATTAN PLOTS                                                              #----
#==============================================================================#----

for (i in seq_len(nrow(files_index))) {
  
  f <- files_index[i, ]
  cat("\n[", i, "/", nrow(files_index), "]", f$run, "|", f$analysis, "|", f$exposure, "\n")
  
  ewas_results <- readRDS(f$path)
  
  # Positions from the annotated file itself (rownames = CpGs, chr, pos)
  anno_pos <- data.frame(chr = ewas_results$chr, pos = ewas_results$pos,
                         row.names = ewas_results$CpG)
  
  # Check: CpGs with position, suggestive and Bonferroni-significant CpGs
  cat("  CpGs with position:", sum(!is.na(anno_pos$pos)), "of", nrow(ewas_results), "\n")
  cat("  Suggestive p <", cutoff_sugg, ":", sum(ewas_results$pvalue < cutoff_sugg, na.rm = TRUE),
      "| Bonferroni < 0.05:", sum(ewas_results$pval_bn < 0.05, na.rm = TRUE), "\n")
  
  # Plot 1. Manhattan plot for one exposure, time point and analysis: genomic
  # position (x) vs -log10(p) (y); suggestive CpGs coloured by direction and
  # labelled "CpG (gene)".
  p_manh <- manhattan_plot(ewas_results, anno_pos,
                           title = paste0(f$exposure, " - ", f$time_point, " (", f$analysis, ")"),
                           suggestive = cutoff_sugg, colour_sugg = TRUE, annot = TRUE)
  
  # Save: same <run>/<analysis>/ structure as the results
  plot_path <- file.path(output_folder, f$run, f$analysis)
  dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
  base_name <- paste0(f$exposure, "_", f$time_point, "_manhplot")
  
  save_plot(p_manh, file.path(plot_path, paste0(base_name, ".png")), width = 24, height = 10)
  saveRDS(p_manh, file.path(plot_path, paste0(base_name, ".rds")))
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----