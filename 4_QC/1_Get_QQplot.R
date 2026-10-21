#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: Genomic inflation (lambda, 95% CI) and QQ plot for every
# cleaned EWAS result file. Each QQ plot is saved as .png and as a ggplot
# object (.rds), so 5_Join_QQ_Plots.R can combine any set of them.
#
# Argument (optional): pattern for the run folders to include (regular
# expression; default: all)
#   Rscript script/4_QQ_Plots/4_QQ_Plots.R "sepages_2tp_t[23]_old"
#
# [INPUT]:  results/3_Clean_Results/files_cleaned_index.csv
#           (from 3_Get_Cleaned_Files.R: path, run, cohort, sample,
#            time_point, run_date, analysis, exposure)
# [OUTPUT]: results/4_QQ_Plots/<run>/<analysis>/<exposure>_<tp>_qqplot.png
#           results/4_QQ_Plots/<run>/<analysis>/<exposure>_<tp>_qqplot.rds
#           results/4_QQ_Plots/qq_index.csv   (index + lambda + plot paths)
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_template_plot.R")
source("script/00_functions.R")

index_file <- "results/3_Clean_Results/files_cleaned_index.csv"
out_dir    <- "results/4_QQ_Plots"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
run_pattern <- if (length(args) >= 1 && args[1] != "") args[1] else "."

# Exposure labels for titles
exp_code   <- c("24d",   "3pba",  "tcpy", "dmp",  "dep",  "dap")
exp_labels <- c("2,4-D", "3-PBA", "TCPy", "ΣDMP", "ΣDEP", "ΣDAP")

# Files to plot
files_index <- read.csv(index_file, stringsAsFactors = FALSE)
files_index <- files_index[grepl(run_pattern, files_index$run), ]
cat("Files to plot:", nrow(files_index), "\n")
stopifnot(nrow(files_index) > 0)
stopifnot(all(file.exists(files_index$path)))
stopifnot(all(files_index$exposure %in% exp_code))

#==============================================================================#----
# LAMBDA & QQ PLOT PER FILE                                                    #----
#==============================================================================#----

qq_index <- files_index
qq_index$n_cpgs   <- NA
qq_index$lambda   <- NA
qq_index$lower_ci <- NA
qq_index$upper_ci <- NA
qq_index$png      <- NA
qq_index$rds      <- NA

for (i in seq_len(nrow(files_index))) {
  
  f <- files_index[i, ]
  cat("\n[", i, "/", nrow(files_index), "]", f$run, "|", f$analysis, "|", f$exposure, "\n")
  
  ewas_results <- readRDS(f$path)
  
  # Lambda and 95% CI
  lambda_ci <- inflation(ewas_results$pvalue)
  qq_index$n_cpgs[i]   <- sum(!is.na(ewas_results$pvalue))
  qq_index$lambda[i]   <- round(lambda_ci[["lambda"]],   3)
  qq_index$lower_ci[i] <- round(lambda_ci[["lower_ci"]], 3)
  qq_index$upper_ci[i] <- round(lambda_ci[["upper_ci"]], 3)
  cat("Lambda:", qq_index$lambda[i], "[", qq_index$lower_ci[i], ",", qq_index$upper_ci[i], "]\n")
  
  # QQ plot: observed vs expected -log10 p-values; lambda [95% CI] in the corner
  exp_label <- exp_labels[match(f$exposure, exp_code)]
  plot_title <- paste0(exp_label, " - ", f$time_point,
                       ifelse(f$analysis == "main", "", paste0(" (", f$analysis, ")")))
  p_qq <- qq_plot(ewas_results$pvalue, title = plot_title)
  
  # Save: same <run>/<analysis>/ structure as the results
  plot_path <- file.path(out_dir, f$run, f$analysis)
  dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
  base_name <- paste0(f$exposure, "_", f$time_point, "_qqplot")
  png_file  <- file.path(plot_path, paste0(base_name, ".png"))
  rds_file  <- file.path(plot_path, paste0(base_name, ".rds"))
  
  save_plot(p_qq, png_file, width = 12, height = 12)
  saveRDS(p_qq, rds_file)
  
  qq_index$png[i] <- png_file
  qq_index$rds[i] <- rds_file
}

#==============================================================================#----
# SAVE INDEX + LAMBDA TABLE                                                    #----
#==============================================================================#----

qq_index[, c("run", "analysis", "exposure", "n_cpgs", "lambda", "lower_ci", "upper_ci")]
write.csv(qq_index, file.path(out_dir, "qq_index.csv"), row.names = FALSE)
cat("Index written:", file.path(out_dir, "qq_index.csv"), "\n")

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----