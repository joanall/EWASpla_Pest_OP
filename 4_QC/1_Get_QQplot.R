#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: Genomic inflation (lambda, 95% CI) and QQ plot for every
# cleaned EWAS result file. Each QQ plot is saved as .png and as a ggplot
# object (.rds), so 2_Join_QQplots.R can combine any set of them.
#
# Arguments:
#   1) cohort: bisc / sepages / meta / bisc_mix / sepages_mix / meta_mix
#   2) pattern for the run folders to include (regular expression; default: all)
#   3) analyses to include, comma-separated (default: all), e.g. "s6"
#   Rscript script/4_QC/1_Get_QQplot.R sepages "sepages_2tp_t[23]_old"
#   Rscript script/4_QC/1_Get_QQplot.R bisc "^bisc_tp_(t2|t3|avg)_20261022$" s6
#   Rscript script/4_QC/1_Get_QQplot.R meta "meta_tp_(t2|t3|avg)_20261106"
# qq_index.csv is updated, not overwritten: rows of the files plotted now are
# replaced, rows of other files are kept.
#
# [INPUT]:  <clean dir>/files_cleaned_index.csv (from 3_Get_Cleaned_Files.R:
#           path, run, cohort, sample, time_point, run_date, analysis, exposure)
#           <clean dir> = results/<cohort>/3_Clean_Results or results/7_Meta/3_Read_GWAMA
# [OUTPUT]: results/<cohort>/4_QC/<run>/<analysis>/<exposure>_<tp>_qqplot.png
#           results/<cohort>/4_QC/<run>/<analysis>/<exposure>_<tp>_qqplot.rds
#           results/<cohort>/4_QC/qq_index.csv   (index + lambda + plot paths)
#           (meta: results/7_Meta/4_QC/...)
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_paths.R")
source("script/00_template_plot.R")
source("script/00_functions.R")

args <- commandArgs(trailingOnly = TRUE)
cohort      <- if (length(args) >= 1) args[1] else "bisc"
run_pattern <- if (length(args) >= 2 && args[2] != "") args[2] else "."
analyses    <- if (length(args) >= 3 && args[3] != "") strsplit(args[3], ",")[[1]] else NULL
check_cohort(cohort)

index_file <- file.path(clean_dir(cohort), "files_cleaned_index.csv")
out_dir    <- results_dir(cohort, "4_QC", create = TRUE)

# Exposure labels for titles
exp_code   <- c("24d",   "3pba",  "tcpy", "dmp",  "dep",  "dap",  "mix")
exp_labels <- c("2,4-D", "3-PBA", "TCPy", "ΣDMP", "ΣDEP", "ΣDAP", "Mixture")

# Files to plot
files_index <- read.csv(index_file, stringsAsFactors = FALSE)
files_index <- files_index[grepl(run_pattern, files_index$run), ]
if (!is.null(analyses)) {
  setdiff(analyses, unique(files_index$analysis))   # should be character(0)
  files_index <- files_index[files_index$analysis %in% analyses, ]
}
table(files_index$run, files_index$analysis)       # check: files per run x analysis
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

# Update the existing index: drop old rows of the files plotted now, add the new ones
qq_index_file <- file.path(out_dir, "qq_index.csv")
if (file.exists(qq_index_file)) {
  qq_index_old <- read.csv(qq_index_file, stringsAsFactors = FALSE)
  n_old <- nrow(qq_index_old)
  qq_index_old <- qq_index_old[!qq_index_old$path %in% qq_index$path, ]
  cat("Index: kept", nrow(qq_index_old), "of", n_old, "old rows | added", nrow(qq_index), "\n")
  qq_index_old <- qq_index_old[, colnames(qq_index_old) %in% colnames(qq_index), drop = FALSE]
  for (col_new in setdiff(colnames(qq_index), colnames(qq_index_old))) qq_index_old[[col_new]] <- NA
  qq_index <- rbind(qq_index_old[, colnames(qq_index)], qq_index)
}
stopifnot(!anyDuplicated(qq_index$path))
table(qq_index$run, qq_index$analysis)   # check: all runs / analyses in the index
# Plain csv on purpose: read back by later scripts (not a results table)
write.csv(qq_index, qq_index_file, row.names = FALSE)
cat("Index written:", qq_index_file, "\n")

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----