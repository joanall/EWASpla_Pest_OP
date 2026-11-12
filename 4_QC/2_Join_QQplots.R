#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: Combine saved QQ plots (from 1_Get_QQplot.R) into multi-panel
# figures. Each figure is defined in FIGURES below: which plots go in it, the
# layout and the output name. Add or edit figures there.
#
# Arguments:
#   1) cohort:   bisc / sepages / meta / bisc_mix / sepages_mix / meta_mix
#   2) sample:   2tp / tp (run folder name: <cohort>_<sample>_<tp>_<date>)
#   3) run date: YYYYMMDD of the run folders
#   Rscript script/4_QC/2_Join_QQplots.R sepages 2tp 20261021
#   Rscript script/4_QC/2_Join_QQplots.R bisc tp 20261022
#   Rscript script/4_QC/2_Join_QQplots.R meta tp 20261106
#
# [INPUT]:  results/<cohort>/4_QC/qq_index.csv
#           results/<cohort>/4_QC/<run>/<analysis>/<exposure>_<tp>_qqplot.rds
# [OUTPUT]: results/<cohort>/4_QC/joined/<figure_name>.png
#           (meta: results/7_Meta/4_QC/...)
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_paths.R")
source("script/00_template_plot.R")
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
cohort     <- if (length(args) >= 1) args[1] else "sepages"
sample_set <- if (length(args) >= 2) args[2] else "2tp"
run_date   <- if (length(args) >= 3) args[3] else "20261021"
check_cohort(cohort)
cohort_label <- cohort_labels[[cohort]]

qq_index <- read.csv(results_dir(cohort, "4_QC", "qq_index.csv"), stringsAsFactors = FALSE)
out_dir  <- results_dir(cohort, "4_QC", "joined", create = TRUE)

exp_order <- c("24d", "3pba", "tcpy", "dmp", "dep", "dap", "mix")

run_t2   <- paste(cohort, sample_set, "t2",  run_date, sep = "_")
run_t3   <- paste(cohort, sample_set, "t3",  run_date, sep = "_")
run_avg  <- paste(cohort, sample_set, "avg", run_date, sep = "_")

# Time points: run folder and label for titles
tp_runs   <- c(t2 = run_t2, t3 = run_t3, avg = run_avg)
tp_labels <- c(t2 = "second trimester", t3 = "third trimester",
               avg = "average of second and third trimester")

# Check: every run folder present in the QQ index (should be character(0))
setdiff(tp_runs, unique(qq_index$run))

# Available plots
table(qq_index$run, qq_index$analysis)

#==============================================================================#----
# FIGURES                                                                      #----
#==============================================================================#----
# One list element per figure. Fields:
#   runs, analyses, exposures : which plots to include (all combinations found)
#   ncol                      : panels per row
#   title                     : figure title
#   width, height             : size in cm

figures <- list(
  
  # 1. Main model, six exposures, t2
  main_t2 = list(runs = run_t2, analyses = "main",
                 exposures = exp_order, ncol = 3,
                 title = paste(cohort_label, "- main models, second trimester"),
                 width = 24, height = 16),
  
  # 2. Main model, six exposures, t3
  main_t3 = list(runs = run_t3, analyses = "main",
                 exposures = exp_order, ncol = 3,
                 title = paste(cohort_label, "- main models, third trimester"),
                 width = 24, height = 16),
  
  # 3. Main model, six exposures, average t2-t3
  main_avg = list(runs = run_avg, analyses = "main",
                  exposures = exp_order, ncol = 3,
                  title = paste(cohort_label, "- main models, average of second and third trimester"),
                  width = 24, height = 16)
)
names(figures) <- paste(cohort, sample_set, names(figures), sep = "_")

#==============================================================================#----
# FIGURES - main and sensitivity analyses, one per exposure and time point      #----
#==============================================================================#----

for (e in exp_order) {
  for (tp in names(tp_runs)) {
    figures[[paste0(cohort, "_", sample_set, "_", e, "_sens_", tp)]] <-
      list(runs = tp_runs[[tp]],
           analyses = c("main", paste0("s", 1:6)),
           exposures = e, ncol = 3,
           title = paste0(cohort_label, " - ", e, ", main and sensitivity analyses, ",
                          tp_labels[[tp]]),
           width = 24, height = 16)
  }
}

# Check: 3 main figures + 6 exposures x 3 time points = 21
cat("Figures defined:", length(figures), "\n")
stopifnot(!any(duplicated(names(figures))))
#==============================================================================#----
# BUILD & SAVE                                                                 #----
#==============================================================================#----

for (fig_name in names(figures)) {
  
  fig <- figures[[fig_name]]
  cat("\nFigure:", fig_name, "\n")
  
  # Select plots, ordered by run, analysis, exposure
  sel <- qq_index[qq_index$run %in% fig$runs &
                    qq_index$analysis %in% fig$analyses &
                    qq_index$exposure %in% fig$exposures, ]
  sel <- sel[order(match(sel$run, fig$runs),
                   match(sel$analysis, fig$analyses),
                   match(sel$exposure, fig$exposures)), ]
  
  # Check: all requested plots found
  n_expected <- length(fig$runs) * length(fig$analyses) * length(fig$exposures)
  cat("Plots found:", nrow(sel), "of", n_expected, "\n")
  if (nrow(sel) == 0) { cat("Nothing to plot, skipped\n"); next }
  if (nrow(sel) < n_expected) warning(fig_name, ": ", n_expected - nrow(sel), " plot(s) missing")
  stopifnot(all(file.exists(sel$rds)))
  
  # Load plot objects and join
  plot_list <- lapply(sel$rds, readRDS)
  p_joined  <- wrap_plots(plot_list, ncol = fig$ncol) +
    plot_annotation(title = fig$title)
  
  save_plot(p_joined, file.path(out_dir, paste0(fig_name, ".png")),
            width = fig$width, height = fig$height)
  cat("Saved:", file.path(out_dir, paste0(fig_name, ".png")), "\n")
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----