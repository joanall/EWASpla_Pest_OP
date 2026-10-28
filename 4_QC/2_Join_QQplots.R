#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: Combine saved QQ plots (from 4_QQ_Plots.R) into multi-panel
# figures. Each figure is defined in FIGURES below: which plots go in it, the
# layout and the output name. Add or edit figures there.
#
# [INPUT]:  results/4_QQ_Plots/qq_index.csv
#           results/4_QQ_Plots/<run>/<analysis>/<exposure>_<tp>_qqplot.rds
# [OUTPUT]: results/4_QQ_Plots/joined/<figure_name>.png
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_template_plot.R")
library(patchwork)

qq_index <- read.csv("results/4_QQ_Plots/qq_index.csv", stringsAsFactors = FALSE)
out_dir  <- "results/4_QQ_Plots/joined"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

exp_order <- c("24d", "3pba", "tcpy", "dmp", "dep", "dap")

run_date <- "20261021"
run_t2   <- paste0("sepages_2tp_t2_", run_date)
run_t3   <- paste0("sepages_2tp_t3_", run_date)

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
  
  # 1. SEPAGES main model, six exposures, t2
  sepages_main_t2 = list(runs = run_t2, analyses = "main",
                         exposures = exp_order, ncol = 3,
                         title = "SEPAGES - main models, second trimester",
                         width = 24, height = 16),
  
  # 2. SEPAGES main model, six exposures, t3
  sepages_main_t3 = list(runs = run_t3, analyses = "main",
                         exposures = exp_order, ncol = 3,
                         title = "SEPAGES - main models, third trimester",
                         width = 24, height = 16)
)

#==============================================================================#----
# FIGURES - main and sensitivity analyses, one per exposure and time point      #----
#==============================================================================#----

for (e in exp_order) {
  for (tp in c("t2", "t3")) {
    run_tp <- if (tp == "t2") run_t2 else run_t3
    figures[[paste0("sepages_", e, "_sens_", tp)]] <-
      list(runs = run_tp,
           analyses = c("main", paste0("s", 1:5)),
           exposures = e, ncol = 3,
           title = paste0("SEPAGES - ", e, ", main and sensitivity analyses, ",
                          toupper(tp)),
           width = 24, height = 16)
  }
}

# Check: 2 main figures + 6 exposures x 2 time points = 14
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