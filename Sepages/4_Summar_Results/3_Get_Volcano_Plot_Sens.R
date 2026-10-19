#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#

# [DESCRIPTION]: Volcano plot for each sensitivity EWAS result. Bonferroni-
# significant CpGs coloured by direction: hypermethylated vs hypomethylated.
#  s1 = Excluding pregnancy complications
#  s2 = Excluding non-Caucasian
#  s3 = Adjusting for specific gravity
#  s4 = Adjusting for fruit and vegetable intake
#  s5 = Without gestational age
#
# [INPUT]:  results/3_Clean_Results/sens_t{2,3}/s{1..5}/{24d,3pba,tcpy,dmp,dep,dap}_t{2,3}.rds
#           data.frame, rows = CpGs: Effect, SE, tvalue, pvalue, CpG, pval_bn
# [OUTPUT]: results/4_Analyze_Results/sens_t{2,3}/s{1..5}/<exposure>_t{2,3}_volcanoplot.png
#           results/4_Analyze_Results/sens_t2_t3/s{1..5}/volcano_<exposure>_t2_t3.png

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_template_plot.R")
source("script/00_functions.R")

input_folder  <- "results/3_Clean_Results/"
output_folder <- "results/4_Analyze_Results/"

time_points <- c("t2", "t3")
sens        <- c("s1", "s2", "s3", "s4", "s5")
sens_labels <- c("Excl. pregnancy complications", "Excl. non-Caucasian",
                 "Adj. specific gravity", "Adj. fruit and vegetables", "Without GA")

exposures  <- c("24d",   "3pba",  "tcpy", "dmp",       "dep",       "dap")
exp_labels <- c("2,4-D", "3-PBA", "TCPy", "DMP (sum)", "DEP (sum)", "DAP (sum)")

# Check: all input files exist
files_in <- NULL
for (tp in time_points) for (s in sens) {
  files_in <- c(files_in, paste0(input_folder, "sens_", tp, "/", s, "/", exposures, "_", tp, ".rds"))
}
length(files_in)                   # 60
files_in[!file.exists(files_in)]   # should be character(0)

#==============================================================================#----
# VOLCANO PLOTS - SINGLE + COMBINED T2 - T3, PER SENSITIVITY                   #----
#==============================================================================#----

for (j in seq_along(sens)) {
  s <- sens[j]
  cat("\n==================== Sensitivity", s, "-", sens_labels[j], "====================\n")
  
  for (i in seq_along(exposures)) {
    
    ewas_both <- NULL
    for (tp in time_points) {
      file_in <- paste0(input_folder, "sens_", tp, "/", s, "/", exposures[i], "_", tp, ".rds")
      cat("Exposure:", exposures[i], "| Time:", tp, "\n")
      ewas_results <- readRDS(file_in)
      
      # Check: Bonferroni-significant CpGs
      cat("  Bonferroni < 0.05:", sum(ewas_results$pval_bn < 0.05, na.rm = TRUE), "\n")
      
      # Plot 1. Volcano plot for one exposure, time point and sensitivity
      # analysis: effect (x) vs -log10(p) (y); Bonferroni CpGs by direction.
      p_volcano <- volcano_plot(ewas_results,
                                title = paste0(exp_labels[i], " - ", tp, " (", sens_labels[j], ")"))
      save_plot(p_volcano,
                paste0(output_folder, "sens_", tp, "/", s, "/", exposures[i], "_", tp, "_volcanoplot.png"),
                width = 14, height = 12)
      
      ewas_results$Time <- tp
      ewas_both <- rbind(ewas_both, ewas_results[, c("CpG", "Effect", "pvalue", "pval_bn", "Time")])
    }
    
    # Plot 2. Volcano plots of one exposure at t2 and t3 side by side for one
    # sensitivity analysis (same axes, to compare time points).
    p_both <- volcano_plot(ewas_both, title = paste0(exp_labels[i], " (", sens_labels[j], ")")) +
      facet_wrap(~ Time, ncol = 2)
    save_plot(p_both,
              paste0(output_folder, "sens_t2_t3/", s, "/volcano_", exposures[i], "_t2_t3.png"),
              width = 24, height = 12)
  }
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----