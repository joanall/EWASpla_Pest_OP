#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#

# [DESCRIPTION]: Volcano plot for each EWAS result. Bonferroni-significant
# CpGs coloured by direction: hypermethylated vs hypomethylated.
#
# [INPUT]:  results/3_Clean_Results/main_t{2,3}/{24d,3pba,tcpy,dmp,dep,dap}_t{2,3}.rds
#           data.frame, rows = CpGs: Effect, SE, tvalue, pvalue, CpG, pval_bn
# [OUTPUT]: results/4_Analyze_Results/main_t{2,3}/<exposure>_t{2,3}_volcanoplot.png
#           results/4_Analyze_Results/volcano_<exposure>_t2_t3.png (t2 next to t3)

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_template_plot.R")
source("script/00_functions.R")

input_folder <- "results/3_Clean_Results/"
output_folder <- "results/4_Analyze_Results/"

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)


exposures  <- c("24d",   "3pba",  "tcpy", "dmp",       "dep",       "dap")
exp_labels <- c("2,4-D", "3-PBA", "TCPy", "DMP (sum)", "DEP (sum)", "DAP (sum)")

# Check: all input files exist
files_in <- c(paste0(input_folder, "main_t2/", exposures, "_t2.rds"),
              paste0(input_folder, "main_t3/", exposures, "_t3.rds"))
files_in[!file.exists(files_in)]   # should be character(0)

#==============================================================================#----
# VOLCANO PLOTS - SINGLE + COMBINED T2 - T3                                    #----
#==============================================================================#----

for (i in seq_along(exposures)) {
  
  ewas_both <- NULL
  for (tp in c("t2", "t3")) {
    file_in <- paste0(input_folder, "main_", tp, "/", exposures[i], "_", tp, ".rds")
    cat("Exposure:", exposures[i], "| Time:", tp, "\n")
    ewas_results <- readRDS(file_in)
    
    # Check: Bonferroni-significant CpGs
    cat("  Bonferroni < 0.05:", sum(ewas_results$pval_bn < 0.05, na.rm = TRUE), "\n")
    
    # Plot 1. Volcano plot for one exposure and time point: effect (x) vs
    # -log10(p) (y); Bonferroni CpGs coloured by direction.
    p_volcano <- volcano_plot(ewas_results, title = paste(exp_labels[i], "-", tp))
    save_plot(p_volcano,
              paste0(output_folder, "main_", tp, "/", exposures[i], "_", tp, "_volcanoplot.png"),
              width = 14, height = 12)
    
    ewas_results$Time <- tp
    ewas_both <- rbind(ewas_both, ewas_results[, c("CpG", "Effect", "pvalue", "pval_bn", "Time")])
  }
  
  # Plot 2. Volcano plots of one exposure at t2 and t3 side by side
  # (same axes, to compare effect sizes between time points).
  p_both <- volcano_plot(ewas_both, title = exp_labels[i]) +
    facet_wrap(~ Time, ncol = 2)
  print(p_both)
  save_plot(p_both, paste0(output_folder, "volcano_", exposures[i], "_t2_t3.png"),
            width = 24, height = 12)
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----