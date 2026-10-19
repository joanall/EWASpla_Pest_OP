#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: Manhattan plot for each sensitivity EWAS result. Bonferroni-
# significant CpGs coloured by direction: hypermethylated vs hypomethylated.
#  s1 = Excluding pregnancy complications
#  s2 = Excluding non-Caucasian
#  s3 = Adjusting for specific gravity
#  s4 = Adjusting for fruit and vegetable intake
#  s5 = Without gestational age
#
# [INPUT]:  results/3_Clean_Results/sens_t{2,3}/s{1..5}/{24d,3pba,tcpy,dmp,dep,dap}_t{2,3}.rds
#           data.frame, rows = CpGs: Effect, SE, tvalue, pvalue, CpG, pval_bn
# [OUTPUT]: results/4_Analyze_Results/sens_t{2,3}/s{1..5}/<exposure>_t{2,3}_manhattanplot.png

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
# Annotation: only chromosome and position needed
annotation <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19))
annotation <- annotation[, c("chr", "pos")]

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
# MANHATTAN PLOTS - SINGLE + COMBINED T2 - T3                                  #----
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
      
      # Check: CpGs found in the annotation
      cat("  CpGs in annotation:", sum(ewas_results$CpG %in% rownames(annotation)),
          "of", nrow(ewas_results), "\n")
      
      # Plot 1. Manhattan plot for one exposure, time point and sensitivity
      # analysis: genomic position (x) vs -log10(p) (y); Bonferroni CpGs by direction.
      p_manh <- manhattan_plot(ewas_results, annotation,
                               title = paste0(exp_labels[i], " - ", tp, " (", sens_labels[j], ")"))
      save_plot(p_manh,
                paste0(output_folder, "sens_", tp, "/", s, "/", exposures[i], "_", tp, "_manhattanplot.png"),
                width = 24, height = 10)
      
      ewas_results$Time <- tp
      ewas_both <- rbind(ewas_both, ewas_results[, c("CpG", "Effect", "pvalue", "pval_bn", "Time")])
    }
    
    # Plot 2. Manhattan plots of one exposure at t2 (top) and t3 (bottom) for
    # one sensitivity analysis, same y-axis to compare time points.
    p_both <- manhattan_plot(ewas_both, annotation,
                             title = paste0(exp_labels[i], " (", sens_labels[j], ")")) +
      facet_wrap(~ Time, ncol = 1)
    save_plot(p_both,
              paste0(output_folder, "sens_t2_t3/", s, "/manhattan_", exposures[i], "_t2_t3.png"),
              width = 24, height = 16)
  }
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----