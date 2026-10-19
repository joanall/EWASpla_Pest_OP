#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#

# [DESCRIPTION]: Manhattan plot for each EWAS result. Bonferroni-significant
# CpGs coloured by direction; dashed line = Bonferroni, dotted = suggestive
# (p < 1e-5).
#
# [INPUT]:  results/3_Clean_Results/main_t{2,3}/{24d,3pba,tcpy,dmp,dep,dap}_t{2,3}.rds
#           data.frame, rows = CpGs: Effect, SE, tvalue, pvalue, CpG, pval_bn
# [OUTPUT]: results/4_Analyze_Results/main_t{2,3}/<exposure>_t{2,3}_manhattan.png
#           results/4_Analyze_Results/manhattan_<exposure>_t2_t3.png (t2 above t3)

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
source("script/00_template_plot.R")
source("script/00_functions.R")

input_folder <- "results/3_Clean_Results/"
output_folder <- "results/4_Analyze_Results/"

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)


exposures  <- c("24d",   "3pba",  "tcpy", "dmp",       "dep",       "dap")
exp_labels <- c("2,4-D", "3-PBA", "TCPy", "DMP (sum)", "DEP (sum)", "DAP (sum)")

# Annotation: only chromosome and position needed
annotation <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19))
annotation <- annotation[, c("chr", "pos")]
dim(annotation)

# Check: all input files exist
files_in <- c(paste0(input_folder, "main_t2/", exposures, "_t2.rds"),
              paste0(input_folder, "main_t3/", exposures, "_t3.rds"))
files_in[!file.exists(files_in)]   # should be character(0)

#==============================================================================#----
# MANHATTAN PLOTS - SINGLE + COMBINED T2 - T3                                  #----
#==============================================================================#----

for (i in seq_along(exposures)) {
  
  ewas_both <- NULL
  for (tp in c("t2", "t3")) {
    file_in <- paste0(input_folder, "main_", tp, "/", exposures[i], "_", tp, ".rds")
    cat("Exposure:", exposures[i], "| Time:", tp, "\n")
    ewas_results <- readRDS(file_in)
    
    # Check: CpGs found in the annotation
    cat("  CpGs in annotation:", sum(ewas_results$CpG %in% rownames(annotation)),
        "of", nrow(ewas_results), "\n")
    
    # Plot 1. Manhattan plot for one exposure and time point: genomic position
    # (x) vs -log10(p) (y); Bonferroni CpGs coloured by direction.
    p_manh <- manhattan_plot(ewas_results, annotation, title = paste(exp_labels[i], "-", tp))
    save_plot(p_manh,
              paste0(output_folder, "main_", tp, "/", exposures[i], "_", tp, "_manhattan.png"),
              width = 24, height = 10)
    
    ewas_results$Time <- tp
    ewas_both <- rbind(ewas_both, ewas_results[, c("CpG", "Effect", "pvalue", "pval_bn", "Time")])
  }
  
  # Plot 2. Manhattan plots of one exposure at t2 (top) and t3 (bottom),
  # same y-axis to compare signals between time points.
  p_both <- manhattan_plot(ewas_both, annotation, title = exp_labels[i]) +
    facet_wrap(~ Time, ncol = 1)
  print(p_both)
  save_plot(p_both, paste0(output_folder, "manhattan_", exposures[i], "_t2_t3.png"),
            width = 24, height = 16)
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----