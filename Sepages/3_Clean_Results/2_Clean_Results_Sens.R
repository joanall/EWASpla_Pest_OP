#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: Remove unreliable probes from the sensitivity EWAS results
# (filter_EPIC list) and compute FDR and Bonferroni adjusted p-values.
#
# [INPUT]:  results/2_EWAS_Results/sens_t{2,3}/s{1..5}/{24d,3pba,tcpy,dmp,dep,dap}_t{2,3}.rds
#           data.frame, rows = CpGs: Effect, SE, tvalue, pvalue, CpG
# [OUTPUT]: results/3_Clean_Results/sens_t{2,3}/s{1..5}/{24d,...,dap}_t{2,3}.rds
#           data.frame, rows = CpGs: Effect, SE, tvalue, pvalue, CpG, pval_fdr, pval_bn
#           + excluded probes and QC descriptives per model (.txt)

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_functions.R")

input_folder  <- "results/2_EWAS_Results/"
output_folder <- "results/3_Clean_Results/"

time_points <- c("t2", "t3")
sens        <- c("s1", "s2", "s3", "s4", "s5")
exposures   <- c("24d", "3pba", "tcpy", "dmp", "dep", "dap")

# Parameters to exclude CpG probes
ethnic  <- "GMAF1p"
array   <- "EPIC"
exclude <- c("control_probes", "noncpg_probes", "Sex", "MASK_mapping",
             "MASK_sub30_copy", "MASK_extBase", "MASK_typeINextBaseSwitch",
             "MASK_snp5_common", "MASK_snp5_ethnic",
             "Unrel_450_EPIC_pla_restrict")
# 'MASK_snp5_ethnic' -> remove SNPs common in European population

# Check: all input files exist
files_in <- NULL
for (tp in time_points) for (s in sens) {
  files_in <- c(files_in, paste0(input_folder, "sens_", tp, "/", s, "/", exposures, "_", tp, ".rds"))
}
length(files_in)                   # 2 x 5 x 6 = 60
files_in[!file.exists(files_in)]   # should be character(0)

#==============================================================================#----
# LOAD, EXCLUDE UNRELIABLE PROBES & ADJUST P-VALUES                            #----
#==============================================================================#----

for (tp in time_points) {
  for (s in sens) {
    
    out_dir <- paste0(output_folder, "sens_", tp, "/", s, "/")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    for (e in exposures) {
      
      file_in  <- paste0(input_folder, "sens_", tp, "/", s, "/", e, "_", tp, ".rds")
      exp_name <- paste0(e, "_", tp)
      cat("\n", tp, "|", s, "|", e, "\n")
      
      ewas_results <- readRDS(file_in)
      cat("Number CpGs before exclusion", nrow(ewas_results), "\n")
      
      ewas_results$probeID <- rownames(ewas_results)
      ewas_results <- exclude_CpGs(ewas_results, "probeID", exclude,
                                   ethnic     = ethnic,
                                   filename   = paste0(out_dir, exp_name, "_excluded_probes.txt"),
                                   fileresume = paste0(out_dir, exp_name, "_cpg_qc_descriptives.txt"),
                                   artype     = array)
      
      # Keep model statistics only, as numeric; CpG as rownames and column
      cpg_ids <- ewas_results$probeID
      ewas_results <- ewas_results[, c("Effect", "SE", "tvalue", "pvalue")]
      ewas_results <- as.data.frame(sapply(ewas_results, as.numeric))
      rownames(ewas_results) <- cpg_ids
      ewas_results$CpG <- cpg_ids
      cat("Number CpGs after exclusion", nrow(ewas_results), "\n")
      
      # Adjusted p-values on the filtered CpGs
      ewas_results$pval_fdr <- p.adjust(ewas_results$pvalue, method = "fdr")
      ewas_results$pval_bn  <- p.adjust(ewas_results$pvalue, method = "bonferroni")
      cat("Bonferroni < 0.05:", sum(ewas_results$pval_bn < 0.05, na.rm = TRUE), "\n")
      
      saveRDS(ewas_results, paste0(out_dir, exp_name, ".rds"))
    }
  }
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----