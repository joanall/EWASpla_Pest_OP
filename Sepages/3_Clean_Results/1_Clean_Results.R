#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: Remove unreliable probes of the EWAS based on the filter_EPIC list
# list based on xxxx papers. and compute pvalue adjsuted for bonferroni 

# [INPUT]:  results/2_EWAS_Results/main_t{2,3}/{24d,3pba, tcpy,dmp,dep,dap}_t{2,3}.rds
#           data.frame, rows = CpGs: Effect, SE, tvalue, pvalue, CpG

# [OUTPUT]: results/3_Clean_Results/main_t{2,3}/{24d,3pba, tcpy,dmp,dep,dap}_t{2,3}.rds
#           data.frame, rows = CpGs: Effect, SE, tvalue, pvalue, CpG, pval_bn


#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_functions.R")
files_t2 <- c("results/2_EWAS_Results/main_t2/24d_t2.rds",
              "results/2_EWAS_Results/main_t2/3pba_t2.rds", 
              "results/2_EWAS_Results/main_t2/tcpy_t2.rds",
              "results/2_EWAS_Results/main_t2/dmp_t2.rds",
              "results/2_EWAS_Results/main_t2/dep_t2.rds",
              "results/2_EWAS_Results/main_t2/dap_t2.rds")

files_t3 <- c("results/2_EWAS_Results/main_t3/24d_t3.rds",
              "results/2_EWAS_Results/main_t3/3pba_t3.rds", 
              "results/2_EWAS_Results/main_t3/tcpy_t3.rds",
              "results/2_EWAS_Results/main_t3/dmp_t3.rds",
              "results/2_EWAS_Results/main_t3/dep_t3.rds",
              "results/2_EWAS_Results/main_t3/dap_t3.rds")


results_folder_t2 <- "results/3_Clean_Results/main_t2/"
results_folder_t3 <- "results/3_Clean_Results/main_t3/"

dir.create(results_folder_t2, recursive = TRUE, showWarnings = FALSE)
dir.create(results_folder_t3, recursive = TRUE, showWarnings = FALSE)

# Define input parameters to exclude CpGs probes
ethnic <- 'GMAF1p'
array <- "EPIC"
exclude <- c('control_probes', 'noncpg_probes', 'Sex','MASK_mapping', 
             'MASK_sub30_copy', 'MASK_extBase', 'MASK_typeINextBaseSwitch', 
             'MASK_snp5_common', 'MASK_snp5_ethnic',
             'Unrel_450_EPIC_pla_restrict')
# 'MASK_snp5_ethnic' -> remove SNPs common in European population 

#==============================================================================#----
# LOAD DATA & EXCLUDE UNRELIABLE PROBES  T2                                    #----
#==============================================================================#----
for (i in seq_along(files_t2)) {
  
  exp_filename  <- tools::file_path_sans_ext(basename(files_t2[i]))
  cat('Exposure', exp_filename, "\n")
  
  # Load EWAS results 
  ewas_results <- readRDS(files_t2[i]) 
  cat('Number CpGs before exclusion', nrow(ewas_results),"\n")
  
  n <- ncol(ewas_results) 
  ewas_results <- cbind(ewas_results, probeID = rownames(ewas_results))
  # alternative: ewas_results$probeID <- ewas_results$CpG
  

  ewas_results <- exclude_CpGs(ewas_results, "probeID", 
                               exclude, 
                               ethnic = ethnic, 
                               filename = paste0(results_folder_t2, exp_filename,'_excluded_probes.txt'), 
                               fileresume = paste0(results_folder_t2, exp_filename,'_cpg_qc_descriptives.txt'), 
                               artype = array )
  
  # Arrange output from exclude_CpGs 
  # - CpGs as rownames
  # - Remove 1st column (probeID) 
  # - Remove last columns with exclusion criteria added by function
  # - methyaltion values as numeric 
  
  cpg_ids <- ewas_results$probeID
  ewas_results <- ewas_results[, c("Effect", "SE", "tvalue", "pvalue")]   # model statistics only
  ewas_results <- as.data.frame(sapply(ewas_results, as.numeric))
  rownames(ewas_results) <- cpg_ids
  ewas_results$CpG <- cpg_ids
  cat("Number CpGs after exclusion", nrow(ewas_results), "\n")
  
  
  # Adjusted p-values on the filtered CpGs
  ewas_results$pval_bn  <- p.adjust(ewas_results$pvalue, method = "bonferroni")
  cat("Bonferroni < 0.05:", sum(ewas_results$pval_bn < 0.05, na.rm = TRUE), "\n")

  saveRDS(ewas_results, paste0(results_folder_t2, basename(files_t2[i])))
  }

#==============================================================================#----
# LOAD DATA & EXCLUDE UNRELIABLE PROBES  T3                                    #----
#==============================================================================#----
for (i in seq_along(files_t3)) {
  
  exp_filename  <- tools::file_path_sans_ext(basename(files_t3[i]))
  cat('Exposure', exp_filename, "\n")
  
  # Load EWAS results 
  ewas_results <- readRDS(files_t3[i]) 
  cat('Number CpGs before exclusion', nrow(ewas_results),"\n")
  
  n <- ncol(ewas_results) 
  ewas_results <- cbind(ewas_results, probeID = rownames(ewas_results))
  # alternative: ewas_results$probeID <- ewas_results$CpG
  
  
  ewas_results <- exclude_CpGs(ewas_results, "probeID", 
                               exclude, 
                               ethnic = ethnic, 
                               filename = paste0(results_folder_t3, exp_filename,'_excluded_probes.txt'), 
                               fileresume = paste0(results_folder_t3, exp_filename,'_cpg_qc_descriptives.txt'), 
                               artype = array )
  
  # Arrange output from exclude_CpGs 
  # - CpGs as rownames
  # - Remove 1st column (probeID) 
  # - Remove last columns with exclusion criteria added by function
  # - methyaltion values as numeric 
  
  cpg_ids <- ewas_results$probeID
  ewas_results <- ewas_results[, c("Effect", "SE", "tvalue", "pvalue")]   # model statistics only
  ewas_results <- as.data.frame(sapply(ewas_results, as.numeric))
  rownames(ewas_results) <- cpg_ids
  ewas_results$CpG <- cpg_ids
  cat("Number CpGs after exclusion", nrow(ewas_results), "\n")
  
  
  # Adjusted p-values on the filtered CpGs
  ewas_results$pval_bn  <- p.adjust(ewas_results$pvalue, method = "bonferroni")
  cat("Bonferroni < 0.05:", sum(ewas_results$pval_bn < 0.05, na.rm = TRUE), "\n")
  
  saveRDS(ewas_results, paste0(results_folder_t3, basename(files_t3[i])))
}

#==============================================================================#----
# END OF SCRIPT.                                                               #----
#==============================================================================#----
