#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: Sensitivity EWAS of pesticide exposure at t3 on placental DNAm
# (SEPAGES). Robust linear regression per CpG, run in parallel by chunks.
# Outer loop: sensitivity analyses. Inner loop: exposures.
#
# Sensitivity analyses (all with main covariates + cell types, except as noted):
#  s1 = Excluding pregnancy complications
#  s2 = Excluding non-Caucasian (ethnicity removed from covariates: constant)
#  s3 = Adjusting for specific gravity
#  s4 = Adjusting for fruit and vegetable intake
#  s5 = Without gestational age
#
# Run all sensitivities, or one per job (e.g. SLURM array 1-5):
#   Rscript script/2_Run_Models/Run_Sens_Models_t3.R        # all
#   Rscript script/2_Run_Models/Run_Sens_Models_t3.R 3      # only s3
#
# Continuous exposures: log2 imputed. 2,4-D: detected (1) vs not detected (0).
#
# [INPUT]:  db/subsets/df_s1_t3.rds ... df_s5_t3.rds
#           SEPAGES_20230615_PreprocessedBetas_nooutliers.RData (beta values,
#           rows = CpGs, columns = Basename)
# [OUTPUT]: results/2_EWAS_Results/sens_t3/s{1..5}/<exposure>_t3.rds
#           data.frame, rows = CpGs: Effect, SE, tvalue, pvalue, CpG
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
library(MASS)
library(sfsmisc)
library(foreach)
library(doParallel)
source("script/00_functions.R")

betas_path     <- "db/plaDNAm_data/SEPAGES_20230615_PreprocessedBetas_nooutliers.RData"
results_folder <- "results/2_EWAS_Results/sens_t3"
dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)

n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 8))

# Sensitivity analyses to run: argument from the terminal, or all
args <- commandArgs(trailingOnly = TRUE)
which_sens <- if (length(args) == 0) 1:5 else as.numeric(args[1])
# which_sens <- 2   # TEST - remove for the real run

# Log file: everything printed goes to the console and to the log
log_folder <- file.path(results_folder, "logs")
dir.create(log_folder, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_folder, paste0("sens_", paste(which_sens, collapse = "-"), "_",
                                         format(Sys.time(), "%Y%m%d_%H%M"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
cat("Start:", format(Sys.time()), "\n")
cat("Sensitivity analyses:", paste0("s", which_sens, collapse = ", "), "\n")

#==============================================================================#----
# PARAMETERS                                                                   #----
#==============================================================================#----

# Exposures (same as main models)
exposures     <- c("mo_pest_24d_bi_t3", "mo_pest_3pba_i_log2_t3", "mo_pest_tcpy_i_log2_t3",
                   "mo_pest_dmp_ms_i_log2_t3", "mo_pest_dep_ms_i_log2_t3",
                   "mo_pest_dap_ms_i_log2_t3")
analysis_name <- c("24d_t3", "3pba_t3", "tcpy_t3", "dmp_t3", "dep_t3", "dap_t3")

# Covariates (same names as main models)
main_cov   <- c("mo_age", "mo_edu", "mo_smoke", "po_gd_lmp_days",
                "ch_ethnicity", "ch_sex", "Meanlog2oddsContamination")
cell_types <- c("Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC")
# Syncytiotrophoblast left out as reference: the 6 proportions sum to 1

# Subset and covariates for each sensitivity analysis
sens_subset <- c(s1 = "db/subsets/df_s1_t3.rds",
                 s2 = "db/subsets/df_s2_t3.rds",
                 s3 = "db/subsets/df_s3_t3.rds",
                 s4 = "db/subsets/df_s4_t3.rds",
                 s5 = "db/subsets/df_s5_t3.rds")

sens_covars <- list(
  s1 = c(main_cov, cell_types),                                  # excl. complications
  s2 = c(setdiff(main_cov, "ch_ethnicity"), cell_types),         # excl. non-Caucasian
  s3 = c(main_cov, cell_types, "mo_pool_sg_T1"),                 # + specific gravity
  s4 = c(main_cov, cell_types, "fv_g1000kcal"),                  # + fruit & vegetables
  s5 = c(setdiff(main_cov, "po_gd_lmp_days"), cell_types)        # without GA
)

#==============================================================================#----
# LOAD BETAS (ONCE)                                                            #----
#==============================================================================#----

betas_obj <- load(betas_path)
meth_all  <- as.matrix(get(betas_obj[1]))
rm(list = betas_obj)
cat("Number of CpGs:", nrow(meth_all), "| Samples:", ncol(meth_all), "\n")
range(meth_all, na.rm = TRUE)

# meth_all <- meth_all[1:50, , drop = FALSE]   # TEST - remove for the real run

cl <- makeCluster(n_cores)
registerDoParallel(cl)

#==============================================================================#----
# LOOP: SENSITIVITY ANALYSES x EXPOSURES                                       #----
#==============================================================================#----

for (s in paste0("s", which_sens)) {
  
  cat("\n==================== Sensitivity", s, "====================\n")
  
  # 1. Subset and covariates
  metadata <- readRDS(sens_subset[[s]])
  covars   <- sens_covars[[s]]
  cat("Subset:", sens_subset[[s]], "| N:", nrow(metadata), "\n")
  cat("Covariates:", paste(covars, collapse = ", "), "\n")
  
  # Check: all variables exist, missingness, no constant covariate
  print(setdiff(c(exposures, covars), colnames(metadata)))
  print(colSums(is.na(metadata[, covars])))
  n_levels <- sapply(metadata[, covars], function(x) length(unique(na.omit(x))))
  if (any(n_levels < 2)) stop("Constant covariate(s) in ", s, ": ",
                              paste(names(n_levels)[n_levels < 2], collapse = ", "))
  
  # 2. Match and order betas with this subset
  metadata <- metadata[metadata$Basename %in% colnames(meth_all), ]
  meth_data <- meth_all[, metadata$Basename, drop = FALSE]
  meth_data <- meth_data[1:50,]
  stopifnot(identical(colnames(meth_data), metadata$Basename))
  cat("Samples in the analysis:", ncol(meth_data), "\n")
  
  # 3. Split CpGs into chunks
  chunk_id   <- cut(seq_len(nrow(meth_data)), breaks = n_cores, labels = FALSE)
  cpg_chunks <- lapply(split(seq_len(nrow(meth_data)), chunk_id),
                       function(idx) meth_data[idx, , drop = FALSE])
  rm(meth_data)
  
  out_folder <- file.path(results_folder, s)
  dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)
  
  # 4. One EWAS per exposure
  for (e in seq_along(exposures)) {
    
    exposure <- exposures[e]
    n_model  <- sum(complete.cases(metadata[, c(exposure, covars)]))
    cat("\n--", s, "|", analysis_name[e], "| N in the model:", n_model, "\n")
    
    start_time <- proc.time()
    chunk_results_list <- foreach(chunk = cpg_chunks,
                                  .packages = c("MASS", "sfsmisc"),
                                  .export   = c("rlm_single_cpg", "exposure", "covars", "metadata")) %dopar% {
                                    rlm_single_cpg(cpg_chunk  = chunk,
                                                   exposures  = exposure,
                                                   covariates = covars,
                                                   metadata   = metadata)
                                  }
    run_time <- proc.time() - start_time
    cat("Run time (min):", round(run_time["elapsed"] / 60, 1), "\n")
    
    results <- do.call(rbind, lapply(chunk_results_list, function(x) x[[exposure]]))
    rownames(results) <- results$CpG
    cat("Failed CpGs:", sum(is.na(results$pvalue)), "\n")
    
    saveRDS(results, file.path(out_folder, paste0(analysis_name[e], ".rds")))
  }
  rm(cpg_chunks)
}

stopCluster(cl)

cat("\nEnd:", format(Sys.time()), "\n")
sessionInfo()
sink(type = "message")
sink()
close(log_con)

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----