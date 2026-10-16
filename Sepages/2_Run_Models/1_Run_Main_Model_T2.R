#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: EWAS of pesticide exposure at t2 on placental DNAm (SEPAGES).
# Robust linear regression per CpG, run in parallel by chunks.
# Run from the terminal with the model number as argument:
#   Rscript script/Run_Main_Model_T2.R 1

# Models (which_model):
#  1 = DNAm ~ 2,4-D (detected) + main cov + CT     7 = same, without CT
#  2 = DNAm ~ 3-PBA + main cov + CT                8 = same, without CT
#  3 = DNAm ~ TCPy + main cov + CT                 9 = same, without CT
#  4 = DNAm ~ DMP (sum) + main cov + CT           10 = same, without CT
#  5 = DNAm ~ DEP (sum) + main cov + CT           11 = same, without CT
#  6 = DNAm ~ DAP (sum) + main cov + CT           12 = same, without CT

# Main cov: maternal age, education, smoking, GA, child ethnicity, child sex,
#           DNA contamination. CT: 5 cell types (-)
# Continuous exposures: log2 imputed. 2,4-D: detected (1) vs not detected (0).

# [INPUT]:  db/subsets/df_main_t2.rds (metadata, complete cases)
#           SEPAGES_20230615_PreprocessedBetas_nooutliers.RData (beta values,
#           rows = CpGs, columns = Basename)

# [OUTPUT]: results/2_EWAS_Results/main_t2/<exposure>_t2[_noCT].rds
#           data.frame, rows = CpGs: Effect, SE, tvalue, pvalue, CpG, 


#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
library(MASS)
library(sfsmisc)
library(foreach)
library(doParallel)
source("script/00_functions.R")

metadata_path  <- "db/subsets/df_main_t2.rds"
betas_path     <- "db/plaDNAm_data/SEPAGES_20230615_PreprocessedBetas_nooutliers.RData"
results_folder <- "results/2_EWAS_Results/main_t2"
dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)

n_cores <- 8

# Model number from the terminal (default 1 when run interactively)
args <- commandArgs(trailingOnly = TRUE)
which_model <- if (length(args) == 0) 1 else as.numeric(args[1])
# which_model <- 2   # TEST - remove for the real run
cat("Model:", which_model, "\n")

# Log file: everything printed goes to the console and to the log
log_folder <- file.path(results_folder, "logs")
dir.create(log_folder, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_folder, paste0("model_", which_model, "_",
                                         format(Sys.time(), "%Y%m%d_%H%M"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)              # printed output: console + log
sink(log_con, type = "message")          # warnings and errors: log only
cat("Start:", format(Sys.time()), "\n")
cat("Model:", which_model, "\n")
#==============================================================================#----
# LOAD DATA                                                                    #----
#==============================================================================#----

metadata <- readRDS(metadata_path)
dim(metadata)

# Beta values: check the object name inside the .RData
betas_obj <- load(betas_path)
meth_data <- as.matrix(get(betas_obj[1]))
rm(list = betas_obj)

cat("Number of CpGs:", nrow(meth_data), "\n")
cat("Number of samples:", ncol(meth_data), "\n")

# Check: beta values between 0 and 1
range(meth_data, na.rm = TRUE)

#==============================================================================#----
# MATCH & ORDER METADATA AND BETAS                                             #----
#==============================================================================#----

sum(metadata$Basename %in% colnames(meth_data))
sum(colnames(meth_data) %in% metadata$Basename)

# Keep participants with both, and order betas columns as metadata rows
metadata  <- metadata[metadata$Basename %in% colnames(meth_data), ]
meth_data <- meth_data[, metadata$Basename]
meth_data <- meth_data[1:50,]
# Check: same samples in the same order (stops the script if not)
stopifnot(identical(colnames(meth_data), metadata$Basename))
cat("Samples in the analysis:", ncol(meth_data), "\n")

#==============================================================================#----
# PARAMETERS                                                                   #----
#==============================================================================#----

# Exposure
if (which_model %in% c(1, 7))  {exposure <- "mo_pest_24d_bi_t2";       analysis_name <- "24d_t2"}
if (which_model %in% c(2, 8))  {exposure <- "mo_pest_3pba_i_log2_t2";   analysis_name <- "3pba_t2"}
if (which_model %in% c(3, 9))  {exposure <- "mo_pest_tcpy_i_log2_t2";   analysis_name <- "tcpy_t2"}
if (which_model %in% c(4, 10)) {exposure <- "mo_pest_dmp_ms_i_log2_t2"; analysis_name <- "dmp_t2"}
if (which_model %in% c(5, 11)) {exposure <- "mo_pest_dep_ms_i_log2_t2"; analysis_name <- "dep_t2"}
if (which_model %in% c(6, 12)) {exposure <- "mo_pest_dap_ms_i_log2_t2"; analysis_name <- "dap_t2"}

# Covariates
main_cov   <- c("mo_age", "mo_edu", "mo_smoke", "po_gd_lmp_days",
                "ch_ethnicity", "ch_sex", "Meanlog2oddsContamination")
cell_types <- c("Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC")
# Syncytiotrophoblast left out as reference: the 6 proportions sum to 1

if (which_model <= 6) {
  covars <- c(main_cov, cell_types)
} else {
  covars <- main_cov
  analysis_name <- paste0(analysis_name, "_noCT")
}

cat("Exposure:", exposure, "| Analysis:", analysis_name, "\n")
cat("Covariates:", paste(covars, collapse = ", "), "\n")

# Check: all variables exist and missingness
setdiff(c(exposure, covars), colnames(metadata))
colSums(is.na(metadata[, c(exposure, covars)]))
n_model <- sum(complete.cases(metadata[, c(exposure, covars)]))
cat("N in the model:", n_model, "\n")

#==============================================================================#----
# SPLIT CpGs INTO CHUNKS                                                       #----
#==============================================================================#----

chunk_id   <- cut(seq_len(nrow(meth_data)), breaks = n_cores, labels = FALSE)
cpg_chunks <- lapply(split(seq_len(nrow(meth_data)), chunk_id),
                     function(idx) meth_data[idx, , drop = FALSE])
rm(meth_data)

cat("Number of chunks:", length(cpg_chunks), "\n")
cat("CpGs per chunk:", nrow(cpg_chunks[[1]]), "\n")

#==============================================================================#----
# RUN MODEL                                                                    #----
#==============================================================================#----

cl <- makeCluster(n_cores)
registerDoParallel(cl)

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
stopCluster(cl)
cat("Run time (min):", round(run_time["elapsed"] / 60, 1), "\n")

#==============================================================================#----
# FORMAT RESULTS & SAVE                                                        #----
#==============================================================================#----

results <- do.call(rbind, lapply(chunk_results_list, function(x) x[[exposure]]))
rownames(results) <- results$CpG

saveRDS(results, file.path(results_folder, paste0(analysis_name, ".rds")))

cat("End:", format(Sys.time()), "\n")
sessionInfo()                            # package versions, for reproducibility
sink(type = "message")
sink()
close(log_con)

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----