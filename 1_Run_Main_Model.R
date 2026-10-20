#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: EWAS of pesticide exposure on placental DNAm (SEPAGES),
# main model, one time point per run (set time_point below: "t2" or "t3").
# Robust linear regression per CpG, run in parallel by chunks.
# Run from the terminal with the model number as argument:
#   Rscript script/Run_Main_Model.R 1

# Models (which_model):
#  1 = DNAm ~ 2,4-D (detected) + main cov 
#  2 = DNAm ~ 3-PBA + main cov 
#  3 = DNAm ~ TCPy + main cov 
#  4 = DNAm ~ ΣDMP + main cov 
#  5 = DNAm ~ ΣDEP + main cov 
#  6 = DNAm ~ ΣDAP + main cov 

# Main cov: maternal age, education, smoking, GA, child ethnicity, child sex,
#           DNA contamination. CT: 5 cell types (Syncytiotrophoblast = reference).
# Continuous exposures: imputed, log2.
# 2,4-D: factor, Not detected (reference) / Detected.

# [INPUT]:  db/subsets/df_main.rds (one dataset, complete cases, exposures at
#           t2 and t3)
#           SEPAGES_20230615_PreprocessedBetas_nooutliers.RData (beta values,
#           rows = CpGs, columns = Basename)

# [OUTPUT]: results/2_EWAS_Results/main_<tp>_<date>/<exposure>_<tp>.rds
#           data.frame, rows = CpGs: Effect, SE, tvalue, pvalue, CpG
#           results/2_EWAS_Results/main_<tp>_<date>/logs/model_<n>_<time>.log


#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
library(MASS)
library(sfsmisc)
library(foreach)
library(doParallel)
source("script/00_functions.R")

metadata_path  <- "db/subsets/df_main.rds"
betas_path     <- "db/plaDNAm_data/SEPAGES_20230615_PreprocessedBetas_nooutliers.RData"

n_cores <- 8

# Model number and time point from the terminal
# (defaults: model 1, t2 when run interactively)
args <- commandArgs(trailingOnly = TRUE)
which_model <- if (length(args) >= 1) as.numeric(args[1]) else 1
time_point  <- if (length(args) >= 2) args[2] else "t2"
stopifnot(which_model %in% 1:6)
stopifnot(time_point %in% c("t2", "t3"))

# One folder per time point and day (all 6 models of a run go together)
results_folder <- file.path("results/2_EWAS_Results",
                            paste0("main_", time_point, "_", format(Sys.Date(), "%Y%m%d")))
dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)

# Log file: everything printed goes to the console and to the log
log_folder <- file.path(results_folder, "logs")
dir.create(log_folder, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_folder, paste0("model_", which_model, "_",
                                         format(Sys.time(), "%Y%m%d_%H%M"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)              # printed output: console + log
sink(log_con, type = "message")          # warnings and errors: log only
cat("Start:", format(Sys.time()), "\n")
cat("Model:", which_model, "| Time point:", time_point, "\n")

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

# Check: same samples in the same order (stops the script if not)
stopifnot(identical(colnames(meth_data), metadata$Basename))
cat("Samples in the analysis:", ncol(meth_data), "\n")

#==============================================================================#----
# PARAMETERS                                                                   #----
#==============================================================================#----

# Exposure (selected for the chosen time point)
stopifnot(which_model %in% 1:6)
if (which_model == 1) {exposure_base <- "pest_24d_bi";       analysis_name <- "24d"}
if (which_model == 2) {exposure_base <- "pest_3pba_ilog2";   analysis_name <- "3pba"}
if (which_model == 3) {exposure_base <- "pest_tcpy_ilog2";   analysis_name <- "tcpy"}
if (which_model == 4) {exposure_base <- "pest_dmp_ms_ilog2"; analysis_name <- "dmp"}
if (which_model == 5) {exposure_base <- "pest_dep_ms_ilog2"; analysis_name <- "dep"}
if (which_model == 6) {exposure_base <- "pest_dap_ms_ilog2"; analysis_name <- "dap"}

exposure      <- paste0(exposure_base, "_", time_point)
analysis_name <- paste0(analysis_name, "_", time_point)

# Covariates
main_cov   <- c("mo_age", "mo_edu", "mo_smoke", "gest_age",
                "ch_eth", "ch_sex", "Meanlog2oddsContamination")
cell_types <- c("Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC")
# Syncytiotrophoblast left out as reference: the 6 proportions sum to 1
covars <- c(main_cov, cell_types)

cat("Exposure:", exposure, "| Analysis:", analysis_name, "\n")
cat("Covariates:", paste(covars, collapse = ", "), "\n")

# Check: all variables exist and missingness
setdiff(c(exposure, covars), colnames(metadata))   # should be character(0)
colSums(is.na(metadata[, c(exposure, covars)]))

# Keep only the variables this model needs (drops the other time point)
metadata <- metadata[, c("mo_id", "Basename", exposure, covars)]
dim(metadata)

n_model <- sum(complete.cases(metadata[, c(exposure, covars)]))
cat("N in the model:", n_model, "\n")

# Name of the exposure coefficient in the model output
# (a factor gets the level appended, e.g. pest_24d_bi_t2Detected)
if (is.factor(metadata[[exposure]])) {
  print(table(metadata[[exposure]], useNA = "ifany"))
  exposure_coef <- paste0(exposure, levels(metadata[[exposure]])[2])
} else {
  exposure_coef <- exposure
}
cat("Exposure coefficient:", exposure_coef, "\n")

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

# Check: which name the function used for the exposure results
names(chunk_results_list[[1]])
results_key <- intersect(c(exposure, exposure_coef), names(chunk_results_list[[1]]))
stopifnot(length(results_key) == 1)       # stops if neither name is found
cat("Results taken from:", results_key, "\n")

results <- do.call(rbind, lapply(chunk_results_list, function(x) x[[results_key]]))
rownames(results) <- results$CpG

# Check: all CpGs returned
cat("CpGs in results:", nrow(results), "\n")
summary(results$pvalue)

saveRDS(results, file.path(results_folder, paste0(analysis_name, ".rds")))

cat("End:", format(Sys.time()), "\n")
sessionInfo()                            # package versions, for reproducibility
sink(type = "message")
sink()
close(log_con)

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----