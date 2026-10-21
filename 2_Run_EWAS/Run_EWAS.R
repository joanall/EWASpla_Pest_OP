#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: EWAS of pesticide exposure on placental DNAm, for BiSC and
# SEPAGES, main model and sensitivity analyses. One time point per run.
# Robust linear regression per CpG, run in parallel by chunks.
# Outer loop: analyses (main, s1-s5). Inner loop: 6 exposures.
#
# Analyses (main covariates + cell types, except as noted):
#  main = Main model
#  s1   = Excluding pregnancy complications
#  s2   = Excluding non-European (ethnicity removed from covariates)
#  s3   = Adjusting for specific gravity (at the same time point)
#  s4   = Adjusting for fruit and vegetable intake
#  s5   = Without gestational age and without cell types
#
# Exposures: 2,4-D (0 = not detected, 1 = detected), 3-PBA, TCPy (imputed,
#            log2), ΣDMP, ΣDEP, ΣDAP (molar sums, log2). All numeric.
#
# Main covariates: maternal age, education, smoking, GA, child ethnicity,
#   child sex, DNA contamination, 5 cell types (Syncytiotrophoblast = ref).
#   BiSC only: + delivery hospital.
#
# Arguments:
#   1) cohort:     sepages / bisc
#   2) analysis:   main / s1 ... s5 / sens (= s1-s5) / all (= main + s1-s5)
#   3) time point: t2 / t3
#   4) sample:     2tp (women with both time points; default)
#                  tp  (all women with this time point; BiSC only)
#   5) run date:   YYYYMMDD (optional; default today)
#
#   Rscript script/2_Run_Models/Run_EWAS.R sepages all t2
#   Rscript script/2_Run_Models/Run_EWAS.R bisc main t3 tp 20261021
#
# [INPUT]:  subsets (see COHORT SETTINGS) and preprocessed betas
#           (rows = CpGs, columns = Basename)
# [OUTPUT]: results/2_EWAS_Results/<cohort>_<sample>_<tp>_<date>/<analysis>/
#           <exposure>_<tp>.rds: data.frame, rows = CpGs:
#           Effect, SE, tvalue, pvalue, CpG
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
library(MASS)
library(sfsmisc)
library(foreach)
library(doParallel)
source("script/00_functions.R")

n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 8))

# Arguments from the terminal (defaults for interactive use)
args <- commandArgs(trailingOnly = TRUE)
cohort       <- if (length(args) >= 1) args[1] else "sepages"
analysis_arg <- if (length(args) >= 2) args[2] else "main"
time_point   <- if (length(args) >= 3) args[3] else "t2"
sample_set   <- if (length(args) >= 4) args[4] else "2tp"
run_date     <- if (length(args) >= 5) args[5] else format(Sys.Date(), "%Y%m%d")
exposure_arg <- if (length(args) >= 6) args[6] else "all"   # e.g. "24d" or "24d,3pba"

stopifnot(cohort %in% c("sepages", "bisc"))
stopifnot(time_point %in% c("t2", "t3"))
stopifnot(sample_set %in% c("2tp", "tp"))
if (cohort == "sepages" && sample_set == "tp") stop("SEPAGES has only the 2tp sample")

# Analyses to run
if (analysis_arg == "all") {
  which_analyses <- c("main", paste0("s", 1:5))
} else if (analysis_arg == "sens") {
  which_analyses <- paste0("s", 1:5)
} else {
  which_analyses <- analysis_arg
}
stopifnot(all(which_analyses %in% c("main", paste0("s", 1:5))))

#==============================================================================#----
# COHORT SETTINGS                                                              #----
#==============================================================================#----

if (cohort == "sepages") {
  betas_path <- "db/plaDNAm_data/SEPAGES_20230615_PreprocessedBetas_nooutliers.RData"
  subset_dir <- "db/subsets/"
  extra_cov  <- character(0)
}

if (cohort == "bisc") {
  betas_path <- "db/methylation_data/bVals_BISC_nosib_noout_nomask.RDS"  
  subset_dir <- if (sample_set == "2tp") "db/subsets_2tp/" else "db/subsets/"
  extra_cov  <- "hospital"
}

# Subset file per analysis: df_main.rds (2tp) or df_main_t2.rds (tp)
file_suffix <- if (sample_set == "2tp") "" else paste0("_", time_point)
subset_path <- setNames(paste0(subset_dir, "df_", which_analyses, file_suffix, ".rds"),
                        which_analyses)
subset_path
stopifnot(all(file.exists(subset_path)))

#==============================================================================#----
# RESULTS FOLDER & LOG                                                         #----
#==============================================================================#----

results_folder <- file.path("results/2_EWAS_Results",
                            paste(cohort, sample_set, time_point, run_date, sep = "_"))
dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)

log_folder <- file.path(results_folder, "logs")
dir.create(log_folder, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_folder, paste0(paste(which_analyses, collapse = "-"), "_",
                                         format(Sys.time(), "%Y%m%d_%H%M"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)              # printed output: console + log
sink(log_con, type = "message")          # warnings and errors: log only
cat("Start:", format(Sys.time()), "\n")
cat("Cohort:", cohort, "| Sample:", sample_set, "| Time point:", time_point,
    "| Analyses:", paste(which_analyses, collapse = ", "), "\n")

#==============================================================================#----
# PARAMETERS                                                                   #----
#==============================================================================#----

# Exposures, selected for the chosen time point
exposures <- paste0(c("pest_24d_bi", "pest_3pba_ilog2", "pest_tcpy_ilog2",
                      "pest_dmp_ms_ilog2", "pest_dep_ms_ilog2", "pest_dap_ms_ilog2"),
                    "_", time_point)
analysis_name <- paste0(c("24d", "3pba", "tcpy", "dmp", "dep", "dap"), "_", time_point)

# Optional: run only some exposures (6th argument), e.g. "24d" or "24d,3pba"
if (exposure_arg != "all") {
  keep <- sub(paste0("_", time_point, "$"), "", analysis_name) %in%
    strsplit(exposure_arg, ",")[[1]]
  if (!any(keep)) stop("No exposure matches: ", exposure_arg)
  exposures     <- exposures[keep]
  analysis_name <- analysis_name[keep]
}
cat("Exposures:", paste(exposures, collapse = ", "), "\n")

# Covariates
main_cov   <- c("mo_age", "mo_edu", "mo_smoke", "gest_age",
                "ch_eth", "ch_sex", "Meanlog2oddsContamination", extra_cov)
cell_types <- c("Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC")
# Syncytiotrophoblast left out as reference: the 6 proportions sum to 1

analysis_covars <- list(
  main = c(main_cov, cell_types),
  s1   = c(main_cov, cell_types),                                   # excl. complications
  s2   = c(setdiff(main_cov, "ch_eth"), cell_types),                # excl. non-European
  s3   = c(main_cov, cell_types, paste0("pool_sg_", time_point)),   # + SG, same time point
  s4   = c(main_cov, cell_types, "mo_fv_g1000kcal"),                # + fruit & vegetables
  s5   = setdiff(main_cov, "gest_age")                              # without GA and cell types
)

#==============================================================================#----
# LOAD BETAS (ONCE)                                                            #----
#==============================================================================#----

# SEPAGES: .RData (load) | BiSC: .RDS (readRDS)
if (grepl("\\.rds$", betas_path, ignore.case = TRUE)) {
  meth_all <- as.matrix(readRDS(betas_path))
} else {
  betas_obj <- load(betas_path)
  meth_all  <- as.matrix(get(betas_obj[1]))
  rm(list = betas_obj)
}

cat("Number of CpGs:", nrow(meth_all), "| Samples:", ncol(meth_all), "\n")
range(meth_all, na.rm = TRUE)

cat("Number of CpGs:", nrow(meth_all), "| Samples:", ncol(meth_all), "\n")
range(meth_all, na.rm = TRUE)

# meth_all <- meth_all[1:50, , drop = FALSE]   # TEST - comment out for the real run

cl <- makeCluster(n_cores)
registerDoParallel(cl)

#==============================================================================#----
# LOOP: ANALYSES x EXPOSURES                                                   #----
#==============================================================================#----

for (a in which_analyses) {
  
  cat("\n==================== Analysis", a, "====================\n")
  
  # 1. Subset and covariates
  metadata <- readRDS(subset_path[[a]])
  covars   <- analysis_covars[[a]]
  cat("Subset:", subset_path[[a]], "| N:", nrow(metadata), "\n")
  cat("Covariates:", paste(covars, collapse = ", "), "\n")
  
  # Check: all variables exist (stops if not), missingness
  missing_vars <- setdiff(c(exposures, covars), colnames(metadata))
  if (length(missing_vars) > 0) stop("Missing variable(s) in ", a, ": ",
                                     paste(missing_vars, collapse = ", "))
  print(colSums(is.na(metadata[, c(exposures, covars)])))
  
  # Keep only the variables needed (drops the other time point)
  metadata <- metadata[, c("mo_id", "Basename", exposures, covars)]
  
  # Check: exposures numeric, no constant covariate
  not_numeric <- exposures[!sapply(metadata[, exposures], is.numeric)]
  if (length(not_numeric) > 0) stop("Exposure(s) not numeric: ",
                                    paste(not_numeric, collapse = ", "))
  n_levels <- sapply(metadata[, covars], function(x) length(unique(na.omit(x))))
  if (any(n_levels < 2)) stop("Constant covariate(s) in ", a, ": ",
                              paste(names(n_levels)[n_levels < 2], collapse = ", "))
  
  # 2. Match and order betas with this subset
  metadata  <- metadata[metadata$Basename %in% colnames(meth_all), ]
  meth_data <- meth_all[, metadata$Basename, drop = FALSE]
  stopifnot(identical(colnames(meth_data), metadata$Basename))
  cat("Samples in the analysis:", ncol(meth_data), "\n")
  
  # 3. Split CpGs into chunks
  chunk_id   <- cut(seq_len(nrow(meth_data)), breaks = n_cores, labels = FALSE)
  cpg_chunks <- lapply(split(seq_len(nrow(meth_data)), chunk_id),
                       function(idx) meth_data[idx, , drop = FALSE])
  rm(meth_data)
  
  out_folder <- file.path(results_folder, a)
  dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)
  
  # 4. One EWAS per exposure
  for (e in seq_along(exposures)) {
    
    exposure <- exposures[e]
    n_model  <- sum(complete.cases(metadata[, c(exposure, covars)]))
    cat("\n--", a, "|", analysis_name[e], "| N in the model:", n_model, "\n")
    
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
    cat("CpGs in results:", nrow(results), "| Failed CpGs:", sum(is.na(results$pvalue)), "\n")
    
    saveRDS(results, file.path(out_folder, paste0(analysis_name[e], ".rds")))
  }
  rm(cpg_chunks)
}

stopCluster(cl)

cat("\nEnd:", format(Sys.time()), "\n")
sessionInfo()                            # package versions, for reproducibility
sink(type = "message")
sink()
close(log_con)

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----