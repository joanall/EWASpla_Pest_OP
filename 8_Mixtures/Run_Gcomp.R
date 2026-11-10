#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the
# BiSC and SEPAGES cohorts.
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: Mixture EWAS with quantile g-computation (qgcomp.noboot),
# one model per CpG.
#
# Mixture (4 compounds): 3-PBA, TCPy (imputed, log2), ΣDMP, ΣDEP (molar sums, log2).
#   2,4-D not included; ΣDAP not included (= ΣDMP + ΣDEP).
#   avg = mean of the log2 values at t2 and t3.
# Exposures are quantized ONCE per analysis (q quantiles, on the analytic
# sample of complete cases), then qgcomp.noboot is fitted per CpG with q = NULL.
#
# Mixture effect (psi): change in methylation (beta) for a one-quantile increase
# in all 4 compounds at the same time.
#
# Analyses (main covariates + cell types, except as noted):
#  main = Main model
#  s1   = Excluding pregnancy complications
#  s2   = Excluding non-European (ethnicity removed from covariates)
#  s3   = Adjusting for specific gravity (at the same time point)
#  s4   = Adjusting for fruit and vegetable intake
#  s5   = Without gestational age and without cell types
#  (s6, winsorized exposures, is NOT run: winsorizing at P1/P99 does not change
#   the quantile categories, so the result would be identical to main)
#
# Main covariates: maternal age, education, smoking, GA, child ethnicity,
#   child sex, DNA contamination, 5 cell types (Syncytiotrophoblast = ref).
#   BiSC only: + delivery hospital.
#
# Arguments:
#   1) cohort:     sepages / bisc
#   2) analysis:   main / s1 ... s5 / sens (= s1-s5) / all (= main + s1-s5)
#   3) time point: t2 / t3 / avg
#   4) sample:     2tp (women with both time points; default)
#                  tp  (all women with this time point; BiSC only)
#   5) run date:   YYYYMMDD (optional; default today)
#   6) cut-offs:   cohort (cohort-specific quantiles; default)
#                  shared (pooled BiSC + SEPAGES cut-offs, from
#                          results/2_Gcomp_Results/quantile_cutoffs/shared_bisc_sepages_q<q>.rds)
#
# Example:
#   Rscript script/8_Mixtures/Run_Gcomp.R sepages main t3 2tp 20261109
#   Rscript script/8_Mixtures/Run_Gcomp.R bisc main t2 tp 20261110 shared
#
# [INPUT]:  subsets (see COHORT SETTINGS) and preprocessed betas
#           (rows = CpGs, columns = Basename)
# [OUTPUT]: results/2_Gcomp_Results/<cohort>_<sample>_<tp>_<date>/<analysis>/
#           mix_<tp>.rds: data.frame, rows = CpGs:
#           Effect (psi), SE, tvalue, pvalue, CpG,
#           beta_3pba, beta_tcpy, beta_dmp, beta_dep (per-compound coefficients;
#           weights = beta / sum of the betas with the same sign)
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

# install.packages("qgcomp")   # once, not inside the script
library(qgcomp)
library(foreach)
library(doParallel)

n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 8))
# One BLAS/OpenMP thread per worker (avoids nested parallelism; workers inherit this)
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

# Arguments from the terminal (defaults for interactive use)
args <- commandArgs(trailingOnly = TRUE)
cohort       <- if (length(args) >= 1) args[1] else "sepages"
analysis_arg <- if (length(args) >= 2) args[2] else "main"
time_point   <- if (length(args) >= 3) args[3] else "t2"
sample_set   <- if (length(args) >= 4) args[4] else "2tp"
run_date     <- if (length(args) >= 5) args[5] else format(Sys.Date(), "%Y%m%d")
cutoff_type  <- if (length(args) >= 6) args[6] else "cohort"   # cohort / shared
stopifnot(cutoff_type %in% c("cohort", "shared"))

stopifnot(cohort %in% c("sepages", "bisc"))
stopifnot(time_point %in% c("t2", "t3", "avg"))
stopifnot(sample_set %in% c("2tp", "tp"))
if (cohort == "sepages" && sample_set == "tp") stop("SEPAGES has only the 2tp sample")

# Analyses to run (no s6: see description)
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
# Separate folder from the single-exposure EWAS, so 3_Clean_Results scripts
# (which expect the 6 single exposures) do not pick these files up

results_folder <- file.path("results/2_Gcomp_Results",
                            paste(cohort, sample_set, time_point, run_date, cutoff_type, sep = "_"))
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
    "| Analyses:", paste(which_analyses, collapse = ", "),
    "| Cut-offs:", cutoff_type, "| Cores:", n_cores, "\n")

#==============================================================================#----
# PARAMETERS                                                                   #----
#==============================================================================#----

# Mixture components, selected for the chosen time point
mix_exposures <- paste0(c("pest_3pba_ilog2", "pest_tcpy_ilog2",
                          "pest_dmp_ms_ilog2", "pest_dep_ms_ilog2"),
                        "_", time_point)
mix_labels    <- c("3pba", "tcpy", "dmp", "dep")
names(mix_labels) <- mix_exposures
analysis_name <- paste0("mix_", time_point)

q_mix <- 3   # number of quantiles (3 = tertiles, 4 = quartiles)

cat("Mixture:", paste(mix_exposures, collapse = ", "), "| q =", q_mix, "\n")

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
# FUNCTION: QGCOMP PER CpG                                                     #----
#==============================================================================#----
# metadata: exposures already quantized (0, 1, ..., q-1), complete cases only

gcomp_single_cpg <- function(cpg_chunk, mix_exposures, covariates, metadata) {
  
  df_exp_cov <- metadata[, c(mix_exposures, covariates), drop = FALSE]
  f <- as.formula(paste("meth ~", paste(c(mix_exposures, covariates), collapse = " + ")))
  
  # Pre-allocate results: psi + per-compound coefficients
  n_cpg <- nrow(cpg_chunk)
  mix_results <- data.frame(Effect = rep(NA_real_, n_cpg), SE = NA_real_,
                            tvalue = NA_real_, pvalue = NA_real_,
                            CpG = rownames(cpg_chunk), stringsAsFactors = FALSE)
  beta_cols <- paste0("beta_", mix_labels[mix_exposures])
  for (b in beta_cols) mix_results[[b]] <- NA_real_
  
  for (i in seq_len(n_cpg)) {
    
    df_for_model <- df_exp_cov
    df_for_model$meth <- as.numeric(cpg_chunk[i, , drop = FALSE])
    
    # A failing CpG gets NA instead of stopping the chunk
    fit <- tryCatch(qgcomp::qgcomp.noboot(f, expnms = mix_exposures,
                                          data = df_for_model,
                                          q = NULL,             # already quantized
                                          family = gaussian()),
                    error = function(e) NULL)
    if (is.null(fit)) next
    
    psi    <- unname(fit$psi)
    se_psi <- sqrt(unname(fit$var.psi))
    mix_results$Effect[i] <- psi
    mix_results$SE[i]     <- se_psi
    mix_results$tvalue[i] <- psi / se_psi
    mix_results$pvalue[i] <- 2 * pt(-abs(psi / se_psi), df = fit$fit$df.residual)
    mix_results[i, beta_cols] <- fit$fit$coefficients[mix_exposures]
  }
  
  return(mix_results)
}

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

# meth_all <- meth_all[1:50, , drop = FALSE]   # TEST - comment out for the real run

cl <- makeCluster(n_cores)
registerDoParallel(cl)

#==============================================================================#----
# LOOP: ANALYSES                                                               #----
#==============================================================================#----

for (a in which_analyses) {
  
  cat("\n==================== Analysis", a, "====================\n")
  
  # 1. Subset and covariates
  metadata <- readRDS(subset_path[[a]])
  covars   <- analysis_covars[[a]]
  cat("Subset:", subset_path[[a]], "| N:", nrow(metadata), "\n")
  cat("Covariates:", paste(covars, collapse = ", "), "\n")
  
  # Check: all variables exist (stops if not), missingness
  missing_vars <- setdiff(c(mix_exposures, covars), colnames(metadata))
  if (length(missing_vars) > 0) stop("Missing variable(s) in ", a, ": ",
                                     paste(missing_vars, collapse = ", "))
  print(colSums(is.na(metadata[, c(mix_exposures, covars)])))
  
  # Keep only the variables needed, complete cases only
  # (all 4 compounds + covariates: quantiles computed on the analytic sample)
  metadata <- metadata[, c("mo_id", "Basename", mix_exposures, covars)]
  metadata <- metadata[complete.cases(metadata[, c(mix_exposures, covars)]), ]
  
  # Check: exposures numeric, no constant covariate
  not_numeric <- mix_exposures[!sapply(metadata[, mix_exposures], is.numeric)]
  if (length(not_numeric) > 0) stop("Exposure(s) not numeric: ",
                                    paste(not_numeric, collapse = ", "))
  n_levels <- sapply(metadata[, covars], function(x) length(unique(na.omit(x))))
  if (any(n_levels < 2)) stop("Constant covariate(s) in ", a, ": ",
                              paste(names(n_levels)[n_levels < 2], collapse = ", "))
  
  # 2. Match and order betas with this subset
  metadata  <- metadata[metadata$Basename %in% colnames(meth_all), ]
  meth_data <- meth_all[, metadata$Basename, drop = FALSE]
  stopifnot(identical(colnames(meth_data), metadata$Basename))
  cat("Samples in the analysis (complete cases):", ncol(meth_data), "\n")
  
  # 3. Quantize the mixture components once
  if (cutoff_type == "cohort") {
    # Cohort-specific quantiles, on this analytic sample
    mix_q <- qgcomp::quantize(metadata, expnms = mix_exposures, q = q_mix)
    metadata_model <- mix_q$data
    cat("Quantile cut points (cohort-specific):\n"); print(mix_q$breaks)
  } else {
    # Shared BiSC + SEPAGES cut-offs (pooled main analytic samples)
    shared_path <- paste0("results/2_Gcomp_Results/quantile_cutoffs/shared_bisc_sepages_q", q_mix, ".rds")
    stopifnot(file.exists(shared_path))
    shared_cutoffs <- readRDS(shared_path)[[time_point]]
    stopifnot(all(mix_exposures %in% names(shared_cutoffs)))
    metadata_model <- metadata
    for (e in mix_exposures) {
      stopifnot(length(shared_cutoffs[[e]]) == q_mix - 1)
      # Same convention as qgcomp::quantize: categories 0, 1, ..., q-1
      metadata_model[[e]] <- cut(metadata[[e]], breaks = c(-Inf, shared_cutoffs[[e]], Inf),
                                 labels = FALSE, include.lowest = TRUE) - 1
    }
    cat("Quantile cut points (shared):\n"); print(shared_cutoffs[mix_exposures])
  }
  
  # Check: number of women per quantile for each compound
  for (e in mix_exposures) {
    cat(" ", mix_labels[e], ":"); print(table(metadata_model[[e]], useNA = "ifany"))
  }
  # Check: tied values (e.g. many values at the imputed LOD) give fewer than q categories
  n_cat <- sapply(metadata_model[, mix_exposures], function(x) length(unique(x)))
  if (any(n_cat < q_mix)) cat("WARNING: fewer than", q_mix, "categories in:",
                              paste(names(n_cat)[n_cat < q_mix], collapse = ", "), "\n")
  
  # Check: correlation between the compounds (raw values)
  print(round(cor(metadata[, mix_exposures], method = "spearman"), 2))
  
  # 4. Split CpGs into chunks
  chunk_id   <- cut(seq_len(nrow(meth_data)), breaks = n_cores, labels = FALSE)
  cpg_chunks <- lapply(split(seq_len(nrow(meth_data)), chunk_id),
                       function(idx) meth_data[idx, , drop = FALSE])
  rm(meth_data)
  
  out_folder <- file.path(results_folder, a)
  dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)
  
  # 5. One qgcomp model per CpG
  cat("\n--", a, "|", analysis_name, "| N in the model:", nrow(metadata_model), "\n")
  
  start_time <- proc.time()
  chunk_results_list <- foreach(chunk = cpg_chunks,
                                .packages = c("qgcomp"),
                                .export   = c("gcomp_single_cpg", "mix_exposures", "mix_labels",
                                              "covars", "metadata_model")) %dopar% {
                                                gcomp_single_cpg(cpg_chunk     = chunk,
                                                                 mix_exposures = mix_exposures,
                                                                 covariates    = covars,
                                                                 metadata      = metadata_model)
                                              }
  run_time <- proc.time() - start_time
  cat("Run time (min):", round(run_time["elapsed"] / 60, 1), "\n")
  
  results <- do.call(rbind, chunk_results_list)
  rownames(results) <- results$CpG
  cat("CpGs in results:", nrow(results), "| Failed CpGs:", sum(is.na(results$pvalue)), "\n")
  
  # Check: psi = sum of the per-compound coefficients (linear model, no interactions)
  beta_cols <- paste0("beta_", mix_labels)
  stopifnot(isTRUE(all.equal(results$Effect, rowSums(results[, beta_cols]),
                             check.attributes = FALSE)))
  
  saveRDS(results, file.path(out_folder, paste0(analysis_name, ".rds")))
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