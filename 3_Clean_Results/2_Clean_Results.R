#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: Remove unreliable probes from the EWAS results (filter_EPIC
# list, exclude_CpGs()) and compute Bonferroni and FDR adjusted p-values on the
# remaining CpGs. Works on any run folder of the new structure:
#   results/2_EWAS_Results/<cohort>_<sample>_<tp>_<date>/<analysis>/<exposure>.rds
#
# Step 1 writes the list of result files to clean to a txt file (one path per
# line). Edit the txt by hand if needed, then step 2 cleans every file in it.
#
# Arguments (optional):
#   1) pattern for the run folders to include (regular expression; default: all)
#      e.g. "sepages_2tp_t2_old" or "^sepages_2tp_t[23]_20261021$"
#   2) "keep" to reuse an existing txt instead of rewriting it
#
#   Rscript script/3_Clean_Results/Clean_Results.R "sepages_2tp_t[23]_old"
#   Rscript script/3_Clean_Results/Clean_Results.R "" keep
#
# [INPUT]:  results/2_EWAS_Results/<run>/<analysis>/<exposure>_<tp>.rds
#           data.frame, rows = CpGs: Effect, SE, tvalue, pvalue, CpG
# [OUTPUT]: results/3_Clean_Results/<run>/<analysis>/<exposure>_<tp>.rds
#           data.frame, rows = CpGs: Effect, SE, tvalue, pvalue, CpG,
#           pval_bn, pval_fdr
#           results/3_Clean_Results/<run>/<analysis>/<exposure>_<tp>_excluded_probes.txt
#           results/3_Clean_Results/<run>/<analysis>/<exposure>_<tp>_cpg_qc_descriptives.txt
#           results/3_Clean_Results/clean_summary_<time>.csv
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_functions.R")

in_dir    <- "results/2_EWAS_Results"
out_dir   <- "results/3_Clean_Results"
list_file <- file.path(out_dir, "files_to_clean.txt")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters to exclude CpG probes
ethnic  <- "GMAF1p"
array   <- "EPIC"
exclude <- c("control_probes", "noncpg_probes", "Sex", "MASK_mapping",
             "MASK_sub30_copy", "MASK_extBase", "MASK_typeINextBaseSwitch",
             "MASK_snp5_common", "MASK_snp5_ethnic",
             "Unrel_450_EPIC_pla_restrict")
# "MASK_snp5_ethnic" -> remove SNPs common in the population given in `ethnic`

#==============================================================================#----
# READ LIST OF FILES TO CLEAN (written by 1_Get_Files_To_Clean.R)             #----
#==============================================================================#----
stopifnot(file.exists(list_file))
files <- readLines(list_file)
files <- files[files != ""]                       # skip empty lines
stopifnot(length(files) > 0)
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) stop("Files in the list not found:\n",
                                    paste(missing_files, collapse = "\n"))
cat("Files to clean:", length(files), "\n")

#==============================================================================#----
# STEP 2: EXCLUDE UNRELIABLE PROBES & ADJUST P-VALUES                          #----
#==============================================================================#----

summary_list <- vector("list", length(files))

for (i in seq_along(files)) {
  
  in_file  <- files[i]
  # Same relative path in 3_Clean_Results: <run>/<analysis>/<exposure>.rds
  out_file <- sub(paste0("^", in_dir), out_dir, in_file)
  out_path <- dirname(out_file)
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
  exp_filename <- tools::file_path_sans_ext(basename(in_file))
  
  cat("\n[", i, "/", length(files), "]", in_file, "\n")
  
  # Load EWAS results
  ewas_results <- readRDS(in_file)
  n_before <- nrow(ewas_results)
  cat("Number CpGs before exclusion:", n_before, "\n")
  
  ewas_results$probeID <- rownames(ewas_results)
  stopifnot(identical(ewas_results$probeID, as.character(ewas_results$CpG)))
  
  # Exclude unreliable probes
  ewas_results <- exclude_CpGs(ewas_results, "probeID",
                               exclude,
                               ethnic     = ethnic,
                               filename   = file.path(out_path, paste0(exp_filename, "_excluded_probes.txt")),
                               fileresume = file.path(out_path, paste0(exp_filename, "_cpg_qc_descriptives.txt")),
                               artype     = array)
  
  # Keep model statistics only, numeric, CpGs as rownames
  cpg_ids      <- ewas_results$probeID
  ewas_results <- ewas_results[, c("Effect", "SE", "tvalue", "pvalue")]
  ewas_results <- as.data.frame(sapply(ewas_results, as.numeric))
  rownames(ewas_results) <- cpg_ids
  ewas_results$CpG <- cpg_ids
  n_after <- nrow(ewas_results)
  cat("Number CpGs after exclusion:", n_after, "| Excluded:", n_before - n_after, "\n")
  
  # Adjusted p-values on the filtered CpGs
  ewas_results$pval_bn  <- p.adjust(ewas_results$pvalue, method = "bonferroni")
  ewas_results$pval_fdr <- p.adjust(ewas_results$pvalue, method = "fdr")
  n_bn  <- sum(ewas_results$pval_bn  < 0.05, na.rm = TRUE)
  n_fdr <- sum(ewas_results$pval_fdr < 0.05, na.rm = TRUE)
  cat("Bonferroni < 0.05:", n_bn, "| FDR < 0.05:", n_fdr,
      "| NA p-values:", sum(is.na(ewas_results$pvalue)), "\n")
  
  saveRDS(ewas_results, out_file)
  
  # Summary row: run / analysis / exposure from the path
  path_parts <- strsplit(in_file, "/")[[1]]
  summary_list[[i]] <- data.frame(
    run        = path_parts[length(path_parts) - 2],
    analysis   = path_parts[length(path_parts) - 1],
    exposure   = exp_filename,
    n_before   = n_before,
    n_after    = n_after,
    n_na_pval  = sum(is.na(ewas_results$pvalue)),
    n_bn_005   = n_bn,
    n_fdr_005  = n_fdr,
    min_pvalue = min(ewas_results$pvalue, na.rm = TRUE)
  )
}

#==============================================================================#----
# SUMMARY                                                                      #----
#==============================================================================#----

clean_summary <- do.call(rbind, summary_list)
print(clean_summary)

# Check: same number of CpGs after exclusion within each run
tapply(clean_summary$n_after, clean_summary$run, function(x) length(unique(x)))   # all 1

write.csv(clean_summary,
          file.path(out_dir, paste0("clean_summary_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv")),
          row.names = FALSE)

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----