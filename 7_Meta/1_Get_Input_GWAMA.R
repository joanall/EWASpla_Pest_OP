#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date:  October - November 2026
#==============================================================================#
# [DESCRIPTION]: Placental pesticide EWAS - BiSC + SEPAGES fixed-effects meta-analysis
# Builds the GWAMA study files and the GWAMA filelists (.in), one per
# meta-analysis (time point / sensitivity) x exposure.

# [INPUT]: results/bisc/3_Clean_Results/<bisc_run>/<analysis>/<exposure>_<tp>.rds
#         results/sepages/3_Clean_Results/<sepages_run>/<analysis>/<exposure>_<tp>.rds
#         <analysis> = main, s1 ... s6 ; <exposure> = 24d, 3pba, tcpy, dmp, dep, dap
#         Each .rds = data.frame, rows = CpGs: Effect, SE, tvalue, pvalue, CpG,
#         pval_bn, pval_fdr
#         Mixtures (argument meta_mix): results/<cohort>/8_Mixtures/3_Clean_Results/
#         <run>/<analysis>/mix_<tp>.rds
# [OUTPUT]: results/7_Meta/1_Input_GWAMA/<meta_id>/   (meta_id = <tp>_<analysis>, e.g. t2_main)
#           (meta_mix: results/8_Mixtures/7_Meta/1_Input_GWAMA/<meta_id>/)
#   Rscript script/7_Meta/1_Get_Input_GWAMA.R            # EWAS
#   Rscript script/7_Meta/1_Get_Input_GWAMA.R meta_mix   # mixture models
#           bisc_gwamain_<exposure>.txt
#           sepages_gwamain_<exposure>.txt
#           gwama_filelist_<exposure>.in
#         GWAMA study file columns: MARKERNAME, EA, NEA, EAF, BETA, SE, P, N
#
# [NOTE]  : GWAMA is allele-based. For an EWAS, every CpG gets the same dummy
#         alleles (EA = "A", NEA = "G") and EAF = 0.5, so no allele flipping
#         happens and BETA is meta-analysed as it is in both cohorts.
#         Both cohorts must therefore share the same exposure scale (log2).
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_paths.R")

# Argument: meta (single-exposure EWAS, default) or meta_mix (mixture models)
args       <- commandArgs(trailingOnly = TRUE)
meta_type  <- if (length(args) >= 1) args[1] else "meta"
check_cohort(meta_type, c("meta", "meta_mix"))

# EWAS: results/<cohort>/3_Clean_Results; mixtures: results/<cohort>/8_Mixtures/3_Clean_Results
cohort_bisc    <- if (meta_type == "meta_mix") "bisc_mix"    else "bisc"
cohort_sepages <- if (meta_type == "meta_mix") "sepages_mix" else "sepages"
dir_bisc    <- clean_dir(cohort_bisc)
dir_sepages <- clean_dir(cohort_sepages)   # copied from the SEPAGES PC
dir_input   <- results_dir(meta_type, "1_Input_GWAMA")
dir_bisc; dir_sepages; dir_input
stopifnot(dir.exists(dir_bisc), dir.exists(dir_sepages))

#------------------------------------------------------------------------------#
# Design: one row per time point
# - *_run = run folder inside 3_Clean_Results (<cohort>_<sample>_<tp>_<date>)
#   Every analysis folder (main, s1 ... s6) found in BOTH runs is meta-analysed
# - n_*   = sample size of the MAIN model in each cohort. GWAMA only uses N for
#   the n_samples column (the fixed-effects estimate uses SE only), so for the
#   sensitivity analyses with fewer women (s1, s2) n_samples will be too high
#------------------------------------------------------------------------------#
design <- data.frame(
  meta_id     = c("t2",                      "t3",                      "avg"),
  bisc_run    = c("bisc_tp_t2_20261022",     "bisc_tp_t3_20261022",     "bisc_tp_avg_20261022"),
  sepages_run = c("sepages_2tp_t2_20261021", "sepages_2tp_t3_20261021", "sepages_2tp_avg_20261021"),
  n_bisc      = c(378,                       437,                       503),
  n_sepages   = c(312,                       312,                       312),
  stringsAsFactors = FALSE
)

# Mixture models (meta_mix): run folders of 8_Mixtures/4_Run_Gcomp and N of the
# mixture analytic sample (complete cases for the 4 compounds) - FILL IN
design_mix <- data.frame(
  meta_id     = c("t2",                                "t3",                                "avg"),
  bisc_run    = c("bisc_tp_t2_YYYYMMDD_cohort",        "bisc_tp_t3_YYYYMMDD_cohort",        "bisc_tp_avg_YYYYMMDD_cohort"),
  sepages_run = c("sepages_2tp_t2_YYYYMMDD_cohort",    "sepages_2tp_t3_YYYYMMDD_cohort",    "sepages_2tp_avg_YYYYMMDD_cohort"),
  n_bisc      = c(NA,                                  NA,                                  NA),
  n_sepages   = c(NA,                                  NA,                                  NA),
  stringsAsFactors = FALSE
)
if (meta_type == "meta_mix") design <- design_mix
design

stopifnot(!anyDuplicated(design$meta_id))
stopifnot(!is.na(design$n_bisc), !is.na(design$n_sepages))   # all N filled in

#==============================================================================#----
# FUNCTIONS                                                                    #----
#==============================================================================#----
# List the cleaned result files of one run folder (files are read later, one
# at a time, to keep memory low). Returns: analysis, exposure, key, path
list_clean_run <- function(dir_clean, run) {
  dir_run <- file.path(dir_clean, run)
  if (!dir.exists(dir_run)) stop("Run folder not found: ", dir_run)
  analyses <- list.dirs(dir_run, recursive = FALSE, full.names = FALSE)
  analyses <- analyses[analyses %in% c("main", paste0("s", 1:7))]   # skips logs/
  stopifnot(length(analyses) > 0)
  
  files_run <- data.frame(analysis = character(0), exposure = character(0),
                          path = character(0), stringsAsFactors = FALSE)
  for (a in analyses) {
    files_rds <- list.files(file.path(dir_run, a), pattern = "\\.rds$", full.names = TRUE)
    if (length(files_rds) == 0) {                         # folder exists but is empty
      cat("   skipped (no .rds files):", file.path(run, a), "\n")
      next
    }
    exposures <- sub("_(t2|t3|avg)$", "", tools::file_path_sans_ext(basename(files_rds)))
    files_run <- rbind(files_run, data.frame(analysis = a, exposure = exposures,
                                             path = files_rds, stringsAsFactors = FALSE))
  }
  files_run$key <- paste(files_run$analysis, files_run$exposure, sep = "_")
  
  # Check: files per analysis (6 exposures expected)
  cat("  ", run, "\n"); print(table(files_run$analysis))
  files_run
}

# Convert one EWAS result table to GWAMA study format
format_gwama <- function(res, n_cohort) {
  names(res)[names(res) == "CpG"]       <- "probeID"
  if (!"probeID" %in% names(res)) res$probeID <- rownames(res)
  names(res)[names(res) == "Effect"]    <- "BETA"
  names(res)[names(res) == "pvalue"]    <- "P"
  names(res)[names(res) == "estimate"]  <- "BETA"
  names(res)[names(res) == "std.error"] <- "SE"
  names(res)[names(res) == "p.value"]   <- "P"
  names(res)[names(res) == "P_VAL"]     <- "P"
  stopifnot(all(c("probeID", "BETA", "SE", "P") %in% names(res)))
  
  gwama_in <- data.frame(MARKERNAME = res$probeID,
                         EA         = "A",
                         NEA        = "G",
                         EAF        = 0.5,
                         BETA       = res$BETA,
                         SE         = res$SE,
                         P          = res$P,
                         N          = n_cohort,
                         stringsAsFactors = FALSE)
  n_before <- nrow(gwama_in)
  gwama_in <- gwama_in[complete.cases(gwama_in) & gwama_in$SE > 0, ]
  if (nrow(gwama_in) < n_before) cat("    removed", n_before - nrow(gwama_in), "CpGs with NA or SE <= 0\n")
  stopifnot(!anyDuplicated(gwama_in$MARKERNAME))
  gwama_in
}

#==============================================================================#----
# BUILD GWAMA INPUTS                                                           #----
#==============================================================================#----
for (i in seq_len(nrow(design))) {
  tp <- design$meta_id[i]
  cat("\n=== Time point:", tp, "===\n")
  
  files_bisc    <- list_clean_run(dir_bisc,    design$bisc_run[i])
  files_sepages <- list_clean_run(dir_sepages, design$sepages_run[i])
  
  # Analysis x exposure combinations present in both cohorts
  keys_shared <- intersect(files_bisc$key, files_sepages$key)
  if (!setequal(files_bisc$key, files_sepages$key)) {
    warning(tp, ": not in both cohorts - BiSC only: ",
            paste(setdiff(files_bisc$key, files_sepages$key), collapse = ", "),
            " | SEPAGES only: ",
            paste(setdiff(files_sepages$key, files_bisc$key), collapse = ", "))
  }
  stopifnot(length(keys_shared) > 0)
  
  for (key in keys_shared) {
    analysis <- files_bisc$analysis[files_bisc$key == key]
    exposure <- files_bisc$exposure[files_bisc$key == key]
    meta_id  <- paste(tp, analysis, sep = "_")                # e.g. t2_main, t2_s1
    
    dir_out <- file.path(dir_input, meta_id)
    dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
    
    res_bisc    <- readRDS(files_bisc$path[files_bisc$key == key])
    res_sepages <- readRDS(files_sepages$path[files_sepages$key == key])
    
    gwama_bisc    <- format_gwama(res_bisc,    design$n_bisc[i])
    gwama_sepages <- format_gwama(res_sepages, design$n_sepages[i])
    
    file_bisc     <- file.path(dir_out, paste0("bisc_gwamain_",    exposure, ".txt"))
    file_sepages  <- file.path(dir_out, paste0("sepages_gwamain_", exposure, ".txt"))
    file_filelist <- file.path(dir_out, paste0("gwama_filelist_",  exposure, ".in"))
    
    write.table(gwama_bisc,    file_bisc,    sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(gwama_sepages, file_sepages, sep = "\t", quote = FALSE, row.names = FALSE)
    writeLines(c(file_bisc, file_sepages), file_filelist)
    
    # Inline check: CpG overlap between cohorts
    n_shared <- length(intersect(gwama_bisc$MARKERNAME, gwama_sepages$MARKERNAME))
    cat(sprintf("  %-8s %-5s BiSC: %7d | SEPAGES: %7d | shared: %7d\n",
                meta_id, exposure, nrow(gwama_bisc), nrow(gwama_sepages), n_shared))
  }
}

# Final check: one .in per meta-analysis x exposure
filelists <- list.files(dir_input, pattern = "\\.in$", recursive = TRUE)
cat("\n", length(filelists), "GWAMA filelists written in", dir_input, "\n")
print(filelists)

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----