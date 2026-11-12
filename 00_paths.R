#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the
# BiSC and SEPAGES cohorts.
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: All input (db/) and output (results/) locations of the
# project, in one place. The same file is used on both computers, so no script
# needs its paths changed when moving between the BiSC and the SEPAGES PC.
# Load at the top of every script with:
# source("script/00_paths.R")
#
# Every script runs from the project root:
#   <project>/script/    this repository
#   <project>/db/        input data (only the cohort available on that PC)
#   <project>/results/   outputs (same layout on both PCs)
#
# results/ layout:
#   results/bisc/<step>/      BiSC single-cohort results
#   results/sepages/<step>/   SEPAGES single-cohort results
#   results/7_Meta/<step>/    BiSC + SEPAGES meta-analysis (all its steps,
#                             including steps 4-6 re-run on the meta results)
#   results/<cohort>/8_Mixtures/<step>/  mixture models (qgcomp) of one cohort,
#                             including steps 3-6 re-run on the mixture results
#   results/8_Mixtures/<step>/  both cohorts' mixture models: 3_Shared_Cutoffs
#                             (pooled cut-offs) and 7_Meta/ (meta-analysis of
#                             the mixture results, same steps as results/7_Meta)
#   db/annotation_data/cpg_annotation_all.rds   CpG annotation (all cohorts)
#
# <step> = name of the script folder that writes it, so results and scripts
# follow the same numbering:
#   1_Get_Descriptives, 2_Run_EWAS, 3_Clean_Results, 4_QC, 5_Annotation,
#   6_Analyze_DMPs, 8_Mixtures/<n>_<script>
#   7_Meta: 1_Input_GWAMA, 2_Run_GWAMA, 3_Read_GWAMA, 4_QC, 5_Annotation,
#           6_Analyze_DMPs, 7_Forest_Meta, 8_Forest_Meta_Sens, 9_Check_CI
#
# "cohort" argument of the scripts (where the results of each analysis live):
#   bisc, sepages        single-exposure EWAS          results/<cohort>/
#   meta                 meta-analysis of the EWAS     results/7_Meta/
#   bisc_mix, sepages_mix  mixture models              results/<cohort>/8_Mixtures/
#   meta_mix             meta-analysis of the mixtures results/8_Mixtures/7_Meta/
#
# Copying the SEPAGES results to the BiSC PC = copying results/sepages/ as is.
#==============================================================================#

#==============================================================================#----
# COHORTS                                                                      #----
#==============================================================================#----
cohorts_all <- c("bisc", "sepages", "meta", "bisc_mix", "sepages_mix", "meta_mix")
cohorts_mix <- c("bisc_mix", "sepages_mix")   # one cohort, mixture models

# Stop with a clear message if the cohort argument is not valid
check_cohort <- function(cohort, allowed = cohorts_all) {
  if (length(cohort) != 1 || !cohort %in% allowed) {
    stop("cohort must be one of: ", paste(allowed, collapse = ", "),
         " (got: ", paste(cohort, collapse = ", "), ")")
  }
  invisible(cohort)
}

cohort_labels <- c(bisc = "BiSC", sepages = "SEPAGES", meta = "Meta-analysis",
                   bisc_mix = "BiSC - mixture", sepages_mix = "SEPAGES - mixture",
                   meta_mix = "Meta-analysis - mixture")

#==============================================================================#----
# RESULTS                                                                      #----
#==============================================================================#----

# Root folder of one cohort's results
results_root <- function(cohort) {
  if (cohort == "meta")     return("results/7_Meta")
  if (cohort == "meta_mix") return("results/8_Mixtures/7_Meta")
  if (cohort %in% cohorts_mix) return(file.path("results", sub("_mix$", "", cohort), "8_Mixtures"))
  file.path("results", cohort)
}

# Raw (uncleaned) per-CpG results: 2_Run_EWAS, or 4_Run_Gcomp for the mixtures
raw_dir <- function(cohort) {
  check_cohort(cohort, c("bisc", "sepages", cohorts_mix))
  if (cohort %in% cohorts_mix) return(results_dir(cohort, "4_Run_Gcomp"))
  results_dir(cohort, "2_Run_EWAS")
}

# Results folder of one step, e.g. results_dir("bisc", "3_Clean_Results")
# -> "results/bisc/3_Clean_Results". create = TRUE also creates it.
results_dir <- function(cohort, ..., create = FALSE) {
  folder <- file.path(results_root(cohort), ...)
  if (create) dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  folder
}

# Cleaned per-CpG results that feed 4_QC, 5_Annotation and 6_Analyze_DMPs:
# 3_Clean_Results for a cohort, 3_Read_GWAMA for the meta-analyses
clean_dir <- function(cohort) {
  check_cohort(cohort)
  if (cohort %in% c("meta", "meta_mix")) return(results_dir(cohort, "3_Read_GWAMA"))
  results_dir(cohort, "3_Clean_Results")
}

# CpG annotation (same for every cohort; written by 5_Annotation/1_Get_Annotation.R)
annotation_file <- "db/annotation_data/cpg_annotation_all.rds"

# Shared BiSC + SEPAGES quantile cut-offs for the mixture models
# (written by 8_Mixtures/3_Join_Grid_Sepages_Bisc.R, read by 4_Run_Gcomp.R)
shared_cutoffs_dir <- "results/8_Mixtures/3_Shared_Cutoffs"

#==============================================================================#----
# INPUT DATA (db/)                                                             #----
#==============================================================================#----
# Only available on the PC of that cohort: scripts that need them check
# file.exists() and stop on the other PC.

# Preprocessed betas (rows = CpGs, columns = Basename)
db_betas <- function(cohort) {
  check_cohort(cohort, c("bisc", "sepages"))
  if (cohort == "sepages") return("db/plaDNAm_data/SEPAGES_20230615_PreprocessedBetas_nooutliers.RData")
  "db/methylation_data/bVals_BISC_nosib_noout_nomask.RDS"
}

# Folder of the analysis subsets
#   SEPAGES: db/subsets/      (only 2tp: women with both time points)
#   BiSC:    db/subsets_2tp/  (2tp) or db/subsets/ (tp: all women with that time point)
db_subset_dir <- function(cohort, sample_set = "2tp") {
  check_cohort(cohort, c("bisc", "sepages"))
  stopifnot(sample_set %in% c("2tp", "tp"))
  if (cohort == "sepages") {
    if (sample_set != "2tp") stop("SEPAGES has only the 2tp sample")
    return("db/subsets/")
  }
  if (sample_set == "2tp") "db/subsets_2tp/" else "db/subsets/"
}

# One subset file: df_<analysis>.rds (2tp) or df_<analysis>_<tp>.rds (tp)
db_subset <- function(cohort, analysis = "main", sample_set = "2tp", time_point = NULL) {
  suffix <- if (sample_set == "2tp") "" else paste0("_", time_point)
  if (sample_set == "tp" && is.null(time_point)) stop("time_point needed for the tp sample")
  paste0(db_subset_dir(cohort, sample_set), "df_", analysis, suffix, ".rds")
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----
