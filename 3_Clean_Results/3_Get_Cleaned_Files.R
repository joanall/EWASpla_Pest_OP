#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: Write the list of cleaned EWAS result files (one path per
# line) and an index table describing each file, for the downstream scripts
# (results tables, plots, meta-analysis).
#
# Also used for the meta-analysis (cohort = meta): it indexes the meta results
# in results/7_Meta/3_Read_GWAMA, so 4_QC, 5_Annotation and 6_Analyze_DMPs can
# be run on them exactly as on a single cohort.
#
# Arguments:
#   1) cohort: bisc / sepages / meta / bisc_mix / sepages_mix / meta_mix
#   2) pattern for the run folders to include (regular expression; default: all),
#      e.g. "sepages_2tp_t[23]_old"
#   Rscript script/3_Clean_Results/3_Get_Cleaned_Files.R sepages "sepages_2tp_avg_20261021"
#   Rscript script/3_Clean_Results/3_Get_Cleaned_Files.R meta "meta_tp_(t2|t3|avg)_20261106"
#
# [INPUT]:  <clean dir>/<run>/<analysis>/<exposure>_<tp>.rds
#           <clean dir> = results/<cohort>/3_Clean_Results (bisc, sepages)
#                         results/7_Meta/3_Read_GWAMA      (meta)
# [OUTPUT]: <clean dir>/files_cleaned.txt
#           <clean dir>/files_cleaned_index.csv
#           (path, run, cohort, sample, time_point, run_date, analysis, exposure)
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_paths.R")

args <- commandArgs(trailingOnly = TRUE)
cohort      <- if (length(args) >= 1) args[1] else "bisc"
run_pattern <- if (length(args) >= 2 && args[2] != "") args[2] else "."
check_cohort(cohort)

input_dir  <- clean_dir(cohort)
list_file  <- file.path(input_dir, "files_cleaned.txt")
index_file <- file.path(input_dir, "files_cleaned_index.csv")

#==============================================================================#----
# LIST CLEANED FILES                                                           #----
#==============================================================================#----

# Run folders matching the pattern
run_folders <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)
run_folders <- run_folders[grepl(run_pattern, basename(run_folders))]
cat("Run folders:\n"); print(basename(run_folders))
stopifnot(length(run_folders) > 0)

# Cleaned result files inside <run>/<analysis>/ (.rds only, not the QC .txt)
all_files <- character(0)
for (rf in run_folders) {
  for (a in c("main", paste0("s", 1:6))) {
    f <- list.files(file.path(rf, a), pattern = "\\.rds$", full.names = TRUE)
    all_files <- c(all_files, f)
  }
}
stopifnot(length(all_files) > 0)

#==============================================================================#----
# INDEX TABLE                                                                  #----
#==============================================================================#----

run_name <- basename(dirname(dirname(all_files)))     # e.g. sepages_2tp_t2_old
run_bits <- strsplit(run_name, "_")                   # cohort, sample, tp, date

files_index <- data.frame(
  path       = all_files,
  run        = run_name,
  cohort     = sapply(run_bits, function(x) x[1]),
  sample     = sapply(run_bits, function(x) x[2]),
  time_point = sapply(run_bits, function(x) x[3]),
  run_date   = sapply(run_bits, function(x) x[4]),
  analysis   = basename(dirname(all_files)),
  exposure   = sub("_(t2|t3|avg)$", "", tools::file_path_sans_ext(basename(all_files))),
  stringsAsFactors = FALSE
)

# Check: exposure codes without the time point suffix (should be character(0))
setdiff(unique(files_index$exposure), c("24d", "3pba", "tcpy", "dmp", "dep", "dap", "mix"))
stopifnot(all(files_index$exposure %in% c("24d", "3pba", "tcpy", "dmp", "dep", "dap", "mix")))

# Check: exposure time point matches the folder time point
stopifnot(all(grepl(paste0("_", files_index$time_point, ".rds$"), files_index$path)))

# Check: files per run and analysis (6 expected each)
table(files_index$run, files_index$analysis)

#==============================================================================#----
# SAVE                                                                         #----
#==============================================================================#----

writeLines(files_index$path, list_file)
# Plain csv on purpose: read back by later scripts (not a results table)
write.csv(files_index, index_file, row.names = FALSE)
cat("List written:", list_file, "|", nrow(files_index), "files\n")
cat("Index written:", index_file, "\n")

