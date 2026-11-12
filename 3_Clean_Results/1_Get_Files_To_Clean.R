#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date:  October - November 2026
#==============================================================================#
# [DESCRIPTION]: Write the list of EWAS result files to clean to a txt file
# (one path per line). The txt can be edited by hand before cleaning.
#
# Arguments:
#   1) cohort: bisc / sepages (EWAS) or bisc_mix / sepages_mix (mixture models)
#   2) pattern for the run folders to include (regular expression; default: all),
#      e.g. "sepages_2tp_t[23]_old"
#   3) analyses to include, comma-separated (default: main, s1 ... s6), e.g. "s6"
#   Rscript script/3_Clean_Results/1_Get_Files_To_Clean.R sepages "sepages_2tp_t[23]_old"
#   Rscript script/3_Clean_Results/1_Get_Files_To_Clean.R bisc "^bisc_tp_(t2|t3|avg)_20261022$" s6
#   Rscript script/3_Clean_Results/1_Get_Files_To_Clean.R sepages_mix "sepages_2tp_(t2|t3|avg)_20261109"
#
# [INPUT]:  results/<cohort>/2_Run_EWAS/<run>/<analysis>/<exposure>_<tp>.rds
#           (mixtures: results/<cohort>/8_Mixtures/4_Run_Gcomp/<run>/<analysis>/mix_<tp>.rds)
# [OUTPUT]: results/<cohort>/3_Clean_Results/files_to_clean.txt
#           (mixtures: results/<cohort>/8_Mixtures/3_Clean_Results/files_to_clean.txt)
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_paths.R")

args <- commandArgs(trailingOnly = TRUE)
cohort      <- if (length(args) >= 1) args[1] else "bisc"
run_pattern <- if (length(args) >= 2 && args[2] != "") args[2] else "."
check_cohort(cohort, c("bisc", "sepages", cohorts_mix))

in_dir    <- raw_dir(cohort)   # 2_Run_EWAS (or 8_Mixtures/4_Run_Gcomp)
out_dir   <- results_dir(cohort, "3_Clean_Results", create = TRUE)
list_file <- file.path(out_dir, "files_to_clean.txt")

# Optional: analyses to include (comma-separated), e.g. "s6" or "main,s6"; default all
analyses_all <- c("main", paste0("s", 1:7))
analyses <- if (length(args) >= 3 && args[3] != "") strsplit(args[3], ",")[[1]] else analyses_all
setdiff(analyses, analyses_all)    # should be character(0)
stopifnot(all(analyses %in% analyses_all))
cat("Analyses:", analyses, "\n")

#==============================================================================#----
# FOLDER TO CLEAN                                                              #----
#==============================================================================#----
# Run folders matching the pattern
run_folders <- list.dirs(in_dir, recursive = FALSE, full.names = TRUE)
run_folders <- run_folders[grepl(run_pattern, basename(run_folders))]
cat("Run folders:\n"); print(basename(run_folders))
stopifnot(length(run_folders) > 0)

# All result files inside <run>/<analysis>/ (logs folders are not included)
all_files <- character(0)
for (rf in run_folders) {
  for (a in analyses) {
    f <- list.files(file.path(rf, a), pattern = "\\.rds$", full.names = TRUE)
    all_files <- c(all_files, f)
  }
}

# Check: files per run and analysis (6 expected each)
table(basename(dirname(dirname(all_files))), basename(dirname(all_files)))

writeLines(all_files, list_file)
cat("List written:", list_file, "|", length(all_files), "files\n")

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----