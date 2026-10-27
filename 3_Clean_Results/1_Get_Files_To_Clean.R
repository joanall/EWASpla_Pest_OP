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
# Argument (optional): pattern for the run folders to include (regular
# expression; default: all), e.g. "sepages_2tp_t[23]_old"
#   Rscript script/3_Clean_Results/1_Get_Files_To_Clean.R "sepages_2tp_t[23]_old"
#
# [INPUT]:  results/2_EWAS_Results/<run>/<analysis>/<exposure>_<tp>.rds
# [OUTPUT]: results/3_Clean_Results/files_to_clean.txt
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
in_dir    <- "results/2_EWAS_Results"
out_dir   <- "results/3_Clean_Results"
list_file <- file.path(out_dir, "files_to_clean.txt")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
run_pattern <- if (length(args) >= 1 && args[1] != "") args[1] else "."

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
  for (a in c("main", paste0("s", 1:5))) {
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