#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: Write the list of annotated EWAS result files (one path per
# line) and an index table describing each file, for the downstream scripts
# (results tables, plots, meta-analysis).
#
# Argument (optional): pattern for the run folders to include (regular
# expression; default: all), e.g. "sepages_2tp_t[23]_old"
#   Rscript script/5_Annotation/3_Get_Annotated_Files.R "sepages_2tp_t[23]_old"
#
# [INPUT]:  results/5_Annotation/<run>/<analysis>/<exposure>_<tp>_annotated.rds
# [OUTPUT]: results/5_Annotation/files_annotated.txt
#           results/5_Annotation/files_annotated_index.csv
#           (path, run, cohort, sample, time_point, run_date, analysis, exposure)

#==============================================================================#


#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

annot_dir  <- "results/5_Annotation"
list_file  <- file.path(annot_dir, "files_annotated.txt")
index_file <- file.path(annot_dir, "files_annotated_index.csv")
args <- commandArgs(trailingOnly = TRUE)
run_pattern <- if (length(args) >= 1 && args[1] != "") args[1] else "."

#==============================================================================#----
# LIST CLEANED FILES                                                           #----
#==============================================================================#----

# Run folders matching the pattern
run_folders <- list.dirs(annot_dir, recursive = FALSE, full.names = TRUE)
run_folders <- run_folders[grepl(run_pattern, basename(run_folders))]
cat("Run folders:\n"); print(basename(run_folders))
stopifnot(length(run_folders) > 0)

# Cleaned result files inside <run>/<analysis>/ (.rds)
all_files <- character(0)
for (rf in run_folders) {
  for (a in c("main", paste0("s", 1:5))) {
    f <- list.files(file.path(rf, a), pattern = "_annotated\\.rds$", full.names = TRUE)
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
  exposure   = sub("_t[23]_annotated$", "", tools::file_path_sans_ext(basename(all_files))),
  stringsAsFactors = FALSE
)

# Check: exposure time point matches the folder time point
stopifnot(all(endsWith(files_index$path,
                       paste0("_", files_index$time_point, "_annotated.rds"))))
# Check: files per run and analysis (6 expected each)
table(files_index$run, files_index$analysis)

#==============================================================================#----
# SAVE                                                                         #----
#==============================================================================#----

writeLines(files_index$path, list_file)
write.csv(files_index, index_file, row.names = FALSE)
cat("List written:", list_file, "|", nrow(files_index), "files\n")
cat("Index written:", index_file, "\n")

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----
