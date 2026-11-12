#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date:  October - November 2026
#==============================================================================#
# [DESCRIPTION]:  Runs GWAMA (inverse-variance fixed effects, quantitative trait) on every
# filelist (.in) created by 1_Get_Input_GWAMA.R

# [INPUT]: results/7_Meta/1_Input_GWAMA/<meta_id>/gwama_filelist_<exposure>.in
# [OUTPUT]: results/7_Meta/2_Run_GWAMA/<meta_id>/meta_<meta_id>_<exposure>.out (+ .log)
#
# [HOW TO INSTALL GWAMA] (once, in bash):
#   conda activate mqtl_env
#   conda install -c bioconda -c conda-forge gwama
#   which GWAMA
# Add it to .Renviron:
#   PATH="/home/isglobal.lan/jllaurado/bin:/home/isglobal.lan/jllaurado/.conda/envs/mqtl_env/bin:${PATH}"
# Then use - Sys.which("GWAMA")
#
# NOTE: Run from the project root (the .in files use relative paths).
#       No genomic control is applied (GWAMA default; -gc not set).
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
path_to_gwama <- "/home/isglobal.lan/jllaurado/.conda/envs/mqtl_env/bin/GWAMA"
stopifnot(file.exists(path_to_gwama))     # GWAMA found in the environment

source("script/00_paths.R")

# Argument: meta (single-exposure EWAS, default) or meta_mix (mixture models)
args       <- commandArgs(trailingOnly = TRUE)
meta_type  <- if (length(args) >= 1) args[1] else "meta"
check_cohort(meta_type, c("meta", "meta_mix"))

dir_input <- results_dir(meta_type, "1_Input_GWAMA")
dir_run   <- results_dir(meta_type, "2_Run_GWAMA")

filelists <- list.files(dir_input, pattern = "\\.in$", recursive = TRUE, full.names = TRUE)

# Analyses to run: "main" only for now; add e.g. "s1", "s2" ... later
which_analyses <- c("main")
analysis_of_file <- sub("^(t2|t3|avg)_", "", basename(dirname(filelists)))   # t2_main -> main
filelists <- filelists[analysis_of_file %in% which_analyses]

# Check: should be t2_main, t3_main, avg_main with 6 exposures each
print(table(basename(dirname(filelists))))
stopifnot(length(filelists) > 0)
cat(length(filelists), "GWAMA runs to do\n")

#==============================================================================#----
# RUN GWAMA                                                                    #----
#==============================================================================#----
runs_failed <- character(0)

for (filelist in filelists) {
  meta_id  <- basename(dirname(filelist))
  exposure <- sub("^gwama_filelist_(.*)\\.in$", "\\1", basename(filelist))
  
  # Inline check: both study files exist
  study_files <- readLines(filelist)
  stopifnot(length(study_files) == 2, file.exists(study_files))
  
  dir_out <- file.path(dir_run, meta_id)
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  prefix_out <- file.path(dir_out, paste0("meta_", meta_id, "_", exposure))
  
  cat("Running", meta_id, "-", exposure, "... ")
  exit_status <- system2(path_to_gwama,
                         args   = c("--filelist", filelist,
                                    "--output",   prefix_out,
                                    "-qt"),
                         stdout = paste0(prefix_out, ".stdout"),
                         stderr = paste0(prefix_out, ".stderr"))
  
  if (exit_status != 0 || !file.exists(paste0(prefix_out, ".out"))) {
    cat("FAILED (see", paste0(prefix_out, ".stderr"), ")\n")
    runs_failed <- c(runs_failed, prefix_out)
  } else {
    cat("done\n")
  }
}

# Final check
if (length(runs_failed) > 0) {
  warning(length(runs_failed), " GWAMA runs failed:\n", paste(runs_failed, collapse = "\n"))
} else {
  cat("\nAll", length(filelists), "GWAMA runs finished\n")
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----