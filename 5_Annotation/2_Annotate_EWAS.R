#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: Functional annotation of all CpGs of each cleaned EWAS model.
# - hg19 annotation: position, relation to CpG island, gene name and group
# - 15 chromatin states in fetal placenta (FP_15_E091)
# - Partially methylated domains in placenta (PMD_placenta)
# - Germline differentially methylated regions / imprinted regions (IR_Placenta)
#
# Arguments:
#   1) cohort: bisc / sepages / meta / bisc_mix / sepages_mix / meta_mix
#   2) pattern for the run folders to include (regular expression; default: all)
#   Rscript script/5_Annotation/2_Annotate_EWAS.R sepages "sepages_2tp_(t2|t3|avg)_20261021"
#   Rscript script/5_Annotation/2_Annotate_EWAS.R meta "meta_tp_(t2|t3|avg)_20261106"
# Then run 3_Get_Annotated_Files.R to (re)write the index of annotated files.
#
# [INPUT]:  <clean dir>/files_cleaned_index.csv (from 3_Get_Cleaned_Files.R)
#           <clean dir> = results/<cohort>/3_Clean_Results or results/7_Meta/3_Read_GWAMA
#           db/annotation_data/cpg_annotation_all.rds (one row per CpG)
# [OUTPUT]: results/<cohort>/5_Annotation/<run>/<analysis>/<exposure>_<tp>_annotated.rds
#           (meta: results/7_Meta/5_Annotation/...)
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

source("script/00_paths.R")

args <- commandArgs(trailingOnly = TRUE)
cohort      <- if (length(args) >= 1) args[1] else "bisc"
run_pattern <- if (length(args) >= 2 && args[2] != "") args[2] else "."
check_cohort(cohort)

clean_folder  <- paste0(clean_dir(cohort), "/")
index_file    <- file.path(clean_dir(cohort), "files_cleaned_index.csv")
annot_file    <- annotation_file
output_folder <- paste0(results_dir(cohort, "5_Annotation", create = TRUE), "/")

# Cleaned EWAS files (main + sensitivity) from the index
files_index <- read.csv(index_file, stringsAsFactors = FALSE)

# Select runs (argument: pattern, default all)
# files_index <- files_index[1, ]   # TEST - remove for the real run
# run_pattern <- "sepages_2tp_t2_old"   # all 36 files of that run (6 analyses × 6 exposures)
files_index <- files_index[grepl(run_pattern, files_index$run), ]
cat("Files to annotate:", nrow(files_index), "\n")
stopifnot(nrow(files_index) > 0)

files_in <- files_index$path
files_in[!file.exists(files_in)]   # should be character(0)
stopifnot(all(file.exists(files_in)))

# Skip files already annotated (set overwrite_annot <- TRUE to redo all)
overwrite_annot <- FALSE
stopifnot(all(startsWith(files_in, clean_folder)))   # index paths match this cohort's clean dir
files_out  <- sub("\\.rds$", "_annotated.rds", sub(clean_folder, output_folder, files_in, fixed = TRUE))
files_todo <- files_in[overwrite_annot | !file.exists(files_out)]
cat("Already annotated:", sum(file.exists(files_out)), "| To annotate:", length(files_todo), "\n")
table(files_index$analysis[files_in %in% files_todo])   # should be only s6

#==============================================================================#----
# LOAD ANNOTATION                                                              #----
#==============================================================================#----

cpg_annot <- readRDS(annot_file)
dim(cpg_annot)
stopifnot("CpG" %in% colnames(cpg_annot))
stopifnot(!any(duplicated(cpg_annot$CpG)))   # one row per CpG, or merge duplicates rows

# Check: no column name clash with the EWAS results (other than CpG)
ewas_cols <- c("Effect", "SE", "tvalue", "pvalue", "CpG", "pval_bn")
clash <- intersect(setdiff(colnames(cpg_annot), "CpG"), ewas_cols)
if (length(clash) > 0) stop("Annotation columns clash with EWAS results: ",
                            paste(clash, collapse = ", "))

#==============================================================================#----
# MERGE ANNOTATION WITH EACH MODEL                                             #----
#==============================================================================#----
# Output mirrors the cleaned results folder structure, e.g.
#   3_Clean_Results/sepages_2tp_t2_old/s1/24d_t2.rds
#   -> 5_Annotation/sepages_2tp_t2_old/s1/24d_t2_annotated.rds

for (f in files_todo) {  
  file_out <- files_out[match(f, files_in)]   # annotated file for this cleaned file
  cat("Annotating:", f, "->", file_out, "\n")
  ewas_results <- readRDS(f)
  
  ewas_annot <- merge(ewas_results, cpg_annot, by = "CpG", all.x = TRUE)
  ewas_annot <- ewas_annot[order(ewas_annot$pvalue), ]
  stopifnot(nrow(ewas_annot) == nrow(ewas_results))   # merge kept all CpGs once
  
  cat("  CpGs without annotation:", sum(!ewas_results$CpG %in% cpg_annot$CpG), "\n")
  rownames(ewas_annot) <- ewas_annot$CpG
  
  dir.create(dirname(file_out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(ewas_annot, file_out)
}

#==============================================================================#----
# CHECK                                                                        #----
#==============================================================================#----

# Check: one annotated file per cleaned file
stopifnot(sum(file.exists(files_out)) == length(files_in))
cat("Annotated files in", output_folder, ":", length(files_out), "\n")
cat("Next: Rscript script/5_Annotation/3_Get_Annotated_Files.R", cohort, "to update the index\n")

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----