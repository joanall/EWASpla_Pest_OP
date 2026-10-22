#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: Functional annotation of all CpGs of each cleaned EWAS model.
# - hg19 annotation: position, relation to CpG island, gene name and group
# - 15 chromatin states in fetal placenta (FP_15_E091)
# - Partially methylated domains in placenta (PMD_placenta)
# - Germline differentially methylated regions / imprinted regions (IR_Placenta)
#
# Argument (optional): pattern for the run folders to include (regular
# expression; default: all)
#   Rscript script/5_Annotation/Annotate_EWAS.R "sepages_2tp_t[23]_old"
#
# [INPUT]:  results/3_Clean_Results/files_cleaned_index.csv (from 3_Get_Cleaned_Files.R)
#           results/5_Annotation/cpg_annotation_all.rds (one row per CpG)
# [OUTPUT]: results/5_Annotation/<run>/<analysis>/<exposure>_<tp>_annotated.rds
#           results/5_Annotation/files_annotated_index.csv

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

index_file    <- "results/3_Clean_Results/files_cleaned_index.csv"
annot_file    <- "results/5_Annotation/cpg_annotation_all.rds"
clean_folder  <- "results/3_Clean_Results/"
output_folder <- "results/5_Annotation/"
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# Cleaned EWAS files (main + sensitivity) from the index
files_index <- read.csv(index_file, stringsAsFactors = FALSE)

# Select runs (argument: pattern, default all)
args <- commandArgs(trailingOnly = TRUE)
run_pattern <- if (length(args) >= 1 && args[1] != "") args[1] else "."
# files_index <- files_index[1, ]   # TEST - remove for the real run
# run_pattern <- "sepages_2tp_t2_old"   # all 36 files of that run (6 analyses × 6 exposures)
files_index <- files_index[grepl(run_pattern, files_index$run), ]
cat("Files to annotate:", nrow(files_index), "\n")
stopifnot(nrow(files_index) > 0)

files_in <- files_index$path
files_in[!file.exists(files_in)]   # should be character(0)
stopifnot(all(file.exists(files_in)))

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

for (f in files_in) {
  
  cat("Annotating:", f, "\n")
  ewas_results <- readRDS(f)
  
  ewas_annot <- merge(ewas_results, cpg_annot, by = "CpG", all.x = TRUE)
  ewas_annot <- ewas_annot[order(ewas_annot$pvalue), ]
  stopifnot(nrow(ewas_annot) == nrow(ewas_results))   # merge kept all CpGs once
  
  cat("  CpGs without annotation:", sum(!ewas_results$CpG %in% cpg_annot$CpG), "\n")
  rownames(ewas_annot) <- ewas_annot$CpG
  
  file_out <- sub(clean_folder, output_folder, f, fixed = TRUE)
  file_out <- sub("\\.rds$", "_annotated.rds", file_out)
  dir.create(dirname(file_out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(ewas_annot, file_out)
}

#==============================================================================#----
# CHECK & INDEX                                                                #----
#==============================================================================#----

# Check: one annotated file per cleaned file
files_out <- sub("\\.rds$", "_annotated.rds", sub(clean_folder, output_folder, files_in, fixed = TRUE))
stopifnot(sum(file.exists(files_out)) == length(files_in))

# Index of annotated files (for steps 6 and 7)
annot_index <- files_index
annot_index$path_annotated <- files_out
write.csv(annot_index, file.path(output_folder, "files_annotated_index.csv"), row.names = FALSE)
cat("Index written:", file.path(output_folder, "files_annotated_index.csv"), "\n")

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----