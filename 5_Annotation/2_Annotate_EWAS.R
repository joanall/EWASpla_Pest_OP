#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: Functional annotation of all CpGs of each EWAS model:
# - hg19 annotation: position, relation to CpG island, gene name and group
# - 15 chromatin states in fetal placenta (FP_15_E091)
# - Partially methylated domains in placenta (PMD_placenta)
# - Germline differentially methylated regions / imprinted regions (IR_Placenta)
# W ahve the annnoation in a single df and we just need to ,etge it to eqch results 
# [INPUT]: results/5_Annotation/cpg_annotation_all.rds
# [OUTPUT]: 









index_file    <- "results/3_Clean_Results/files_cleaned_index.csv"
clean_folder  <- "results/3_Clean_Results/"
output_folder <- "results/5_Annotation/"

# Cleaned EWAS files (main + sensitivity) from the index
files_index <- read.csv(index_file)
head(files_index)                  # check: column with the file paths
index_col <- "path"                

files_in <- files_index[[index_col]]
length(files_in)
files_in[!file.exists(files_in)]   # should be character(0)


==============================================================================#----
# MERGE ANNOTATION WITH EACH MODEL                                             #----
#==============================================================================#----
# Output mirrors the cleaned results folder structure, e.g.
# 3_Clean_Results/sens_t2/s1/24d_t2.rds -> 5_Annotation/sens_t2/s1/24d_t2_annotated.rds

for (f in files_in) {
  
  cat("Annotating:", f, "\n")
  ewas_results <- readRDS(f)
  
  ewas_annot <- merge(ewas_results, cpg_annot, by = "CpG", all.x = TRUE)
  ewas_annot <- ewas_annot[order(ewas_annot$pvalue), ]
  stopifnot(nrow(ewas_annot) == nrow(ewas_results))   # merge kept all CpGs once
  
  file_out <- sub(clean_folder, output_folder, f, fixed = TRUE)
  file_out <- sub("\\.rds$", "_annotated.rds", file_out)
  dir.create(dirname(file_out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(ewas_annot, file_out)
}

# Check: one annotated file per cleaned file
files_out <- sub("\\.rds$", "_annotated.rds", sub(clean_folder, output_folder, files_in, fixed = TRUE))
sum(file.exists(files_out)) == length(files_in)   # TRUE