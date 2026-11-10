#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: Collect the suggestive CpGs (p < 1e-5) of every annotated EWAS
# result (main + sensitivity) into one long table: one row per CpG per model,
# with the model described by the index columns.
# Arguments (optional): 1) pattern for the runs (default: all)
#                       2) p-value cutoff (default: 1e-5)
#   Rscript script/6_Analyze_Results/5_Get_Suggestive_CpGs_Anno.R "sepages_2tp_t[23]_20261021" 1e-5
#
# [INPUT]:  results/5_Annotation/files_annotated_index.csv
#           -> results/5_Annotation/<run>/<analysis>/<exposure>_<tp>_annotated.rds
# [OUTPUT]: results/6_Analyze_Results/summary_tables/all_suggestive_<pattern>.csv / .rds
#           results/6_Analyze_Results/summary_tables/n_suggestive_<pattern>.csv
#==============================================================================#-

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

index_file <- "results/5_Annotation/files_annotated_index.csv"
outdir     <- "results/6_Analyze_Results/summary_tables"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
run_pattern <- if (length(args) >= 1 && args[1] != "") args[1] else "."
cutoff_sugg <- if (length(args) >= 2) as.numeric(args[2]) else 1e-5
# run_pattern <- "sepages_2tp_t[23]_20261021"   # TEST

files_index <- read.csv(index_file, stringsAsFactors = FALSE)
files_index <- files_index[grepl(run_pattern, files_index$run), ]
cat("Files:", nrow(files_index), "| cutoff:", cutoff_sugg, "\n")
stopifnot(nrow(files_index) > 0, all(file.exists(files_index$path)))

# Model information (from the index) and CpG information (from the results)
index_cols <- c("run", "cohort", "sample", "time_point", "run_date", "analysis", "exposure")
# Required: the script stops if any of these is missing
cpg_cols_req <- c("CpG", "chr", "pos", "gene_name", "gene_group",
                  "Effect", "SE", "pvalue", "pval_bn")

# Annotation: kept when present, filled with NA when a run predates them
cpg_cols_ann <- c("Relation_to_Island", "States15_FP", "PMD", "gDMR")

cpg_cols <- c(cpg_cols_req, cpg_cols_ann)

# Name for the output files ("all" when no pattern)
tag <- if (run_pattern == ".") "all" else gsub("[^A-Za-z0-9_]", "", run_pattern)

#==============================================================================#----
# COLLECT SUGGESTIVE CpGs                                                      #----
#==============================================================================#----

all_sugg  <- NULL
n_per_mod <- NULL

for (i in seq_len(nrow(files_index))) {
  
  f   <- files_index[i, ]
  res <- readRDS(f$path)
  
  # Check: required columns present; annotation columns added as NA if absent
  missing_req <- setdiff(cpg_cols_req, colnames(res))
  if (length(missing_req) > 0) stop("Missing in ", f$path, ": ", paste(missing_req, collapse = ", "))
  
  missing_ann <- setdiff(cpg_cols_ann, colnames(res))
  if (length(missing_ann) > 0) {
    warning("No annotation columns in ", f$path, ": ", paste(missing_ann, collapse = ", "))
    res[, missing_ann] <- NA_character_
  }
  
  sugg <- res[!is.na(res$pvalue) & res$pvalue < cutoff_sugg, cpg_cols]
  
  # Binary flags: CpG inside a partially methylated domain / imprinted region
  sugg$in_pmd  <- as.integer(!is.na(sugg$PMD)  & sugg$PMD  != "")
  sugg$in_gdmr <- as.integer(!is.na(sugg$gDMR) & sugg$gDMR != "")
  
  n_per_mod <- rbind(n_per_mod, data.frame(f[, index_cols],
                                           n_suggestive = nrow(sugg),
                                           n_bonferroni = sum(sugg$pval_bn < 0.05, na.rm = TRUE),
                                           n_pmd        = sum(sugg$in_pmd),
                                           n_gdmr       = sum(sugg$in_gdmr)))
  cat("[", i, "/", nrow(files_index), "]", f$run, f$analysis, f$exposure, "|", nrow(sugg), "\n")
  
  if (nrow(sugg) == 0) next
  all_sugg <- rbind(all_sugg, data.frame(f[rep(1, nrow(sugg)), index_cols], sugg, row.names = NULL))
  
  # Check: annotation coverage among the suggestive CpGs
  round(100 * colMeans(!is.na(all_sugg[, cpg_cols_ann])), 1)
  table(all_sugg$in_pmd, all_sugg$in_gdmr)
  
}

#==============================================================================#----
# SAVE                                                                         #----
#==============================================================================#----

all_sugg <- all_sugg[order(all_sugg$run, all_sugg$analysis, all_sugg$exposure, all_sugg$pvalue), ]

# Check: rows = sum of suggestive per model
nrow(all_sugg) == sum(n_per_mod$n_suggestive)   # TRUE
head(all_sugg)

write.csv(all_sugg, file.path(outdir, paste0("all_suggestive_", tag, ".csv")), row.names = FALSE)
saveRDS(all_sugg,   file.path(outdir, paste0("all_suggestive_", tag, ".rds")))
write.csv(n_per_mod, file.path(outdir, paste0("n_suggestive_", tag, ".csv")), row.names = FALSE)

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----