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
#   Rscript script/6_Analyze_Results/5_Get_All_Suggestive.R "sepages_2tp_t[23]_20261021" 1e-5
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
cpg_cols   <- c("CpG", "chr", "pos", "gene_name", "gene_group",
                "Effect", "SE", "pvalue", "pval_bn")

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
  
  # Check: all needed columns present
  missing_cols <- setdiff(cpg_cols, colnames(res))
  if (length(missing_cols) > 0) stop("Missing in ", f$path, ": ", paste(missing_cols, collapse = ", "))
  
  sugg <- res[!is.na(res$pvalue) & res$pvalue < cutoff_sugg, cpg_cols]
  n_per_mod <- rbind(n_per_mod, data.frame(f[, index_cols], n_suggestive = nrow(sugg),
                                           n_bonferroni = sum(sugg$pval_bn < 0.05, na.rm = TRUE)))
  cat("[", i, "/", nrow(files_index), "]", f$run, f$analysis, f$exposure, "|", nrow(sugg), "\n")
  
  if (nrow(sugg) == 0) next
  all_sugg <- rbind(all_sugg, data.frame(f[rep(1, nrow(sugg)), index_cols], sugg, row.names = NULL))
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