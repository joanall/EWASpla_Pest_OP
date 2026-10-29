#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: OBJECTIVE 1 - Association between prenatal pesticide exposure
# and placental DNA methylation. Summary tables of suggestive CpGs (p < 1e-5)
# from the main models, with their look-up in the sensitivity analyses.
# - Table 1-2: suggestive CpGs at t2 / t3 + effect and p in s1-s5 (look-up),
#              to check whether each association remains after the
#              sensitivity analyses
# Consistency across time points and exposures (Objective 2) is in
#   3b_Get_Consistency_Summary.R, which reads Tables 1-2 produced here.
# Argument: pattern for the runs (one cohort, sample and run date), e.g.
#   Rscript script/6_Analyze_Results/3a_Get_Suggestive_Summary.R "sepages_2tp_t[23]_20261021"
#
# [INPUT]:  results/5_Annotation/files_annotated_index.csv
#           -> results/5_Annotation/<run>/<analysis>/<exposure>_<tp>_annotated.rds
# [OUTPUT]: results/6_Analyze_Results/summary_tables/<cohort>_<sample>_<date>/
#           table1_suggestive_t2.xlsx/.rds, table2_suggestive_t3.xlsx/.rds
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_functions.R")
library(writexl)   # .xlsx output: no separator/decimal issues when opening in Excel

index_file   <- "results/5_Annotation/files_annotated_index.csv"
cutoff_sugg  <- 1e-5    # suggestive threshold
cutoff_nom   <- 0.05    # nominal threshold for look-up replication

args <- commandArgs(trailingOnly = TRUE)
run_pattern <- if (length(args) >= 1 && args[1] != "") args[1] else "."
# run_pattern <- "sepages_2tp_t[23]_20261021"   # TEST

files_index <- read.csv(index_file, stringsAsFactors = FALSE)
files_index <- files_index[grepl(run_pattern, files_index$run), ]
stopifnot(all(file.exists(files_index$path)))

# Check: one cohort, sample and run date (comparisons must not mix runs)
stopifnot(length(unique(files_index$cohort)) == 1,
          length(unique(files_index$sample)) == 1,
          length(unique(files_index$run_date)) == 1)
table(files_index$time_point, files_index$analysis)   # 6 per cell

main_index <- files_index[files_index$analysis == "main", ]
sens_index <- files_index[files_index$analysis != "main", ]

tag    <- paste(main_index$cohort[1], main_index$sample[1], main_index$run_date[1], sep = "_")
outdir <- file.path("results/6_Analyze_Results/summary_tables", tag)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

exposures <- unique(main_index$exposure)
keep_cols <- c("CpG", "chr", "pos", "gene_name", "Relation_to_Island", "States15_FP",
               "Effect", "SE", "pvalue", "pval_bn")

#==============================================================================#----
# 1. SUGGESTIVE CpGs - MAIN MODELS                                             #----
#==============================================================================#----

top_main <- NULL

for (i in seq_len(nrow(main_index))) {
  f   <- main_index[i, ]
  key <- paste(f$exposure, f$time_point, sep = "_")
  res <- readRDS(f$path)
  
  top <- res[!is.na(res$pvalue) & res$pvalue < cutoff_sugg, intersect(keep_cols, colnames(res))]
  cat(key, "| suggestive:", nrow(top), "\n")
  if (nrow(top) == 0) next
  top_main <- rbind(top_main, data.frame(exposure = f$exposure, time_point = f$time_point, top))
}
nrow(top_main)

#==============================================================================#----
# 2. LOOK-UP IN SENSITIVITY ANALYSES                                           #----
#==============================================================================#----

sens_names <- sort(unique(sens_index$analysis))

for (s in sens_names) {
  top_main[[paste0("Effect_", s)]] <- NA
  top_main[[paste0("SE_", s)]]     <- NA
  top_main[[paste0("pvalue_", s)]] <- NA
}

for (i in seq_len(nrow(sens_index))) {
  f   <- sens_index[i, ]
  idx <- which(top_main$exposure == f$exposure & top_main$time_point == f$time_point)
  if (length(idx) == 0) next
  res <- readRDS(f$path)
  m   <- match(top_main$CpG[idx], res$CpG)
  top_main[idx, paste0("Effect_", f$analysis)] <- res$Effect[m]
  top_main[idx, paste0("SE_", f$analysis)]     <- res$SE[m]
  top_main[idx, paste0("pvalue_", f$analysis)] <- res$pvalue[m]
}

# Agreement with the main model: same direction, and p < 0.05
eff_sens <- as.matrix(top_main[, paste0("Effect_", sens_names)])
p_sens   <- as.matrix(top_main[, paste0("pvalue_", sens_names)])
top_main$n_sens_same_dir <- rowSums(sign(eff_sens) == sign(top_main$Effect), na.rm = TRUE)
top_main$n_sens_p005     <- rowSums(p_sens < cutoff_nom & sign(eff_sens) == sign(top_main$Effect), na.rm = TRUE)

#==============================================================================#----
# TABLES 1-2: SUGGESTIVE CpGs PER TIME POINT                                   #----
#==============================================================================#----

for (tp in c("t2", "t3")) {
  tab <- top_main[top_main$time_point == tp, ]
  tab <- tab[order(tab$exposure, tab$pvalue), ]
  n_tab <- ifelse(tp == "t2", 1, 2)
  file_base <- file.path(outdir, paste0("table", n_tab, "_suggestive_", tp))
  write_xlsx(tab, paste0(file_base, ".xlsx"))
  saveRDS(tab, paste0(file_base, ".rds"))
  cat("Table", n_tab, "-", tp, ":", nrow(tab), "suggestive CpGs\n")
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----