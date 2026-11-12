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
# - Table 1-2: suggestive CpGs at t2 / t3 + effect and p in s1-s6 (look-up),
#              to check whether each association remains after the
#              sensitivity analyses
# - Table 7:   same for the pregnancy-average model (only if the avg run is
#              included in the pattern; tables 3-6 are taken by 3b / 3c)
# Consistency across time points and exposures (Objective 2) is in
#   3b_Get_Consistency_Summary.R, which reads Tables 1-2 produced here.
# Arguments: 1) cohort (bisc / sepages / meta / bisc_mix / sepages_mix / meta_mix), 2) pattern for the runs (one cohort, sample and run date), e.g.
#   Rscript script/6_Analyze_DMPs/3a_Get_Suggestive_Summary.R sepages "sepages_2tp_t[23]_20261021"
#   Rscript script/6_Analyze_DMPs/3a_Get_Suggestive_Summary.R bisc "bisc_2tp_(t2|t3|avg)_20261022"
#
# (<cohort> = bisc / sepages; meta = results/7_Meta/)
# [INPUT]:  results/<cohort>/5_Annotation/files_annotated_index.csv
#           -> results/<cohort>/5_Annotation/<run>/<analysis>/<exposure>_<tp>_annotated.rds
# [OUTPUT]: results/<cohort>/6_Analyze_DMPs/summary_tables/<cohort>_<sample>_<date>/
#           table1_suggestive_t2.xlsx/.rds, table2_suggestive_t3.xlsx/.rds,
#           table7_suggestive_avg.xlsx/.rds
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_paths.R")

args <- commandArgs(trailingOnly = TRUE)
cohort      <- if (length(args) >= 1) args[1] else "bisc"
run_pattern <- if (length(args) >= 2 && args[2] != "") args[2] else "."
check_cohort(cohort)

source("script/00_functions.R")

index_file   <- results_dir(cohort, "5_Annotation", "files_annotated_index.csv")
cutoff_sugg  <- 1e-5    # suggestive threshold
cutoff_nom   <- 0.05    # nominal threshold for look-up replication

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
outdir <- file.path(results_dir(cohort, "6_Analyze_DMPs", "summary_tables"), tag)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

exposures <- unique(main_index$exposure)
keep_cols <- c("CpG", "chr", "pos", "gene_name", "UCSC_RefGene_Group","States15_FP",
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
# TABLES 1-2 AND 7: SUGGESTIVE CpGs PER TIME POINT                             #----
#==============================================================================#----

# Table number per time point (3-6 are used by 3b / 3c)
table_number <- c(t2 = 1, t3 = 2, avg = 7)
tps_run <- intersect(names(table_number), unique(main_index$time_point))
tps_run   # time points in the selected runs

for (tp in tps_run) {
  tab <- top_main[top_main$time_point == tp, ]
  tab <- tab[order(tab$exposure, tab$pvalue), ]
  n_tab <- table_number[[tp]]
  file_base <- file.path(outdir, paste0("table", n_tab, "_suggestive_", tp))
  save_table(tab, file_base, png = FALSE)
  saveRDS(tab, paste0(file_base, ".rds"))
  cat("Table", n_tab, "-", tp, ":", nrow(tab), "suggestive CpGs\n")
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----