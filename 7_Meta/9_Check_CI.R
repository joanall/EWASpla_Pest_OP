#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the
# BiSC and SEPAGES cohorts.
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: Check cohort contributions and large SEPAGES CIs in the meta-analysis
# All exposures x all time points (T2, T3, pregnancy average) - main analysis
# Note: the forest plots (BiSC, SEPAGES, meta) show much wider CIs in SEPAGES
# at some CpGs. This script
#  (1) quantifies each cohort's contribution to the pooled fixed-effects
#      estimate (inverse-variance weights), and
#  (2) flags top hits where the SEPAGES SE is unusually inflated
#      relative to the typical SEPAGES/BiSC SE ratio.
# All top hits (with a flag for large SEPAGES CIs) are saved to txt, to be
# checked on the PC holding the SEPAGES methylation data
# (9_1_Check_CI_Sepages.R: CpG-specific diagnostics; plots for all top CpGs. 
# Copy results/7_Meta/9_Check_CI/tophits_ci_detail.txt from BiSC PC to
# the same path on the SEPAGES PC.
#
# [INPUT]: 
# [OUTPUT]: results/7_Meta/9_Check_CI/tophits_ci_detail.xlsx (+ .txt, read by
#         9_1_Check_CI_Sepages.R), tophits_ci_unique.txt
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_paths.R")
source("script/00_template_plot.R")   # save_table()

exposures <- c("24d", "3pba", "tcpy", "dmp", "dep")   # check: file prefixes as in the .rds names
tps       <- c("t2", "t3", "avg")

# Folder per time point (adjust dated folders for T2 / T3)
# Run folders per time point: <cohort>_<sample>_<tp>_<date>
run_meta    <- paste0("meta_tp_",     tps, "_20261106")
run_bisc    <- paste0("bisc_tp_",     tps, "_20261022")
run_sepages <- paste0("sepages_2tp_", tps, "_20261021")
dir_meta    <- setNames(paste0(file.path(clean_dir("meta"),    run_meta,    "main"), "/"), tps)
dir_bisc    <- setNames(paste0(file.path(clean_dir("bisc"),    run_bisc,    "main"), "/"), tps)
dir_sepages <- setNames(paste0(file.path(clean_dir("sepages"), run_sepages, "main"), "/"), tps)
dir_meta; dir_bisc; dir_sepages
stopifnot(dir.exists(dir_meta), dir.exists(dir_bisc), dir.exists(dir_sepages))

# Check: file names in each folder
list.files(dir_meta["avg"])

p_top          <- 1e-5   # top-hit threshold (meta p-value)
se_ratio_cut   <- 1.5    # flag if SEPAGES/BiSC SE ratio > 1.5x its median
dir_out        <- paste0(results_dir("meta", "9_Check_CI", create = TRUE), "/")


#==============================================================================#----
# DEFINE FUNCTION TO CHECK CI                                                  #----
#==============================================================================#----
## ---- 1-3. Function: load, compute weights and SE ratios, flag top hits ----
check_ci <- function(exp_i, tp_i) {
  file_i <- paste0(exp_i, "_", tp_i, ".rds")
  
  # 1. Load meta-analysis + per-cohort results
  meta        <- readRDS(paste0(dir_meta[tp_i],    file_i))
  res_bisc    <- readRDS(paste0(dir_bisc[tp_i],    file_i))
  res_sepages <- readRDS(paste0(dir_sepages[tp_i], file_i))
  
  res_cohorts <- merge(res_bisc[, c("CpG", "Effect", "SE")],
                       res_sepages[, c("CpG", "Effect", "SE")],
                       by = "CpG", suffixes = c("_bisc", "_sepages"))
  meta_full   <- merge(meta[, c("CpG", "Effect", "SE", "pvalue", "i2")],
                       res_cohorts, by = "CpG")
  
  # Check: merge did not lose CpGs
  if (nrow(meta_full) != nrow(meta))
    warning(exp_i, " ", tp_i, ": ", nrow(meta) - nrow(meta_full), " CpGs lost in merge")
  
  # 2. Inverse-variance weights and SE ratio
  w_bisc    <- 1 / meta_full$SE_bisc^2
  w_sepages <- 1 / meta_full$SE_sepages^2
  meta_full$pct_w_bisc    <- round(100 * w_bisc    / (w_bisc + w_sepages), 1)
  meta_full$pct_w_sepages <- round(100 * w_sepages / (w_bisc + w_sepages), 1)
  meta_full$se_ratio      <- meta_full$SE_sepages / meta_full$SE_bisc
  meta_full$se_ratio_rel  <- meta_full$se_ratio / median(meta_full$se_ratio)
  
  # Check: pooled SE matches fixed-effects SE (should be ~0)
  se_diff_max <- max(abs(meta_full$SE - sqrt(1 / (w_bisc + w_sepages))))
  if (se_diff_max > 1e-4)
    warning(exp_i, " ", tp_i, ": pooled SE differs from FE SE (max diff ", se_diff_max, ")")
  
  # 3. Top hits and flagged CpGs
  meta_top <- meta_full[meta_full$pvalue < p_top, ]
  meta_top$flag_large_ci <- meta_top$se_ratio_rel > se_ratio_cut
  
  cat(exp_i, tp_i, "| top hits:", nrow(meta_top),
      "| flagged:", sum(meta_top$flag_large_ci),
      "| median SE ratio:", round(median(meta_full$se_ratio), 2), "\n")
  
  if (nrow(meta_top) == 0) return(NULL)
  data.frame(exposure = exp_i, tp = tp_i,
             meta_top[, c("CpG", "pvalue", "Effect", "SE",
                          "Effect_bisc", "SE_bisc", "Effect_sepages", "SE_sepages",
                          "pct_w_sepages", "se_ratio", "se_ratio_rel", "i2", "flag_large_ci")])
}

#==============================================================================#----
# GET CpGS CI THQT NEEDS CHECK                                                 #----
#==============================================================================#----
tophits_list <- list()
for (tp_i in tps) {
  for (exp_i in exposures) {
    tophits_list[[paste(exp_i, tp_i, sep = "_")]] <- check_ci(exp_i, tp_i)
  }
}

tophits_all <- do.call(rbind, tophits_list)
rownames(tophits_all) <- NULL
# Check: top hits and flagged CpGs per exposure x time point
table(tophits_all$exposure, tophits_all$tp)
table(tophits_all$exposure[tophits_all$flag_large_ci], tophits_all$tp[tophits_all$flag_large_ci])

# Save 
save_table(tophits_all, paste0(dir_out, "tophits_ci_detail"), png = FALSE)   # to read
# Plain-text copy: input of 9_1_Check_CI_Sepages.R on the SEPAGES PC
write.table(tophits_all,
            paste0(dir_out, "tophits_ci_detail.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cpg_tophits_unique <- unique(tophits_all$CpG)
length(cpg_tophits_unique)
writeLines(cpg_tophits_unique, paste0(dir_out, "tophits_ci_unique.txt"))


#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----
