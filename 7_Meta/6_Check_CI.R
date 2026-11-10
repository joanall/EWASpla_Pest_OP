################################################################################
# Check cohort contributions and large SEPAGES CIs in the meta-analysis
# All exposures x all time points (T2, T3, pregnancy average) - main analysis
#
# Aim: the forest plots (BiSC, SEPAGES, meta) show much wider CIs in SEPAGES
#      at some CpGs. This script
#        (1) quantifies each cohort's contribution to the pooled fixed-effects
#            estimate (inverse-variance weights), and
#        (2) flags top hits where the SEPAGES SE is unusually inflated
#            relative to the typical SEPAGES/BiSC SE ratio.
#      Flagged CpGs are saved to txt, to be checked on the PC holding the
#      SEPAGES methylation data (section 4: CpG-specific diagnostics).
################################################################################

## ---- 0. Settings ----
exposures <- c("24d", "3pba", "tcpy", "dmp", "dep")   # check: file prefixes as in the .rds names
tps       <- c("t2", "t3", "avg")

# Folder per time point (adjust dated folders for T2 / T3)
dir_meta    <- c(t2  = "results/7_Meta/3_Read_GWAMA/meta_tp_t2_20261106/main/",
                 t3  = "results/7_Meta/3_Read_GWAMA/meta_tp_t3_20261106/main/",
                 avg = "results/7_Meta/3_Read_GWAMA/meta_tp_avg_20261106/main/")
dir_bisc    <- c(t2  = "results/bisc/3_Clean_Results/bisc_tp_t2_20261022/main/",
                 t3  = "results/bisc/3_Clean_Results/bisc_tp_t3_20261022/main/",
                 avg = "results/bisc/3_Clean_Results/bisc_tp_avg_20261022/main/")
dir_sepages <- c(t2  = "results/sepages/3_Clean_Results/sepages_2tp_t2_20261021/main/",
                 t3  = "results/sepages/3_Clean_Results/sepages_2tp_t3_20261021/main/",
                 avg = "results/sepages/3_Clean_Results/sepages_2tp_avg_20261021/main/")

# Check: file names in each folder
list.files(dir_meta["avg"])

p_top          <- 1e-5   # top-hit threshold (meta p-value)
se_ratio_cut   <- 1.5    # flag if SEPAGES/BiSC SE ratio > 1.5x its median
dir_out        <- "results/7_Meta/4_Check_CI/"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

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
  meta_top  <- meta_full[meta_full$pvalue < p_top, ]
  meta_flag <- meta_top[meta_top$se_ratio_rel > se_ratio_cut, ]
  
  cat(exp_i, tp_i, "| top hits:", nrow(meta_top),
      "| flagged:", nrow(meta_flag),
      "| median SE ratio:", round(median(meta_full$se_ratio), 2), "\n")
  
  if (nrow(meta_flag) == 0) return(NULL)
  data.frame(exposure = exp_i, tp = tp_i,
             meta_flag[, c("CpG", "pvalue", "Effect", "SE",
                           "Effect_bisc", "SE_bisc", "Effect_sepages", "SE_sepages",
                           "pct_w_sepages", "se_ratio", "se_ratio_rel", "i2")])
}

## ---- Run for all exposures x time points ----
flag_list <- list()
for (tp_i in tps) {
  for (exp_i in exposures) {
    flag_list[[paste(exp_i, tp_i, sep = "_")]] <- check_ci(exp_i, tp_i)
  }
}

flagged_all <- do.call(rbind, flag_list)
rownames(flagged_all) <- NULL

# Check: counts per exposure x time point
table(flagged_all$exposure, flagged_all$tp)

## ---- Save ----
write.table(flagged_all,
            paste0(dir_out, "flagged_cpgs_large_ci_detail.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cpg_flagged_unique <- unique(flagged_all$CpG)
length(cpg_flagged_unique)
writeLines(cpg_flagged_unique,
           paste0(dir_out, "flagged_cpgs_large_ci_unique.txt"))