################################################################################
# Section 4 - CpG-specific diagnostics for large SEPAGES CIs
# Run on the PC holding the SEPAGES methylation data
#
# Input : flagged_cpgs_large_ci_detail.txt (from sections 1-3)
# Aim   : for each flagged exposure x time point x CpG, refit the EWAS model
#         in SEPAGES and identify why the SE is large:
#           - fewer samples (probe filtering)
#           - low / high methylation variability (vs all CpGs)
#           - outliers or many down-weighted samples in the robust regression
#           - non-convergence
# Output: table of diagnostics + PDF with per-CpG plots
################################################################################

library(MASS)

## ---- 4.0 Settings and data ----
dir_out <- "results/sepages/4_Check_CI/"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

flagged_all <- read.table("flagged_cpgs_large_ci_detail.txt",
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cpg_large_ci <- unique(flagged_all$CpG)

# Methylation: CpG x sample, same scale as EWAS (e.g. M-values)  <-- adjust
meth_sepages  <- readRDS("path/to/sepages_methylation.rds")
# Phenotype: one row per sample, sample_id matching colnames(meth_sepages)  <-- adjust
pheno_sepages <- readRDS("path/to/sepages_pheno_final.rds")

# Check: CpGs and samples present
sum(cpg_large_ci %in% rownames(meth_sepages)); length(cpg_large_ci)
sum(pheno_sepages$sample_id %in% colnames(meth_sepages)); nrow(pheno_sepages)

# Exposure variable per exposure x time point  <-- adjust to real names
exp_vars <- c(
  "24d_t2"  = "pest_24d_det_t2",       "24d_t3"  = "pest_24d_det_t3",       "24d_avg"  = "pest_24d_det_avg",
  "3pba_t2" = "pest_3pba_ilog2_t2",    "3pba_t3" = "pest_3pba_ilog2_t3",    "3pba_avg" = "pest_3pba_ilog2_avg",
  "tcpy_t2" = "pest_tcpy_ilog2_t2",    "tcpy_t3" = "pest_tcpy_ilog2_t3",    "tcpy_avg" = "pest_tcpy_ilog2_avg",
  "dmp_t2"  = "pest_dmp_ms_ilog2_t2",  "dmp_t3"  = "pest_dmp_ms_ilog2_t3",  "dmp_avg"  = "pest_dmp_ms_ilog2_avg",
  "dep_t2"  = "pest_dep_ms_ilog2_t2",  "dep_t3"  = "pest_dep_ms_ilog2_t3",  "dep_avg"  = "pest_dep_ms_ilog2_avg"
)

# Covariates (main model, as in Run_EWAS)  <-- adjust
covars <- c("covar1", "covar2", "covar3")

# Check: all variables exist
setdiff(c(exp_vars, covars), colnames(pheno_sepages))

## ---- 4.1 Genome-wide reference: methylation SD across all CpGs ----
sd_meth_all <- apply(meth_sepages, 1, sd, na.rm = TRUE)
summary(sd_meth_all)

## ---- 4.2 Function: refit one CpG for one exposure x time point ----
diag_cpg <- function(cpg, exp_var, covars, meth, pheno, sd_meth_all) {
  pheno_cc <- pheno[complete.cases(pheno[, c(exp_var, covars)]), ]
  df_i     <- data.frame(meth = meth[cpg, pheno_cc$sample_id], pheno_cc)
  df_i     <- df_i[!is.na(df_i$meth), ]
  
  form_i <- as.formula(paste("meth ~", exp_var, "+", paste(covars, collapse = " + ")))
  fit_i  <- rlm(form_i, data = df_i, maxit = 200)   # <-- same settings as Run_EWAS
  
  # Covariate-adjusted exposure SD
  form_x  <- as.formula(paste(exp_var, "~", paste(covars, collapse = " + ")))
  sd_x_adj <- sd(resid(lm(form_x, data = df_i)))
  
  # Outliers (> 4 MAD from median)
  meth_mad <- mad(df_i$meth)
  exp_mad  <- mad(df_i[[exp_var]])
  is_meth_out <- abs(df_i$meth - median(df_i$meth)) > 4 * meth_mad
  is_exp_out  <- if (exp_mad > 0) abs(df_i[[exp_var]] - median(df_i[[exp_var]])) > 4 * exp_mad else rep(FALSE, nrow(df_i))
  is_low_w    <- fit_i$w < 0.5
  
  list(
    stats = c(
      n_cc_pheno       = nrow(pheno_cc),
      n_used           = nrow(df_i),
      sd_meth          = sd(df_i$meth),
      sd_meth_pctl     = round(100 * mean(sd_meth_all <= sd(df_i$meth)), 1),
      s_resid          = fit_i$s,
      sd_x_adj         = sd_x_adj,
      n_low_w          = sum(is_low_w),
      n_meth_out       = sum(is_meth_out),
      n_exp_out        = sum(is_exp_out),
      n_low_w_exp_out  = sum(is_low_w & is_exp_out),
      converged        = as.numeric(fit_i$converged),
      beta_refit       = coef(fit_i)[exp_var],
      se_refit         = summary(fit_i)$coefficients[exp_var, "Std. Error"]
    ),
    df  = df_i,
    w   = fit_i$w
  )
}

## ---- 4.3 Run for all flagged rows + plots ----
diag_list <- list()
pdf(paste0(dir_out, "flagged_cpgs_large_ci_plots.pdf"), width = 10, height = 4.5)

for (i in seq_len(nrow(flagged_all))) {
  cpg_i   <- flagged_all$CpG[i]
  key_i   <- paste(flagged_all$exposure[i], flagged_all$tp[i], sep = "_")
  exp_var <- exp_vars[key_i]
  
  if (!cpg_i %in% rownames(meth_sepages)) {
    warning(cpg_i, " not in SEPAGES methylation matrix"); next
  }
  
  res_i <- diag_cpg(cpg_i, exp_var, covars, meth_sepages, pheno_sepages, sd_meth_all)
  diag_list[[i]] <- data.frame(flagged_all[i, c("exposure", "tp", "CpG", "pvalue",
                                                "SE_bisc", "SE_sepages", "se_ratio_rel",
                                                "pct_w_sepages")],
                               t(res_i$stats), row.names = NULL)
  
  # Plots: methylation distribution + methylation vs exposure (red = weight < 0.5)
  par(mfrow = c(1, 2))
  hist(res_i$df$meth, breaks = 40, main = paste(cpg_i, "-", key_i),
       xlab = "Methylation (EWAS scale)", col = "grey80", border = "white")
  plot(res_i$df[[exp_var]], res_i$df$meth, pch = 19, cex = 0.7,
       col = ifelse(res_i$w < 0.5, "red", "grey40"),
       xlab = exp_var, ylab = cpg_i,
       main = paste0("n low weight = ", sum(res_i$w < 0.5),
                     " | SE ratio rel = ", round(flagged_all$se_ratio_rel[i], 2)))
}
dev.off()

diag_large_ci <- do.call(rbind, diag_list)
names(diag_large_ci) <- sub("\\..*$", "", names(diag_large_ci))   # drop "beta_refit.<exp_var>" suffix

# Check: refitted SE reproduces EWAS SE (ratio should be ~1)
summary(diag_large_ci$se_refit / diag_large_ci$SE_sepages)

## ---- 4.4 Save ----
write.table(diag_large_ci, paste0(dir_out, "flagged_cpgs_large_ci_diagnostics.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
diag_large_ci