#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the
# BiSC and SEPAGES cohorts.
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: CpG-specific diagnostics for large SEPAGES CIs (follows 9_Check_CI.R)
# Run on the PC holding the SEPAGES methylation data
# For each top-hit exposure x time point x CpG, refit the EWAS model
#         in SEPAGES and identify why the SE is large:
#           - fewer samples (probe filtering)
#           - low / high methylation variability (vs all CpGs)
#           - outliers or many down-weighted samples in the robust regression
#           - non-convergence

# [INPUT]: results/7_Meta/9_Check_CI/tophits_ci_detail.txt (from 9_Check_CI.R,
#         copied from the BiSC PC to the same path)
# [OUTPUT]: results/sepages/7_Meta/9_Check_CI/ - table of diagnostics (all top
#         hits) + PNG/PDF plots. Kept under results/sepages/ so that copying
#         results/sepages/ to the BiSC PC brings it along.
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
library(MASS)
source("script/00_paths.R")
source("script/00_template_plot.R")   # palettes, save_table()


dir_out <- paste0(results_dir("sepages", "7_Meta", "9_Check_CI", create = TRUE), "/")

#==============================================================================#----
# Get Hit CpGs to be inspected                                                 #----
#==============================================================================#----
flagged_all <- read.table(results_dir("meta", "9_Check_CI", "tophits_ci_detail.txt"),
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Check: top hits and flagged rows
nrow(flagged_all); sum(flagged_all$flag_large_ci)
cpg_large_ci <- unique(flagged_all$CpG)

# Methylation: CpG x sample, same scale as EWAS (e.g. M-values)  <-- adjust
load(db_betas("sepages"))
meth_sepages <- betafinal.nooutlier

# Phenotype: one row per sample, sample_id matching colnames(meth_sepages)  <-- adjust
pheno_sepages <- readRDS(db_subset("sepages", "main", "2tp"))   # same for all the time points

# Check: CpGs and samples present
sum(cpg_large_ci %in% rownames(meth_sepages)); length(cpg_large_ci)
sum(pheno_sepages$Basename %in% colnames(meth_sepages)); nrow(pheno_sepages)

# Exposure variable per exposure x time point  <-- adjust to real names
exp_vars <- c(
  "24d_t2"  = "pest_24d_bi_t2",       "24d_t3"  = "pest_24d_bi_t3",       "24d_avg"  = "pest_24d_any_avg",
  "3pba_t2" = "pest_3pba_ilog2_t2",    "3pba_t3" = "pest_3pba_ilog2_t3",    "3pba_avg" = "pest_3pba_ilog2_avg",
  "tcpy_t2" = "pest_tcpy_ilog2_t2",    "tcpy_t3" = "pest_tcpy_ilog2_t3",    "tcpy_avg" = "pest_tcpy_ilog2_avg",
  "dmp_t2"  = "pest_dmp_ms_ilog2_t2",  "dmp_t3"  = "pest_dmp_ms_ilog2_t3",  "dmp_avg"  = "pest_dmp_ms_ilog2_avg",
  "dep_t2"  = "pest_dep_ms_ilog2_t2",  "dep_t3"  = "pest_dep_ms_ilog2_t3",  "dep_avg"  = "pest_dep_ms_ilog2_avg"
)

# Covariates (main model, as in Run_EWAS)  <-- adjust
main_cov   <- c("mo_age", "mo_edu", "mo_smoke", "gest_age",
                "ch_eth", "ch_sex", "Meanlog2oddsContamination")
cell_types <- c("Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC")
covars <- c(main_cov, cell_types)

# Check: all variables exist
setdiff(c(exp_vars, covars), colnames(pheno_sepages))

#==============================================================================#----
# Genome-wide reference: methylation SD across all CpGs                        #----
#==============================================================================#----
meth_sub    <- meth_sepages[, pheno_sepages$Basename]   # samples used in the EWAS
sd_meth_all <- apply(meth_sub, 1, sd, na.rm = TRUE)
summary(sd_meth_all)

#==============================================================================#----
# Function: refit one CpG for one exposure x time point                        #----
#==============================================================================#----
diag_cpg <- function(cpg, exp_var, covars, meth, pheno, sd_meth_all) {
  pheno_cc <- pheno[complete.cases(pheno[, c(exp_var, covars)]), ]
  df_i     <- data.frame(meth = meth[cpg, pheno_cc$Basename], pheno_cc)
  df_i     <- df_i[!is.na(df_i$meth), ]
  
  form_i <- as.formula(paste("meth ~", exp_var, "+", paste(covars, collapse = " + ")))
  fit_i  <- MASS::rlm(form_i, data = df_i, maxit = 200)   # <-- same settings as Run_EWAS
  
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
      beta_refit       = unname(coef(fit_i)[exp_var]),
      se_refit         = summary(fit_i)$coefficients[exp_var, "Std. Error"]
    ),
    df  = df_i,
    w   = fit_i$w
  )
}

#==============================================================================#----
# Run for all top hits + plots                                                 #----
#==============================================================================#----

dir_plots_flag   <- paste0(dir_out, "plots/flagged/")
dir_plots_noflag <- paste0(dir_out, "plots/not_flagged/")
dir.create(dir_plots_flag,   recursive = TRUE, showWarnings = FALSE)
dir.create(dir_plots_noflag, recursive = TRUE, showWarnings = FALSE)


# Plot function: methylation distribution + methylation vs exposure (red = weight < 0.5)
plot_diag <- function(plot_data) {
  par(mfrow = c(1, 2))
  hist(plot_data$meth, breaks = 40,
       main = paste(plot_data$cpg, "-", plot_data$key),
       xlab = "Methylation (EWAS scale)", col = pal_bright[["grey"]], border = "white")
  plot(plot_data$exposure, plot_data$meth, pch = 19, cex = 0.7,
       col = ifelse(plot_data$w < 0.5, pal_bright[["red"]], "grey40"),
       xlab = plot_data$exp_var, ylab = plot_data$cpg,
       main = paste0("n low weight = ", sum(plot_data$w < 0.5),
                     " | SE ratio rel = ", round(plot_data$se_ratio_rel, 2),
                     ifelse(plot_data$flag, " | FLAGGED", "")))
}

diag_list      <- list()
plot_data_list <- list()
plot_rec_list  <- list()

for (i in seq_len(nrow(flagged_all))) {
  cpg_i   <- flagged_all$CpG[i]
  key_i   <- paste(flagged_all$exposure[i], flagged_all$tp[i], sep = "_")
  exp_var <- exp_vars[key_i]
  name_i  <- paste(cpg_i, key_i, sep = "_")
  
  if (!cpg_i %in% rownames(meth_sepages)) {
    warning(cpg_i, " not in SEPAGES methylation matrix"); next
  }
  
  res_i <- diag_cpg(cpg_i, exp_var, covars, meth_sepages, pheno_sepages, sd_meth_all)
  diag_list[[name_i]] <- data.frame(flagged_all[i, c("exposure", "tp", "CpG", "pvalue",
                                                     "SE_bisc", "SE_sepages", "se_ratio_rel",
                                                     "pct_w_sepages", "flag_large_ci")],
                                    t(res_i$stats), row.names = NULL)
  
  # Plot data (robust object to replot later)
  plot_data_i <- list(cpg = cpg_i, key = key_i, exp_var = unname(exp_var),
                      meth = res_i$df$meth, exposure = res_i$df[[exp_var]],
                      w = res_i$w, se_ratio_rel = flagged_all$se_ratio_rel[i],
                      flag = flagged_all$flag_large_ci[i],
                      sample_id = res_i$df$Basename)
  plot_data_list[[name_i]] <- plot_data_i
  
  # PNG + recorded plot
  dir_png_i <- ifelse(flagged_all$flag_large_ci[i], dir_plots_flag, dir_plots_noflag)
  png(paste0(dir_png_i, name_i, ".png"), width = 10, height = 4.5, units = "in", res = 300)
  dev.control(displaylist = "enable")
  plot_diag(plot_data_i)
  if (flagged_all$flag_large_ci[i]) plot_rec_list[[name_i]] <- recordPlot()
  dev.off()
}

# Check: one PNG per top hit
length(list.files(dir_plots_flag,   pattern = "\\.png$")); sum(flagged_all$flag_large_ci)
length(list.files(dir_plots_noflag, pattern = "\\.png$")); sum(!flagged_all$flag_large_ci)

# Combined PDFs (one page per plot)
name_flag <- names(plot_data_list)[sapply(plot_data_list, function(x) x$flag)]

pdf(paste0(dir_out, "tophits_ci_plots_flagged.pdf"), width = 10, height = 4.5)
for (name_i in name_flag) plot_diag(plot_data_list[[name_i]])
dev.off()

pdf(paste0(dir_out, "tophits_ci_plots_all.pdf"), width = 10, height = 4.5)
for (name_i in names(plot_data_list)) plot_diag(plot_data_list[[name_i]])
dev.off()

# Save plot objects
saveRDS(plot_data_list, paste0(dir_out, "tophits_ci_plot_data.rds"))
saveRDS(plot_rec_list,  paste0(dir_out, "tophits_ci_plot_recorded.rds"))

diag_large_ci <- do.call(rbind, diag_list)

# Check: refitted SE reproduces EWAS SE (ratio should be ~1)
summary(diag_large_ci$se_refit / diag_large_ci$SE_sepages)

# Flagged vs non-flagged top hits
aggregate(cbind(n_used, sd_meth_pctl, s_resid, n_low_w, n_meth_out, n_exp_out) ~ flag_large_ci,
          data = diag_large_ci, FUN = median)
table(diag_large_ci$flag_large_ci, diag_large_ci$converged)

save_table(diag_large_ci, paste0(dir_out, "tophits_ci_diagnostics"), png = FALSE)

#==============================================================================#----
# Genome-wide CpG metrics (exposure-independent)                               #----
#==============================================================================#----
n_nona_all  <- rowSums(!is.na(meth_sub))
med_all     <- apply(meth_sub, 1, median, na.rm = TRUE)
mad_all     <- apply(meth_sub, 1, mad, na.rm = TRUE)
n_out_all   <- rowSums(abs(meth_sub - med_all) > 4 * mad_all, na.rm = TRUE)

cpg_metrics_sepages <- data.frame(
  CpG          = rownames(meth_sub),
  n_nona       = n_nona_all,
  sd_meth      = sd_meth_all[rownames(meth_sub)],
  sd_meth_pctl = round(100 * rank(sd_meth_all[rownames(meth_sub)]) / nrow(meth_sub), 1),
  n_meth_out   = n_out_all
)

saveRDS(cpg_metrics_sepages, paste0(dir_out, "cpg_metrics_sepages_all.rds"))
#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----