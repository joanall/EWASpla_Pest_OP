#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: Descriptive statistics of the chemical compounds, per cohort,
# in the analytical sample:
# 1) Table: distribution of each compound at t2 and t3 (raw, ug/L)
#    N, % > LOD, % > LOQ, P5, P50, P95 (percentiles < LOD shown as "<LOD")
# 2) Boxplots of raw concentrations at t2 and t3 (values >= LOD, log10 axis)
# 3) Spearman correlations between compounds at t2 and t3 (raw)
# 4) ICC between t2 and t3 for each continuous exposure (log2 imputed)
# Argument: cohort (sepages / bisc)
#   Rscript script/1_Descriptives/3_Get_Descr_Chem.R sepages
#
# [INPUT]:  SEPAGES: db/subsets/df_main.rds (both time points)
#           BiSC:    db/subsets/df_main_t2.rds, df_main_t3.rds
# [OUTPUT]: results/1_Descriptives/descr_chem/<cohort>/
#           table_pest_distribution.csv/.png, boxplot_exposures.png,
#           cor_compounds_t2.png, cor_compounds_t3.png, icc_t2_t3.csv/.png
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_template_plot.R")

args   <- commandArgs(trailingOnly = TRUE)
cohort <- if (length(args) >= 1) args[1] else "sepages"
stopifnot(cohort %in% c("sepages", "bisc"))

outdir <- file.path("results/1_Descriptives/descr_chem", cohort, "")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cohort_label <- if (cohort == "sepages") "SEPAGES" else "BiSC"

# Compounds: same order in all vectors
compounds <- c("2,4-D", "3-PBA", "TCPy", "DMP", "DMTP", "DMDTP", "DEP", "DETP", "SG")
chem_raw_t2 <- c("pest_24d_raw_t2", "pest_3pba_raw_t2", "pest_tcpy_raw_t2",
                 "pest_dmp_raw_t2", "pest_dmtp_raw_t2", "pest_dmdtp_raw_t2",
                 "pest_dep_raw_t2", "pest_detp_raw_t2", "pool_sg_t2")
chem_raw_t3 <- sub("_t2$", "_t3", chem_raw_t2)

#==============================================================================#----
# COHORT SETTINGS                                                              #----
#==============================================================================#----
# CHECK the LOD/LOQ values (ug/L) against each lab's documentation

if (cohort == "sepages") {
  df_t2 <- readRDS("db/subsets/df_main.rds")    # same women at t2 and t3
  df_t3 <- df_t2
  lod <- c(0.15, 0.20, 0.10, 1.20, 0.20, 0.21, 0.38, 0.33, NA)
  loq <- c(0.50, 0.60, 0.30, 4.00, 0.66, 0.63, 1.30, 1.00, NA)
}

if (cohort == "bisc") {
  df_t2 <- readRDS("db/subsets/df_main_t2.rds")
  df_t3 <- readRDS("db/subsets/df_main_t3.rds")
  lod <- c(0.15, 0.20, 0.10, 1.20, 0.20, 0.21, 0.38, 0.33, NA)   
  loq <- c(0.50, 0.60, 0.30, 4.00, 0.66, 0.63, 1.30, 1.00, NA)
}
nrow(df_t2); nrow(df_t3)
data.frame(compounds, lod, loq)   # check alignment

# Keep only compounds available at both time points (e.g. DMDTP may be missing)
keep <- chem_raw_t2 %in% colnames(df_t2) & chem_raw_t3 %in% colnames(df_t3)
compounds[!keep]                  # dropped compounds
compounds   <- compounds[keep]
chem_raw_t2 <- chem_raw_t2[keep]
chem_raw_t3 <- chem_raw_t3[keep]
lod <- lod[keep]
loq <- loq[keep]

# Raw values as numeric
for (v in chem_raw_t2) df_t2[[v]] <- as.numeric(df_t2[[v]])
for (v in chem_raw_t3) df_t3[[v]] <- as.numeric(df_t3[[v]])

#==============================================================================#----
# 1. TABLE - DISTRIBUTION AT T2 AND T3                                         #----
#==============================================================================#----

# N, % > LOD, % > LOQ, P5, P50, P95 (percentiles below LOD -> "<LOD")
desc_chem <- function(x, lod, loq) {
  q <- quantile(x, c(0.05, 0.50, 0.95), na.rm = TRUE)
  q_txt <- ifelse(!is.na(lod) & q < lod, "<LOD", sprintf("%.2f", q))
  c(N       = sum(!is.na(x)),
    pct_lod = round(100 * mean(x >= lod, na.rm = TRUE), 1),
    pct_loq = round(100 * mean(x >= loq, na.rm = TRUE), 1),
    P5 = q_txt[1], P50 = q_txt[2], P95 = q_txt[3])
}

tab_t2 <- t(mapply(desc_chem, df_t2[chem_raw_t2], lod, loq))
tab_t3 <- t(mapply(desc_chem, df_t3[chem_raw_t3], lod, loq))

tab_dist <- cbind(tab_t2, tab_t3)
rownames(tab_dist) <- compounds
colnames(tab_dist) <- paste(rep(c("N", "% > LOD", "% > LOQ", "P5", "P50", "P95"), 2),
                            rep(c("t2", "t3"), each = 6))
tab_dist

# Table. Distribution of pesticide metabolites (ug/L) and specific gravity at
# t2 and t3: N, % above LOD and LOQ, P5, P50, P95.
save_table(tab_dist, paste0(outdir, "table_pest_distribution"), base_size = 9)

#==============================================================================#----
# 2. BOXPLOTS - RAW CONCENTRATIONS                                             #----
#==============================================================================#----
# Values below LOD excluded (works for any non-detect coding); log10 y-axis.

df_box <- NULL
for (i in seq_along(compounds)) {
  df_box <- rbind(df_box,
                  data.frame(compound = compounds[i], time = "t2", lod = lod[i],
                             value = df_t2[[chem_raw_t2[i]]]),
                  data.frame(compound = compounds[i], time = "t3", lod = lod[i],
                             value = df_t3[[chem_raw_t3[i]]]))
}
df_box$compound <- factor(df_box$compound, levels = compounds)
df_box <- df_box[!is.na(df_box$value) & df_box$value > 0 &
                   (is.na(df_box$lod) | df_box$value >= df_box$lod), ]

# Check: N >= LOD per compound and time point
table(df_box$compound, df_box$time)

# Plot 1. Boxplots of raw concentration (y, log10 scale) per time point
# (x, fill), one panel per compound; values >= LOD only.
p1 <- ggplot(df_box, aes(x = time, y = value, fill = time)) +
  geom_boxplot(outlier.size = 0.6, width = 0.6) +
  facet_wrap(~ compound, scales = "free_y", nrow = 2) +
  scale_y_log10() +
  scale_fill_manual(values = pal_time) +
  labs(title = "Pesticide exposure during pregnancy",
       subtitle = paste0(cohort_label, ", raw concentrations (values >= LOD, log10 scale)"),
       x = "Time point", y = "Concentration (ug/L)", fill = "Time point") +
  theme_ewas() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
save_plot(p1, paste0(outdir, "boxplot_exposures.png"), width = 20, height = 14)

#==============================================================================#----
# 3. CORRELATIONS - RAW CONCENTRATIONS                                         #----
#==============================================================================#----
# Spearman: rank-based, so non-detects (a constant value) tie at the lowest rank.

for (tp in c("t2", "t3")) {
  df_tp   <- if (tp == "t2") df_t2 else df_t3
  vars_tp <- if (tp == "t2") chem_raw_t2 else chem_raw_t3
  
  cor_mat <- cor(df_tp[, vars_tp], method = "spearman", use = "pairwise.complete.obs")
  dimnames(cor_mat) <- list(compounds, compounds)
  print(round(cor_mat, 2))
  
  # Long format for ggplot
  df_cor <- expand.grid(x = compounds, y = compounds)
  df_cor$rho <- as.vector(cor_mat)
  df_cor$x <- factor(df_cor$x, levels = compounds)
  df_cor$y <- factor(df_cor$y, levels = rev(compounds))
  
  # Plot 2. Spearman correlation heatmap between compounds at one time point;
  # fill and text = rho.
  p2 <- ggplot(df_cor, aes(x = x, y = y, fill = rho)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
    fill_diverging +
    labs(title = paste("Correlation between compounds -", tp),
         subtitle = paste0(cohort_label, ", Spearman rho (raw concentrations)"),
         x = NULL, y = NULL, fill = "rho") +
    coord_fixed() +
    theme_ewas()
  save_plot(p2, paste0(outdir, "cor_compounds_", tp, ".png"), width = 16, height = 14)
}

#==============================================================================#----
# 4. REPRODUCIBILITY T2 - T3 - ICC                                             #----
#==============================================================================#----
# ICC = between-woman variance / (between + within), one-way ANOVA estimator
# on log2 imputed values. SEPAGES: same women at both time points.
# BiSC: women from both subsets (unbalanced; n0 corrects for it).

cont_exp <- data.frame(
  compound = c("3-PBA",           "TCPy",            "DMP (sum)",
               "DEP (sum)",       "DAP (sum)"),
  var_t2   = c("pest_3pba_ilog2_t2",   "pest_tcpy_ilog2_t2",   "pest_dmp_ms_ilog2_t2",
               "pest_dep_ms_ilog2_t2", "pest_dap_ms_ilog2_t2")
)
cont_exp$var_t3 <- sub("_t2$", "_t3", cont_exp$var_t2)

# Check: log2 variables exist
setdiff(cont_exp$var_t2, colnames(df_t2))
setdiff(cont_exp$var_t3, colnames(df_t3))

icc_tab <- data.frame(Compound = cont_exp$compound, N_women = NA, N_both = NA, ICC = NA)
for (i in seq_len(nrow(cont_exp))) {
  
  # Long format: one row per woman and time point
  df_long <- rbind(data.frame(mo_id = df_t2$mo_id, value = df_t2[[cont_exp$var_t2[i]]]),
                   data.frame(mo_id = df_t3$mo_id, value = df_t3[[cont_exp$var_t3[i]]]))
  df_long <- df_long[!is.na(df_long$value), ]
  df_long$mo_id <- factor(df_long$mo_id)
  
  # One-way ANOVA: between-woman (MSB) and within-woman (MSW) mean squares
  aov_tab <- anova(lm(value ~ mo_id, data = df_long))
  msb <- aov_tab["mo_id", "Mean Sq"]
  msw <- aov_tab["Residuals", "Mean Sq"]
  
  # n0: average samples per woman, corrected for unbalanced data
  n_i <- table(df_long$mo_id)
  k   <- length(n_i)
  N   <- sum(n_i)
  n0  <- (N - sum(n_i^2) / N) / (k - 1)
  
  var_between <- max((msb - msw) / n0, 0)          # truncate at 0
  icc_tab$N_women[i] <- k
  icc_tab$N_both[i]  <- sum(n_i == 2)
  icc_tab$ICC[i]     <- round(var_between / (var_between + msw), 2)
}
colnames(icc_tab) <- c("Compound", "N women", "N with both samples", "ICC")
icc_tab

# Table. Intraclass correlation coefficient (ICC) between t2 and t3 for each
# continuous exposure (log2 imputed concentrations).
save_table(icc_tab, paste0(outdir, "icc_t2_t3"))

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----