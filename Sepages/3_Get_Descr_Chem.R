#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#

# [DESCRIPTION]: Descriptive statistics of the chemical compounds in SEPAGES
# 1) Table 1: LOD and LOQ of each compound (ug/L)
# 2) Table 2: distribution of each compound at t2 and t3 (raw, ug/L)
#    - N, % > LOD, % > LOQ, P5, P50, P95 (percentiles < LOD shown as "<LOD")
# 3) Boxplots of the continuous exposures (log2, imputed) at t2 and t3
# 3) Boxplots of raw concentrations at t2 and t3 (detected values, log10 axis)
# 4) Spearman correlations between compounds at t2 and t3 (raw)
# 5) ICC between t2 and t3 for each continuous exposure (log2 imputed)
#
#
# [INPUT]:  db/subsets/df_main_t2.rds, df_main_t3.rds
# [OUTPUT]: results/descriptives_chem/table1_lod_loq.csv / .png
#           results/descriptives_chem/table2_distribution.csv / .png
#           results/descriptives_chem/boxplot_exposures.png
#           results/descriptives_chem/cor_compounds_t2.png, cor_compounds_t3.png
#           results/descriptives_chem/icc_t2_t3.csv / .png

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_template_plot.R")

outdir <- "results/descriptives_chem/"

df_main_t2 <- readRDS("db/subsets/df_main_t2.rds")
df_main_t3 <- readRDS("db/subsets/df_main_t3.rds")
nrow(df_main_t2); nrow(df_main_t3)

lod_placeholder <- 0.00001   # raw values < LOD are stored as 0.00001


# Raw compounds: same order in all three vectors
compounds <- c("2,4-D", "3-PBA", "TCPy", "DMP", "DMTP", "DMDTP", "DEP", "DETP", "SG")

chem_raw_t2 <- c("mo_pest_24d_raw_t2", "mo_pest_3pba_raw_t2", "mo_pest_tcpy_raw_t2",
                 "mo_pest_dmp_raw_t2", "mo_pest_dmtp_raw_t2", "mo_pest_dmdtp_raw_t2",
                 "mo_pest_dep_raw_t2", "mo_pest_detp_raw_t2", "mo_pool_sg_T1")
chem_raw_t3 <- c("mo_pest_24d_raw_t3", "mo_pest_3pba_raw_t3", "mo_pest_tcpy_raw_t3",
                 "mo_pest_dmp_raw_t3", "mo_pest_dmtp_raw_t3", "mo_pest_dmdtp_raw_t3",
                 "mo_pest_dep_raw_t3", "mo_pest_detp_raw_t3", "mo_pool_sg_T3")

#==============================================================================#----
# 1. TABLE 1 - LOD AND LOQ                                                     #----
#==============================================================================#----
# OP metabolites from the lab documentation (ug/L).
# PENDING: LOD/LOQ of 2,4-D, 3-PBA and TCPy (fill in from documentation)
# Note: DEDTP was measured by the lab but is not in our data.

lod_tab <- data.frame(
  compound = c("DMP", "DMTP", "DMDTP", "DEP", "DETP", "DEDTP", "2,4-D", "3-PBA", "TCPy"),
  stem     = c("dmp", "dmtp", "dmdtp", "dep", "detp", "dedtp", "24d",   "3pba",  "tcpy"),
  LOD      = c(1.20,  0.20,   0.21,    0.38,  0.33,   0.15,    NA,      NA,      NA),
  LOQ      = c(4.00,  0.66,   0.63,    1.30,  1.00,   0.50,    NA,      NA,      NA)
)
lod_tab

# Table 1. LOD and LOQ (ug/L) of each pesticide metabolite.
save_table(lod_tab[, c("compound", "LOD", "LOQ")], paste0(outdir, "table1_lod_loq"))

#==============================================================================#----
# 2. TABLE 2 - DISTRIBUTION AT T2 AND T3                                       #----
#==============================================================================#----
# % > LOD and % > LOQ from raw values vs LOD/LOQ.
# If LOD is not known yet: % > LOD from raw > placeholder, % > LOQ = NA.
# ug/L; NA = not known yet (2,4-D, 3-PBA, TCPy) or not applicable (SG)
lod <- c(NA, NA, NA, 1.20, 0.20, 0.21, 0.38, 0.33, NA)
loq <- c(NA, NA, NA, 4.00, 0.66, 0.63, 1.30, 1.00, NA)

# Check: all variables exist
setdiff(chem_raw_t2, colnames(df_main_t2))
setdiff(chem_raw_t3, colnames(df_main_t3))

# N, % > LOD, % > LOQ, P5, P50, P95 (percentiles below LOD -> "<LOD")
desc_chem <- function(x, lod, loq) {
  q <- quantile(x, c(0.05, 0.50, 0.95), na.rm = TRUE)
  q_txt <- ifelse(!is.na(lod) & q < lod, "<LOD", sprintf("%.2f", q))
  c(N         = sum(!is.na(x)),
    pct_lod   = round(100 * mean(x >= lod, na.rm = TRUE), 1),
    pct_loq   = round(100 * mean(x >= loq, na.rm = TRUE), 1),
    P5 = q_txt[1], P50 = q_txt[2], P95 = q_txt[3])
}

tab2_t2 <- t(mapply(desc_chem, df_main_t2[chem_raw_t2], lod, loq))
tab2_t3 <- t(mapply(desc_chem, df_main_t3[chem_raw_t3], lod, loq))
rownames(tab2_t2) <- chem_raw_t2 
rownames(tab2_t3) <- chem_raw_t3
colnames(tab2_t2) <- c("N", "% > LOD", "% > LOQ", "P5", "P50", "P95")
colnames(tab2_t3) <- colnames(tab2_t2)

# t2 and t3 side by side
tab2 <- cbind(tab2_t2, tab2_t3)
colnames(tab2) <- paste(colnames(tab2), rep(c("t2", "t3"), each = 6))
tab2

# Table 2. Distribution of pesticide metabolites (ug/L) and specific gravity
# at t2 and t3: N, % above LOD and LOQ, P5, P50, P95.
save_table(tab2, paste0(outdir, "table2_distribution"))

#==============================================================================#----
# 3. BOXPLOTS - RAW CONCENTRATIONS                                             #----
#==============================================================================#----
# Non-detects (stored as 0.00001) excluded; y-axis on log10 scale.
df_box <- NULL
for (i in seq_along(compounds)) {
  df_box <- rbind(df_box,
                  data.frame(compound = compounds[i], time = "t2",
                             value = df_main_t2[[chem_raw_t2[i]]]),
                  data.frame(compound = compounds[i], time = "t3",
                             value = df_main_t3[[chem_raw_t3[i]]]))
}
df_box$compound <- factor(df_box$compound, levels = compounds)
df_box <- df_box[!is.na(df_box$value) & df_box$value > lod_placeholder, ]

# Check: N detected per compound and time point
table(df_box$compound, df_box$time)

# Plot 1. Boxplots of raw concentration (y, log10 scale) per time point
# (x, fill), one panel per compound; detected values only.
p1 <- ggplot(df_box, aes(x = time, y = value, fill = time)) +
  geom_boxplot(outlier.size = 0.6, width = 0.6) +
  facet_wrap(~ compound, scales = "free_y", nrow = 2) +
  scale_y_log10() +
  scale_fill_manual(values = pal_time) + 
  labs(title = "Pesticide exposure during pregnancy",
       subtitle = "SEPAGES, raw concentrations (detected values only, log10 scale)",
       x = "Time point", y = "Concentration (ug/L)", fill = "Time point") +
  theme_ewas() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
p1
save_plot(p1, paste0(outdir, "boxplot_exposures.png"), width = 20, height = 14)

#==============================================================================#----
# 4. CORRELATIONS - RAW CONCENTRATIONS                                         #----
#==============================================================================#----
# Spearman: rank-based, so non-detects (0.00001) tie at the lowest rank.

cor_t2 <- cor(df_main_t2[, chem_raw_t2], method = "spearman", use = "pairwise.complete.obs")
cor_t3 <- cor(df_main_t3[, chem_raw_t3], method = "spearman", use = "pairwise.complete.obs")
dimnames(cor_t2) <- list(compounds, compounds)
dimnames(cor_t3) <- list(compounds, compounds)
round(cor_t2, 2); round(cor_t3, 2)

for (tp in c("t2", "t3")) {
  cor_mat <- if (tp == "t2") cor_t2 else cor_t3
  
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
         subtitle = "SEPAGES, Spearman rho (raw concentrations)",
         x = NULL, y = NULL, fill = "rho") +
    coord_fixed() +
    theme_ewas()
  print(p2)
  save_plot(p2, paste0(outdir, "cor_compounds_", tp, ".png"), width = 16, height = 14)
}


#==============================================================================#----
# 5. REPRODUCIBILITY T2 - T3 - ICC                                             #----
#==============================================================================#----
# ICC = between-woman variance / (between + within), random-intercept model
# on log2 imputed values (raw values with the placeholder would distort it)

cont_exp <- data.frame(
  compound = c("3-PBA", "TCPy", "DMP (sum)", "DEP (sum)", "DAP (sum)"),
  stem     = c("3pba_i", "tcpy_i", "dmp_ms_i", "dep_ms_i", "dap_ms_i")
)

# Check: log2 variables exist in the subsets
setdiff(paste0("mo_pest_", cont_exp$stem, "_log2_t2"), colnames(df_main_t2))
setdiff(paste0("mo_pest_", cont_exp$stem, "_log2_t3"), colnames(df_main_t3))

icc_tab <- data.frame(Compound = cont_exp$compound, N_women = NA, N_both = NA, ICC = NA)
for (i in seq_len(nrow(cont_exp))) {
  v_t2 <- paste0("mo_pest_", cont_exp$stem[i], "_log2_t2")
  v_t3 <- paste0("mo_pest_", cont_exp$stem[i], "_log2_t3")
  
  # Long format: one row per woman and time point
  df_long <- rbind(data.frame(ident = df_main_t2$ident, value = df_main_t2[[v_t2]]),
                   data.frame(ident = df_main_t3$ident, value = df_main_t3[[v_t3]]))
  df_long <- df_long[!is.na(df_long$value), ]
  df_long$ident <- factor(df_long$ident)
  
  # One-way ANOVA: between-woman (MSB) and within-woman (MSW) mean squares
  aov_tab <- anova(lm(value ~ ident, data = df_long))
  msb <- aov_tab["ident", "Mean Sq"]
  msw <- aov_tab["Residuals", "Mean Sq"]
  
  # n0: average samples per woman, corrected for unbalanced data
  n_i <- table(df_long$ident)
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

# Table 3. Intraclass correlation coefficient (ICC) between t2 and t3 for each
# continuous exposure (log2 imputed concentrations).
save_table(icc_tab, paste0(outdir, "icc_t2_t3"))
#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----