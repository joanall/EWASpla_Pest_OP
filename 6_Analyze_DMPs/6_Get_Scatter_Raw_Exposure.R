#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the
# BiSC and SEPAGES cohorts.
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org
# Date: November 2026
#==============================================================================#
# [DESCRIPTION]: Outlier check of the top CpGs of one exposure (default TCPy).
# For the top CpGs of the main model (p < 1e-5, max n_top_plot), scatter plot of
# the RAW urinary concentration (x, ng/mL) against placental DNAm (y, beta
# values), one panel per CpG. Points above P99 (or below P1) of the exposure in
# the analytic sample are highlighted (same cut-offs as s6 winsorizing).
# Optional (add_fit_lines = TRUE): two crude (unadjusted) lm lines per panel,
# all women vs excluding outliers, plus Table 1 with both slopes.
#
# Plot 1: x = raw concentration (linear scale)
# Plot 2: x = raw concentration on a log2 axis (closer to what the model saw)
# Each plot: multi-page PDF (panels_per_page CpGs per page) + one PNG per CpG
# Table 1 (only if add_fit_lines): per CpG, crude slope with / without outliers
#
# Arguments:
#   1) cohort:     sepages / bisc
#   2) time point: t2 / t3 / avg
#   3) sample:     2tp / tp (BiSC only)
#   4) run date:   YYYYMMDD of the EWAS run
#   5) exposure:   24d excluded (binary); 3pba / tcpy / dmp / dep (default tcpy)
# Example:
#   Rscript script/6_Analyze_DMPs/6_Get_Scatter_Raw_Exposure.R bisc avg tp 20261022 tcpy
#
# [INPUT]:  results/5_Annotation/files_annotated_index.csv
#           -> results/5_Annotation/<run>/main/<exposure>_<tp>_annotated.rds
#           subsets (same as Run_EWAS.R) and preprocessed betas
# [OUTPUT]: results/6_Analyze_Results/scatter_raw/<run>/
#           scatter_raw_<exposure>_<tp>.pdf, scatter_log2_<exposure>_<tp>.pdf
#           solo/scatter_raw_<exposure>_<tp>_<CpG>.png (and _log2_)
#           slopes_outliers_<exposure>_<tp>.csv / .png (only if add_fit_lines)
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_template_plot.R")

index_file   <- "results/bisc/5_Annotation/files_annotated_index.csv"
cutoff_sugg  <- 1e-5          # top CpGs: suggestive in the main model
n_top_plot      <- 50            # max CpGs to plot (ordered by p-value)
panels_per_page <- 6             # CpGs per PDF page (2 rows x 3 columns)
outlier_prob    <- c(0.01, 0.99) # same percentiles as s6
add_fit_lines   <- FALSE         # TRUE: crude lm lines (all / without outliers) + Table 1

args <- commandArgs(trailingOnly = TRUE)
cohort     <- if (length(args) >= 1) args[1] else "bisc"
time_point <- if (length(args) >= 2) args[2] else "t2"
sample_set <- if (length(args) >= 3) args[3] else "tp"
run_date   <- if (length(args) >= 4) args[4] else "20261021"
exp_code   <- if (length(args) >= 5) args[5] else "tcpy"

stopifnot(cohort %in% c("sepages", "bisc"))
stopifnot(time_point %in% c("t2", "t3", "avg"))
stopifnot(exp_code %in% c("3pba", "tcpy", "dmp", "dep"))
# ΣDMP / ΣDEP have no single raw variable: the raw molar sums (pest_*_ms_i_*) are used

exp_labels <- c("3pba" = "3-PBA", "tcpy" = "TCPy", "dmp" = "ΣDMP", "dep" = "ΣDEP")
exp_label  <- exp_labels[[exp_code]]
exp_unit   <- if (exp_code %in% c("dmp", "dep")) "µmol/L" else "ng/mL"

run_name <- paste(cohort, sample_set, time_point, run_date, sep = "_")
outdir   <- file.path("results/6_Analyze_Results/scatter_raw", run_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

#==============================================================================#----
# COHORT SETTINGS (same as Run_EWAS.R)                                         #----
#==============================================================================#----

if (cohort == "sepages") {
  betas_path  <- "db/plaDNAm_data/SEPAGES_20230615_PreprocessedBetas_nooutliers.RData"
  subset_path <- "db/subsets/df_main.rds"
}
if (cohort == "bisc") {
  betas_path  <- "db/methylation_data/bVals_BISC_nosib_noout_nomask.RDS"
  subset_path <- if (sample_set == "2tp") "db/subsets_2tp/df_main.rds" else
    paste0("db/subsets/df_main_", time_point, ".rds")
}
stopifnot(file.exists(subset_path))

#==============================================================================#----
# TOP CpGs OF THE MAIN MODEL                                                   #----
#==============================================================================#----

files_index <- read.csv(index_file, stringsAsFactors = FALSE)
res_path <- files_index$path[files_index$run == run_name &
                               files_index$analysis == "main" &
                               files_index$exposure == exp_code]
res_path
stopifnot(length(res_path) == 1, file.exists(res_path))

res_main <- readRDS(res_path)
res_main <- res_main[order(res_main$pvalue), ]
top_cpgs <- head(res_main[!is.na(res_main$pvalue) & res_main$pvalue < cutoff_sugg, ], n_top_plot)
nrow(top_cpgs)    # N CpGs to plot

# If nothing is suggestive, take the top n_top_plot by p-value anyway
if (nrow(top_cpgs) == 0) {
  cat("No CpG with p <", cutoff_sugg, "- plotting the top", n_top_plot, "by p-value\n")
  top_cpgs <- head(res_main, n_top_plot)
}
top_cpgs[, c("CpG", "gene_name", "Effect", "pvalue")]

#==============================================================================#----
# EXPOSURE (RAW) AND METHYLATION                                               #----
#==============================================================================#----

metadata <- readRDS(subset_path)
nrow(metadata)

# Raw exposure variable for this time point
if (exp_code %in% c("3pba", "tcpy")) {
  raw_t2 <- paste0("pest_", exp_code, "_raw_t2")
  raw_t3 <- paste0("pest_", exp_code, "_raw_t3")
} else {
  raw_t2 <- paste0("pest_", exp_code, "_ms_i_t2")   # molar sums, natural scale
  raw_t3 <- paste0("pest_", exp_code, "_ms_i_t3")
}
setdiff(c(raw_t2, raw_t3), colnames(metadata))      # should be character(0)

if (time_point == "t2")  metadata$exp_raw <- as.numeric(metadata[[raw_t2]])
if (time_point == "t3")  metadata$exp_raw <- as.numeric(metadata[[raw_t3]])
if (time_point == "avg") {
  # Arithmetic mean of the raw values at the available time points
  metadata$exp_raw <- rowMeans(cbind(as.numeric(metadata[[raw_t2]]),
                                     as.numeric(metadata[[raw_t3]])), na.rm = TRUE)
  metadata$exp_raw[is.nan(metadata$exp_raw)] <- NA
}
summary(metadata$exp_raw)    # all >= 0
sum(is.na(metadata$exp_raw)) # women without raw value (dropped from the plot)

# Betas: load once, keep only the top CpGs
if (grepl("\\.rds$", betas_path, ignore.case = TRUE)) {
  meth_all <- readRDS(betas_path)
} else {
  betas_obj <- load(betas_path)
  meth_all  <- get(betas_obj[1])
  rm(list = betas_obj)
}
meth_top <- as.matrix(meth_all[rownames(meth_all) %in% top_cpgs$CpG, , drop = FALSE])
rm(meth_all); gc()
setdiff(top_cpgs$CpG, rownames(meth_top))   # should be character(0)

# Match samples
metadata <- metadata[!is.na(metadata$exp_raw) & metadata$Basename %in% colnames(meth_top), ]
meth_top <- meth_top[, metadata$Basename, drop = FALSE]
stopifnot(identical(colnames(meth_top), metadata$Basename))
nrow(metadata)   # N in the plot (should be close to N in the model)

# Outliers: below P1 or above P99 of the exposure in this sample
outlier_limits <- quantile(metadata$exp_raw, probs = outlier_prob, na.rm = TRUE, names = FALSE)
outlier_limits
metadata$outlier <- metadata$exp_raw < outlier_limits[1] | metadata$exp_raw > outlier_limits[2]
table(metadata$outlier)

#==============================================================================#----
# LONG DATA: ONE ROW PER WOMAN PER CpG                                         #----
#==============================================================================#----

make_label <- function(cpg, gene) {
  gene <- as.character(gene)
  gene[is.na(gene)] <- ""
  gene <- sapply(strsplit(gene, ";"), function(g) paste(unique(g[g != ""]), collapse = ";"))
  ifelse(gene == "", cpg, paste0(cpg, " (", gene, ")"))
}
cpg_label <- make_label(top_cpgs$CpG, top_cpgs$gene_name)

df_long <- NULL
for (i in seq_len(nrow(top_cpgs))) {
  cpg <- top_cpgs$CpG[i]
  df_long <- rbind(df_long, data.frame(
    CpG     = cpg,
    label   = cpg_label[i],
    exp_raw = metadata$exp_raw,
    beta    = as.numeric(meth_top[cpg, ]),
    outlier = metadata$outlier))
}
df_long$label   <- factor(df_long$label, levels = cpg_label)   # ordered by p-value
df_long$outlier <- factor(ifelse(df_long$outlier, "Outlier (<P1 or >P99)", "Rest"),
                          levels = c("Rest", "Outlier (<P1 or >P99)"))
df_long <- df_long[!is.na(df_long$beta), ]
table(df_long$label, df_long$outlier)   # same N per CpG (except NA betas)

# Data without outliers, for the second fit line
df_long_no_out <- df_long[df_long$outlier == "Rest", ]

#==============================================================================#----
# TABLE 1: CRUDE SLOPES WITH / WITHOUT OUTLIERS                                #----
#==============================================================================#----
# Crude lm(beta ~ exposure), raw scale; not the EWAS model (no covariates, not rlm)
# Only if add_fit_lines = TRUE

if (add_fit_lines) {
  slopes_tab <- NULL
  for (i in seq_len(nrow(top_cpgs))) {
    cpg <- top_cpgs$CpG[i]
    d_all    <- df_long[df_long$CpG == cpg, ]
    d_no_out <- d_all[d_all$outlier == "Rest", ]
    fit_all    <- summary(lm(beta ~ exp_raw, data = d_all))$coefficients
    fit_no_out <- summary(lm(beta ~ exp_raw, data = d_no_out))$coefficients
    slopes_tab <- rbind(slopes_tab, data.frame(
      CpG              = cpg_label[i],
      ewas_effect      = signif(top_cpgs$Effect[i], 3),
      ewas_p           = signif(top_cpgs$pvalue[i], 3),
      n_all            = nrow(d_all),
      n_outliers       = nrow(d_all) - nrow(d_no_out),
      slope_all        = signif(fit_all["exp_raw", "Estimate"], 3),
      p_all            = signif(fit_all["exp_raw", "Pr(>|t|)"], 3),
      slope_no_out     = signif(fit_no_out["exp_raw", "Estimate"], 3),
      p_no_out         = signif(fit_no_out["exp_raw", "Pr(>|t|)"], 3)))
  }
  slopes_tab$same_sign <- sign(slopes_tab$slope_all) == sign(slopes_tab$slope_no_out)
  slopes_tab
  save_table(slopes_tab, file.path(outdir, paste0("slopes_outliers_", exp_code, "_", time_point)))
}

#==============================================================================#----
# PLOTS                                                                        #----
#==============================================================================#----

# Plot 2 (log2 axis) uses values > 0 only
df_long_pos <- df_long[df_long$exp_raw > 0, ]
nrow(df_long) - nrow(df_long_pos)   # rows with raw value = 0, dropped in Plot 2

subtitle_txt <- paste0(run_name, ", main model | N = ", nrow(metadata),
                       "\nOutliers (red): < ", signif(outlier_limits[1], 3),
                       " or > ", signif(outlier_limits[2], 3), " ", exp_unit)
if (add_fit_lines) subtitle_txt <- paste0(subtitle_txt,
                                          "\nLines: crude lm, all women (solid) vs without outliers (dashed)")

# Scatter of exposure (x) vs DNAm beta (y), one panel per CpG in df_plot;
# red = exposure outliers; log2_axis = TRUE puts x on a log2 axis.
make_scatter <- function(df_plot, log2_axis = FALSE, n_col = 3, n_row = NULL) {
  x_lab <- paste0(exp_label, " ", time_point, " (", exp_unit,
                  if (log2_axis) ", log2 axis" else "", ")")
  p <- ggplot(df_plot, aes(x = exp_raw, y = beta)) +
    geom_point(aes(colour = outlier), size = 1.1, alpha = 0.6) +
    scale_colour_manual(values = c("Rest" = "grey55",
                                   "Outlier (<P1 or >P99)" = pal_bright[["red"]]),
                        drop = FALSE) +
    facet_wrap(~ label, scales = "free", ncol = n_col, nrow = n_row, drop = FALSE) +
    labs(title = paste0(exp_label, if (log2_axis) " (raw, log2 axis)" else " (raw)",
                        " vs placental DNAm"),
         subtitle = subtitle_txt, x = x_lab, y = "DNAm (beta value)", colour = NULL) +
    theme_ewas() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 0.5))
  if (log2_axis) p <- p + scale_x_continuous(trans = "log2")
  if (add_fit_lines) {
    df_no_out <- df_plot[df_plot$outlier == "Rest", ]
    p <- p +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                  colour = "black", linewidth = 0.6) +
      geom_smooth(data = df_no_out, method = "lm", formula = y ~ x, se = FALSE,
                  colour = pal_bright[["blue"]], linetype = "dashed", linewidth = 0.6)
  }
  p
}

# Pages: CpGs in order of p-value, panels_per_page per page
page_id <- ceiling(seq_along(cpg_label) / panels_per_page)
n_pages <- max(page_id)
n_pages

# Plot 1 and Plot 2. Multi-page PDF, panels_per_page CpGs per page (2 x 3).
for (log2_axis in c(FALSE, TRUE)) {
  df_scale  <- if (log2_axis) df_long_pos else df_long
  scale_tag <- if (log2_axis) "log2" else "raw"
  pdf_file  <- file.path(outdir, paste0("scatter_", scale_tag, "_", exp_code, "_", time_point, ".pdf"))
  
  pdf(pdf_file, width = 26 / 2.54, height = 20 / 2.54)   # cm -> inches
  for (pg in seq_len(n_pages)) {
    labels_page <- cpg_label[page_id == pg]
    df_page <- df_scale[df_scale$label %in% labels_page, ]
    # Last page: pad with empty slots (blank labels) so panels keep the same size
    n_empty <- panels_per_page - length(labels_page)
    pad_levels <- if (n_empty > 0) strrep(" ", seq_len(n_empty)) else character(0)
    df_page$label <- factor(df_page$label, levels = c(labels_page, pad_levels))
    # Fixed 2 x 3 grid: same panel size on the last (incomplete) page
    p_page <- make_scatter(df_page, log2_axis = log2_axis, n_col = 3, n_row = 2) +
      labs(caption = paste0("Page ", pg, " of ", n_pages, " | CpGs ordered by EWAS p-value"))
    print(p_page)
  }
  dev.off()
  cat("Saved:", pdf_file, "\n")
}

# Solo plots: one PNG per CpG and scale, in solo/
solo_dir <- file.path(outdir, "solo")
dir.create(solo_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(cpg_label)) {
  cpg <- top_cpgs$CpG[i]
  for (log2_axis in c(FALSE, TRUE)) {
    df_scale  <- if (log2_axis) df_long_pos else df_long
    scale_tag <- if (log2_axis) "log2" else "raw"
    df_cpg <- df_scale[df_scale$CpG == cpg, ]
    df_cpg$label <- droplevels(df_cpg$label)   # one panel only
    p_solo <- make_scatter(df_cpg, log2_axis = log2_axis, n_col = 1)
    save_plot(p_solo, file.path(solo_dir, paste0("scatter_", scale_tag, "_", exp_code, "_",
                                                 time_point, "_", sprintf("%02d", i), "_", cpg, ".png")),
              width = 14, height = 13)
  }
}
length(list.files(solo_dir, pattern = paste0("_", exp_code, "_", time_point, "_")))  # 2 x N CpGs

cat("Saved in:", outdir, "\n")

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----