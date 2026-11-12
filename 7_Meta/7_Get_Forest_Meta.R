#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the
# BiSC and SEPAGES cohorts.
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: Forest plots of the top CpGs of the BiSC + SEPAGES meta-analysis:
# for each CpG, effect (95% CI) in the meta-analysis, BiSC and SEPAGES, to see
# whether one cohort drives the meta-analysis result.
# - Top CpGs: meta p < cutoff_sugg (max n_top_plot); if none, top n_top_plot by p
# - Cohort effects: the GWAMA input files (same BETA / SE that went into GWAMA)
# - Right-hand text: weight of BiSC in the meta-analysis (inverse variance,
#   = 1/SE^2 share) and heterogeneity (I2, Cochran's Q p-value)
# - Filled point: p < cutoff_sugg in that estimate
# Label: CpG (gene) - genes de-duplicated; CpG only if intergenic.
#
# Arguments:
#   1) meta_id:  <tp>_<analysis>, e.g. t2_main, t3_main, avg_main, t2_s6
#   2) exposure: 24d / 3pba / tcpy / dmp / dep / dap / all (default all)
# Example:
#   Rscript script/7_Meta/7_Get_Forest_Meta.R t2_main tcpy
#
# [INPUT]:  results/7_Meta/3_Read_GWAMA/meta_<meta_id>.rds (list by exposure)
#           results/7_Meta/1_Input_GWAMA/<meta_id>/bisc_gwamain_<exposure>.txt,
#           sepages_gwamain_<exposure>.txt
#           db/annotation_data/cpg_annotation_all.rds (gene names)
# [OUTPUT]: results/7_Meta/7_Forest_Meta/<meta_id>/forest_meta_<exposure>.png
#           results/7_Meta/7_Forest_Meta/<meta_id>/forest_meta_<exposure>.xlsx
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_paths.R")
source("script/00_template_plot.R")

dir_input   <- results_dir("meta", "1_Input_GWAMA")
dir_read    <- results_dir("meta", "3_Read_GWAMA")
annot_file  <- annotation_file
cutoff_sugg <- 1e-5    # top CpGs and filled points
n_top_plot  <- 20      # max CpGs per plot

args <- commandArgs(trailingOnly = TRUE)
meta_id      <- if (length(args) >= 1) args[1] else "t2_main"
exposure_arg <- if (length(args) >= 2) args[2] else "all"

outdir <- results_dir("meta", "7_Forest_Meta", meta_id)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Meta-analysis results: named list by exposure
meta_file <- file.path(dir_read, paste0("meta_", meta_id, ".rds"))
stopifnot(file.exists(meta_file))
meta_list <- readRDS(meta_file)
names(meta_list)    # exposures available

exposures <- if (exposure_arg == "all") names(meta_list) else strsplit(exposure_arg, ",")[[1]]
setdiff(exposures, names(meta_list))    # should be character(0)
exposures <- intersect(exposures, names(meta_list))
stopifnot(length(exposures) > 0)

# Gene names
annotation <- readRDS(annot_file)
stopifnot(all(c("CpG", "gene_name") %in% names(annotation)))

# Labels and colours
exp_code   <- c("24d",   "3pba",  "tcpy", "dmp",  "dep",  "dap",  "mix")
exp_labels <- c("2,4-D", "3-PBA", "TCPy", "ΣDMP", "ΣDEP", "ΣDAP", "Mixture")
shape_source <- c("Meta-analysis" = 18, "BiSC" = 16, "SEPAGES" = 16)   # diamond for meta

make_label <- function(cpg, gene) {
  gene <- as.character(gene)
  gene[is.na(gene)] <- ""
  gene <- sapply(strsplit(gene, ";"), function(g) paste(unique(g[g != ""]), collapse = ";"))
  ifelse(gene == "", cpg, paste0(cpg, " (", gene, ")"))
}

#==============================================================================#----
# FOREST PLOTS - META vs COHORTS                                               #----
#==============================================================================#----

for (e in exposures) {
  
  exp_label <- exp_labels[match(e, exp_code)]
  
  # 1. Top CpGs of the meta-analysis (already ordered by p_value in 3_Read_GWAMA)
  meta_e <- meta_list[[e]]
  meta_e <- meta_e[order(meta_e$p_value), ]
  top_meta <- head(meta_e[meta_e$p_value < cutoff_sugg, ], n_top_plot)
  if (nrow(top_meta) == 0) {
    cat(e, ": no CpG with p <", cutoff_sugg, "- plotting the top", n_top_plot, "by p-value\n")
    top_meta <- head(meta_e, n_top_plot)
  }
  cpgs <- top_meta$probeID
  cat("\n--", meta_id, "|", e, "| CpGs plotted:", length(cpgs), "\n")
  
  # 2. Cohort estimates: the GWAMA input files
  gwama_bisc    <- read.delim(file.path(dir_input, meta_id, paste0("bisc_gwamain_",    e, ".txt")))
  gwama_sepages <- read.delim(file.path(dir_input, meta_id, paste0("sepages_gwamain_", e, ".txt")))
  m_bisc    <- match(cpgs, gwama_bisc$MARKERNAME)
  m_sepages <- match(cpgs, gwama_sepages$MARKERNAME)
  stopifnot(!anyNA(m_bisc), !anyNA(m_sepages))    # n_studies == 2, so all present
  
  beta_bisc <- gwama_bisc$BETA[m_bisc];       se_bisc <- gwama_bisc$SE[m_bisc]
  beta_sep  <- gwama_sepages$BETA[m_sepages]; se_sep  <- gwama_sepages$SE[m_sepages]
  
  # Check: inverse-variance weighted beta reproduces the GWAMA beta
  w_bisc <- 1 / se_bisc^2
  w_sep  <- 1 / se_sep^2
  beta_ivw <- (beta_bisc * w_bisc + beta_sep * w_sep) / (w_bisc + w_sep)
  stopifnot(max(abs(beta_ivw - top_meta$beta)) < 1e-4)
  
  # 3. Labels and per-CpG summary (also saved as table)
  gene  <- annotation$gene_name[match(cpgs, annotation$CpG)]
  label <- make_label(cpgs, gene)
  
  forest_tab <- data.frame(
    CpG          = cpgs,
    label        = label,
    beta_meta    = top_meta$beta,     se_meta    = top_meta$se,    p_meta    = top_meta$p_value,
    beta_bisc    = beta_bisc,         se_bisc    = se_bisc,        p_bisc    = gwama_bisc$P[m_bisc],
    beta_sepages = beta_sep,          se_sepages = se_sep,         p_sepages = gwama_sepages$P[m_sepages],
    weight_bisc  = round(100 * w_bisc / (w_bisc + w_sep), 1),     # % of the meta weight
    i2           = round(top_meta$i2, 2),
    q_p          = signif(top_meta$q_p_value, 2),
    direction    = top_meta$effects,                              # GWAMA, e.g. "++"
    stringsAsFactors = FALSE)
  forest_tab$same_direction <- sign(forest_tab$beta_bisc) == sign(forest_tab$beta_sepages)
  print(forest_tab[, c("label", "p_meta", "p_bisc", "p_sepages", "weight_bisc", "i2", "same_direction")])
  save_table(forest_tab, file.path(outdir, paste0("forest_meta_", e)), png = FALSE)
  
  # 4. Long data: one row per CpG per estimate
  df_forest <- rbind(
    data.frame(label = label, source = "Meta-analysis", Effect = forest_tab$beta_meta,
               se = forest_tab$se_meta, p = forest_tab$p_meta),
    data.frame(label = label, source = "BiSC", Effect = forest_tab$beta_bisc,
               se = forest_tab$se_bisc, p = forest_tab$p_bisc),
    data.frame(label = label, source = "SEPAGES", Effect = forest_tab$beta_sepages,
               se = forest_tab$se_sepages, p = forest_tab$p_sepages))
  df_forest$lower      <- df_forest$Effect - 1.96 * df_forest$se
  df_forest$upper      <- df_forest$Effect + 1.96 * df_forest$se
  df_forest$suggestive <- df_forest$p < cutoff_sugg
  df_forest$label  <- factor(df_forest$label, levels = rev(label))         # top CpG at the top
  df_forest$source <- factor(df_forest$source, levels = c("Meta-analysis", "BiSC", "SEPAGES"))
  table(df_forest$source, df_forest$suggestive)
  
  # Right-hand text: BiSC weight and heterogeneity, placed past the widest CI
  x_range <- range(c(df_forest$lower, df_forest$upper))
  x_text  <- x_range[2] + 0.05 * diff(x_range)
  df_text <- data.frame(label = factor(label, levels = rev(label)),
                        text  = paste0("BiSC wt ", forest_tab$weight_bisc, "% | I2 ",
                                       forest_tab$i2, " | Q p ", forest_tab$q_p))
  
  # Plot 1. Meta vs cohorts: top meta CpGs (y), effect with 95% CI (x) in the
  # meta-analysis (black diamond), BiSC and SEPAGES (colour); filled = p < cutoff;
  # text = BiSC share of the inverse-variance weight and heterogeneity.
  p_forest <- ggplot(df_forest, aes(x = Effect, y = label, colour = source)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0,
                   position = position_dodge(width = 0.7), linewidth = 0.5) +
    geom_point(aes(shape = source, alpha = suggestive), size = 2.2,
               position = position_dodge(width = 0.7)) +
    geom_text(data = df_text, aes(x = x_text, y = label, label = text),
              inherit.aes = FALSE, hjust = 0, size = 2.8, colour = "grey30") +
    scale_colour_manual(values = pal_source) +
    scale_shape_manual(values = shape_source) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.45),
                       labels = c(`TRUE` = paste("p <", cutoff_sugg), `FALSE` = "not suggestive")) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.6))) +   # room for the text
    labs(title = paste0(exp_label, ": meta-analysis vs cohorts (", meta_id, ")"),
         subtitle = paste0("Top ", length(cpgs), " CpGs of the BiSC + SEPAGES fixed-effects meta-analysis"),
         x = "Effect (95% CI)", y = NULL, colour = NULL, shape = NULL, alpha = NULL) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, shape = unname(shape_source))),
           shape  = "none") +
    theme_ewas() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "bottom")
  
  save_plot(p_forest, file.path(outdir, paste0("forest_meta_", e, ".png")),
            width = 22, height = 5 + 0.9 * length(cpgs))
  cat("Saved forest_meta:", e, "\n")
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----