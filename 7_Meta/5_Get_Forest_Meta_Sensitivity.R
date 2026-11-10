#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the
# BiSC and SEPAGES cohorts.
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org
# Date: November 2026
#==============================================================================#
# [DESCRIPTION]: Meta-analysis counterpart of 6_Analyze_DMPs/4a_Get_Forest_Sensitivity.R.
# Forest plots of the top suggestive CpGs of the main BiSC + SEPAGES
# meta-analysis (effect with 95% CI) compared with the meta-analyses of the
# sensitivity analyses (S1-S6): one plot per CpG, models on the y-axis.
# - Top CpGs: meta main p < cutoff_sugg (max n_top_plot per exposure x tp)
# - Filled point: the association holds in that model (p < 0.05 and same
#   direction as the main meta-analysis; same criterion as n_sens_p005 in 3a).
#   The main model is always filled (all CpGs are suggestive there).
# - Sensitivity meta-analyses not run (no <sX>/ folder) are skipped;
#   a CpG missing in one sensitivity meta-analysis is left empty in that row.
# Label: CpG (gene) - genes de-duplicated; CpG only if intergenic.
#
# Arguments:
#   1) run date of the meta-analysis folders (default: 20261106)
#   2) time points, comma-separated (default: t2,t3,avg)
#   3) exposures, comma-separated (default: all in the main meta-analysis)
# Example:
#   Rscript script/7_Meta/5_Get_Forest_Meta_Sensitivity.R 20261106 t2,t3,avg tcpy
#
# [INPUT]:  results/7_Meta/3_Read_GWAMA/meta_tp_<tp>_<date>/<analysis>/<exposure>_<tp>.rds
#           results/5_Annotation/cpg_annotation_all.rds (gene names)
# [OUTPUT]: results/7_Meta/5_Forest_Meta_Sens/<date>/
#           table_meta_suggestive_<tp>.xlsx/.rds (look-up, as Tables 1-2 / 7)
#           <exposure>_<tp>/forest_meta_sens_<exposure>_<tp>_<CpG>.png
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_template_plot.R")
library(writexl)

dir_read    <- "results/7_Meta/3_Read_GWAMA"
annot_file  <- "results/5_Annotation/cpg_annotation_all.rds"
outdir      <- "results/7_Meta/5_Forest_Meta_Sens"
cutoff_sugg <- 1e-5    # top CpGs in the main meta-analysis
cutoff_nom  <- 0.05    # the association holds in a sensitivity model: p < 0.05, same direction
n_top_plot  <- 20      # max CpGs per exposure x time point
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
run_date <- if (length(args) >= 1 && args[1] != "") args[1] else "20261106"
tps_arg  <- if (length(args) >= 2 && args[2] != "") strsplit(args[2], ",")[[1]] else c("t2", "t3", "avg")
exp_arg  <- if (length(args) >= 3 && args[3] != "") strsplit(args[3], ",")[[1]] else NULL
outdir   <- file.path(outdir, run_date)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
stopifnot(all(tps_arg %in% c("t2", "t3", "avg")))

sens_all <- paste0("s", 1:6)

pal_models <- c(main = "grey20",
                s1 = pal_bright[["blue"]],  s2 = pal_bright[["cyan"]],
                s3 = pal_bright[["green"]], s4 = pal_bright[["yellow"]],
                s5 = pal_bright[["purple"]], s6 = pal_bright[["red"]])

# Model descriptions (y-axis labels)
model_labels <- c(main = "Main model",
                  s1   = "S1: Excluding pregnancy complications",
                  s2   = "S2: Europeans only",
                  s3   = "S3: Adjusted for specific gravity",
                  s4   = "S4: Adjusted for fruit and vegetable intake",
                  s5   = "S5: Not adjusted for gestational age and cell type",
                  s6   = "S6: Exposure outliers winsorized (P1-P99)")

# Exposure labels for titles
exp_code   <- c("24d",   "3pba",  "tcpy", "dmp",  "dep",  "dap")
exp_labels <- c("2,4-D", "3-PBA", "TCPy", "ΣDMP", "ΣDEP", "ΣDAP")

# Label: "CpG (gene)", genes de-duplicated ("A;A;B" -> "A;B"); CpG only if no gene
make_label <- function(cpg, gene) {
  gene <- as.character(gene)
  gene[is.na(gene)] <- ""
  gene <- sapply(strsplit(gene, ";"), function(g) paste(unique(g[g != ""]), collapse = ";"))
  ifelse(gene == "", cpg, paste0(cpg, " (", gene, ")"))
}

# Gene names
annotation <- readRDS(annot_file)
stopifnot(all(c("CpG", "gene_name") %in% names(annotation)))

#==============================================================================#----
# FOREST PLOTS - MAIN vs SENSITIVITY META-ANALYSES                             #----
#==============================================================================#----

for (tp in tps_arg) {
  
  # Folder of this time point: meta_tp_<tp>_<date>/<analysis>/<exposure>_<tp>.rds
  dir_tp <- file.path(dir_read, paste0("meta_tp_", tp, "_", run_date))
  if (!dir.exists(dir_tp)) { cat("No folder", dir_tp, "- skipped\n"); next }
  
  # Reads all <exposure>_<tp>.rds of one analysis into a list named by exposure
  read_meta_analysis <- function(analysis) {
    files_a <- list.files(file.path(dir_tp, analysis), pattern = paste0("_", tp, "\\.rds$"), full.names = TRUE)
    if (length(files_a) == 0) return(NULL)
    meta_a <- lapply(files_a, readRDS)
    names(meta_a) <- sub(paste0("_", tp, "\\.rds$"), "", basename(files_a))
    meta_a
  }
  
  # Main meta-analysis of this time point
  meta_main <- read_meta_analysis("main")
  if (is.null(meta_main)) { cat("No main meta-analysis in", dir_tp, "- skipped\n"); next }
  cat("\n==", tp, "| main exposures:", names(meta_main), "\n")
  
  # Sensitivity meta-analyses available for this time point
  meta_sens <- list()
  for (s in sens_all) {
    meta_s <- read_meta_analysis(s)
    if (!is.null(meta_s)) meta_sens[[s]] <- meta_s
  }
  sens_names <- names(meta_sens)
  cat("   sensitivity meta-analyses found:", sens_names, "\n")
  sens_missing <- setdiff(sens_all, sens_names)
  if (length(sens_missing) > 0) cat("   not found (skipped):", sens_missing, "\n")
  models <- c("main", sens_names)
  
  exposures <- if (is.null(exp_arg)) names(meta_main) else intersect(exp_arg, names(meta_main))
  stopifnot(length(exposures) > 0)
  stopifnot(all(c("CpG", "Effect", "SE", "pvalue") %in% colnames(meta_main[[exposures[1]]])))
  meta_main[[exposures[1]]][1:2, c("CpG", "Effect", "SE", "pvalue")]   # check: columns read
  
  # 1. Look-up table: suggestive CpGs of the main meta-analysis + sensitivity estimates
  top_meta <- NULL
  for (e in exposures) {
    res_e <- meta_main[[e]]
    res_e <- res_e[!is.na(res_e$pvalue), ]
    res_e <- res_e[order(res_e$pvalue), ]
    top_e <- head(res_e[res_e$pvalue < cutoff_sugg, ], n_top_plot)
    cat(e, "| suggestive:", sum(res_e$pvalue < cutoff_sugg), "| plotted:", nrow(top_e), "\n")
    if (nrow(top_e) == 0) next
    
    tab_e <- data.frame(exposure = e, time_point = tp, CpG = top_e$CpG,
                        Effect = top_e$Effect, SE = top_e$SE, pvalue = top_e$pvalue,
                        stringsAsFactors = FALSE)
    for (s in sens_names) {
      res_s <- meta_sens[[s]][[e]]
      if (is.null(res_s)) {
        tab_e[[paste0("Effect_", s)]] <- NA
        tab_e[[paste0("SE_", s)]]     <- NA
        tab_e[[paste0("pvalue_", s)]] <- NA
        next
      }
      m <- match(tab_e$CpG, res_s$CpG)
      tab_e[[paste0("Effect_", s)]] <- res_s$Effect[m]
      tab_e[[paste0("SE_", s)]]     <- res_s$SE[m]
      tab_e[[paste0("pvalue_", s)]] <- res_s$pvalue[m]
    }
    cat("   look-up columns:", ncol(tab_e), "| CpGs:", nrow(tab_e), "\n")   # check: 6 + 3 x n sens
    top_meta <- rbind(top_meta, tab_e)
  }
  if (is.null(top_meta)) { cat("No suggestive CpGs at", tp, "\n"); next }
  stopifnot(is.data.frame(top_meta))
  dim(top_meta)   # check: CpGs x columns
  
  top_meta$gene_name <- annotation$gene_name[match(top_meta$CpG, annotation$CpG)]
  top_meta$label     <- make_label(top_meta$CpG, top_meta$gene_name)
  
  # Agreement with the main meta-analysis: same direction, and p < 0.05
  if (length(sens_names) > 0) {
    eff_sens <- as.matrix(top_meta[, paste0("Effect_", sens_names), drop = FALSE])
    p_sens   <- as.matrix(top_meta[, paste0("pvalue_", sens_names), drop = FALSE])
    top_meta$n_sens_same_dir <- rowSums(sign(eff_sens) == sign(top_meta$Effect), na.rm = TRUE)
    top_meta$n_sens_p005     <- rowSums(p_sens < cutoff_nom & sign(eff_sens) == sign(top_meta$Effect), na.rm = TRUE)
  }
  table(top_meta$exposure)   # CpGs per exposure
  
  file_base <- file.path(outdir, paste0("table_meta_suggestive_", tp))
  write_xlsx(top_meta, paste0(file_base, ".xlsx"))
  saveRDS(top_meta, paste0(file_base, ".rds"))
  
  # 2. One forest plot per CpG
  for (e in unique(top_meta$exposure)) {
    exp_label <- exp_labels[match(e, exp_code)]
    tab <- top_meta[top_meta$exposure == e, ]
    
    plot_dir <- file.path(outdir, paste0(e, "_", tp))
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    
    for (k in seq_len(nrow(tab))) {
      cpg_row <- tab[k, ]
      
      df_fp <- NULL
      for (mdl in models) {
        if (mdl == "main") {
          eff   <- cpg_row$Effect
          se    <- cpg_row$SE
          holds <- TRUE
        } else {
          eff   <- cpg_row[[paste0("Effect_", mdl)]]
          se    <- cpg_row[[paste0("SE_", mdl)]]
          p_mdl <- cpg_row[[paste0("pvalue_", mdl)]]
          holds <- !is.na(p_mdl) & p_mdl < cutoff_nom & sign(eff) == sign(cpg_row$Effect)
        }
        df_fp <- rbind(df_fp, data.frame(model = mdl, Effect = eff,
                                         lower = eff - 1.96 * se, upper = eff + 1.96 * se,
                                         holds = holds))
      }
      # Main model on top, then S1-S6
      df_fp$model_label <- factor(model_labels[df_fp$model], levels = rev(model_labels[models]))
      
      # Check: sensitivity models where the association holds (should equal n_sens_p005)
      n_holds <- sum(df_fp$holds[df_fp$model != "main"])
      cat("  ", cpg_row$CpG, "| holds in", n_holds, "of", length(sens_names),
          "sensitivity meta-analyses (n_sens_p005 =", cpg_row$n_sens_p005, ")\n")
      df_fp$holds <- factor(df_fp$holds, levels = c(TRUE, FALSE))
      df_fp_plot <- df_fp[!is.na(df_fp$Effect), ]   # CpG missing in a sensitivity meta: empty row
      
      # Plot 1. One suggestive CpG of the main meta-analysis: effect with 95% CI
      # (x) in the main and each sensitivity meta-analysis (y); filled = p < 0.05
      # and same direction as the main meta-analysis.
      p_fp <- ggplot(df_fp_plot, aes(x = Effect, y = model_label, colour = model, shape = holds)) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
        geom_pointrange(aes(xmin = lower, xmax = upper), size = 0.4) +
        scale_y_discrete(drop = FALSE) +
        scale_colour_manual(values = pal_models, guide = "none") +
        scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), drop = FALSE,
                           labels = c(`TRUE`  = paste("p <", cutoff_nom, "& same direction"),
                                      `FALSE` = "not replicated")) +
        labs(title = paste0(exp_label, " - ", tp, ": ", cpg_row$label),
             subtitle = paste0("BiSC + SEPAGES meta-analysis, main model p = ", signif(cpg_row$pvalue, 2)),
             x = "Effect (95% CI)", y = NULL, shape = NULL) +
        theme_ewas() +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
              legend.position = "bottom")
      
      save_plot(p_fp, file.path(plot_dir, paste0("forest_meta_sens_", e, "_", tp, "_", cpg_row$CpG, ".png")),
                width = 18, height = 3 + 1.2 * length(models))
    }
    cat("Saved forest_meta_sens:", e, tp, "-", nrow(tab), "CpGs\n")
  }
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----