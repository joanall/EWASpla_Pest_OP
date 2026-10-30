#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the
# BiSC and SEPAGES cohorts.
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: OBJECTIVE 1 - Forest plots of the top suggestive CpGs of the
# main models (effect with 95% CI) compared with the sensitivity analyses
# (S1-S5): one plot per CpG, models on the y-axis. Visual counterpart of
# Tables 1-2.
# - Filled point: the association holds in that model (p < 0.05 and same
#   direction as the main model; same criterion as n_sens_p005 in Tables 1-2).
#   The main model is always filled (all CpGs are suggestive there).
# Label: CpG (gene) - genes de-duplicated; CpG only if intergenic.
# Uses Tables 1-2 from 3a_Get_Suggestive_Summary.R (run it first).
# Argument: pattern for the runs (one cohort, sample and run date), e.g.
#   Rscript script/6_Analyze_Results/4a_Get_Forest_Sensitivity.R "sepages_2tp_t[23]_20261021"
#
# [INPUT]:  results/5_Annotation/files_annotated_index.csv (run selection)
#           results/6_Analyze_Results/summary_tables/<tag>/table1_suggestive_t2.rds,
#           table2_suggestive_t3.rds (with Effect_/SE_/pvalue_ s1-s5)
# [OUTPUT]: results/6_Analyze_Results/summary_tables/<tag>/forest_sens/<exposure>_<tp>/
#           forest_sens_<exposure>_<tp>_<CpG>.png
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_template_plot.R")
source("script/00_functions.R")

index_file  <- "results/5_Annotation/files_annotated_index.csv"
cutoff_nom  <- 0.05    # the association holds in a sensitivity model: p < 0.05, same direction
n_top_plot  <- 20      # CpGs per plot

args <- commandArgs(trailingOnly = TRUE)
run_pattern <- if (length(args) >= 1 && args[1] != "") args[1] else "."
# run_pattern <- "sepages_2tp_t[23]_20261021"   # TEST

files_index <- read.csv(index_file, stringsAsFactors = FALSE)
files_index <- files_index[grepl(run_pattern, files_index$run), ]
stopifnot(length(unique(files_index$cohort)) == 1,
          length(unique(files_index$sample)) == 1,
          length(unique(files_index$run_date)) == 1)

main_index <- files_index[files_index$analysis == "main", ]
tag    <- paste(main_index$cohort[1], main_index$sample[1], main_index$run_date[1], sep = "_")
outdir <- file.path("results/6_Analyze_Results/summary_tables", tag)

exposures  <- unique(main_index$exposure)
sens_names <- sort(unique(files_index$analysis[files_index$analysis != "main"]))
models     <- c("main", sens_names)

pal_models <- c(main = "grey20",
                s1 = pal_bright[["blue"]],  s2 = pal_bright[["cyan"]],
                s3 = pal_bright[["green"]], s4 = pal_bright[["yellow"]],
                s5 = pal_bright[["purple"]])

# Model descriptions (y-axis labels)
model_labels <- c(main = "Main model",
                  s1   = "S1: Excluding pregnancy complications",
                  s2   = "S2: Europeans only",
                  s3   = "S3: Adjusted for specific gravity",
                  s4   = "S4: Adjusted for fruit and vegetable intake",
                  s5   = "S5: Not adjusted for gestational age and cell type")
stopifnot(all(models %in% names(model_labels)))   # every model has a description

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

# Suggestive CpGs of the main models with sensitivity look-up (3a)
file_t1 <- file.path(outdir, "table1_suggestive_t2.rds")
file_t2 <- file.path(outdir, "table2_suggestive_t3.rds")
stopifnot(file.exists(file_t1), file.exists(file_t2))   # run 3a first
top_main <- rbind(readRDS(file_t1), readRDS(file_t2))

# Check: sensitivity Effect, SE and pvalue columns present (should be character(0))
setdiff(c(paste0("Effect_", sens_names), paste0("SE_", sens_names), paste0("pvalue_", sens_names)),
        colnames(top_main))

#==============================================================================#----
# FOREST PLOTS - MAIN vs SENSITIVITY                                           #----
#==============================================================================#----

for (e in exposures) {
  exp_label <- exp_labels[match(e, exp_code)]
  
  for (tp in c("t2", "t3")) {
    
    tab <- top_main[top_main$exposure == e & top_main$time_point == tp, ]
    if (nrow(tab) == 0) next
    tab <- head(tab[order(tab$pvalue), ], n_top_plot)
    tab$label <- make_label(tab$CpG, tab$gene_name)
    
    # One folder per exposure x time point, one plot per CpG
    plot_dir <- file.path(outdir, "forest_sens", paste0(e, "_", tp))
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
      # Main model on top, then S1-S5
      df_fp$model_label <- factor(model_labels[df_fp$model], levels = rev(model_labels[models]))
      
      # Check: sensitivity models where the association holds (should equal n_sens_p005)
      n_holds <- sum(df_fp$holds[df_fp$model != "main"])
      cat("  ", cpg_row$CpG, "| holds in", n_holds, "of", length(sens_names),
          "sensitivity models (n_sens_p005 =", cpg_row$n_sens_p005, ")\n")
      df_fp$holds <- factor(df_fp$holds, levels = c(TRUE, FALSE))
      
      # Plot 1. One suggestive CpG: effect with 95% CI (x) in the main model and
      # each sensitivity analysis (y); filled = p < 0.05 and same direction as
      # the main model.
      p_fp <- ggplot(df_fp, aes(x = Effect, y = model_label, colour = model, shape = holds)) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
        geom_pointrange(aes(xmin = lower, xmax = upper), size = 0.4) +
        scale_colour_manual(values = pal_models, guide = "none") +
        scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), drop = FALSE,
                           labels = c(`TRUE`  = paste("p <", cutoff_nom, "& same direction"),
                                      `FALSE` = "not replicated")) +
        labs(title = paste0(exp_label, " - ", tp, ": ", cpg_row$label),
             subtitle = paste0(tag, ", main model p = ", signif(cpg_row$pvalue, 2)),
             x = "Effect (95% CI)", y = NULL, shape = NULL) +
        theme_ewas() +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
              legend.position = "bottom")
      
      save_plot(p_fp, file.path(plot_dir, paste0("forest_sens_", e, "_", tp, "_", cpg_row$CpG, ".png")),
                width = 18, height = 8)
    }
    cat("Saved forest_sens:", e, tp, "-", nrow(tab), "CpGs\n")
  }
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----