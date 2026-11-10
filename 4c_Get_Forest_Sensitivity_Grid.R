#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the
# BiSC and SEPAGES cohorts.
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: OBJECTIVE 1 - Overview grid of forest plots, one figure per
# time point: every suggestive CpG of the main models (all exposures) in its
# own small panel, effect with 95% CI in the main model and S1-S5.
# Same content as 4a_Get_Forest_Sensitivity.R, but all CpGs at once.
# - Filled point: p < 0.05 and same direction as the main model (as
#   n_sens_p005 in Tables 1-2). The main model is always filled.
# - Panels ordered by exposure, then main-model p-value; free x-axis per panel.
# - Short model names on the y-axis (S1-S5), described in the caption.
# Uses Tables 1-2 from 3a_Get_Suggestive_Summary.R (run it first).
# Argument: pattern for the runs (one cohort, sample and run date), e.g.
#   Rscript script/6_Analyze_Results/4c_Get_Forest_Sensitivity_Grid.R "sepages_2tp_t[23]_20261021"
#
# [INPUT]:  results/5_Annotation/files_annotated_index.csv (run selection)
#           results/6_Analyze_Results/summary_tables/<tag>/table1_suggestive_t2.rds,
#           table2_suggestive_t3.rds (with Effect_/SE_/pvalue_ s1-s5)
# [OUTPUT]: results/6_Analyze_Results/summary_tables/<tag>/forest_sens_grid_<tp>.png
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_template_plot.R")
source("script/00_functions.R")

index_file  <- "results/5_Annotation/files_annotated_index.csv"
cutoff_nom  <- 0.05    # the association holds in a sensitivity model: p < 0.05, same direction
n_top_plot  <- 20      # max CpGs per exposure (as in 4a)
n_col_grid  <- 6       # panels per row

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

sens_names <- sort(unique(files_index$analysis[files_index$analysis != "main"]))
models     <- c("main", sens_names)

pal_models <- c(main = "grey20",
                s1 = pal_bright[["blue"]],  s2 = pal_bright[["cyan"]],
                s3 = pal_bright[["green"]], s4 = pal_bright[["yellow"]],
                s5 = pal_bright[["purple"]])

# Short names (y-axis) and descriptions (caption)
model_short  <- c(main = "Main", s1 = "S1", s2 = "S2", s3 = "S3", s4 = "S4", s5 = "S5")
model_labels <- c(s1 = "S1: Excluding pregnancy complications",
                  s2 = "S2: Europeans only",
                  s3 = "S3: Adjusted for specific gravity",
                  s4 = "S4: Adjusted for fruit and vegetable intake",
                  s5 = "S5: Not adjusted for gestational age and cell type")
stopifnot(all(models %in% names(model_short)),
          all(sens_names %in% names(model_labels)))   # every model has a name and description

# Exposure labels and order
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
# FOREST PLOT GRID - ALL SUGGESTIVE CpGs PER TIME POINT                        #----
#==============================================================================#----

for (tp in c("t2", "t3")) {
  
  tab <- top_main[top_main$time_point == tp, ]
  if (nrow(tab) == 0) next
  
  # Order: exposure (fixed order), then main-model p-value; top n per exposure
  tab <- tab[order(match(tab$exposure, exp_code), tab$pvalue), ]
  tab <- do.call(rbind, lapply(split(tab, factor(tab$exposure, levels = unique(tab$exposure))),
                               function(x) head(x, n_top_plot)))
  tab$panel <- paste0(exp_labels[match(tab$exposure, exp_code)], ": ",
                      make_label(tab$CpG, tab$gene_name))
  table(tab$exposure)   # CpGs per exposure in the grid
  
  # Long format: one row per CpG x model
  df_grid <- NULL
  for (mdl in models) {
    if (mdl == "main") {
      eff   <- tab$Effect
      se    <- tab$SE
      holds <- rep(TRUE, nrow(tab))
    } else {
      eff   <- tab[[paste0("Effect_", mdl)]]
      se    <- tab[[paste0("SE_", mdl)]]
      p_mdl <- tab[[paste0("pvalue_", mdl)]]
      holds <- !is.na(p_mdl) & p_mdl < cutoff_nom & sign(eff) == sign(tab$Effect)
    }
    df_grid <- rbind(df_grid, data.frame(panel = tab$panel, model = mdl, Effect = eff,
                                         lower = eff - 1.96 * se, upper = eff + 1.96 * se,
                                         holds = holds))
  }
  df_grid$panel       <- factor(df_grid$panel, levels = tab$panel)
  df_grid$model_short <- factor(model_short[df_grid$model], levels = rev(model_short[models]))
  df_grid$model       <- factor(df_grid$model, levels = models)
  
  # Check: sensitivity models where the association holds, per CpG (= n_sens_p005)
  n_holds <- tapply(df_grid$holds[df_grid$model != "main"],
                    df_grid$panel[df_grid$model != "main"], sum)
  stopifnot(all(n_holds[tab$panel] == tab$n_sens_p005))
  df_grid$holds <- factor(df_grid$holds, levels = c(TRUE, FALSE))
  
  n_rows_grid <- ceiling(nrow(tab) / n_col_grid)
  cat(tp, ":", nrow(tab), "CpGs in", n_rows_grid, "rows\n")
  
  # Plot 1. Grid of forest plots, one panel per suggestive CpG at this time
  # point (all exposures): effect with 95% CI (x) in the main model and S1-S5
  # (y); filled = p < 0.05 and same direction as the main model.
  p_grid <- ggplot(df_grid, aes(x = Effect, y = model_short, colour = model, shape = holds)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_pointrange(aes(xmin = lower, xmax = upper), size = 0.25) +
    facet_wrap(~ panel, ncol = n_col_grid, scales = "free_x",
               labeller = label_wrap_gen(width = 28)) +
    scale_colour_manual(values = pal_models, guide = "none") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), drop = FALSE,
                       labels = c(`TRUE`  = paste("p <", cutoff_nom, "& same direction"),
                                  `FALSE` = "not replicated")) +
    labs(title = paste0("Suggestive CpGs at ", tp, ": main vs sensitivity analyses"),
         subtitle = paste0(tag, ", ", nrow(tab), " suggestive CpGs (main model, p < 1e-5)"),
         caption = paste(model_labels[sens_names], collapse = "\n"),
         x = "Effect (95% CI)", y = NULL, shape = NULL) +
    theme_ewas() +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 6),
          strip.text   = element_text(size = 7),
          plot.caption = element_text(hjust = 0, size = 8),
          legend.position = "bottom")
  
  save_plot(p_grid, file.path(outdir, paste0("forest_sens_grid_", tp, ".png")),
            width = 40, height = 6 + 6 * n_rows_grid)
  cat("Saved forest_sens_grid:", tp, "\n")
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----