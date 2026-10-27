#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: Forest plots of the top suggestive CpGs (effect with 95% CI):
# - Main vs sensitivity analyses (s1-s5), per exposure and time point
# - t2 vs t3, per exposure (CpGs suggestive at either time point)
# Label: first gene - CpG (CpG only if intergenic).
# Argument: pattern for the runs (one cohort, sample and run date), e.g.
#   Rscript script/6_Analyze_Results/6_Get_Forest_Plots.R "sepages_2tp_t[23]_20261021"
#
# [INPUT]:  results/6_Analyze_Results/summary_tables/<tag>/table1_suggestive_t2.rds,
#           table2_suggestive_t3.rds (with Effect_/SE_/pvalue_ s1-s5)
#           results/5_Annotation/files_annotated_index.csv (main results, t2 vs t3)
# [OUTPUT]: results/6_Analyze_Results/summary_tables/<tag>/forest_sens_<exposure>_<tp>.png
#           results/6_Analyze_Results/summary_tables/<tag>/forest_t2_t3_<exposure>.png
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_template_plot.R")
source("script/00_functions.R")

index_file  <- "results/5_Annotation/files_annotated_index.csv"
cutoff_sugg <- 1e-5
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

# Label: first gene - CpG (CpG only if no gene)
make_label <- function(cpg, gene) {
  gene1 <- sub(";.*", "", gene)
  ifelse(is.na(gene1) | gene1 == "", cpg, paste0(gene1, " - ", cpg))
}

# Suggestive CpGs of the main models with sensitivity look-up (summary script)
top_main <- rbind(readRDS(file.path(outdir, "table1_suggestive_t2.rds")),
                  readRDS(file.path(outdir, "table2_suggestive_t3.rds")))

# Check: sensitivity effect and SE columns present
setdiff(c(paste0("Effect_", sens_names), paste0("SE_", sens_names)), colnames(top_main))

#==============================================================================#----
# 1. FOREST PLOTS - MAIN vs SENSITIVITY                                        #----
#==============================================================================#----

for (e in exposures) {
  for (tp in c("t2", "t3")) {
    
    tab <- top_main[top_main$exposure == e & top_main$time_point == tp, ]
    if (nrow(tab) == 0) next
    tab <- head(tab[order(tab$pvalue), ], n_top_plot)
    tab$label <- make_label(tab$CpG, tab$gene_name)
    
    df_fp <- NULL
    for (mdl in models) {
      eff <- if (mdl == "main") tab$Effect else tab[[paste0("Effect_", mdl)]]
      se  <- if (mdl == "main") tab$SE     else tab[[paste0("SE_", mdl)]]
      df_fp <- rbind(df_fp, data.frame(label = tab$label, model = mdl, Effect = eff,
                                       lower = eff - 1.96 * se, upper = eff + 1.96 * se))
    }
    df_fp$label <- factor(df_fp$label, levels = rev(tab$label))
    df_fp$model <- factor(df_fp$model, levels = models)
    
    # Plot 1. Main vs sensitivity analyses: one row per suggestive CpG (y),
    # effect with 95% CI (x), colour = model.
    p_fp <- ggplot(df_fp, aes(x = Effect, y = label, colour = model)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_pointrange(aes(xmin = lower, xmax = upper),
                      position = position_dodge(width = 0.7), size = 0.2) +
      scale_colour_manual(values = pal_models) +
      labs(title = paste0(e, " - ", tp, ": main vs sensitivity analyses"),
           subtitle = paste0(tag, ", top ", nrow(tab), " suggestive CpGs (main model)"),
           x = "Effect (95% CI)", y = NULL, colour = "Model") +
      theme_ewas() +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
    
    save_plot(p_fp, file.path(outdir, paste0("forest_sens_", e, "_", tp, ".png")),
              width = 18, height = 4 + 0.7 * nrow(tab))
    cat("Saved forest_sens:", e, tp, "\n")
  }
}

#==============================================================================#----
# 2. FOREST PLOTS - T2 vs T3                                                   #----
#==============================================================================#----

for (e in exposures) {
  
  cpgs <- unique(top_main$CpG[top_main$exposure == e])
  if (length(cpgs) == 0) next
  
  # Main results of this exposure at both time points
  res_t2 <- readRDS(main_index$path[main_index$exposure == e & main_index$time_point == "t2"])
  res_t3 <- readRDS(main_index$path[main_index$exposure == e & main_index$time_point == "t3"])
  m2 <- match(cpgs, res_t2$CpG)
  m3 <- match(cpgs, res_t3$CpG)
  
  # Top CpGs by their best p-value across time points
  best_p <- pmin(res_t2$pvalue[m2], res_t3$pvalue[m3], na.rm = TRUE)
  ord    <- head(order(best_p), n_top_plot)
  cpgs <- cpgs[ord]; m2 <- m2[ord]; m3 <- m3[ord]
  label <- make_label(cpgs, res_t2$gene_name[m2])
  
  df_tp <- rbind(
    data.frame(label = label, time = "t2", Effect = res_t2$Effect[m2], se = res_t2$SE[m2],
               suggestive = res_t2$pvalue[m2] < cutoff_sugg),
    data.frame(label = label, time = "t3", Effect = res_t3$Effect[m3], se = res_t3$SE[m3],
               suggestive = res_t3$pvalue[m3] < cutoff_sugg))
  df_tp$lower <- df_tp$Effect - 1.96 * df_tp$se
  df_tp$upper <- df_tp$Effect + 1.96 * df_tp$se
  df_tp$label <- factor(df_tp$label, levels = rev(label))
  
  # Plot 2. t2 vs t3: CpGs suggestive at either time point (y), effect with
  # 95% CI at each time point (x, colour); filled = suggestive at that time point.
  p_tp <- ggplot(df_tp, aes(x = Effect, y = label, colour = time, shape = suggestive)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_pointrange(aes(xmin = lower, xmax = upper),
                    position = position_dodge(width = 0.6), size = 0.3) +
    scale_colour_manual(values = pal_time) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                       labels = c(`TRUE` = paste("p <", cutoff_sugg), `FALSE` = "not suggestive")) +
    labs(title = paste0(e, ": t2 vs t3"),
         subtitle = paste0(tag, ", CpGs suggestive at t2 or t3 (top ", length(cpgs), ")"),
         x = "Effect (95% CI)", y = NULL, colour = "Time point", shape = NULL) +
    theme_ewas() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  save_plot(p_tp, file.path(outdir, paste0("forest_t2_t3_", e, ".png")),
            width = 18, height = 4 + 0.6 * length(cpgs))
  cat("Saved forest_t2_t3:", e, "\n")
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----