#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the
# BiSC and SEPAGES cohorts.
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: OBJECTIVE 2 - Forest plots comparing t2 and t3 (main models),
# per exposure: CpGs suggestive at either time point, effect with 95% CI at
# each time point. Visual counterpart of Table 3.
# - Filled point: suggestive (p < 1e-5) at that time point.
# Label: CpG (gene) - genes de-duplicated; CpG only if intergenic.
# Uses Tables 1-2 from 3a_Get_Suggestive_Summary.R (run it first).
# Argument: pattern for the runs (one cohort, sample and run date), e.g.
#   Rscript script/6_Analyze_Results/4b_Get_Forest_Time_Points.R "sepages_2tp_t[23]_20261021"
#
# [INPUT]:  results/5_Annotation/files_annotated_index.csv
#           -> results/5_Annotation/<run>/main/<exposure>_<tp>_annotated.rds
#           results/6_Analyze_Results/summary_tables/<tag>/table1_suggestive_t2.rds,
#           table2_suggestive_t3.rds
# [OUTPUT]: results/6_Analyze_Results/summary_tables/<tag>/forest_t2_t3_<exposure>.png
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_template_plot.R")
source("script/00_functions.R")

index_file  <- "results/5_Annotation/files_annotated_index.csv"
cutoff_sugg <- 1e-5    # suggestive threshold (filled points)
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
stopifnot(all(file.exists(main_index$path)))
table(main_index$exposure, main_index$time_point)   # 1 per cell

tag    <- paste(main_index$cohort[1], main_index$sample[1], main_index$run_date[1], sep = "_")
outdir <- file.path("results/6_Analyze_Results/summary_tables", tag)

exposures <- unique(main_index$exposure)

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

# Suggestive CpGs of the main models (3a)
file_t1 <- file.path(outdir, "table1_suggestive_t2.rds")
file_t2 <- file.path(outdir, "table2_suggestive_t3.rds")
stopifnot(file.exists(file_t1), file.exists(file_t2))   # run 3a first
top_main <- rbind(readRDS(file_t1), readRDS(file_t2))
table(top_main$exposure, top_main$time_point)           # n suggestive per model

#==============================================================================#----
# FOREST PLOTS - T2 vs T3                                                      #----
#==============================================================================#----

for (e in exposures) {
  
  exp_label <- exp_labels[match(e, exp_code)]
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
  
  # Gene from t2, or t3 if the CpG is missing at t2
  gene  <- ifelse(is.na(m2), res_t3$gene_name[m3], res_t2$gene_name[m2])
  label <- make_label(cpgs, gene)
  
  df_tp <- rbind(
    data.frame(label = label, time = "t2", Effect = res_t2$Effect[m2], se = res_t2$SE[m2],
               suggestive = !is.na(res_t2$pvalue[m2]) & res_t2$pvalue[m2] < cutoff_sugg),
    data.frame(label = label, time = "t3", Effect = res_t3$Effect[m3], se = res_t3$SE[m3],
               suggestive = !is.na(res_t3$pvalue[m3]) & res_t3$pvalue[m3] < cutoff_sugg))
  df_tp$lower <- df_tp$Effect - 1.96 * df_tp$se
  df_tp$upper <- df_tp$Effect + 1.96 * df_tp$se
  df_tp$label <- factor(df_tp$label, levels = rev(label))
  
  # Check: suggestive CpGs per time point among the plotted ones
  print(table(df_tp$time, df_tp$suggestive))
  
  # Plot 1. t2 vs t3: CpGs suggestive at either time point (y), effect with
  # 95% CI at each time point (x, colour); filled = suggestive at that time point.
  p_tp <- ggplot(df_tp, aes(x = Effect, y = label, colour = time, shape = suggestive)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_pointrange(aes(xmin = lower, xmax = upper),
                    position = position_dodge(width = 0.6), size = 0.3) +
    scale_colour_manual(values = pal_time) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                       labels = c(`TRUE` = paste("p <", cutoff_sugg), `FALSE` = "not suggestive")) +
    labs(title = paste0(exp_label, ": t2 vs t3"),
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