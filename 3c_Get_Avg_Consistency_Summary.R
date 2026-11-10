#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the
# BiSC and SEPAGES cohorts.
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: OBJECTIVE 2 - Consistency between the pregnancy-average model
# and the trimester-specific models (t2, t3), main models only, per exposure.
# Question: is each average signal consistent across trimesters, or driven by
# one of them?
# - Table 5: CpGs suggestive (p < 1e-5) in the average model, with their
#            effect and p at t2 and t3; holds = p < 0.05 and same direction as
#            the average; driven_by = both / t2 only / t3 only / neither.
# - Table 6: summary per exposure: n suggestive per time point, overlaps,
#            avg CpGs holding at t2 / t3, t2 and t3 CpGs holding in the avg,
#            correlation of effects (all CpGs) avg vs t2 and avg vs t3.
# - Plot:    forest plot per exposure, CpGs suggestive at any time point,
#            effect with 95% CI at t2, t3 and avg; filled = suggestive there.
# Note: the average is the mean of the log2 t2 and t3 values, so it is not
# independent of the trimester models; overlap is expected and is not a
# replication. In BiSC, use the 2tp sample so all models have the same women.
# Self-contained: suggestive CpGs are taken from the annotated main results
# (does not need 3a).
# Argument: pattern for the runs (one cohort, sample and run date, with t2, t3
# and avg), e.g.
#   Rscript script/6_Analyze_Results/3c_Get_Avg_Consistency_Summary.R "sepages_2tp_(t2|t3|avg)_20261021"
#
# [INPUT]:  results/5_Annotation/files_annotated_index.csv
#           -> results/5_Annotation/<run>/main/<exposure>_<tp>_annotated.rds
# [OUTPUT]: results/6_Analyze_Results/summary_tables/<cohort>_<sample>_<date>/
#           table5_avg_vs_trimesters.xlsx/.rds, table6_avg_summary.xlsx,
#           forest_t2_t3_avg_<exposure>.png
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_template_plot.R")
source("script/00_functions.R")
library(writexl)   # .xlsx output: no separator/decimal issues when opening in Excel

index_file  <- "results/5_Annotation/files_annotated_index.csv"
cutoff_sugg <- 1e-5    # suggestive threshold
cutoff_nom  <- 0.05    # the association holds at the other time point: p < 0.05, same direction
n_top_plot  <- 20      # CpGs per forest plot

args <- commandArgs(trailingOnly = TRUE)
run_pattern <- if (length(args) >= 1 && args[1] != "") args[1] else "."
run_pattern <- "sepages_2tp_(t2|t3|avg)_20261021"   # TEST

files_index <- read.csv(index_file, stringsAsFactors = FALSE)
files_index <- files_index[grepl(run_pattern, files_index$run), ]
stopifnot(all(file.exists(files_index$path)))

# Check: one cohort, sample and run date (comparisons must not mix runs)
stopifnot(length(unique(files_index$cohort)) == 1,
          length(unique(files_index$sample)) == 1,
          length(unique(files_index$run_date)) == 1)

main_index <- files_index[files_index$analysis == "main", ]
time_points <- c("t2", "t3", "avg")

# Check: the three time points present, one file per exposure x time point
table(main_index$exposure, main_index$time_point)
stopifnot(all(time_points %in% main_index$time_point))

tag    <- paste(main_index$cohort[1], main_index$sample[1], main_index$run_date[1], sep = "_")
outdir <- file.path("results/6_Analyze_Results/summary_tables", tag)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Exposure labels and order
exp_code   <- c("24d",   "3pba",  "tcpy", "dmp",  "dep",  "dap")
exp_labels <- c("2,4-D", "3-PBA", "TCPy", "ΣDMP", "ΣDEP", "ΣDAP")
exposures  <- exp_code[exp_code %in% main_index$exposure]

# Colours per time point (avg not in pal_time)
pal_tp <- c(pal_time[c("t2", "t3")], avg = "grey20")

# Label: "CpG (gene)", genes de-duplicated ("A;A;B" -> "A;B"); CpG only if no gene
make_label <- function(cpg, gene) {
  gene <- as.character(gene)
  gene[is.na(gene)] <- ""
  gene <- sapply(strsplit(gene, ";"), function(g) paste(unique(g[g != ""]), collapse = ";"))
  ifelse(gene == "", cpg, paste0(cpg, " (", gene, ")"))
}

#==============================================================================#----
# 1. MAIN RESULTS AT t2, t3 AND avg                                            #----
#==============================================================================#----

main_res <- list()
for (i in seq_len(nrow(main_index))) {
  f   <- main_index[i, ]
  key <- paste(f$exposure, f$time_point, sep = "_")
  res <- readRDS(f$path)
  main_res[[key]] <- res[, intersect(c("CpG", "gene_name", "Effect", "SE", "pvalue"), colnames(res))]
  cat(key, "| CpGs:", nrow(res), "| suggestive:", sum(res$pvalue < cutoff_sugg, na.rm = TRUE), "\n")
}

#==============================================================================#----
# TABLE 5: AVG SUGGESTIVE CpGs WITH LOOK-UP AT t2 AND t3                       #----
#==============================================================================#----

table5 <- NULL
for (e in exposures) {
  res_avg <- main_res[[paste(e, "avg", sep = "_")]]
  res_t2  <- main_res[[paste(e, "t2",  sep = "_")]]
  res_t3  <- main_res[[paste(e, "t3",  sep = "_")]]
  
  sugg_avg <- res_avg[!is.na(res_avg$pvalue) & res_avg$pvalue < cutoff_sugg, ]
  cat("\n", e, "| suggestive in avg:", nrow(sugg_avg), "\n")
  if (nrow(sugg_avg) == 0) next
  
  m_t2 <- match(sugg_avg$CpG, res_t2$CpG)
  m_t3 <- match(sugg_avg$CpG, res_t3$CpG)
  
  tab <- data.frame(exposure   = e,
                    CpG        = sugg_avg$CpG,
                    gene_name  = sugg_avg$gene_name,
                    Effect_avg = sugg_avg$Effect,
                    SE_avg     = sugg_avg$SE,
                    pvalue_avg = sugg_avg$pvalue,
                    Effect_t2  = res_t2$Effect[m_t2],
                    pvalue_t2  = res_t2$pvalue[m_t2],
                    Effect_t3  = res_t3$Effect[m_t3],
                    pvalue_t3  = res_t3$pvalue[m_t3],
                    stringsAsFactors = FALSE)
  
  # Holds at the trimester: p < 0.05 and same direction as the average
  tab$holds_t2 <- !is.na(tab$pvalue_t2) & tab$pvalue_t2 < cutoff_nom & sign(tab$Effect_t2) == sign(tab$Effect_avg)
  tab$holds_t3 <- !is.na(tab$pvalue_t3) & tab$pvalue_t3 < cutoff_nom & sign(tab$Effect_t3) == sign(tab$Effect_avg)
  tab$sugg_t2  <- !is.na(tab$pvalue_t2) & tab$pvalue_t2 < cutoff_sugg
  tab$sugg_t3  <- !is.na(tab$pvalue_t3) & tab$pvalue_t3 < cutoff_sugg
  tab$driven_by <- ifelse(tab$holds_t2 & tab$holds_t3, "both",
                          ifelse(tab$holds_t2, "t2 only",
                                 ifelse(tab$holds_t3, "t3 only", "neither")))
  
  table5 <- rbind(table5, tab[order(tab$pvalue_avg), ])
}

# Check: CpGs per exposure and driven_by
if (!is.null(table5)) {
  print(table(table5$exposure, factor(table5$driven_by, levels = c("both", "t2 only", "t3 only", "neither"))))
  write_xlsx(table5, file.path(outdir, "table5_avg_vs_trimesters.xlsx"))
  saveRDS(table5, file.path(outdir, "table5_avg_vs_trimesters.rds"))
} else {
  cat("No suggestive CpGs in the average models: Table 5 not written\n")
}

#==============================================================================#----
# TABLE 6: SUMMARY PER EXPOSURE                                                #----
#==============================================================================#----

table6 <- NULL
for (e in exposures) {
  res_avg <- main_res[[paste(e, "avg", sep = "_")]]
  res_t2  <- main_res[[paste(e, "t2",  sep = "_")]]
  res_t3  <- main_res[[paste(e, "t3",  sep = "_")]]
  
  cpg_avg <- res_avg$CpG[!is.na(res_avg$pvalue) & res_avg$pvalue < cutoff_sugg]
  cpg_t2  <- res_t2$CpG[!is.na(res_t2$pvalue) & res_t2$pvalue < cutoff_sugg]
  cpg_t3  <- res_t3$CpG[!is.na(res_t3$pvalue) & res_t3$pvalue < cutoff_sugg]
  
  # Reverse look-up: t2 / t3 suggestive CpGs that hold in the avg model
  m_t2_avg <- match(cpg_t2, res_avg$CpG)
  m_t3_avg <- match(cpg_t3, res_avg$CpG)
  holds_t2_in_avg <- res_avg$pvalue[m_t2_avg] < cutoff_nom &
    sign(res_avg$Effect[m_t2_avg]) == sign(res_t2$Effect[match(cpg_t2, res_t2$CpG)])
  holds_t3_in_avg <- res_avg$pvalue[m_t3_avg] < cutoff_nom &
    sign(res_avg$Effect[m_t3_avg]) == sign(res_t3$Effect[match(cpg_t3, res_t3$CpG)])
  
  # Correlation of effects across all CpGs
  rho_avg_t2 <- cor(res_avg$Effect, res_t2$Effect[match(res_avg$CpG, res_t2$CpG)],
                    method = "spearman", use = "complete.obs")
  rho_avg_t3 <- cor(res_avg$Effect, res_t3$Effect[match(res_avg$CpG, res_t3$CpG)],
                    method = "spearman", use = "complete.obs")
  
  # avg CpGs by driven_by (from Table 5)
  tab5_e <- if (is.null(table5)) NULL else table5[table5$exposure == e, ]
  n_driven <- function(x) if (is.null(tab5_e)) 0 else sum(tab5_e$driven_by == x)
  
  table6 <- rbind(table6, data.frame(
    exposure              = e,
    n_sugg_avg            = length(cpg_avg),
    n_sugg_t2             = length(cpg_t2),
    n_sugg_t3             = length(cpg_t3),
    n_sugg_avg_and_t2     = length(intersect(cpg_avg, cpg_t2)),
    n_sugg_avg_and_t3     = length(intersect(cpg_avg, cpg_t3)),
    n_sugg_all_three      = length(Reduce(intersect, list(cpg_avg, cpg_t2, cpg_t3))),
    n_avg_holds_both      = n_driven("both"),
    n_avg_holds_t2_only   = n_driven("t2 only"),
    n_avg_holds_t3_only   = n_driven("t3 only"),
    n_avg_holds_neither   = n_driven("neither"),
    n_t2_holds_in_avg     = sum(holds_t2_in_avg, na.rm = TRUE),
    n_t3_holds_in_avg     = sum(holds_t3_in_avg, na.rm = TRUE),
    rho_effects_avg_t2    = round(rho_avg_t2, 3),
    rho_effects_avg_t3    = round(rho_avg_t3, 3)
  ))
}

# Check: avg CpGs split by driven_by add up to the avg suggestive CpGs
stopifnot(all(table6$n_avg_holds_both + table6$n_avg_holds_t2_only +
                table6$n_avg_holds_t3_only + table6$n_avg_holds_neither == table6$n_sugg_avg))
table6
write_xlsx(table6, file.path(outdir, "table6_avg_summary.xlsx"))

#==============================================================================#----
# FOREST PLOTS - t2, t3 AND avg                                                #----
#==============================================================================#----

for (e in exposures) {
  
  exp_label <- exp_labels[match(e, exp_code)]
  res_tp <- list(t2  = main_res[[paste(e, "t2",  sep = "_")]],
                 t3  = main_res[[paste(e, "t3",  sep = "_")]],
                 avg = main_res[[paste(e, "avg", sep = "_")]])
  
  # CpGs suggestive at any time point, top by best p-value
  cpgs <- unique(unlist(lapply(res_tp, function(r) r$CpG[!is.na(r$pvalue) & r$pvalue < cutoff_sugg])))
  if (length(cpgs) == 0) { cat(e, ": no suggestive CpGs, no forest plot\n"); next }
  p_mat  <- sapply(res_tp, function(r) r$pvalue[match(cpgs, r$CpG)])
  p_mat  <- matrix(p_mat, ncol = length(time_points))
  best_p <- apply(p_mat, 1, min, na.rm = TRUE)
  cpgs   <- head(cpgs[order(best_p)], n_top_plot)
  
  # Gene from the avg results (t2 / t3 if missing there)
  gene <- res_tp$avg$gene_name[match(cpgs, res_tp$avg$CpG)]
  gene[is.na(gene)] <- res_tp$t2$gene_name[match(cpgs[is.na(gene)], res_tp$t2$CpG)]
  label <- make_label(cpgs, gene)
  
  df_fp <- NULL
  for (tp in time_points) {
    m <- match(cpgs, res_tp[[tp]]$CpG)
    df_fp <- rbind(df_fp, data.frame(label  = label, time = tp,
                                     Effect = res_tp[[tp]]$Effect[m],
                                     se     = res_tp[[tp]]$SE[m],
                                     suggestive = !is.na(res_tp[[tp]]$pvalue[m]) &
                                       res_tp[[tp]]$pvalue[m] < cutoff_sugg))
  }
  df_fp$lower <- df_fp$Effect - 1.96 * df_fp$se
  df_fp$upper <- df_fp$Effect + 1.96 * df_fp$se
  df_fp$label <- factor(df_fp$label, levels = rev(label))
  df_fp$time  <- factor(df_fp$time, levels = time_points)
  
  # Check: suggestive CpGs per time point among the plotted ones
  print(table(df_fp$time, df_fp$suggestive))
  
  # Plot 1. t2, t3 and avg: CpGs suggestive at any time point (y), effect with
  # 95% CI at each time point (x, colour); filled = suggestive at that time point.
  p_fp <- ggplot(df_fp, aes(x = Effect, y = label, colour = time, shape = suggestive)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_pointrange(aes(xmin = lower, xmax = upper),
                    position = position_dodge(width = 0.7), size = 0.3) +
    scale_colour_manual(values = pal_tp, labels = c(t2 = "t2", t3 = "t3", avg = "Average")) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                       labels = c(`TRUE` = paste("p <", cutoff_sugg), `FALSE` = "not suggestive")) +
    labs(title = paste0(exp_label, ": t2, t3 and average"),
         subtitle = paste0(tag, ", CpGs suggestive at any time point (top ", length(cpgs), ")"),
         x = "Effect (95% CI)", y = NULL, colour = "Time point", shape = NULL) +
    theme_ewas() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  save_plot(p_fp, file.path(outdir, paste0("forest_t2_t3_avg_", e, ".png")),
            width = 18, height = 4 + 0.7 * length(cpgs))
  cat("Saved forest_t2_t3_avg:", e, "\n")
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----