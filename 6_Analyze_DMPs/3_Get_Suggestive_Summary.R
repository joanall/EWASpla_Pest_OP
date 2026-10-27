#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: Summary tables of suggestive CpGs (p < 1e-5) from the main
# models, and their consistency across sensitivity analyses, time points and
# exposures.
# - Table 1-2: suggestive CpGs at t2 / t3 + effect and p in s1-s5 (look-up)
# - Table 3:   t2 vs t3 per exposure: intersection, look-up replication
#              (same direction & p < 0.05), correlation of effects (all CpGs)
# - Table 4:   CpGs suggestive for >= 2 exposures at the same time point
# - Plot:      shared suggestive CpGs between all exposure x time point pairs
# Argument: pattern for the runs (one cohort, sample and run date), e.g.
#   Rscript script/6_Analyze_Results/4_Get_Summary_Tables.R "sepages_2tp_t[23]_20261021"
#
# [INPUT]:  results/5_Annotation/files_annotated_index.csv
#           -> results/5_Annotation/<run>/<analysis>/<exposure>_<tp>_annotated.rds
# [OUTPUT]: results/6_Analyze_Results/summary_tables/<cohort>_<sample>_<date>/
#           table1_suggestive_t2.csv/.rds, table2_suggestive_t3.csv/.rds,
#           table3_t2_vs_t3.csv, table4_shared_exposures.csv/.rds,
#           overlap_matrix.csv, overlap_heatmap.png
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_template_plot.R")
source("script/00_functions.R")

index_file   <- "results/5_Annotation/files_annotated_index.csv"
cutoff_sugg  <- 1e-5    # suggestive threshold
cutoff_nom   <- 0.05    # nominal threshold for look-up replication

args <- commandArgs(trailingOnly = TRUE)
run_pattern <- if (length(args) >= 1 && args[1] != "") args[1] else "."
# run_pattern <- "sepages_2tp_t[23]_20261021"   # TEST

files_index <- read.csv(index_file, stringsAsFactors = FALSE)
files_index <- files_index[grepl(run_pattern, files_index$run), ]
stopifnot(all(file.exists(files_index$path)))

# Check: one cohort, sample and run date (comparisons must not mix runs)
stopifnot(length(unique(files_index$cohort)) == 1,
          length(unique(files_index$sample)) == 1,
          length(unique(files_index$run_date)) == 1)
table(files_index$time_point, files_index$analysis)   # 6 per cell

main_index <- files_index[files_index$analysis == "main", ]
sens_index <- files_index[files_index$analysis != "main", ]

tag    <- paste(main_index$cohort[1], main_index$sample[1], main_index$run_date[1], sep = "_")
outdir <- file.path("results/6_Analyze_Results/summary_tables", tag)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

exposures <- unique(main_index$exposure)
keep_cols <- c("CpG", "chr", "pos", "gene_name", "Relation_to_Island", "States15_FP",
               "Effect", "SE", "pvalue", "pval_bn")

#==============================================================================#----
# 1. SUGGESTIVE CpGs - MAIN MODELS                                             #----
#==============================================================================#----
# Also keep Effect and p of all CpGs of each main model for the t2-t3 look-up

top_main    <- NULL
main_lookup <- list()

for (i in seq_len(nrow(main_index))) {
  f   <- main_index[i, ]
  key <- paste(f$exposure, f$time_point, sep = "_")
  res <- readRDS(f$path)
  
  main_lookup[[key]] <- res[, c("CpG", "Effect", "pvalue")]
  
  top <- res[!is.na(res$pvalue) & res$pvalue < cutoff_sugg, intersect(keep_cols, colnames(res))]
  cat(key, "| suggestive:", nrow(top), "\n")
  if (nrow(top) == 0) next
  top_main <- rbind(top_main, data.frame(exposure = f$exposure, time_point = f$time_point, top))
}
nrow(top_main)

#==============================================================================#----
# 2. LOOK-UP IN SENSITIVITY ANALYSES                                           #----
#==============================================================================#----

sens_names <- sort(unique(sens_index$analysis))

for (s in sens_names) {
  top_main[[paste0("Effect_", s)]] <- NA
  top_main[[paste0("SE_", s)]]     <- NA
  top_main[[paste0("pvalue_", s)]] <- NA
  }


for (i in seq_len(nrow(sens_index))) {
  f   <- sens_index[i, ]
  idx <- which(top_main$exposure == f$exposure & top_main$time_point == f$time_point)
  if (length(idx) == 0) next
  res <- readRDS(f$path)
  m   <- match(top_main$CpG[idx], res$CpG)
  top_main[idx, paste0("Effect_", f$analysis)] <- res$Effect[m]
  top_main[idx, paste0("SE_", f$analysis)]     <- res$SE[m]
  top_main[idx, paste0("pvalue_", f$analysis)] <- res$pvalue[m]
}

# Agreement with the main model: same direction, and p < 0.05
eff_sens <- as.matrix(top_main[, paste0("Effect_", sens_names)])
p_sens   <- as.matrix(top_main[, paste0("pvalue_", sens_names)])
top_main$n_sens_same_dir <- rowSums(sign(eff_sens) == sign(top_main$Effect), na.rm = TRUE)
top_main$n_sens_p005     <- rowSums(p_sens < cutoff_nom & sign(eff_sens) == sign(top_main$Effect), na.rm = TRUE)

#==============================================================================#----
# 3. LOOK-UP IN THE OTHER TIME POINT                                           #----
#==============================================================================#----

top_main$other_tp        <- ifelse(top_main$time_point == "t2", "t3", "t2")
top_main$Effect_other_tp <- NA
top_main$pvalue_other_tp <- NA

for (e in exposures) {
  for (tp in c("t2", "t3")) {
    idx <- which(top_main$exposure == e & top_main$time_point == tp)
    lk  <- main_lookup[[paste(e, ifelse(tp == "t2", "t3", "t2"), sep = "_")]]
    if (length(idx) == 0 || is.null(lk)) next
    m <- match(top_main$CpG[idx], lk$CpG)
    top_main$Effect_other_tp[idx] <- lk$Effect[m]
    top_main$pvalue_other_tp[idx] <- lk$pvalue[m]
  }
}
top_main$replicated_other_tp <- top_main$pvalue_other_tp < cutoff_nom &
  sign(top_main$Effect_other_tp) == sign(top_main$Effect)

#==============================================================================#----
# TABLES 1-2: SUGGESTIVE CpGs PER TIME POINT                                   #----
#==============================================================================#----

for (tp in c("t2", "t3")) {
  tab <- top_main[top_main$time_point == tp, ]
  tab <- tab[order(tab$exposure, tab$pvalue), ]
  n_tab <- ifelse(tp == "t2", 1, 2)
  file_base <- file.path(outdir, paste0("table", n_tab, "_suggestive_", tp))
  write.csv(tab, paste0(file_base, ".csv"), row.names = FALSE)
  saveRDS(tab, paste0(file_base, ".rds"))
  cat("Table", n_tab, "-", tp, ":", nrow(tab), "suggestive CpGs\n")
}

#==============================================================================#----
# TABLE 3: T2 vs T3 PER EXPOSURE                                               #----
#==============================================================================#----

table3 <- NULL
for (e in exposures) {
  cpg_t2 <- top_main$CpG[top_main$exposure == e & top_main$time_point == "t2"]
  cpg_t3 <- top_main$CpG[top_main$exposure == e & top_main$time_point == "t3"]
  
  # Correlation of effects across all CpGs (overall concordance t2 - t3)
  lk_t2 <- main_lookup[[paste(e, "t2", sep = "_")]]
  lk_t3 <- main_lookup[[paste(e, "t3", sep = "_")]]
  eff_t3_matched <- lk_t3$Effect[match(lk_t2$CpG, lk_t3$CpG)]
  rho_all <- cor(lk_t2$Effect, eff_t3_matched, method = "spearman", use = "complete.obs")
  
  table3 <- rbind(table3, data.frame(
    exposure          = e,
    n_sugg_t2         = length(cpg_t2),
    n_sugg_t3         = length(cpg_t3),
    n_sugg_both       = length(intersect(cpg_t2, cpg_t3)),
    n_t2_replic_in_t3 = sum(top_main$replicated_other_tp[top_main$exposure == e & top_main$time_point == "t2"], na.rm = TRUE),
    n_t3_replic_in_t2 = sum(top_main$replicated_other_tp[top_main$exposure == e & top_main$time_point == "t3"], na.rm = TRUE),
    rho_effects_all   = round(rho_all, 3)
  ))
}
table3
write.csv(table3, file.path(outdir, "table3_t2_vs_t3.csv"), row.names = FALSE)

#==============================================================================#----
# TABLE 4: CpGs SHARED BETWEEN EXPOSURES (SAME TIME POINT)                     #----
#==============================================================================#----

shared <- aggregate(exposure ~ CpG + time_point, data = top_main,
                    FUN = function(x) paste(sort(unique(x)), collapse = ";"))
shared$n_exposures <- lengths(strsplit(shared$exposure, ";"))
shared <- shared[shared$n_exposures >= 2, ]
shared$gene_name <- top_main$gene_name[match(shared$CpG, top_main$CpG)]
shared <- shared[order(shared$time_point, -shared$n_exposures), ]
shared
write.csv(shared, file.path(outdir, "table4_shared_exposures.csv"), row.names = FALSE)
saveRDS(shared, file.path(outdir, "table4_shared_exposures.rds"))

#==============================================================================#----
# OVERLAP MATRIX & HEATMAP (ALL EXPOSURE x TIME POINT PAIRS)                   #----
#==============================================================================#----

top_main$combo <- paste(top_main$exposure, top_main$time_point, sep = "_")
combo_order    <- as.vector(outer(exposures, c("t2", "t3"), paste, sep = "_"))
combo_order    <- combo_order[combo_order %in% top_main$combo]

# Incidence (CpG x model) -> shared counts; diagonal = n suggestive per model
incidence  <- unclass(table(top_main$CpG, factor(top_main$combo, levels = combo_order))) > 0
overlap_mx <- crossprod(incidence * 1)
overlap_mx
write.csv(overlap_mx, file.path(outdir, "overlap_matrix.csv"))

df_ov   <- expand.grid(x = combo_order, y = combo_order)
df_ov$n <- as.vector(overlap_mx)
df_ov$x <- factor(df_ov$x, levels = combo_order)
df_ov$y <- factor(df_ov$y, levels = rev(combo_order))

# Plot 1. Number of shared suggestive CpGs between each pair of main models
# (exposure x time point); diagonal = suggestive CpGs of each model.
p_ov <- ggplot(df_ov, aes(x = x, y = y, fill = n)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = n), size = 3) +
  fill_heat +
  labs(title = "Shared suggestive CpGs between models",
       subtitle = paste0(tag, ", main models, p < ", cutoff_sugg),
       x = NULL, y = NULL, fill = "n CpGs") +
  coord_fixed() +
  theme_ewas()
save_plot(p_ov, file.path(outdir, "overlap_heatmap.png"), width = 18, height = 16)

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----