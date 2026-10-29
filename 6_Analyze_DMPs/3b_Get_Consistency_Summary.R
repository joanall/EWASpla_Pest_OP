#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: October - November 2026
#==============================================================================#
# [DESCRIPTION]: OBJECTIVE 2 - Consistency of the associations across time
# points and exposures. Starts from the suggestive CpGs (Tables 1-2) produced
# by 3a_Get_Suggestive_Summary.R (run it first with the same pattern).
# - Table 3:   t2 vs t3 per exposure: intersection, look-up replication
#              (same direction & p < 0.05), correlation of effects (all CpGs)
# - Table 4:   CpGs suggestive for >= 2 exposures at the same time point
# - Plot:      shared suggestive CpGs between all exposure x time point pairs
# Argument: pattern for the runs (one cohort, sample and run date), e.g.
#   Rscript script/6_Analyze_Results/3b_Get_Consistency_Summary.R "sepages_2tp_t[23]_20261021"
#
# [INPUT]:  results/5_Annotation/files_annotated_index.csv
#           -> results/5_Annotation/<run>/main/<exposure>_<tp>_annotated.rds
#           results/6_Analyze_Results/summary_tables/<cohort>_<sample>_<date>/
#           table1_suggestive_t2.rds, table2_suggestive_t3.rds
# [OUTPUT]: results/6_Analyze_Results/summary_tables/<cohort>_<sample>_<date>/
#           table3_t2_vs_t3.csv, table3_lookup_other_tp.csv/.rds,
#           table4_shared_exposures.csv/.rds,
#           overlap_matrix.csv, overlap_heatmap.png
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_template_plot.R")
source("script/00_functions.R")

index_file   <- "results/5_Annotation/files_annotated_index.csv"
cutoff_sugg  <- 1e-5    # suggestive threshold (as in 3a; plot subtitle only)
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

main_index <- files_index[files_index$analysis == "main", ]
table(main_index$exposure, main_index$time_point)   # 1 per cell

tag    <- paste(main_index$cohort[1], main_index$sample[1], main_index$run_date[1], sep = "_")
outdir <- file.path("results/6_Analyze_Results/summary_tables", tag)

exposures <- unique(main_index$exposure)

#==============================================================================#----
# 1. LOAD SUGGESTIVE CpGs (TABLES 1-2 FROM 3a)                                 #----
#==============================================================================#----

file_t1 <- file.path(outdir, "table1_suggestive_t2.rds")
file_t2 <- file.path(outdir, "table2_suggestive_t3.rds")
stopifnot(file.exists(file_t1), file.exists(file_t2))   # run 3a first

top_main <- rbind(readRDS(file_t1), readRDS(file_t2))
table(top_main$exposure, top_main$time_point)           # n suggestive per model

#==============================================================================#----
# 2. EFFECT AND p OF ALL CpGs OF EACH MAIN MODEL                               #----
#==============================================================================#----

main_lookup <- list()
for (i in seq_len(nrow(main_index))) {
  f   <- main_index[i, ]
  key <- paste(f$exposure, f$time_point, sep = "_")
  res <- readRDS(f$path)
  main_lookup[[key]] <- res[, c("CpG", "Effect", "pvalue")]
}
names(main_lookup)

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
table(top_main$time_point, top_main$replicated_other_tp, useNA = "ifany")

# CpG-level look-up (supplementary detail for Table 3)
lookup_cols <- c("exposure", "time_point", "CpG", "gene_name", "Effect", "pvalue",
                 "other_tp", "Effect_other_tp", "pvalue_other_tp", "replicated_other_tp")
lookup_tab  <- top_main[order(top_main$exposure, top_main$time_point, top_main$pvalue),
                        intersect(lookup_cols, colnames(top_main))]
write.csv(lookup_tab, file.path(outdir, "table3_lookup_other_tp.csv"), row.names = FALSE)
saveRDS(lookup_tab, file.path(outdir, "table3_lookup_other_tp.rds"))

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
stopifnot(all(diag(overlap_mx) == table(top_main$combo)[combo_order]))   # diagonal = Tables 1-2
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