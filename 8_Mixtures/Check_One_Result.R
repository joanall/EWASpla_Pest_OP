

res_mix <- readRDS("results/2_Gcomp_Results/sepages_2tp_t2_20261109/main/mix_t2.rds")
beta_cols <- c("beta_3pba", "beta_tcpy", "beta_dmp", "beta_dep")

sum(is.na(res_mix$pvalue))                       # failed CpGs (should be ~0)
lambda_mix <- median(qchisq(1 - res_mix$pvalue, 1), na.rm = TRUE) / qchisq(0.5, 1)
lambda_mix                                       # genomic inflation

# Table of hits
res_mix$bonf <- p.adjust(res_mix$pvalue, method = "bonferroni")
cat("BF < 0.05:", sum(res_mix$bonf < 0.05, na.rm = TRUE),
    "| p < 1e-5:", sum(res_mix$pvalue < 1e-5, na.rm = TRUE), "\n")

# Weights from the per-compound betas (within each direction)
beta_mat <- as.matrix(res_mix[, beta_cols])
beta_pos <- pmax(beta_mat, 0); beta_neg <- pmin(beta_mat, 0)
w_pos <- beta_pos / rowSums(beta_pos); w_neg <- beta_neg / rowSums(beta_neg)
colnames(w_pos) <- sub("beta_", "wpos_", beta_cols)
colnames(w_neg) <- sub("beta_", "wneg_", beta_cols)
res_mix <- cbind(res_mix, round(w_pos, 2), round(w_neg, 2))

# Top hits + annotation
cpg_annot <- readRDS("results/5_Annotation/cpg_annotation_all.rds")
hits_mix  <- res_mix[!is.na(res_mix$pvalue) & res_mix$pvalue < 1e-5, ]
hits_mix  <- merge(hits_mix, cpg_annot, by = "CpG", all.x = TRUE)
hits_mix  <- hits_mix[order(hits_mix$pvalue), ]
nrow(hits_mix)

# Volcano and Manhattan plots: reuse your 6_Analyze_DMPs 

# Weighted plot 
top_n <- min(20, nrow(hits_mix))
beta_top <- t(as.matrix(hits_mix[1:top_n, beta_cols]))
colnames(beta_top) <- hits_mix$CpG[1:top_n]
barplot(beta_top, beside = FALSE, las = 2, cex.names = 0.6,
        col = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"),
        legend.text = sub("beta_", "", beta_cols), ylab = "Coefficient per tertile")
abline(h = 0)


# dominant compound 
dominant_mix <- sub("beta_", "", beta_cols)[apply(abs(as.matrix(hits_mix[, beta_cols])), 1, which.max)]
table(dominant_mix)




Keep cohort-specific tertiles, and describe the result as “per cohort-specific 
tertile increase”. This is common in consortium papers and is acceptable.
Use shared cut points: compute the tertile breaks once on both cohorts’ exposures 
pooled together, then pass them to quantize() in each cohort. Only the cut points 
need to be shared, not the individual-level data. This is cleaner if you want psi 
to mean the same exposure contrast in both cohorts.



res_mix_bisc    <- readRDS("results/2_Gcomp_Results/bisc_2tp_t2_<date>/main/mix_t2.rds")
res_mix_sepages <- readRDS("results/2_Gcomp_Results/sepages_2tp_t2_<date>/main/mix_t2.rds")

cpg_common <- intersect(res_mix_bisc$CpG, res_mix_sepages$CpG)
cat("CpGs in both cohorts:", length(cpg_common), "\n")
b_bisc <- res_mix_bisc[cpg_common, ];    b_sepages <- res_mix_sepages[cpg_common, ]
stopifnot(identical(b_bisc$CpG, b_sepages$CpG))

w_bisc    <- 1 / b_bisc$SE^2
w_sepages <- 1 / b_sepages$SE^2
meta_mix <- data.frame(CpG = cpg_common,
                       Effect = (b_bisc$Effect * w_bisc + b_sepages$Effect * w_sepages) / (w_bisc + w_sepages),
                       SE     = sqrt(1 / (w_bisc + w_sepages)))
meta_mix$pvalue <- 2 * pnorm(-abs(meta_mix$Effect / meta_mix$SE))

# Heterogeneity (Cochran's Q, I2) and direction agreement
q_het <- w_bisc * (b_bisc$Effect - meta_mix$Effect)^2 + w_sepages * (b_sepages$Effect - meta_mix$Effect)^2
meta_mix$p_het <- pchisq(q_het, df = 1, lower.tail = FALSE)
meta_mix$I2    <- pmax(0, (q_het - 1) / q_het)
meta_mix$same_direction <- sign(b_bisc$Effect) == sign(b_sepages$Effect)
meta_mix$fdr   <- p.adjust(meta_mix$pvalue, method = "BH")