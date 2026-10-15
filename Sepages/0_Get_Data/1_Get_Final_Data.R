#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]:Filter data based on information available at each time point
# and for each sensitivity analysis. Main subsets require:
# 1) Placental DNAm (outcome)
# 2) Pesticide data at the time point (exposure)
# 3) Main covariates: maternal age, education, smoking, GA, child ethnicity,
#    child sex (+ cell types and DNA contamination, complete in DNAm samples)
#
# Subsets, for t2 and t3:
#  - Main
#  - S1: Excluding pregnancy complications
#  - S2: Excluding non-Caucasian
#  - S3: Adjusting for specific gravity (complete SG)
#  - S4: Adjusting for fruit and vegetable intake (complete FFQ)
#  - S5: Without GA (main covariates minus GA)
#
# [INPUT]:  db/db_pest_pladnam_incomplete_cases.rds
# [OUTPUT]: db/subsets/df_main_t2.rds, df_main_t3.rds, df_s1_t2.rds ... df_s5_t3.rds
#==============================================================================#
#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

infile <- "db/db_pest_pladnam_incomplete_cases.rds"
outdir <- "db/subsets/"
dir.create(outdir, showWarnings = FALSE)

df <- readRDS(infile)
nrow(df) # 484

#==============================================================================#----
# DNAm + EXPOSURE                                                              #----
#==============================================================================#----

# 1. Participants with placental DNAm
with_dnam <- df[!is.na(df$Meanlog2oddsContamination), ]
nrow(with_dnam) # 389

# 2. Participants with pesticide data at each time point
with_dnam_chem_t2 <- with_dnam[!is.na(with_dnam$mo_pest_24d_cat_t2), ] # 389
with_dnam_chem_t3 <- with_dnam[!is.na(with_dnam$mo_pest_24d_cat_t3), ] # 386

#==============================================================================#----
# MAIN COVARIATES                                                              #----
#==============================================================================#----

main_cov    <- c("mo_age", "mo_edu", "mo_smoke", "po_gd_lmp_days","ch_ethnicity", "ch_sex")
main_cov_s5 <- c("mo_age", "mo_edu", "mo_smoke","ch_ethnicity", "ch_sex")   

# Check: all covariates exist 
setdiff(main_cov, colnames(df))

# Check: missingness per covariate
colSums(is.na(with_dnam_chem_t2[, intersect(main_cov, colnames(df))]))
# mo_age         mo_edu       mo_smoke po_gd_lmp_days   ch_ethnicity 
#.     0              2             35              0             42 

colSums(is.na(with_dnam_chem_t3[, intersect(main_cov, colnames(df))]))
# mo_age         mo_edu       mo_smoke po_gd_lmp_days   ch_ethnicity 
#      0              2             35              0             44 

# Main subsets
df_main_t2 <- with_dnam_chem_t2[complete.cases(with_dnam_chem_t2[, main_cov]), ] # 314
df_main_t3 <- with_dnam_chem_t3[complete.cases(with_dnam_chem_t3[, main_cov]), ] # 313

#==============================================================================#----
# SENSITIVITY SUBSETS                                                          #----
#==============================================================================#----

# S1 - Excluding pregnancy complications
df_s1_t2 <- df_main_t2[!is.na(df_main_t2$preg_compl) & df_main_t2$preg_compl == "No", ] # 279
df_s1_t3 <- df_main_t3[!is.na(df_main_t3$preg_compl) & df_main_t3$preg_compl == "No", ] # 278

# S2 - Excluding non-Caucasian
df_s2_t2 <- df_main_t2[df_main_t2$ch_ethnicity == "Caucasian", ] # 292
df_s2_t3 <- df_main_t3[df_main_t3$ch_ethnicity == "Caucasian", ] # 291

# S3 - Adjusting for specific gravity (complete SG at the matching time point)
df_s3_t2 <- df_main_t2[!is.na(df_main_t2$mo_pool_sg_T1), ] # 314
df_s3_t3 <- df_main_t3[!is.na(df_main_t3$mo_pool_sg_T3), ] # 313

# S4 - Adjusting for fruit and vegetable intake
df_s4_t2 <- df_main_t2[!is.na(df_main_t2$fv_g1000kcal), ] # 281
df_s4_t3 <- df_main_t3[!is.na(df_main_t3$fv_g1000kcal), ] # 280

# S5 - Without GA (starts from DNAm + chem, not from main)
df_s5_t2 <- with_dnam_chem_t2[complete.cases(with_dnam_chem_t2[, main_cov_s5]), ] # 314
df_s5_t3 <- with_dnam_chem_t3[complete.cases(with_dnam_chem_t3[, main_cov_s5]), ] # 313

# Check: N per subset
n_subsets <- data.frame(
  subset = c("main", "s1_no_compl", "s2_caucasian", "s3_sg", "s4_fruitveg", "s5_no_ga"),
  t2 = c(nrow(df_main_t2), nrow(df_s1_t2), nrow(df_s2_t2),
         nrow(df_s3_t2),   nrow(df_s4_t2), nrow(df_s5_t2)),
  t3 = c(nrow(df_main_t3), nrow(df_s1_t3), nrow(df_s2_t3),
         nrow(df_s3_t3),   nrow(df_s4_t3), nrow(df_s5_t3))
)
n_subsets

#==============================================================================#----
# SAVE                                                                         #----
#==============================================================================#----

saveRDS(df_main_t2, paste0(outdir, "df_main_t2.rds"))
saveRDS(df_main_t3, paste0(outdir, "df_main_t3.rds"))
saveRDS(df_s1_t2,   paste0(outdir, "df_s1_t2.rds"))
saveRDS(df_s1_t3,   paste0(outdir, "df_s1_t3.rds"))
saveRDS(df_s2_t2,   paste0(outdir, "df_s2_t2.rds"))
saveRDS(df_s2_t3,   paste0(outdir, "df_s2_t3.rds"))
saveRDS(df_s3_t2,   paste0(outdir, "df_s3_t2.rds"))
saveRDS(df_s3_t3,   paste0(outdir, "df_s3_t3.rds"))
saveRDS(df_s4_t2,   paste0(outdir, "df_s4_t2.rds"))
saveRDS(df_s4_t3,   paste0(outdir, "df_s4_t3.rds"))
saveRDS(df_s5_t2,   paste0(outdir, "df_s5_t2.rds"))
saveRDS(df_s5_t3,   paste0(outdir, "df_s5_t3.rds"))
write.csv(n_subsets, paste0(outdir, "n_subsets.csv"), row.names = FALSE)

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----
