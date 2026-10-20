#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: Recode the data based on the analysis plan and built
# the subsets for each model. 
# All subsets require:
# 1) Placental DNAm (outcome)
# 2) Pesticide data at the time point (exposure)
# Main model: maternal age, education, smoking, GA, child ethnicity,
#             child sex, cell types and DNA contamination
#
# Subsets, for t2 and t3:
#  - Main
#  - S1: Excluding pregnancy complications
#  - S2: Excluding non-Caucasian
#  - S3: Adjusting for specific gravity (complete SG)
#  - S4: Adjusting for fruit and vegetable intake (complete FFQ)
#  - S5: Without GA (main covariates minus GA and Cell Types)
#
# [INPUT]:  db/db_pest_pladnam_incomplete_cases.rds
# [OUTPUT]: db/subsets/df_main_t2.rds, df_main_t3.rds, df_s1_t2.rds ... df_s5_t3.rds

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
library(dplyr)

infile <- "db/db_pest_pladnam_incomplete_cases.rds"
outdir <- "db/subsets/"
dir.create(outdir, showWarnings = FALSE)

df <- readRDS(infile)
nrow(df) # 484


#==============================================================================#----
# RENAME                                                                       #----
#==============================================================================#----
# Safety coy 
df_cc <- df

# Old -> new names
# Only variables whose name changes are listed; all others keep their name
# (mo_age, mo_par, fa_ethnicity, mo_ethnicity, ch_sex, Basename,
#  cell types, Meanlog2oddsContamination).
rename_map <- matrix(c(
  # Covariates
  "ident",                 "mo_id",
  "mo_dipl",               "mo_edu",
  "mo_tob_gr_anyt_yn1_n2", "mo_smoke",
  "mo_bmi_bepr",           "mo_bmi",
  "fv_g1000kcal",          "mo_fv_g1000kcal",
  "po_gd_lmp",             "gest_age",        # still in WEEKS: x7 at recode
  "mt1haa1_q40",           "preg_compl",
  
  # Exposures for models (before recoding / log2)
  "mo_pest_24d_cat_t2",    "pest_24d_cat_t2",
  "mo_pest_24d_cat_t3",    "pest_24d_cat_t3",
  "mo_pest_3pba_i_t2",     "pest_3pba_i_t2",
  "mo_pest_3pba_i_t3",     "pest_3pba_i_t3",
  "mo_pest_tcpy_i_t2",     "pest_tcpy_i_t2",
  "mo_pest_tcpy_i_t3",     "pest_tcpy_i_t3",
  "mo_pest_dep_ms_i_t2",   "pest_dep_ms_i_t2",
  "mo_pest_dep_ms_i_t3",   "pest_dep_ms_i_t3",
  "mo_pest_dmp_ms_i_t2",   "pest_dmp_ms_i_t2",
  "mo_pest_dmp_ms_i_t3",   "pest_dmp_ms_i_t3",
  "mo_pest_dap_ms_i_t2",   "pest_dap_ms_i_t2",
  "mo_pest_dap_ms_i_t3",   "pest_dap_ms_i_t3",
  
  # Raw concentrations (descriptives)
  "mo_pest_24d_raw_t2",    "pest_24d_raw_t2",
  "mo_pest_24d_raw_t3",    "pest_24d_raw_t3",
  "mo_pest_3pba_raw_t2",   "pest_3pba_raw_t2",
  "mo_pest_3pba_raw_t3",   "pest_3pba_raw_t3",
  "mo_pest_tcpy_raw_t2",   "pest_tcpy_raw_t2",
  "mo_pest_tcpy_raw_t3",   "pest_tcpy_raw_t3",
  "mo_pest_dmp_raw_t2",    "pest_dmp_raw_t2",
  "mo_pest_dmp_raw_t3",    "pest_dmp_raw_t3",
  "mo_pest_dmtp_raw_t2",   "pest_dmtp_raw_t2",
  "mo_pest_dmtp_raw_t3",   "pest_dmtp_raw_t3",
  "mo_pest_dmdtp_raw_t2",  "pest_dmdtp_raw_t2",
  "mo_pest_dmdtp_raw_t3",  "pest_dmdtp_raw_t3",
  "mo_pest_dep_raw_t2",    "pest_dep_raw_t2",
  "mo_pest_dep_raw_t3",    "pest_dep_raw_t3",
  "mo_pest_detp_raw_t2",   "pest_detp_raw_t2",
  "mo_pest_detp_raw_t3",   "pest_detp_raw_t3",
  
  # Specific gravity (check: T1 in the source name = 2nd-trimester pool?)
  "mo_pool_sg_T1",         "pool_sg_t2",
  "mo_pool_sg_T3",         "pool_sg_t3"
), ncol = 2, byrow = TRUE)

old_names <- rename_map[, 1]
new_names <- rename_map[, 2]

# Rename
names(df)[match(old_names, names(df))] <- new_names


#==============================================================================#----
# RECODE                                                                       #----
#==============================================================================#----

#------------------------------------------------------------------------------#
# Covariates                                                                   #
#------------------------------------------------------------------------------#

# Maternal age - Continuous - years
df$mo_age <- as.numeric(df$mo_age)
summary(df$mo_age)

# Maternal education - 2 cat - University / Non-university
# Raw: 2 = 29 | 3 = 54 | 4 = 127 | 5 = 272 | NA = 2
df$mo_edu <- ifelse(df$mo_edu >= 3, "University", "Non-university")
df$mo_edu <- factor(df$mo_edu, levels = c("Non-university", "University"))

# Maternal smoking - 2 cat - Smoker / Non-smoker
# Raw: 0 = 408 | 1 = 20 | NA = 56
df$mo_smoke <- ifelse(df$mo_smoke == 1, "Smoker",
                      ifelse(df$mo_smoke == 0, "Non-smoker", NA))
df$mo_smoke <- factor(df$mo_smoke, levels = c("Non-smoker", "Smoker"))

# Pre-pregnancy BMI - Continuous - kg/m2
df$mo_bmi <- as.numeric(df$mo_bmi)
summary(df$mo_bmi)

# Parity - 2 cat - Nulliparous / Multiparous
# Raw: 0 = 22 | 1 = 214 | 2 = 48
df$mo_parity <- ifelse(df$mo_par == 0, "Nulliparous",
                       ifelse(df$mo_par > 0,  "Multiparous", NA))
df$mo_parity <- factor(df$mo_parity, levels = c("Nulliparous", "Multiparous"))

# Fruit and vegetable intake - Continuous - g/1000kcal/day
df$mo_fv_g1000kcal <- as.numeric(df$mo_fv_g1000kcal)
summary(df$mo_fv_g1000kcal)

# Child ethnicity - 2 cat - European / Non-European (both parents Caucasian = 1)
df$ch_eth <- ifelse(df$mo_ethnicity == 1 & df$fa_ethnicity == 1,
                    "European", "Non-European")
df$ch_eth <- factor(df$ch_eth, levels = c("European", "Non-European"))
table(df$mo_ethnicity, df$fa_ethnicity, useNA = "ifany")
table(df$ch_eth, useNA = "ifany")

# Child sex - 2 cat - Male / Female (1 = Male, 2 = Female)
df$ch_sex <- factor(df$ch_sex, levels = c(1, 2), labels = c("Male", "Female"))
table(df_cc$ch_sex, df$ch_sex, useNA = "ifany")

# Gestational age at birth - Continuous - weeks -> days
# Raw (weeks): min 28.57 | median 40.00 | max 43.43 | NA = 4
df$gest_age <- as.numeric(df$gest_age) * 7
summary(df$gest_age)

# Pregnancy complications - 2 cat - Yes / No (1 = Yes; 0, 2 = No)
df$preg_compl <- ifelse(df$preg_compl == 1, "Yes",
                        ifelse(df$preg_compl %in% c(0, 2), "No", NA))
df$preg_compl <- factor(df$preg_compl, levels = c("No", "Yes"))
table(df_cc$mt1haa1_q40, df$preg_compl, useNA = "ifany")

#------------------------------------------------------------------------------#
# Methylation-related and design covariates                                    #
#------------------------------------------------------------------------------#
ct_vars <- c("Trophoblasts", "Stromal", "Hofbauer", "Endothelial",
             "nRBC", "Syncytiotrophoblast")
for (v in ct_vars) df[[v]] <- as.numeric(df[[v]])
summary(df[, ct_vars])

df$Meanlog2oddsContamination <- as.numeric(df$Meanlog2oddsContamination)
summary(df$Meanlog2oddsContamination)

df$pool_sg_t2 <- as.numeric(df$pool_sg_t2)
df$pool_sg_t3 <- as.numeric(df$pool_sg_t3)
summary(df[, c("pool_sg_t2", "pool_sg_t3")])

#------------------------------------------------------------------------------#
# Pesticides - raw concentrations (descriptives) - Continuous - ng/mL          #
#------------------------------------------------------------------------------#
raw_vars <- grep("^pest_.*_raw_t[23]$", names(df), value = TRUE)
raw_vars                                              # 16 variables
for (v in raw_vars) df[[v]] <- as.numeric(df[[v]])
summary(df[, raw_vars])                               # all >= 0

#------------------------------------------------------------------------------#
# Pesticides - exposures for models                                            #
#------------------------------------------------------------------------------#

# 2,4-D - 2 cat - Detected / Not detected
# 1 = <LOD -> Not detected; 2 = LOD-LOQ, 3 = >LOQ -> Detected
df$pest_24d_bi_t2 <- ifelse(df$pest_24d_cat_t2 >= 2, "Detected", "Not detected")
df$pest_24d_bi_t2 <- factor(df$pest_24d_bi_t2, levels = c("Not detected", "Detected"))
table(df$pest_24d_cat_t2, df$pest_24d_bi_t2, useNA = "ifany")

df$pest_24d_bi_t3 <- ifelse(df$pest_24d_cat_t3 >= 2, "Detected", "Not detected")
df$pest_24d_bi_t3 <- factor(df$pest_24d_bi_t3, levels = c("Not detected", "Detected"))
table(df$pest_24d_cat_t3, df$pest_24d_bi_t3, useNA = "ifany")

# Imputed compounds and molar sums - Continuous - log2
# Inputs are on the natural scale (ng/mL or umol/L): all must be > 0 before log2
i_vars <- c("pest_3pba_i_t2",   "pest_3pba_i_t3",
            "pest_tcpy_i_t2",   "pest_tcpy_i_t3",
            "pest_dep_ms_i_t2", "pest_dep_ms_i_t3",
            "pest_dmp_ms_i_t2", "pest_dmp_ms_i_t3",
            "pest_dap_ms_i_t2", "pest_dap_ms_i_t3")
sapply(df[, i_vars], function(x) sum(x <= 0, na.rm = TRUE))   # all 0

df$pest_3pba_ilog2_t2   <- log2(as.numeric(df$pest_3pba_i_t2))
df$pest_3pba_ilog2_t3   <- log2(as.numeric(df$pest_3pba_i_t3))
df$pest_tcpy_ilog2_t2   <- log2(as.numeric(df$pest_tcpy_i_t2))
df$pest_tcpy_ilog2_t3   <- log2(as.numeric(df$pest_tcpy_i_t3))
df$pest_dep_ms_ilog2_t2 <- log2(as.numeric(df$pest_dep_ms_i_t2))
df$pest_dep_ms_ilog2_t3 <- log2(as.numeric(df$pest_dep_ms_i_t3))
df$pest_dmp_ms_ilog2_t2 <- log2(as.numeric(df$pest_dmp_ms_i_t2))
df$pest_dmp_ms_ilog2_t3 <- log2(as.numeric(df$pest_dmp_ms_i_t3))
df$pest_dap_ms_ilog2_t2 <- log2(as.numeric(df$pest_dap_ms_i_t2))
df$pest_dap_ms_ilog2_t3 <- log2(as.numeric(df$pest_dap_ms_i_t3))

log2_vars <- grep("_ilog2_t[23]$", names(df), value = TRUE)
summary(df[, log2_vars])
sapply(df[, log2_vars], function(x) sum(is.infinite(x)))   # all 0

#------------------------------------------------------------------------------#
# Final check: class of every analysis variable                                #
#------------------------------------------------------------------------------#
sapply(df[, c("mo_age", "mo_edu", "mo_smoke", "mo_bmi", "mo_parity",
              "mo_fv_g1000kcal", "ch_eth", "ch_sex", "gest_age", "preg_compl",
              "pest_24d_bi_t2", "pest_24d_bi_t3", log2_vars)], class)



#==============================================================================#----
# MAIN DATA                                                                    #----
#==============================================================================#----
# All subsets require placental DNAm (outcome) and pesticide data (exposure)

# 1. Participants with placental DNAm
with_dnam <- df[!is.na(df$Meanlog2oddsContamination), ]
nrow(with_dnam) # 389

# 2. Participants with all model exposures at both T2 and T3
expo_vars <- c("pest_24d_bi_t2",       "pest_24d_bi_t3",
               "pest_3pba_ilog2_t2",   "pest_3pba_ilog2_t3",
               "pest_tcpy_ilog2_t2",   "pest_tcpy_ilog2_t3",
               "pest_dep_ms_ilog2_t2", "pest_dep_ms_ilog2_t3",
               "pest_dmp_ms_ilog2_t2", "pest_dmp_ms_ilog2_t3",
               "pest_dap_ms_ilog2_t2", "pest_dap_ms_ilog2_t3")
colSums(is.na(with_dnam[, expo_vars]))

with_dnam_chem <- with_dnam[complete.cases(with_dnam[, expo_vars]), ] # 385
nrow(with_dnam_chem)

#==============================================================================#----
# COVARIATE SETS                                                               #----
#==============================================================================#----
ct_vars <- c("Trophoblasts", "Stromal", "Hofbauer", "Endothelial",
             "nRBC", "Syncytiotrophoblast")

# Main model
main_cov <- c("mo_age", "mo_edu", "mo_smoke", "gest_age", "ch_eth", "ch_sex",
              ct_vars, "Meanlog2oddsContamination")

# S5: main minus GA and cell types
main_cov_s5 <- c("mo_age", "mo_edu", "mo_smoke", "ch_eth", "ch_sex",
                 "Meanlog2oddsContamination")

colSums(is.na(with_dnam_chem[, main_cov]))
# Fro, 385 ,other with plaNDAm and Chem at T2 and T3:
# - 2 mising in education
# - 42 ,issing in ethnicity 
# - 35 missing in smoking 
#==============================================================================#----
# MAIN SUBSET                                                                  #----
#==============================================================================#----
df_main <- with_dnam_chem[complete.cases(with_dnam_chem[, main_cov]), ] # 312 ,others
nrow(df_main)

#==============================================================================#----
# SENSITIVITY SUBSETS                                                          #----
#==============================================================================#----

# S1 - Excluding pregnancy complications (NA in preg_compl also excluded)
df_s1 <- df_main[df_main$preg_compl %in% "No", ]
table(df_main$preg_compl, useNA = "ifany")
nrow(df_s1)

# S2 - Excluding non-European
df_s2 <- df_main[df_main$ch_eth == "European", ]
nrow(df_s2)

# S3 - Adjusting for specific gravity (complete SG at both time points)
df_s3 <- df_main[!is.na(df_main$pool_sg_t2) & !is.na(df_main$pool_sg_t3), ]
nrow(df_s3)

# S4 - Adjusting for fruit and vegetable intake (complete FFQ)
df_s4 <- df_main[!is.na(df_main$mo_fv_g1000kcal), ]
nrow(df_s4)

# S5 - Without GA and cell types (starts from DNAm + chem, not from main)
df_s5 <- with_dnam_chem[complete.cases(with_dnam_chem[, main_cov_s5]), ]
nrow(df_s5)

#==============================================================================#----
# SAVE                                                                         #----
#==============================================================================#----

# Check: N per subset
n_subsets <- data.frame(
  subset = c("main", "s1_no_compl", "s2_european", "s3_sg", "s4_fruitveg", "s5_no_ga_ct"),
  n      = c(nrow(df_main), nrow(df_s1), nrow(df_s2),
             nrow(df_s3),   nrow(df_s4), nrow(df_s5))
)
n_subsets
#subset   n
#1        main 312
#2 s1_no_compl 277
#3 s2_european 290
#4       s3_sg 312
#5 s4_fruitveg 279
#6 s5_no_ga_ct 312

saveRDS(df_main, paste0(outdir, "df_main.rds"))
saveRDS(df_s1,   paste0(outdir, "df_s1.rds"))
saveRDS(df_s2,   paste0(outdir, "df_s2.rds"))
saveRDS(df_s3,   paste0(outdir, "df_s3.rds"))
saveRDS(df_s4,   paste0(outdir, "df_s4.rds"))
saveRDS(df_s5,   paste0(outdir, "df_s5.rds"))
write.csv(n_subsets, paste0(outdir, "n_subsets.csv"), row.names = FALSE)

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----
