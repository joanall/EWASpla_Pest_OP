#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: April - October 2026                                                                              
#==============================================================================#
# [DESCRIPTION]: Rename and recode BiSC data following the analysis plan, and
# build the analysis subsets. All subsets require:
# 1) Placental DNAm (outcome)
# 2) Pesticide data at the time point(s) (exposure)
# 3) Main covariates: maternal age, education, smoking, GA, child ethnicity,
#    child sex, delivery hospital, cell types and DNA contamination
#
# BiSC samples differ a lot between time points, so three sample sets are built:
#  - T2 only : women with exposure at t2 (12w pool)   -> db/subsets/*_t2.rds
#  - T3 only : women with exposure at t3 (32w pool)   -> db/subsets/*_t3.rds
#  - Both    : women with exposure at t2 AND t3       -> db/subsets_2tp/*.rds
#              (same file names as SEPAGES, so the same model scripts apply)
#
# Subsets in each set:
#  - Main
#  - S1: Excluding pregnancy complications
#  - S2: Excluding non-Europeans
#  - S3: Adjusting for specific gravity (complete SG at the time point)
#  - S4: Adjusting for fruit and vegetable intake (complete FFQ)
#  - S5: Without GA and cell types
#
# Naming: BiSC 12w pool = t2, 32w pool = t3 (harmonized with SEPAGES).
#
# [INPUT]:  db/0_Get_Data/final_df_incomplete_cases.rds
# [OUTPUT]: db/subsets/df_main_t2.rds ... df_s5_t3.rds
#           db/subsets_2tp/df_main.rds ... df_s5.rds
#           db/subsets/n_subsets.csv

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

infile     <- "db/0_Get_Data/final_df_incomplete_cases.rds"
outdir     <- "db/subsets/"
outdir_2tp <- "db/subsets_2tp/"
dir.create(outdir,     recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_2tp, recursive = TRUE, showWarnings = FALSE)

df <- readRDS(infile)
dim(df) # 1080 x 57

#==============================================================================#----
# RENAME                                                                       #----
#==============================================================================#----
# Untouched copy, to check recodes later with table(df_orig$x, df$y)
df_orig <- df

# Old -> new names, one pair per line.
# Variables not listed keep their name (Basename, cell types,
# Meanlog2oddsContamination).
rename_map <- matrix(c(
  # Covariates
  "ID_mother",             "mo_id",
  "age_12w_m",             "mo_age",
  "educ_level_m_3cat",     "mo_edu",
  "smoke_sust_m",          "mo_smoke",
  "bmi_m_12w",             "mo_bmi",
  "fruit_veg_intake",      "mo_fv_g1000kcal",
  "gestage_0y_c_days",     "gest_age",          # already in DAYS in BiSC
  "any_compl_m_0y_2c",     "preg_compl",
  "parity_m_2cat",         "mo_parity",         # already 2 categories
  "sex_0y_c",              "ch_sex",
  "ethnicity_c_4cat",      "ch_eth",
  "hospital_delivery",     "hospital",
  "dna_contamination",    "Meanlog2oddsContamination",
  
  # Exposures for models
  "cc.2,4-D_m_12w_dic",    "pest_24d_bi_t2",    # already detected / not detected
  "cc.2,4-D_m_32w_dic",    "pest_24d_bi_t3",
  "cc.3-PBA_m_12w_log2i",  "pest_3pba_ilog2_t2", # already imputed + log2
  "cc.3-PBA_m_32w_log2i",  "pest_3pba_ilog2_t3",
  "cc.TCPy_m_12w_log2i",   "pest_tcpy_ilog2_t2",
  "cc.TCPy_m_32w_log2i",   "pest_tcpy_ilog2_t3",
  "sum_DEP_12w",           "pest_dep_ms_i_t2",  # molar sums, natural scale
  "sum_DEP_32w",           "pest_dep_ms_i_t3",
  "sum_DMP_12w",           "pest_dmp_ms_i_t2",
  "sum_DMP_32w",           "pest_dmp_ms_i_t3",
  "sum_DAP_12w",           "pest_dap_ms_i_t2",
  "sum_DAP_32w",           "pest_dap_ms_i_t3",
  
  # Raw concentrations (descriptives)
  "cc.2,4-D_m_12w_raw",    "pest_24d_raw_t2",
  "cc.2,4-D_m_32w_raw",    "pest_24d_raw_t3",
  "cc.3-PBA_m_12w_raw",    "pest_3pba_raw_t2",
  "cc.3-PBA_m_32w_raw",    "pest_3pba_raw_t3",
  "cc.TCPy_m_12w_raw",     "pest_tcpy_raw_t2",
  "cc.TCPy_m_32w_raw",     "pest_tcpy_raw_t3",
  "cc.DMP_m_12w_raw",      "pest_dmp_raw_t2",
  "cc.DMP_m_32w_raw",      "pest_dmp_raw_t3",
  "cc.DMTP_m_12w_raw",     "pest_dmtp_raw_t2",
  "cc.DMTP_m_32w_raw",     "pest_dmtp_raw_t3",
  "cc.DMDTP_m_12w_raw",    "pest_dmdtp_raw_t2", # may not exist in BiSC
  "cc.DMDTP_m_32w_raw",    "pest_dmdtp_raw_t3",
  "cc.DEP_m_12w_raw",      "pest_dep_raw_t2",
  "cc.DEP_m_32w_raw",      "pest_dep_raw_t3",
  "cc.DETP_m_12w_raw",     "pest_detp_raw_t2",
  "cc.DETP_m_32w_raw",     "pest_detp_raw_t3",
  
  # Specific gravity
  "ur_pool_grav_m_12w",    "pool_sg_t2",
  "ur_pool_grav_m_32w",    "pool_sg_t3"
), ncol = 2, byrow = TRUE)

# Check: old names not found in the data (only DMDTP expected)
setdiff(rename_map[, 1], names(df))
rename_map <- rename_map[rename_map[, 1] %in% names(df), ]

old_names <- rename_map[, 1]
new_names <- rename_map[, 2]

stopifnot(!any(duplicated(old_names)))
stopifnot(!any(duplicated(new_names)))
stopifnot(!any(new_names %in% names(df)))

# Rename
names(df)[match(old_names, names(df))] <- new_names

# Check
setdiff(names(df_orig), names(df))   # old names that changed
setdiff(names(df), names(df_orig))   # their new names

#==============================================================================#----
# RECODE                                                                       #----
#==============================================================================#----
# After each recode: same number of NA as before (no label silently lost)

#------------------------------------------------------------------------------#
# Covariates                                                                   #
#------------------------------------------------------------------------------#

# Maternal age - Continuous - years
df$mo_age <- as.numeric(df$mo_age)
summary(df$mo_age)

# Maternal education - 2 cat - University / Non-university
# CHECK: Secondary grouped with Non-university (as in the analysis plan)
table(df$mo_edu, useNA = "ifany")
df$mo_edu <- as.character(df$mo_edu)
df$mo_edu <- ifelse(df$mo_edu == "University", "University",
                    ifelse(df$mo_edu %in% c("Secondary", "Primary or less"), "Non-university", NA))
df$mo_edu <- factor(df$mo_edu, levels = c("Non-university", "University"))
table(df_orig$educ_level_m_3cat, df$mo_edu, useNA = "ifany")
stopifnot(sum(is.na(df$mo_edu)) == sum(is.na(df_orig$educ_level_m_3cat)))

# Maternal smoking - 2 cat - Smoker / Non-smoker
# sustained smoker -> Smoker; non-sustained smoker, non-smoker -> Non-smoker
table(df$mo_smoke, useNA = "ifany")
df$mo_smoke <- as.character(df$mo_smoke)
df$mo_smoke <- ifelse(df$mo_smoke %in% c("sustained smoker", "smoker"), "Smoker",
                      ifelse(df$mo_smoke %in% c("non-sustained smoker", "non-smoker"), "Non-smoker", NA))
df$mo_smoke <- factor(df$mo_smoke, levels = c("Non-smoker", "Smoker"))
table(df_orig$smoke_sust_m, df$mo_smoke, useNA = "ifany")
stopifnot(sum(is.na(df$mo_smoke)) == sum(is.na(df_orig$smoke_sust_m)))

# Pre-pregnancy BMI - Continuous - kg/m2
df$mo_bmi <- as.numeric(df$mo_bmi)
summary(df$mo_bmi)

# Parity - 2 cat - Nulliparous / Multiparous
# CHECK the labels printed below and adjust the two strings if needed
table(df$mo_parity, useNA = "ifany")
df$mo_parity <- as.character(df$mo_parity)
df$mo_parity <- ifelse(tolower(df$mo_parity) %in% c("nulliparous", "0"), "Nulliparous",
                       ifelse(tolower(df$mo_parity) %in% c("multiparous", ">=1", "1"), "Multiparous", NA))
df$mo_parity <- factor(df$mo_parity, levels = c("Nulliparous", "Multiparous"))
table(df_orig$parity_m_2cat, df$mo_parity, useNA = "ifany")
stopifnot(sum(is.na(df$mo_parity)) == sum(is.na(df_orig$parity_m_2cat)))

# Fruit and vegetable intake - Continuous - g/1000kcal/day
df$mo_fv_g1000kcal <- as.numeric(df$mo_fv_g1000kcal)
summary(df$mo_fv_g1000kcal)

# Child ethnicity - 2 cat - European / Non-European
# eur -> European; eur_lat_amr, lat_amr, other -> Non-European
table(df$ch_eth, useNA = "ifany")
df$ch_eth <- as.character(df$ch_eth)
df$ch_eth <- ifelse(df$ch_eth == "eur", "European",
                    ifelse(df$ch_eth %in% c("eur_lat_amr", "lat_amr", "other"), "Non-European", NA))
df$ch_eth <- factor(df$ch_eth, levels = c("European", "Non-European"))
table(df_orig$ethnicity_c_4cat, df$ch_eth, useNA = "ifany")
stopifnot(sum(is.na(df$ch_eth)) == sum(is.na(df_orig$ethnicity_c_4cat)))

# Child sex - 2 cat - Male / Female
# BiSC labels: "Boy" -> Male, "Girl" -> Female
table(df$ch_sex, useNA = "ifany")
df$ch_sex <- as.character(df$ch_sex)
df$ch_sex <- ifelse(df$ch_sex == "Boy",  "Male",
                    ifelse(df$ch_sex == "Girl", "Female", NA))
df$ch_sex <- factor(df$ch_sex, levels = c("Male", "Female"))
table(df_orig$sex_0y_c, df$ch_sex, useNA = "ifany")
stopifnot(sum(is.na(df$ch_sex)) == sum(is.na(df_orig$sex_0y_c)))


# Gestational age at birth - Continuous - days (already in days in BiSC)
df$gest_age <- as.numeric(df$gest_age)
summary(df$gest_age)            # expected range ~ 200-300 days

# Pregnancy complications - 2 cat - Yes / No
table(df$preg_compl, useNA = "ifany")
df$preg_compl <- as.character(df$preg_compl)
df$preg_compl <- ifelse(tolower(df$preg_compl) == "yes", "Yes",
                        ifelse(tolower(df$preg_compl) == "no",  "No", NA))
df$preg_compl <- factor(df$preg_compl, levels = c("No", "Yes"))
table(df_orig$any_compl_m_0y_2c, df$preg_compl, useNA = "ifany")
stopifnot(sum(is.na(df$preg_compl)) == sum(is.na(df_orig$any_compl_m_0y_2c)))

# Delivery hospital - 4 cat - HSJD / Maternitat / HSPau / Other
# Home and Other births grouped as Other
table(df$hospital, useNA = "ifany")
df$hospital <- as.character(df$hospital)
df$hospital <- ifelse(df$hospital %in% c("Home", "Other"), "Other", df$hospital)
df$hospital <- factor(df$hospital, levels = c("HSjD", "Maternitat", "HSPau", "Other"))
table(df_orig$hospital_delivery, df$hospital, useNA = "ifany")
stopifnot(sum(is.na(df$hospital)) == sum(is.na(df_orig$hospital_delivery)))

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
raw_vars                                              # 14 (16 if DMDTP exists)
for (v in raw_vars) df[[v]] <- as.numeric(df[[v]])
summary(df[, raw_vars])                               # all >= 0

#------------------------------------------------------------------------------#
# Pesticides - exposures for models                                            #
#------------------------------------------------------------------------------#

# 2,4-D - Binary numeric - 0 = Not detected (reference), 1 = Detected
# BiSC labels: "undetected" -> 0, "detected" -> 1
# Kept numeric (not factor) so the model coefficient keeps the variable name
table(df$pest_24d_bi_t2, useNA = "ifany")
#   detected undetected       <NA>
#        266        429        385

df$pest_24d_bi_t2 <- as.character(df$pest_24d_bi_t2)
df$pest_24d_bi_t2 <- ifelse(df$pest_24d_bi_t2 == "detected", 1,
                            ifelse(df$pest_24d_bi_t2 == "undetected", 0, NA))
table(df_orig[["cc.2,4-D_m_12w_dic"]], df$pest_24d_bi_t2, useNA = "ifany")
stopifnot(sum(is.na(df$pest_24d_bi_t2)) == sum(is.na(df_orig[["cc.2,4-D_m_12w_dic"]])))

table(df$pest_24d_bi_t3, useNA = "ifany")
df$pest_24d_bi_t3 <- as.character(df$pest_24d_bi_t3)
df$pest_24d_bi_t3 <- ifelse(df$pest_24d_bi_t3 == "detected", 1,
                            ifelse(df$pest_24d_bi_t3 == "undetected", 0, NA))
table(df_orig[["cc.2,4-D_m_32w_dic"]], df$pest_24d_bi_t3, useNA = "ifany")
stopifnot(sum(is.na(df$pest_24d_bi_t3)) == sum(is.na(df_orig[["cc.2,4-D_m_32w_dic"]])))

class(df$pest_24d_bi_t2)   # "numeric"
class(df$pest_24d_bi_t3)   # "numeric"

# 3-PBA and TCPy - already imputed and log2-transformed in BiSC
df$pest_3pba_ilog2_t2 <- as.numeric(df$pest_3pba_ilog2_t2)
df$pest_3pba_ilog2_t3 <- as.numeric(df$pest_3pba_ilog2_t3)
df$pest_tcpy_ilog2_t2 <- as.numeric(df$pest_tcpy_ilog2_t2)
df$pest_tcpy_ilog2_t3 <- as.numeric(df$pest_tcpy_ilog2_t3)
summary(df[, c("pest_3pba_ilog2_t2", "pest_3pba_ilog2_t3",
               "pest_tcpy_ilog2_t2", "pest_tcpy_ilog2_t3")])   # negatives expected

# Molar sums (umol/L) - log2-transform here
ms_vars <- c("pest_dep_ms_i_t2", "pest_dep_ms_i_t3",
             "pest_dmp_ms_i_t2", "pest_dmp_ms_i_t3",
             "pest_dap_ms_i_t2", "pest_dap_ms_i_t3")
for (v in ms_vars) df[[v]] <- as.numeric(df[[v]])
summary(df[, ms_vars])
# Check: all > 0 before log2 (a 0 would give -Inf)
n_nonpos <- sapply(df[, ms_vars], function(x) sum(x <= 0, na.rm = TRUE))
n_nonpos
if (any(n_nonpos > 0)) stop("Molar sums with values <= 0: ",
                            paste(names(n_nonpos)[n_nonpos > 0], collapse = ", "))

df$pest_dep_ms_ilog2_t2 <- log2(df$pest_dep_ms_i_t2)
df$pest_dep_ms_ilog2_t3 <- log2(df$pest_dep_ms_i_t3)
df$pest_dmp_ms_ilog2_t2 <- log2(df$pest_dmp_ms_i_t2)
df$pest_dmp_ms_ilog2_t3 <- log2(df$pest_dmp_ms_i_t3)
df$pest_dap_ms_ilog2_t2 <- log2(df$pest_dap_ms_i_t2)
df$pest_dap_ms_ilog2_t3 <- log2(df$pest_dap_ms_i_t3)

#------------------------------------------------------------------------------#
# Final check: class of every analysis variable                                #
#------------------------------------------------------------------------------#
expo_t2 <- c("pest_24d_bi_t2",       "pest_3pba_ilog2_t2",   "pest_tcpy_ilog2_t2",
             "pest_dep_ms_ilog2_t2", "pest_dmp_ms_ilog2_t2", "pest_dap_ms_ilog2_t2")
expo_t3 <- c("pest_24d_bi_t3",       "pest_3pba_ilog2_t3",   "pest_tcpy_ilog2_t3",
             "pest_dep_ms_ilog2_t3", "pest_dmp_ms_ilog2_t3", "pest_dap_ms_ilog2_t3")

sapply(df[, c("mo_age", "mo_edu", "mo_smoke", "mo_bmi", "mo_parity",
              "mo_fv_g1000kcal", "ch_eth", "ch_sex", "gest_age", "preg_compl",
              "hospital", expo_t2, expo_t3)], class)
# Exposures must all be numeric (models fail on factors)
stopifnot(all(sapply(df[, c(expo_t2, expo_t3)], is.numeric)))
sapply(df[, c(expo_t2, expo_t3)], function(x) sum(is.infinite(x)))   # all 0

#==============================================================================#----
# COVARIATE SETS                                                               #----
#==============================================================================#----

# Main model
main_cov <- c("mo_age", "mo_edu", "mo_smoke", "gest_age", "ch_eth", "ch_sex",
              "hospital", ct_vars, "Meanlog2oddsContamination")

# S5: main minus GA and cell types
main_cov_s5 <- c("mo_age", "mo_edu", "mo_smoke", "ch_eth", "ch_sex",
                 "hospital", "Meanlog2oddsContamination")

setdiff(c(main_cov, main_cov_s5), colnames(df))   # should be character(0)

#==============================================================================#----
# AVERAGES                                                                     #----
#==============================================================================#----
# Pregnancy-average exposure (main subset: all women have both time points)
# - Continuous: mean of the two log2 values (= log2 of the geometric mean)
# - 2,4-D: number of time points detected (0, 1, 2)

# 1. TCPy, 3-PBA - log2
df$pest_tcpy_ilog2_avg <- rowMeans(df[, c("pest_tcpy_ilog2_t2", "pest_tcpy_ilog2_t3")], na.rm = TRUE)
df$pest_3pba_ilog2_avg <- rowMeans(df[, c("pest_3pba_ilog2_t2", "pest_3pba_ilog2_t3")], na.rm = TRUE)

# 2. OP molar sums - log2
df$pest_dep_ms_ilog2_avg <- rowMeans(df[, c("pest_dep_ms_ilog2_t2", "pest_dep_ms_ilog2_t3")], na.rm = TRUE)
df$pest_dmp_ms_ilog2_avg <- rowMeans(df[, c("pest_dmp_ms_ilog2_t2", "pest_dmp_ms_ilog2_t3")], na.rm = TRUE)
df$pest_dap_ms_ilog2_avg <- rowMeans(df[, c("pest_dap_ms_ilog2_t2", "pest_dap_ms_ilog2_t3")], na.rm = TRUE)

# 2,4-D average: detected at >= 1 available time point (0/1; NA if both missing)
df$pest_24d_any_avg <- ifelse(is.na(df$pest_24d_bi_t2) & is.na(df$pest_24d_bi_t3), NA,
                              as.numeric(rowSums(df[, c("pest_24d_bi_t2", "pest_24d_bi_t3")],
                                                 na.rm = TRUE) >= 1))

# 4. Specific gravity - mean of both pools
df$pool_sg_avg <- rowMeans(df[, c("pool_sg_t2", "pool_sg_t3")], na.rm = TRUE)

# rowMeans of two NA gives NaN -> NA
avg_vars <- c("pest_tcpy_ilog2_avg", "pest_3pba_ilog2_avg", "pest_dep_ms_ilog2_avg",
              "pest_dmp_ms_ilog2_avg", "pest_dap_ms_ilog2_avg", "pool_sg_avg")
for (v in avg_vars) df[[v]][is.nan(df[[v]])] <- NA

# Check: averages lie between the two time points
summary(df[, c(avg_vars, "pest_24d_any_avg")])

#==============================================================================#----
# DNAm + EXPOSURE                                                              #----
#==============================================================================#----

# 1. Participants with placental DNAm
with_dnam <- df[!is.na(df$Meanlog2oddsContamination), ]
nrow(with_dnam) # 565

# 2. Participants with all model exposures at t2, at t3, and at both
with_dnam_chem_t2  <- with_dnam[complete.cases(with_dnam[, expo_t2]), ]
with_dnam_chem_t3  <- with_dnam[complete.cases(with_dnam[, expo_t3]), ]
with_dnam_chem_2tp <- with_dnam[complete.cases(with_dnam[, c(expo_t2, expo_t3)]), ]

nrow(with_dnam_chem_t2) # 384
nrow(with_dnam_chem_t3) # 442
nrow(with_dnam_chem_2tp) # 316

colSums(is.na(with_dnam_chem_t2[, main_cov]))
# We have: 6 Na in smoking
colSums(is.na(with_dnam_chem_t3[, main_cov]))
# We have: 5 Na in smoking

#==============================================================================#----
# SUBSETS - T2 ONLY                                                            #----
#==============================================================================#----
df_main_t2 <- with_dnam_chem_t2[complete.cases(with_dnam_chem_t2[, main_cov]), ] # 378
df_s1_t2   <- df_main_t2[df_main_t2$preg_compl %in% "No", ]
df_s2_t2   <- df_main_t2[df_main_t2$ch_eth %in% "European", ]
df_s3_t2   <- df_main_t2[!is.na(df_main_t2$pool_sg_t2), ]
df_s4_t2   <- df_main_t2[!is.na(df_main_t2$mo_fv_g1000kcal), ]
df_s5_t2   <- with_dnam_chem_t2[complete.cases(with_dnam_chem_t2[, main_cov_s5]), ]

#==============================================================================#----
# SUBSETS - T3 ONLY                                                            #----
#==============================================================================#----
df_main_t3 <- with_dnam_chem_t3[complete.cases(with_dnam_chem_t3[, main_cov]), ]
df_s1_t3   <- df_main_t3[df_main_t3$preg_compl %in% "No", ]
df_s2_t3   <- df_main_t3[df_main_t3$ch_eth %in% "European", ]
df_s3_t3   <- df_main_t3[!is.na(df_main_t3$pool_sg_t3), ]
df_s4_t3   <- df_main_t3[!is.na(df_main_t3$mo_fv_g1000kcal), ]
df_s5_t3   <- with_dnam_chem_t3[complete.cases(with_dnam_chem_t3[, main_cov_s5]), ]

#==============================================================================#----
# SUBSETS - BOTH TIME POINTS                                                   #----
#==============================================================================#----
df_main_2tp <- with_dnam_chem_2tp[complete.cases(with_dnam_chem_2tp[, main_cov]), ]
df_s1_2tp   <- df_main_2tp[df_main_2tp$preg_compl %in% "No", ]
df_s2_2tp   <- df_main_2tp[df_main_2tp$ch_eth %in% "European", ]
df_s3_2tp   <- df_main_2tp[!is.na(df_main_2tp$pool_sg_t2) & !is.na(df_main_2tp$pool_sg_t3), ]
df_s4_2tp   <- df_main_2tp[!is.na(df_main_2tp$mo_fv_g1000kcal), ]
df_s5_2tp   <- with_dnam_chem_2tp[complete.cases(with_dnam_chem_2tp[, main_cov_s5]), ]

# Check: both-time-point main = women in both single-time-point mains
identical(sort(df_main_2tp$mo_id),
          sort(intersect(df_main_t2$mo_id, df_main_t3$mo_id)))   # TRUE

# Check: hospital levels per subset (sparse "Other" can be a problem in S2)
table(df_main_t2$hospital); table(df_main_t3$hospital); table(df_main_2tp$hospital)

#==============================================================================#----
# SUBSETS - AVERAGES                                                           #----
#==============================================================================#----

# Women with all 6 average exposures available (>= 1 time point each)
expo_avg <- c("pest_24d_any_avg", "pest_3pba_ilog2_avg", "pest_tcpy_ilog2_avg",
              "pest_dep_ms_ilog2_avg", "pest_dmp_ms_ilog2_avg", "pest_dap_ms_ilog2_avg")
stopifnot(all(sapply(df[, expo_avg], is.numeric)))

with_dnam_chem_avg <- with_dnam[complete.cases(with_dnam[, expo_avg]), ]
nrow(with_dnam_chem_avg)

# Number of urine samples per woman (to report)
with_dnam_chem_avg$n_samples_pest <- rowSums(!is.na(with_dnam_chem_avg[, c("pest_tcpy_ilog2_t2",
                                                                           "pest_tcpy_ilog2_t3")]))
table(with_dnam_chem_avg$n_samples_pest)

df_main_avg <- with_dnam_chem_avg[complete.cases(with_dnam_chem_avg[, main_cov]), ]
df_s1_avg   <- df_main_avg[df_main_avg$preg_compl %in% "No", ]
df_s2_avg   <- df_main_avg[df_main_avg$ch_eth %in% "European", ]
df_s3_avg   <- df_main_avg[!is.na(df_main_avg$pool_sg_avg), ]
df_s4_avg   <- df_main_avg[!is.na(df_main_avg$mo_fv_g1000kcal), ]
df_s5_avg   <- with_dnam_chem_avg[complete.cases(with_dnam_chem_avg[, main_cov_s5]), ]

# Check: average main contains the 2tp main, plus women with one sample
all(df_main_2tp$mo_id %in% df_main_avg$mo_id)   # TRUE
table(df_main_avg$n_samples_pest)


#==============================================================================#----
# SAVE                                                                         #----
#==============================================================================#----

# Check: N per subset
n_subsets <- data.frame(
  subset = c("main", "s1_no_compl", "s2_european", "s3_sg", "s4_fruitveg", "s5_no_ga_ct"),
  n_t2   = c(nrow(df_main_t2), nrow(df_s1_t2), nrow(df_s2_t2),
             nrow(df_s3_t2),   nrow(df_s4_t2), nrow(df_s5_t2)),
  n_t3   = c(nrow(df_main_t3), nrow(df_s1_t3), nrow(df_s2_t3),
             nrow(df_s3_t3),   nrow(df_s4_t3), nrow(df_s5_t3)),
  n_2tp  = c(nrow(df_main_2tp), nrow(df_s1_2tp), nrow(df_s2_2tp),
             nrow(df_s3_2tp),   nrow(df_s4_2tp), nrow(df_s5_2tp)),
  n_avg  = c(nrow(df_main_avg), nrow(df_s1_avg), nrow(df_s2_avg),
             nrow(df_s3_avg),   nrow(df_s4_avg), nrow(df_s5_avg))
)
n_subsets

# T2 only and T3 only
saveRDS(df_main_t2, paste0(outdir, "df_main_t2.rds"))
saveRDS(df_s1_t2,   paste0(outdir, "df_s1_t2.rds"))
saveRDS(df_s2_t2,   paste0(outdir, "df_s2_t2.rds"))
saveRDS(df_s3_t2,   paste0(outdir, "df_s3_t2.rds"))
saveRDS(df_s4_t2,   paste0(outdir, "df_s4_t2.rds"))
saveRDS(df_s5_t2,   paste0(outdir, "df_s5_t2.rds"))

saveRDS(df_main_t3, paste0(outdir, "df_main_t3.rds"))
saveRDS(df_s1_t3,   paste0(outdir, "df_s1_t3.rds"))
saveRDS(df_s2_t3,   paste0(outdir, "df_s2_t3.rds"))
saveRDS(df_s3_t3,   paste0(outdir, "df_s3_t3.rds"))
saveRDS(df_s4_t3,   paste0(outdir, "df_s4_t3.rds"))
saveRDS(df_s5_t3,   paste0(outdir, "df_s5_t3.rds"))

# Both time points (same file names as SEPAGES)
saveRDS(df_main_2tp, paste0(outdir_2tp, "df_main.rds"))
saveRDS(df_s1_2tp,   paste0(outdir_2tp, "df_s1.rds"))
saveRDS(df_s2_2tp,   paste0(outdir_2tp, "df_s2.rds"))
saveRDS(df_s3_2tp,   paste0(outdir_2tp, "df_s3.rds"))
saveRDS(df_s4_2tp,   paste0(outdir_2tp, "df_s4.rds"))
saveRDS(df_s5_2tp,   paste0(outdir_2tp, "df_s5.rds"))

# Averages, all available data (bisc ... avg tp)
saveRDS(df_main_avg, paste0(outdir, "df_main_avg.rds"))
saveRDS(df_s1_avg,   paste0(outdir, "df_s1_avg.rds"))
saveRDS(df_s2_avg,   paste0(outdir, "df_s2_avg.rds"))
saveRDS(df_s3_avg,   paste0(outdir, "df_s3_avg.rds"))
saveRDS(df_s4_avg,   paste0(outdir, "df_s4_avg.rds"))
saveRDS(df_s5_avg,   paste0(outdir, "df_s5_avg.rds"))

write.csv(n_subsets, paste0(outdir, "n_subsets.csv"), row.names = FALSE)

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----