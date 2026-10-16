#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: Get from data needed for the analysis.

# [EXTENDED INFO] --------------------------------------------------------------
# List of covariables to have in the analysis(descriptive or models)
# 6 Maternal-Related: 
# - Maternal age at enrollment - Continuous - years
# - Maternal educational level - Categorical - University / Non-University
# - Maternal smoking during pregnancy - Categorical - Smoker / Non-smoker
# - Pre-pregnancy BMI - Continuous - kg/m2
# - Parity - Categorical - Nuliiparous / Multiparous 
# - Fruit and Vegetable Intake - Continutous - g/1000kcal/day

# 4 Child-Related: 
# - Child Ethnicity - Categorical - European / Latin American / Other
# - Child Sex - Categorical - Male / Female
# - Gestational Age at birth - Continuous - days 
# - Pregnancy Complications - Categorical -  Yes/No 

# 7 Methylation-Related: (From DNM data - 2_Get_DNAm_covars)
# - Estimate 6 Cell type Proportions - Continuous - Proportion % 
# - DNA contamination - Continuous - Mean log2 odds

# 1 Design / Study - Related: 
# - Specific Gravity - Continuous - mg/L

# CHEMICAL COMPUNDS DATA: 8 compounds
# 2 pesticide biomarkers:
# - 2,4-Dichlorophenoxyacetic acid (2,4-D) — herbicide  
# - Pyrethorids - specific metabolite 3-Phenoxybenzoic acid (3-PBA)  
# 6 organophosphate pesticide (OP) metabolites
# - 1 specific biomarker: TCPy (metabolite of chlorpyrifos) 
# - 5 non-specific OP metabolites: DMP, DMTP, DMDTP, DEP, DETP

# Values: 
# - 2,4-D — dichotomized
# - 3-PBA — continuous and dichotomized
# - TCPy — continuous
# - molar sums - 
# to check is in Values: Raw micromililiter/mol and 
# - Compounds have been log transformed to follow normal distribution. 
#-------------------------------------------------------------------------------


# [INPUT]: 
# - data.frame object nwith Rows = participants and Columns = Variables. N = xxx 

# [OUTPUT]: 

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
# Sociodemographic + Pest + Diet
infile1 <- "db/pest_cov_data/data_sepages_joana_v2_261013.rds" 

# Betas methylation
infile2 <- "db/plaDNAm_data/SEPAGES_ids+basenames.rds"
infile3 <- "db/plaDNAm_data/SEPAGES_20230615_Preprocessed.RData"
infile4 <- "db/plaDNAm_data/SEPAGES_20230615_Recommended_Samples_to_Remove.csv"

outfile <- "db/db_pest_pladnam_incomplete_cases.rds"

# Load data (variables + pesticides)
df <- readRDS(infile1) # 484 x 1490

# Get all available vars listed and explained                         
variable_labels <- lapply(df, function(x) {
  lbl <- attr(x, "label")
  if (is.null(lbl)) "No label" else paste(lbl, collapse = " ")
})
labels_df <- data.frame(Variable = names(variable_labels),
                        Label    = unlist(variable_labels),
                        stringsAsFactors = FALSE)
write.table(labels_df, file = "db/variable_labels_output.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Remove value-label classes (haven_labelled) -> plain numeric/character.
# Keeps the variable label (attr "label"), so no need to load haven.
for (v in names(df)) {
  if (any(grepl("labelled", class(df[[v]])))) {
    lbl <- attr(df[[v]], "label")
    df[[v]] <- as.vector(unclass(df[[v]]))
    attr(df[[v]], "label") <- lbl
  }
}

# Check: no labelled variables left
sum(sapply(df, function(x) any(grepl("labelled", class(x)))))   # 0


# Load data relates to methylation 
load(infile3)
ids <- readRDS(infile2) # 480 x 2
samples_to_rm <- read.csv(infile4) # 395 x 10 

#==============================================================================#----
# Sociodemographic                                                             #----
#==============================================================================#----

df_cov <- df[ ,c("ident",
                  "mo_age",      # maternal age
                  "mo_dipl",     # maternal educational level
                  "mo_tob_gr_anyt_yn1_n2", # maternal smoking during pregnancy 
                  "mo_bmi_bepr", # maternal pre pregnancy BMI kg/m2
                  "mo_par",      # maternal parity 
                  
                  "fa_ethnicity", # father ethnicity
                  "mo_ethnicity", # maternal ethnicity
                  "ch_sex",       # child sex
                  "po_gd_lmp",   # gestational duration
                  "mt1haa1_q40"  # pregnancy complications
                  )] 

#==============================================================================#----
# Diet - related.                                                              #----
#==============================================================================#----
# Keep women with a completed FFQ 
df_new <- df[!is.na(df$mo_alim_tot_sfalc_v2_t2), ]
nrow(df_new) # 387 - OK

# Vegetable item columns (Blandine's VegetablesrMEDpT2)
veg_cols <- c("mo_alim_qttj7_t2",  "mo_alim_qttj8_t2",
              "mo_alim_qttj9_t2",  "mo_alim_qttj10_t2",
              "mo_alim_qttj11_t2", "mo_alim_qttj12_t2",
              "mo_alim_qttj13_t2", "mo_alim_qttj14_t2",
              "mo_alim_qttj15_t2", "mo_alim_qttj16_t2",
              "mo_alim_qttj17_t2", "mo_alim_qttj18_t2",
              "mo_alim_qttj19_t2", "mo_alim_qttj20_t2",
              "mo_alim_qttj22_t2", "mo_alim_qttj23_t2",
              "mo_alim_qttj35_t2")

# Fruit item columns (Blandine's FruitrMEDpT2, includes qttj95 and qttj109)
fruit_cols <- c("mo_alim_qttj85_t2", "mo_alim_qttj86_t2",
                "mo_alim_qttj87_t2", "mo_alim_qttj88_t2",
                "mo_alim_qttj89_t2", "mo_alim_qttj90_t2",
                "mo_alim_qttj91_t2", "mo_alim_qttj92_t2",
                "mo_alim_qttj93_t2", "mo_alim_qttj94_t2",
                "mo_alim_qttj95_t2", "mo_alim_qttj109_t2")

# Convert to numeric (Blandine applies as.numeric to everything)
for (v in c(veg_cols, fruit_cols, "mo_alim_totv2_t2")) {
  df_new[[v]] <- as.numeric(df_new[[v]])
}

# Check: energy intake looks plausible (kcal/day) and has no zeros
summary(df_new$mo_alim_totv2_t2)
sum(df_new$mo_alim_totv2_t2 == 0, na.rm = TRUE)

# Raw daily intake (g/day)
df_new$veg_gday   <- rowSums(df_new[, veg_cols],   na.rm = TRUE)
df_new$fruit_gday <- rowSums(df_new[, fruit_cols], na.rm = TRUE)

# Energy-adjusted intake (g per 1000 kcal), as in Blandine's code
df_new$veg_g1000kcal   <- round(df_new$veg_gday   * 1000 / df_new$mo_alim_totv2_t2, 1)
df_new$fruit_g1000kcal <- round(df_new$fruit_gday * 1000 / df_new$mo_alim_totv2_t2, 1)

# Combined fruit + vegetable intake (g/1000kcal/day)
df_new$fv_g1000kcal <- df_new$veg_g1000kcal + df_new$fruit_g1000kcal
summary(df_new$fv_g1000kcal)

df_diet <- df_new[, c("ident", "veg_g1000kcal", "fruit_g1000kcal", "fv_g1000kcal")]
rm(df_new)

#==============================================================================#----
# Pesticides                                                                   #----
#==============================================================================#----

df_pest <- df[ ,c("ident",
                   "mo_pest_24d_cat_t2",   "mo_pest_24d_cat_t3",
                   "mo_pest_3pba_i_t2",    "mo_pest_3pba_i_t3",
                   "mo_pest_tcpy_i_t2",    "mo_pest_tcpy_i_t3",
                   
                   "mo_pest_dep_ms_i_t2",  "mo_pest_dep_ms_i_t3",
                   "mo_pest_dmp_ms_i_t2",  "mo_pest_dmp_ms_i_t3",
                   "mo_pest_dap_ms_i_t2",  "mo_pest_dap_ms_i_t3",
                   
                   "mo_pest_24d_raw_t2",   "mo_pest_24d_raw_t3",
                   "mo_pest_3pba_raw_t2",  "mo_pest_3pba_raw_t3",
                   "mo_pest_tcpy_raw_t2",  "mo_pest_tcpy_raw_t3",
                   "mo_pest_dmp_raw_t2",   "mo_pest_dmp_raw_t3",
                   "mo_pest_dmtp_raw_t2",  "mo_pest_dmtp_raw_t3",
                   "mo_pest_dmdtp_raw_t2", "mo_pest_dmdtp_raw_t3",
                   "mo_pest_dep_raw_t2",   "mo_pest_dep_raw_t3",
                   "mo_pest_detp_raw_t2",  "mo_pest_detp_raw_t3", 
                   "mo_pool_sg_T1",        "mo_pool_sg_T3"
)]

#==============================================================================#----
# Methylation - related                                                        #----
#==============================================================================#----

# 1. Get CT 
cells <- as.data.frame(processedOut$Omega) # 389 samples x 6 CT
cells$Basename <- rownames(cells)
rm(processedOut)

# 2. Get DNA contamination 
samples_to_rm <- read.csv(infile4) # N = 395 x 10 

# 3. Merge cells + Contamination  
meth_cov <- merge(cells, 
                  samples_to_rm[c("Basename","Meanlog2oddsContamination")], 
                  by = "Basename")

# 4. Get IDs & it 
meth_cov$ident <- ids$id[match(meth_cov$Basename, ids$Basename)]

rm(samples_to_rm,cells, ids)

#==============================================================================#----
# GET N & MERGE                                                                #----
#==============================================================================#----
# 6. Get all data together (keep incomplete cases)

# Check: one row per ident in each block
sum(duplicated(df_cov$ident)); sum(duplicated(df_diet$ident))
sum(duplicated(df_pest$ident)); sum(duplicated(meth_cov$ident))
sum(is.na(meth_cov$ident))  # should be 0

df_sel <- merge(df_cov, df_diet,  by = "ident", all.x = TRUE)
df_sel <- merge(df_sel, df_pest,  by = "ident", all.x = TRUE)
df_sel <- merge(df_sel, meth_cov, by = "ident", all.x = TRUE)
nrow(df_sel) # should be 484

# N with placental DNAm
sum(!is.na(df_sel$Basename)) # 389 

#==============================================================================#----
# CHECK & RECODE IF NEEDED                                                     #----
#==============================================================================#----

# Check Maternal age - Continuous - years - OK 
summary(df_sel$mo_age)                       

# Check Maternal education - 2 cat - University / Non-university
table(df_sel$mo_dipl, useNA = "ifany")       
#   2    3    4    5 <NA> 
#  29   54  127  272    2 
df_sel$mo_edu <- ifelse(df_sel$mo_dipl >= 3, "University", "Non-university")

# Check Maternal smoking - 2 cat - Smoker / Non-smoker - OK 
table(df_sel$mo_tob_gr_anyt_yn1_n2, useNA = "ifany")
#.   0    1 <NA> 
#  408   20   56 
df_sel$mo_smoke <- ifelse(df_sel$mo_tob_gr_anyt_yn1_n2 == 1, "Smoker",
                          ifelse(df_sel$mo_tob_gr_anyt_yn1_n2 == 0, "Non-smoker", NA))

# Check Pre-pregnancy BMI - Continuous - kg/m2 - OK 
summary(df_sel$mo_bmi_bepr)                 

# Check Parity - 2 cat - Nulliparous / Multiparous - OK 
table(df_sel$mo_par, useNA = "ifany")
df_sel$mo_parity <- ifelse(df_sel$mo_par == 0, "Nulliparous",
                           ifelse(df_sel$mo_par > 0, "Multiparous", NA))
table(df_sel$mo_parity, df_sel$mo_parity, useNA = "ifany")

# Check Fruit and vegetable intake - Continuous - g/1000kcal/day - OK 
summary(df_sel$fv_g1000kcal)

# Check Child ethnicity - 2 cat - Caucasian / Others (1 = Caucasian) - OK 
table(df_sel$mo_ethnicity, df_sel$fa_ethnicity, useNA = "ifany")
df_sel$ch_ethnicity <- ifelse(df_sel$mo_ethnicity == 1 & df_sel$fa_ethnicity == 1,
                              "Caucasian", "Others")

# Check Child sex - 2 cat - Male / Female (1 = Male, 2 = Female) - OK 
table(df_sel$ch_sex, useNA = "ifany")

# Check Gestational age at birth - Continuous - days
summary(df_sel$po_gd_lmp)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 28.57   39.00   40.00   39.71   40.71   43.43       4 
df_sel$po_gd_lmp_days <- df_sel$po_gd_lmp*7

# Check Pregnancy complications - 2 cat - Yes / No - OK 
table(df_sel$mt1haa1_q40, useNA = "ifany")
df_sel$preg_compl <- ifelse(df_sel$mt1haa1_q40 == 1, "Yes",
                            ifelse(df_sel$mt1haa1_q40 %in% c(0, 2), "No", NA))

# Check: pesticides log2-transformed
# Transformed (_i_, _ms_i_) should have negative values; raw should be all >= 0
summary(df_sel[, c("mo_pest_tcpy_raw_t2","mo_pest_tcpy_i_t2")])
# Min.   : 0.3583     Min.   : 0.3583  
# 1st Qu.: 2.1280     1st Qu.: 2.1280  
# Median : 3.2199     Median : 3.2199  
# Mean   : 4.2788     Mean   : 4.2788  
# 3rd Qu.: 4.7526     3rd Qu.: 4.7526  
# Max.   :52.6753     Max.   :52.6753  
# NA's   :5           NA's   :5   

# Transform pest to log2
df_sel$mo_pest_3pba_i_log2_t2 <- log2(df_sel$mo_pest_3pba_i_t2)
df_sel$mo_pest_3pba_i_log2_t3 <- log2(df_sel$mo_pest_3pba_i_t3)
df_sel$mo_pest_tcpy_i_log2_t2 <- log2(df_sel$mo_pest_tcpy_i_t2)
df_sel$mo_pest_tcpy_i_log2_t3 <- log2(df_sel$mo_pest_tcpy_i_t3)

# Transform 2,4-D to detected / not detected
# 1 = <LOD -> 0 (not detected); 2 = LOD-LOQ, 3 = >LOQ -> 1 (detected)
df_sel$mo_pest_24d_bi_t2 <- ifelse(df_sel$mo_pest_24d_cat_t2 >= 2, 1, 0)
df_sel$mo_pest_24d_bi_t3 <- ifelse(df_sel$mo_pest_24d_cat_t3 >= 2, 1, 0)

# Check: molar sums in umol/L - OK 
summary(df_sel$mo_pest_dep_ms_i_t2)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.00142 0.01254 0.02001 0.03002 0.03033 1.75660       6 

# Transform to log2
df_sel$mo_pest_dmp_ms_i_log2_t2 <- log2(df_sel$mo_pest_dmp_ms_i_t2)
df_sel$mo_pest_dmp_ms_i_log2_t3 <- log2(df_sel$mo_pest_dmp_ms_i_t3)
df_sel$mo_pest_dep_ms_i_log2_t2 <- log2(df_sel$mo_pest_dep_ms_i_t2)
df_sel$mo_pest_dep_ms_i_log2_t3 <- log2(df_sel$mo_pest_dep_ms_i_t3)
df_sel$mo_pest_dap_ms_i_log2_t2 <- log2(df_sel$mo_pest_dap_ms_i_t2)
df_sel$mo_pest_dap_ms_i_log2_t3 <- log2(df_sel$mo_pest_dap_ms_i_t3)

saveRDS(df_sel, outfile)

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----
