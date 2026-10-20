#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: Join data sources to obtain all data needed for the analysis in 
# a single db. 
# Mother-related (6): Age at enrolment, educational level, smoking during pregnancy,
#                       pre-pregnancy BMI, parity, fruit and vegetables intake.
# Child-related (4): Ethnicity, Sex, Gestational Age at Birth, Pregnancy Complications.
# Methylation-Related (7): stimate 6 Cell type Proportions, DNA contamination
# Design-Related (1): Specific Gravity

# [INPUT]: 
# - infile1: data.frame with Rows = participants and Columns = Variables. Dim: # 484 x 1490 
# - infile2: data.frame with Rows = participants with methylation, Col = Ids 480 x 2
# - infile3: processedOut named object. GRset Object with methylation metadata. 
# - infile4: data.frame data/frame with tehcnincal related variables # 395 x 10 

# [OUTPUT]: 
# "db/db_pest_pladnam_incomplete_cases.rds"

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


# Load data related to methylation 
load(infile3)
ids <- readRDS(infile2) # 480 x 2
samples_to_rm <- read.csv(infile4) # 395 x 10 

#==============================================================================#----
# Get Fruit Veg Intake                                                         #----
#==============================================================================#----
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


# Raw daily intake (g/day)
df$veg_gday   <- rowSums(df[, veg_cols],   na.rm = TRUE)
df$fruit_gday <- rowSums(df[, fruit_cols], na.rm = TRUE)

# Energy-adjusted intake (g per 1000 kcal), as in Blandine's code
df$veg_g1000kcal   <- round(df$veg_gday   * 1000 / df$mo_alim_totv2_t2, 1)
df$fruit_g1000kcal <- round(df$fruit_gday * 1000 / df$mo_alim_totv2_t2, 1)

# Combined fruit + vegetable intake (g/1000kcal/day)
df$fv_g1000kcal <- df$veg_g1000kcal + df$fruit_g1000kcal


#==============================================================================#----
# Sociodemographic                                                             #----
#==============================================================================#----

df_cov <- df[ ,c("ident",
                  "mo_age",      # maternal age
                  "mo_dipl",     # maternal educational level
                  "mo_tob_gr_anyt_yn1_n2", # maternal smoking during pregnancy 
                  "mo_bmi_bepr", # maternal pre pregnancy BMI kg/m2
                  "mo_par",      # maternal parity 
                 "fv_g1000kcal", # maternal fruit and vegetables intake
                  
                  "fa_ethnicity", # father ethnicity
                  "mo_ethnicity", # maternal ethnicity
                  "ch_sex",       # child sex
                  "po_gd_lmp",   # gestational duration
                  "mt1haa1_q40"  # pregnancy complications
                  )] 
summary(is.na(df_cov))

# From the 484 mothers; we have missing data on:
# - 2 maternal education
# - 56 tabacco smoke
# - 4 bmi 
# - 97 fruit & vegetables intake
# - 74 father ethnicity 
# - 21 mother ethnicity 
# - 4 child sex
# - 4 gestational age 
# - 8 pregnancy complications 
   


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

# From 484 recruited mothers, the one's that have pesticide exposure data:
# 2,4-D; 3-BPA; TCPy -->  479 for T2 and 457 T3
# DEP; DMP; DAP --> 478 for T2 and 455 T3
# 2,4-D; 3-BPA; TCPy | T2: 5 | T3: 27
# DEP; DMP; DAP | T2: 6 | T3: 25
# SG | T2: 5 | T3: 25


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
# Get all data together (keep incomplete cases)
df_sel <- merge(df_cov, df_pest,  by = "ident", all.x = TRUE)
df_sel <- merge(df_sel, meth_cov, by = "ident", all.x = TRUE)
nrow(df_sel) # should be 484
saveRDS(df_sel, outfile)


#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----
