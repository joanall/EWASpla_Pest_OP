#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: April - October 2026                                                                              
#==============================================================================#

# [DESCRIPTION]: Get from "general" df covariables needed for this analysis, either 
# for descriptives or for modelling. 

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

# 2 Design / Study - Realted: 
# - Study Center - Cateogrcal - HSJ / HC / HSCSP
# - Specific Gravity - Continuous - mg/L (From Chemicals db - 1_1_Get_OP_Data_BiSC & From 1_1_Get_Pest_Data_BiS)

# Question: do we have the season?
# [INPUT]: .rdata object nwith Rows = participants and Columns = Variables. N = 1080 

# [OUTPUT]: 1.rds file, Rows = Participants, Cols = Vars - Related to Analysis N = 1080


#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

# Load data
load("../db/general_covars_data/1126_jllaurado_mQTLs.RData") # df2 1080 160
load("../db/general_covars_data/2023_10_05_Camille_BiSC_23_018.RData") # camille_bic_23_018 1080 174 
diet_df <- readRDS("../db/general_covars_data/diet_data_08012025.rds") # 809 25 

outfile <- "db/0_Get_Data/all_covars.rds"

#==============================================================================#----
# SELECT VARS OF INTEREST                                                      #----
#==============================================================================#----


# Get from each df covars of interest 
cov_df1 <- camille_bisc_23_018[ , c("id_mother", 
                                   "smoke_sust_m"),
                                ]

cov_df2 <- df2[ , c("id_bisc", 
                  "age_12w_m", 
                  "educ_level_m_3cat", 
                  "bmi_m_12w",
                  "parity_m_2cat", 
                  "ethnicity_c_4cat", 
                  "sex_0y_c",
                  "gestage_0y_c_days", 
                  "any_compl_m_0y_2c", 
                  "hospital_delivery")
                ]

cov_df3 <- diet_df[, c('id_mother', 
                     'fruit', 
                     'vegetables', 
                     'energy_intake')
                   ]


# Join Covariables in 1 df
# id_bisc in df2 == id_mother in camille_df & diet 

final_df <- merge(cov_df2, cov_df1,by.x = 'id_bisc', 'id_mother')
final_df <- merge(final_df, cov_df3, by.x = 'id_bisc', by.y = 'id_mother', all = TRUE) # 1080 14 
final_df$fruit_veg_intake <- 1000*((final_df$fruit + final_df$vegetables)/final_df$energy_intake)

rm(camille_bisc_23_018, df2, diet_df, cov_df1, cov_df2, cov_df3)

# Save df
saveRDS(final_df, outfile)


#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----