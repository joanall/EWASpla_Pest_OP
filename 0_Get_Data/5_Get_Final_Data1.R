#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: April - October 2026                                                                              
#==============================================================================#

# [DESCRIPTION]: Join all DBs needed to perfom EWAS analsysis on Pesticides - DNAm.
# Join exposure data (pest + op), covariates (general) + methylation covars
# [INPUT]: 

# [OUTPUT]: do we want cc or not? do we want different sample sizes in each? 



#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

infile <- "db/final_df_incomplete_cases.rds"
df <- readRDS(infile) # 1080 57 


#==============================================================================#----
# APPROACH 1 - RLM                                                             #----
#==============================================================================#----

# 1) We need all participants to have outcome and exposure 
# 2) All participants DNAm 
# 3) All participants Chem  
# 4) All participants main covars - mat age, ed, smoke, ga, chil eth, child sex, hospital (+ ct, dna cont)

# Remove all participants without placental DNAm data 
with_dnam <- df[!is.na(df$dna_contamination),] # 565 

# Remove all participants without chem data at each time point 
with_dnam_chem_12w <- with_dnam[!is.na(with_dnam$Diclor_12w),] # 384 
with_dnam_chem_32w <- with_dnam[!is.na(with_dnam$Diclor_32w),] # 442
with_dnam_chem_avg <- with_dnam[!is.na(with_dnam$Diclor_avg),] # 510
with_dnam_chem_12w_32w <- with_dnam[!is.na(with_dnam$Diclor_12w) & is.na(with_dnam$Diclor_32w),] # 68 ? 

# Remove all participants without main covariables
main_cov <- c('age_12w_m', 'educ_level_m_3cat', 'smoke_sust_m', 'gestage_0y_c_days', 
              'ethnicity_c_4cat', 'sex_0y_c', 'hospital_delivery')

with_dnam_chem_12w_cov <- with_dnam_chem_12w[complete.cases(with_dnam_chem_12w[, main_cov]), ] # 378 
with_dnam_chem_32w_cov <- with_dnam_chem_32w[complete.cases(with_dnam_chem_32w[, main_cov]), ] # 437 
with_dnam_chem_avg_cov <- with_dnam_chem_avg[complete.cases(with_dnam_chem_avg[, main_cov]), ] # 503


# Sensitivity models - Only in AVG? 
# S1 - Excl. Pregnancy Complications  N  = XXX
# S2 - Excl. Non Europeans N == XXX
# S3 - Adjusting SG
# S4 - Adjusting FruitVeg
# S5 - Wo GA

s1 <- with_dnam_chem_avg_cov[with_dnam_chem_avg_cov$any_compl_m_0y_2c == "no",] # 395 
s2 <- with_dnam_chem_avg_cov[with_dnam_chem_avg_cov$ethnicity_c_4cat == "eur",] # 351 
s3 <- with_dnam_chem_avg_cov[complete.cases(with_dnam_chem_avg_cov$SG_avg),] # 503 
s4 <- with_dnam_chem_avg_cov[complete.cases(with_dnam_chem_avg_cov$fruit_veg_intake),] # 422 
s5_cov <- c('age_12w_m', 'educ_level_m_3cat', 'smoke_sust_m','ethnicity_c_4cat', 'sex_0y_c', 'hospital_delivery')
s5 <- with_dnam_chem_avg_cov[complete.cases(s5_cov),] # 503 

# Save -- decide if ok diff sample sizes



#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----
