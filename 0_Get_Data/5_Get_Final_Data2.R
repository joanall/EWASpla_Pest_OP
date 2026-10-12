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
# APPROACH 2 - RLMER                                                           #----
#==============================================================================#----

# Wide Format - one column indicate time and another one indicated the what 
# Example: 
# Idmother     value pest 1 time 


# Get Chem x each week 
df_12w <- data.frame(df[,c('ID_mother', 'Diclor_12w', 'TCPy_12w', 'PBA_12w', 
                           'sum_DEP_12w', 'sum_DMP_12w', 'sum_DAP_12w')])

df_32w <- data.frame(df[,c('ID_mother', 'Diclor_32w', 'TCPy_32w', 'PBA_32w', 
                           'sum_DEP_32w', 'sum_DMP_32w', 'sum_DAP_32w')])
# Add time variable
df_12w$week <- 12
df_32w$week <- 32

# Standardize names: remove _12w / _32w so columns match
names(df_12w) <- gsub("_12w", "", names(df_12w), fixed = TRUE)
names(df_32w) <- gsub("_32w", "", names(df_32w), fixed = TRUE)

# Check order IDs and stack Join dfs
df_long <- rbind(df_12w, df_32w)

# Get covars and attach 
main_cov <- c('ID_mother',
              'age_12w_m', 'educ_level_m_3cat', 'smoke_sust_m', 'gestage_0y_c_days', 
              'ethnicity_c_4cat', 'sex_0y_c', 'hospital_delivery', 
              "Trophoblasts", "Stromal", "Hofbauer", "Endothelial",   
              "nRBC", "Syncytiotrophoblast", "dna_contamination")

df_cov <- df[, main_cov]
df_long_cov <- merge(df_long, df_cov, by = 'ID_mother') # 2160 

df_long_cov_cc <- df_long_cov[complete.cases(df_long_cov[,main_cov]),] # 1106 (553) 

# Save ? 



#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----
