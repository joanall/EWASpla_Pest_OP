#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#

# [DESCRIPTION]: Descriptive statistics for BISC for main the 12w and 32 and for bt tp 
# Table 1: participant characteristics - all (N = ) vs main 12w / 32w / vs 12and32w

# [INPUT]:
# - db/final_df_incomplete_cases.rds
# - db/subsets/df_main_12w.rds, df_main_32w.rds 
# - db/subsets_2tp/df_main_12w.rds, df_main_32w.rds 
# [OUTPUT]: results/descriptives/table1.csv / .png
#           results/descriptives/table2_exposures.csv / .png

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
library(tableone)
source("script/00_template_plot.R")

outdir <- "results/descriptives/"

df_all     <- readRDS("db/final_df_incomplete_cases.rds")   # 484
df_main_12w <- readRDS("db/subsets/df_main_12w.rds")
df_main_32w <- readRDS("db/subsets/df_main_32w.rds")
df_main_12w_2tp <- readRDS("db/subsets_2tp/df_main_12w.rds")
df_main_32w_2tp <- readRDS("db/subsets_2tp/df_main_32w.rds")

nrow(df_all); nrow(df_main_12w); nrow(df_main_32w); nrow(df_main_12w_2tp); nrow(df_main_32w_2tp)
# 

#==============================================================================#----
# TABLE 1 - PARTICIPANT CHARACTERISTICS                                        #----
#==============================================================================#----

t1_vars   <- c("age_12w_m","educ_level_m_3cat", "smoke_sust_m", "bmi_m_12w", "parity_m_2cat",
               "fruit_veg_intake", "ethnicity_c_4cat", "sex_0y_c", "any_compl_m_0y_2c", "hospital_delivery", 
               "gestage_0y_c_days")  

t1_cat    <- c("educ_level_m_3cat","smoke_sust_m", "parity_m_2cat", "ethnicity_c_4cat", "sex_0y_c", "any_compl_m_0y_2c")
t1_nonnor <- c("fruit_veg_intake")

# Check: all variables exist
setdiff(t1_vars, colnames(df_all))

tab1_all <- CreateTableOne(vars = t1_vars, factorVars = t1_cat, data = df_all)
tab1_12w  <- CreateTableOne(vars = t1_vars, factorVars = t1_cat, data = df_main_12w)
tab1_32w  <- CreateTableOne(vars = t1_vars, factorVars = t1_cat, data = df_main_32w)
tab1_12w_2tp  <- CreateTableOne(vars = t1_vars, factorVars = t1_cat, data = df_main_12w_2tp)
tab1_32w_2tp  <- CreateTableOne(vars = t1_vars, factorVars = t1_cat, data = df_main_32w_2tp)

tab1_all_out <- print(tab1_all, nonnormal = t1_nonnor, showAllLevels = TRUE, printToggle = FALSE)
tab1_12w_out  <- print(tab1_12w,  nonnormal = t1_nonnor, showAllLevels = TRUE, printToggle = FALSE)
tab1_32w_out  <- print(tab1_32w,  nonnormal = t1_nonnor, showAllLevels = TRUE, printToggle = FALSE)
tab1_12w_2tp_out  <- print(tab1_12w_2tp,  nonnormal = t1_nonnor, showAllLevels = TRUE, printToggle = FALSE)
tab1_32w_2tp_out  <- print(tab1_32w_2tp,  nonnormal = t1_nonnor, showAllLevels = TRUE, printToggle = FALSE)
#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#

# [DESCRIPTION]: Descriptive statistics for BiSC
# Table 1: participant characteristics - all vs main 12w / main 32w /
# women with both time points (12w and 32w)
#
# [INPUT]:  db/final_df_incomplete_cases.rds
#           db/subsets/df_main_12w.rds, df_main_32w.rds
#           db/subsets_2tp/df_main_12w.rds, df_main_32w.rds
# [OUTPUT]: results/descriptives/table1.csv / .png

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
library(tableone)
source("script/00_template_plot.R")

outdir <- "results/descriptives/"

df_all          <- readRDS("db/0_Get_Data/final_df_incomplete_cases.rds")      # 1080
df_main_12w     <- readRDS("db/subsets/df_main_12w.rds")
df_main_32w     <- readRDS("db/subsets/df_main_32w.rds")
df_main_12w_2tp <- readRDS("db/subsets_2tp/df_main_12w.rds")
df_main_32w_2tp <- readRDS("db/subsets_2tp/df_main_32w.rds")

nrow(df_all); nrow(df_main_12w); nrow(df_main_32w); nrow(df_main_12w_2tp); nrow(df_main_32w_2tp)

# Check: 2tp subsets have the same women -> one column in Table 1
identical(sort(df_main_12w_2tp$ID_mother), sort(df_main_32w_2tp$ID_mother))

#==============================================================================#----
# TABLE 1 - PARTICIPANT CHARACTERISTICS                                        #----
#==============================================================================#----

t1_vars   <- c("age_12w_m", "educ_level_m_3cat", "smoke_sust_m", "bmi_m_12w", "parity_m_2cat",
               "fruit_veg_intake", "ethnicity_c_4cat", "sex_0y_c", "any_compl_m_0y_2c",
               "hospital_delivery", "gestage_0y_c_days")
t1_cat    <- c("educ_level_m_3cat", "smoke_sust_m", "parity_m_2cat", "ethnicity_c_4cat",
               "sex_0y_c", "any_compl_m_0y_2c", "hospital_delivery")
t1_nonnor <- c("fruit_veg_intake")

# Check: all variables exist
setdiff(t1_vars, colnames(df_all))

tab1_all     <- CreateTableOne(vars = t1_vars, factorVars = t1_cat, data = df_all)
tab1_12w     <- CreateTableOne(vars = t1_vars, factorVars = t1_cat, data = df_main_12w)
tab1_32w     <- CreateTableOne(vars = t1_vars, factorVars = t1_cat, data = df_main_32w)
tab1_2tp     <- CreateTableOne(vars = t1_vars, factorVars = t1_cat, data = df_main_12w_2tp)

tab1_all_out <- print(tab1_all, nonnormal = t1_nonnor, showAllLevels = TRUE, printToggle = FALSE)
tab1_12w_out <- print(tab1_12w, nonnormal = t1_nonnor, showAllLevels = TRUE, printToggle = FALSE)
tab1_32w_out <- print(tab1_32w, nonnormal = t1_nonnor, showAllLevels = TRUE, printToggle = FALSE)
tab1_2tp_out <- print(tab1_2tp, nonnormal = t1_nonnor, showAllLevels = TRUE, printToggle = FALSE)

# Check: same rows in the same order before combining
identical(rownames(tab1_all_out), rownames(tab1_12w_out))
identical(rownames(tab1_all_out), rownames(tab1_32w_out))
identical(rownames(tab1_all_out), rownames(tab1_2tp_out))

tab1 <- cbind(tab1_all_out, tab1_12w_out[, "Overall"], tab1_32w_out[, "Overall"],
              tab1_2tp_out[, "Overall"])
colnames(tab1) <- c("Level", "All", "Main 12w", "Main 32w", "Both time points")
tab1

# Table 1. Participant characteristics: rows = covariates (mean (SD),
# median [IQR] or n (%)), columns = all participants, main subsets and women
# with both time points.
save_table(tab1, paste0(outdir, "table1"))

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----
# Check: same rows in the same order before combining
identical(rownames(tab1_all_out), rownames(tab1_12w_out))
identical(rownames(tab1_all_out), rownames(tab1_32w_out))

tab1 <- cbind(tab1_all_out, tab1_12w_out[, "Overall"], tab1_32w_out[, "Overall"])
colnames(tab1) <- c("Level", "All", "Main 12w", "Main 32w", "Main 12w 2tp", "Main 32w 2tp")
tab1

# Table 1. Participant characteristics: rows = covariates (mean (SD),
# median [IQR] or n (%)), columns = all participants and main subsets.
save_table(tab1, paste0(outdir, "table1"))

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----