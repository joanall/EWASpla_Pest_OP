#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#

# [DESCRIPTION]: Descriptive statistics for SEPAGES
# Table 1: participant characteristics - all (N = 484) vs main t2 / main t3

# [INPUT]:  db/db_pest_pladnam_incomplete_cases.rds
#           db/subsets/df_main_t2.rds, df_main_t3.rds
# [OUTPUT]: results/1_Descriptives/table1.csv 

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
library(haven)
library(tableone)
source("script/00_template_plot.R")

outdir <- "results/descriptives/"

df_all     <- readRDS("db/db_pest_pladnam_incomplete_cases.rds")   # 484
df_main_t2 <- readRDS("db/subsets/df_main_t2.rds")
df_main_t3 <- readRDS("db/subsets/df_main_t3.rds")
nrow(df_all); nrow(df_main_t2); nrow(df_main_t3)
# 484 314 313 

#==============================================================================#----
# TABLE 1 - PARTICIPANT CHARACTERISTICS                                        #----
#==============================================================================#----

t1_vars   <- c("mo_age", "mo_edu", "mo_smoke", "mo_bmi_bepr", "mo_parity",
               "fv_g1000kcal", "ch_ethnicity", "ch_sex", "preg_compl")  # + "ga_days"
t1_cat    <- c("mo_edu", "mo_smoke", "mo_parity", "ch_ethnicity", "ch_sex", "preg_compl")
t1_nonnor <- c("fv_g1000kcal")

# Check: all variables exist
setdiff(t1_vars, colnames(df_all))

tab1_all <- CreateTableOne(vars = t1_vars, factorVars = t1_cat, data = df_all)
tab1_t2  <- CreateTableOne(vars = t1_vars, factorVars = t1_cat, data = df_main_t2)
tab1_t3  <- CreateTableOne(vars = t1_vars, factorVars = t1_cat, data = df_main_t3)

tab1_all_out <- print(tab1_all, nonnormal = t1_nonnor, showAllLevels = TRUE, printToggle = FALSE)
tab1_t2_out  <- print(tab1_t2,  nonnormal = t1_nonnor, showAllLevels = TRUE, printToggle = FALSE)
tab1_t3_out  <- print(tab1_t3,  nonnormal = t1_nonnor, showAllLevels = TRUE, printToggle = FALSE)

# Check: same rows in the same order before combining
identical(rownames(tab1_all_out), rownames(tab1_t2_out))
identical(rownames(tab1_all_out), rownames(tab1_t3_out))

tab1 <- cbind(tab1_all_out, tab1_t2_out[, "Overall"], tab1_t3_out[, "Overall"])
colnames(tab1) <- c("Level", "All", "Main t2", "Main t3")
tab1

# Table 1. Participant characteristics: rows = covariates (mean (SD),
# median [IQR] or n (%)), columns = all participants and main subsets.
save_table(tab1, paste0(outdir, "table1"))

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----