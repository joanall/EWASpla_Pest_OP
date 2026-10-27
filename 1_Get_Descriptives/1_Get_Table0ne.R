#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: Table 1 - participant characteristics, per cohort.
# - SEPAGES: all participants vs main subset (women with both time points)
# - BiSC:    all participants vs main t2 / main t3 / both time points /
#            average (>= 1 time point)
#
# [INPUT SEPAGES]: db/subsets/df_all_recoded.rds, db/subsets/df_main.rds
# [INPUT BiSC]:    db/subsets/df_all_recoded.rds, df_main_t2.rds, df_main_t3.rds,
#                  df_main_avg.rds; db/subsets_2tp/df_main.rds
# [OUTPUT]: results/1_Descriptives/table1_sepages.csv / .png
#           results/1_Descriptives/table1_bisc.csv / .png
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
library(tableone)
source("script/00_template_plot.R")

outdir <- "results/1_Descriptives/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Variables (same names in both cohorts; hospital only in BiSC)
t1_vars   <- c("mo_age", "mo_edu", "mo_smoke", "mo_bmi", "mo_parity",
               "mo_fv_g1000kcal", "ch_eth", "ch_sex", "preg_compl", "gest_age")
t1_cat    <- c("mo_edu", "mo_smoke", "mo_parity", "ch_eth", "ch_sex", "preg_compl")
t1_nonnor <- c("mo_fv_g1000kcal")

# One tableone column (as character matrix) for a data frame
table1_column <- function(df, vars, cat_vars, nonnormal) {
  tab <- CreateTableOne(vars = vars, factorVars = cat_vars, data = df)
  print(tab, nonnormal = nonnormal, showAllLevels = TRUE, printToggle = FALSE)
}

#==============================================================================#----
# SEPAGES                                                                      #----
#==============================================================================#----
# Run this section on the SEPAGES project (paths relative to its root)

sep_all  <- readRDS("db/subsets/df_main_recoded.rds")
sep_main <- readRDS("db/subsets/df_main.rds")
nrow(sep_all); nrow(sep_main)   # 484, 312

# Check: all variables exist
setdiff(t1_vars, colnames(sep_all))
setdiff(t1_vars, colnames(sep_main))

sep_all_out  <- table1_column(sep_all,  t1_vars, t1_cat, t1_nonnor)
sep_main_out <- table1_column(sep_main, t1_vars, t1_cat, t1_nonnor)

# Check: same rows in the same order before combining
identical(rownames(sep_all_out), rownames(sep_main_out))

tab1_sepages <- cbind(sep_all_out, sep_main_out[, "Overall"])
colnames(tab1_sepages) <- c("Level", "All", "Main (both time points)")
tab1_sepages

# Table 1 SEPAGES. Participant characteristics: rows = covariates (mean (SD),
# median [IQR] or n (%)), columns = all participants and main subset.
save_table(tab1_sepages, paste0(outdir, "table1_sepages"))

#==============================================================================#----
# BiSC                                                                         #----
#==============================================================================#----
# Run this section on the BiSC project (paths relative to its root)

t1_vars_bisc <- c(t1_vars, "hospital")
t1_cat_bisc  <- c(t1_cat,  "hospital")

bisc_all  <- readRDS("db/subsets/df_all_recoded.rds")
bisc_t2   <- readRDS("db/subsets/df_main_t2.rds")
bisc_t3   <- readRDS("db/subsets/df_main_t3.rds")
bisc_2tp  <- readRDS("db/subsets_2tp/df_main.rds")
bisc_avg  <- readRDS("db/subsets/df_main_avg.rds")
nrow(bisc_all); nrow(bisc_t2); nrow(bisc_t3); nrow(bisc_2tp); nrow(bisc_avg)

# Check: all variables exist
setdiff(t1_vars_bisc, colnames(bisc_all))

bisc_all_out <- table1_column(bisc_all, t1_vars_bisc, t1_cat_bisc, t1_nonnor)
bisc_t2_out  <- table1_column(bisc_t2,  t1_vars_bisc, t1_cat_bisc, t1_nonnor)
bisc_t3_out  <- table1_column(bisc_t3,  t1_vars_bisc, t1_cat_bisc, t1_nonnor)
bisc_2tp_out <- table1_column(bisc_2tp, t1_vars_bisc, t1_cat_bisc, t1_nonnor)
bisc_avg_out <- table1_column(bisc_avg, t1_vars_bisc, t1_cat_bisc, t1_nonnor)

# Check: same rows in the same order before combining
identical(rownames(bisc_all_out), rownames(bisc_t2_out))
identical(rownames(bisc_all_out), rownames(bisc_t3_out))
identical(rownames(bisc_all_out), rownames(bisc_2tp_out))
identical(rownames(bisc_all_out), rownames(bisc_avg_out))

tab1_bisc <- cbind(bisc_all_out, bisc_t2_out[, "Overall"], bisc_t3_out[, "Overall"],
                   bisc_2tp_out[, "Overall"], bisc_avg_out[, "Overall"])
colnames(tab1_bisc) <- c("Level", "All", "Main t2", "Main t3",
                         "Both time points", "Average (>= 1 time point)")
tab1_bisc

# Table 1 BiSC. Participant characteristics: rows = covariates (mean (SD),
# median [IQR] or n (%)), columns = all participants, main subsets per time
# point, women with both time points, and the average-exposure subset.
save_table(tab1_bisc, paste0(outdir, "table1_bisc"), base_size = 9)

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----