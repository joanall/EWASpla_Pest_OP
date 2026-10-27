#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: April - October 2026                                                                              
#==============================================================================#

# [DESCRIPTION]: Get Pesticides Exposure data from BiSC cohort from 12w (N = 695) 
# and 32w (N = 755) and built a unique data frame. 

# [INPUT]: 2 .rds files, data.frame with Rows = Participants IDs, Cols = Exposures_week
# Both files N = 920 Participant, some with missing. 

# [OUTPUT]: 1.rds file, All pesticides and op metabolites data in single data.frame. N = XXXX


#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

# Load libraries

# Define input path files 
file_pest_12w_32w <- "db/0_Get_Data/all_pest_12w_32w.rds"
file_op_12w_32w <- "db/0_Get_Data/all_op_12w_32w.rds"

# Define output path files 
output_file <- "db/0_Get_Data/all_pest_op.rds"


#==============================================================================#----
# JOIN PEST + OP                                                               #----
#==============================================================================#----

# 1. Load both dfs
pest_all <- readRDS(file_pest_12w_32w) # 920 x 15 vars (3 pest x 2 tp + 1 ID)
op_all <- readRDS(file_op_12w_32w) # 920 x 15 vars (5 op x 2 tp + 3 molar sum + 1 ID)

# 2. Join them by Participant (mother) ID
pest_op <- merge(pest_all, 
                 op_all, 
                 by = 'ID_mother', 
                 all = TRUE) # 920 x 21

# 3. Save Data 
saveRDS(pest_op, output_file)

#==============================================================================#----
# CHECK MISSIGNESS                                                             #----
#==============================================================================#----

# Overall Missing Data ---------------------------------------------------------
summary(is.na(pest_op))

# N with Chemical Data 
# Note: Data Chemicals 12w = 695 yes // 225 no data 
# Note: Data Chemicals 32w = 755 yes // 165 no 

# Check N with missing in each imepoint ----------------------------------------

# 1. Define cols chemicals each timepoints 
cols_12w <- c("cc.DMP_m_12w_log2i", "cc.DMTP_m_12w_log2i", "cc.DMDTP_m_12w_dic",
              "cc.DEP_m_12w_log2i", "cc.DETP_m_12w_log2i",
              "sum_DMP_12w","sum_DEP_12w","sum_DAP_12w")
cols_32w <- c("cc.DMP_m_32w_log2i", "cc.DMTP_m_32w_log2i", "cc.DMDTP_m_32w_dic",
              "cc.DEP_m_32w_log2i", "cc.DETP_m_32w_log2i",
              "sum_DMP_32w","sum_DEP_32w","sum_DAP_32w")

# 2. Flag each participant as missing at a timepoint if ANY chemicals at that timepoint are NA
missing_12w <- apply(pest_op[, cols_12w], 1, function(row) any(is.na(row))) # factor true/false
missing_32w <- apply(pest_op[, cols_32w], 1, function(row) any(is.na(row)))

# 3. Count missings 
n_12w <- sum(!missing_12w) # 695 — have 12w data
n_12w_only <- sum(!missing_12w & missing_32w) # 165 — 12w only, no 32w

n_32w <- sum(!missing_32w) # 755 — have 32w data
n_32w_only <- sum(missing_12w & !missing_32w) # 225 — 32w only, no 12w

n_both <- sum(!missing_12w & !missing_32w) # 530 — have both
n_total <- nrow(pest_op) # 920

# Sanity check (530 + 165 + 225 = 920)
n_both + n_12w_only + n_32w_only == n_total # TRUE

#==============================================================================#----
# END OF CODE                                                                  #----
#==============================================================================#----