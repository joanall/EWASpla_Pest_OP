#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: April - October 2026                                                                              
#==============================================================================#

# [DESCRIPTION]: Get OP metabolite exposure data from BiSC at 12w (N = 695)
# and 32w (N = 755) and build a single data frame.

# [INPUT]:  BISC_db_ops_m_12w_20240219_temp_molarsum.rdata,
#           BISC_db_ops_m_32w_20240219_temp_molarsum.rdata (object dataAll)
# [OUTPUT]: db/all_op_12w_32w.rds - one row per mother, 12w and 32w variables

# OP METABOLITES: DMP, DMTP, DMDTP, DEP, DETP
# Variables kept per metabolite and time point:
# - _raw   : raw concentration (ug/L), for descriptives
# - _log2i : log2 imputed (renamed without suffix)
# Molar sums (umol/L, not log-transformed), for the models:
# - sum_DMP (DMP + DMTP + DMDTP), sum_DEP (DEP + DETP), sum_DAP (all)

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

file_op_12w <- "../db/chemicals/organophosphorous_urine/BISC_db_ops_m_12w_20240219_temp_molarsum.rdata"
file_op_32w <- "../db/chemicals/organophosphorous_urine/BISC_db_ops_m_32w_20240219_temp_molarsum.rdata"
output_file <- "db/0_Get_Data/all_op_12w_32w.rds"

#==============================================================================#----
# OP METABOLITES 12w (T2)                                                      #----
#==============================================================================#----

load(file_op_12w) # dataAll 695 x 39


# Old names and new names in the same order
cov_12w <- c("ID_mother",
             "cc.DMP_m_12w_raw",   "cc.DMTP_m_12w_raw",   "cc.DMDTP_m_12w_raw",
             "cc.DEP_m_12w_raw",   "cc.DETP_m_12w_raw",
             "cc.DMP_m_12w_log2i", "cc.DMTP_m_12w_log2i", "cc.DMDTP_m_12w_dic",
             "cc.DEP_m_12w_log2i", "cc.DETP_m_12w_log2i",
             "sum_DMP_12w",        "sum_DEP_12w",         "sum_DAP_12w")


# Check: all old names exist
op_12w <- dataAll[, cov_12w]

rm(dataAll)

#==============================================================================#----
# OP METABOLITES 32w (T3)                                                      #----
#==============================================================================#----

load(file_op_32w) # dataAll 755 x 39
cov_32w <- c("ID_mother",
             "cc.DMP_m_32w_raw",   "cc.DMTP_m_32w_raw",   "cc.DMDTP_m_32w_raw",
             "cc.DEP_m_32w_raw",   "cc.DETP_m_32w_raw",
             "cc.DMP_m_32w_log2i", "cc.DMTP_m_32w_log2i", "cc.DMDTP_m_32w_dic",
             "cc.DEP_m_32w_log2i", "cc.DETP_m_32w_log2i",
             "sum_DMP_32w",        "sum_DEP_32w",         "sum_DAP_32w")

op_32w <- dataAll[, cov_32w]

rm(dataAll)

#==============================================================================#----
# JOIN OP METABOLITES 12w + 32w                                                #----
#==============================================================================#----

op_all <- merge(op_12w, op_32w, by = "ID_mother", all = TRUE)
nrow(op_all) # 920

# Mothers with both time points
length(intersect(op_12w$ID_mother, op_32w$ID_mother)) # 530
# - 530 with both, 165 only 12w, 225 only 32w

saveRDS(op_all, output_file)

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----