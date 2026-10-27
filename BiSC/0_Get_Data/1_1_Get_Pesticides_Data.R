#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: April - October 2026                                                                              
#==============================================================================#

# [DESCRIPTION]: Get pesticide exposure data from BiSC at 12w (N = 695) and
# 32w (N = 755) and build a single data frame.

# [INPUT]:  BISC_db_pest_m_12w_20240729.rdata, BISC_db_pest_m_32w_20240729.rdata
#           (object dataAll, rows = mothers)
# [OUTPUT]: db/all_pest_12w_32w.rds - one row per mother, 12w and 32w variables

# CHEMICAL COMPOUNDS: 2,4-D, TCPy, 3-PBA
# Variables kept per compound and time point:
# - _raw   : raw concentration (ug/L), for descriptives
# - _log2i : log2 imputed, for the models (TCPy, 3-PBA)
# - _dic   : detected / undetected, for the models (2,4-D)
# + urine specific gravity (SG), for the sensitivity analysis

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

file_pest_12w <- "../db/chemicals/pesticides_pah_urine/BISC_db_pest_m_12w_20240729.rdata"
file_pest_32w <- "../db/chemicals/pesticides_pah_urine/BISC_db_pest_m_32w_20240729.rdata"
output_file   <- "db/0_Get_Data/all_pest_12w_32w.rds"

#==============================================================================#----
# PESTICIDES 12w (T2)                                                          #----
#==============================================================================#----

load(file_pest_12w) # dataAll 695 x 56

# Check: all pesticide variables available (raw, log2i, dic, OP metabolites?)
grep("^cc\\.", colnames(dataAll), value = TRUE)

# Old names and new names in the same order
cov_12w <- c("ID_mother",
             "cc.2,4-D_m_12w_raw", "cc.TCPy_m_12w_raw",   "cc.3-PBA_m_12w_raw",
             "cc.2,4-D_m_12w_dic", "cc.TCPy_m_12w_log2i", "cc.3-PBA_m_12w_log2i",
             "ur_pool_grav_m_12w")

pest_12w <- dataAll[, cov_12w]

rm(dataAll)

#==============================================================================#----
# PESTICIDES 32w (T3)                                                          #----
#==============================================================================#----

load(file_pest_32w) # dataAll 755 x 37

cov_32w <- c("ID_mother",
             "cc.2,4-D_m_32w_raw", "cc.TCPy_m_32w_raw",   "cc.3-PBA_m_32w_raw",
             "cc.2,4-D_m_32w_dic", "cc.TCPy_m_32w_log2i", "cc.3-PBA_m_32w_log2i",
             "ur_pool_grav_m_32w")

pest_32w <- dataAll[, cov_32w]

rm(dataAll)

#==============================================================================#----
# JOIN PESTICIDES 12w + 32w                                                    #----
#==============================================================================#----

pest_all <- merge(pest_12w, pest_32w, by = "ID_mother", all = TRUE)
nrow(pest_all) # 920

# Mothers with both time points
length(intersect(pest_12w$ID_mother, pest_32w$ID_mother)) # 530
# - 530 with both, 165 only 12w, 225 only 32w

saveRDS(pest_all, output_file)

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----