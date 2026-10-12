#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: April - October 2026                                                                              
#==============================================================================#

# [DESCRIPTION]: OP metabolites can have multiple parent pesticides and it's dificutl
# to know from whihc ompunds they come from. Therefore, it make sense to sutdy them in grups
# fro that we'll compute the molar sum of OP metabolites in both time points 

# [INPUT]: 1 rds file with  Rows = Participants IDs, Cols = Exposures_week
# Both files N = 920 Participant, some with missing. 

# [OUTPUT]: 1.rds file, All pesticides  in single data.frame. N = XXXX
# with averages between time points. 

# CHEMICAL COMPUNDS DATA: 8 compounds
# 2 pesticide biomarkers:
# - 2,4-Dichlorophenoxyacetic acid (2,4-D) — herbicide  
# - Pyrethorids - specific metabolite 3-Phenoxybenzoic acid (3-PBA) 
# 6 organophosphate pesticide (OP) metabolites
# - 1 specific biomarker: TCPy (metabolite of chlorpyrifos)
# - 5 non-specific OP metabolites: DMP, DMTP, DMDTP, DEP, DETP

# [CONSIDERATIONS]: 
# Pregnancy-average exposure: 
# - Continuous variables:
#     - 2 time points: geometric mean of both timepoints (log2 scale)
#     - 1 time point: single value available 
# - Dichotomous variables: 3-level categorical
#  - 0 not detected at either
#  - 1 detected at one time point 
#  - 2 detected at both 

# Molar Sum of OP metabolites: 
# ΣDEP (diethylphosphate metabolites):  DEP/154.1 + DETP/170.17
# ΣDMP (dimethylphosphate metabolites): DMP/126.05 + DMTP/142.12
# ΣDAP (total dialkylphosphate metabolites): ΣDEP + ΣDMP

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

# Load libraries

# Define input path files 
in_file_op_12w <- "../db/chemicals/organophosphorous_urine/BISC_db_ops_m_12w_20240219_temp.rdata"
in_file_op_32w <- "../db/chemicals/organophosphorous_urine/BISC_db_ops_m_32w_20240219_temp.rdata"

# Define output path files 
out_file_op_12w <- "../db/chemicals/organophosphorous_urine/BISC_db_ops_m_12w_20240219_temp_molarsum.rdata"
out_file_op_32w <- "../db/chemicals/organophosphorous_urine/BISC_db_ops_m_32w_20240219_temp_molarsum.rdata"

#==============================================================================#----
# OP metabolites 12w (T2)                                                      #----
#==============================================================================#----

# 1. Load op metabolites at 12 weeks (T3)
load(in_file_op_12w) # dataAll 695 x 39

# 2. Back-transform imputed log2 values to raw scale (ng/mL)
DEP_12w <- 2^as.numeric(dataAll$cc.DEP_m_12w_log2i)
DETP_12w <- 2^as.numeric(dataAll$cc.DETP_m_12w_log2i)
DMP_12w <- 2^as.numeric(dataAll$cc.DMP_m_12w_log2i)
DMTP_12w <- 2^as.numeric(dataAll$cc.DMTP_m_12w_log2i)

# 3. Molar sums (umol/L) - ng/mL divided by molecular weight (g/mol) gives umol/L
sum_DEP_12w <- DEP_12w/154.1  + DETP_12w/170.17
sum_DMP_12w <- DMP_12w/126.05 + DMTP_12w/142.12
sum_DAP_12w <- sum_DEP_12w + sum_DMP_12w

# 4. Add molar sum to db and save
dataAll$sum_DEP_12w <- sum_DEP_12w
dataAll$sum_DMP_12w <- sum_DMP_12w
dataAll$sum_DAP_12w <- sum_DAP_12w

save(dataAll, file = out_file_op_12w)

#==============================================================================#----
# OP metabolites 32w (T3)                                                      #----
#==============================================================================#----
# 1. Load op metabolites at 32 weeks (T3)
load(in_file_op_32w) # dataAll 755 x 39

# 2. Back-transform imputed log2 values to raw scale (ng/mL)
DEP_32w <- 2^as.numeric(dataAll$cc.DEP_m_32w_log2i)
DETP_32w <- 2^as.numeric(dataAll$cc.DETP_m_32w_log2i)
DMP_32w <- 2^as.numeric(dataAll$cc.DMP_m_32w_log2i)
DMTP_32w <- 2^as.numeric(dataAll$cc.DMTP_m_32w_log2i)

# 3. Molar sums (umol/L) - ng/mL divided by molecular weight (g/mol) gives umol/L
sum_DEP_32w <- DEP_32w/154.1  + DETP_32w/170.17
sum_DMP_32w <- DMP_32w/126.05 + DMTP_32w/142.12
sum_DAP_32w <- sum_DEP_32w + sum_DMP_32w

# 4. Add molar sum to db and save
dataAll$sum_DEP_32w <- sum_DEP_32w
dataAll$sum_DMP_32w <- sum_DMP_32w
dataAll$sum_DAP_32w <- sum_DAP_32w

save(dataAll, file = out_file_op_32w)

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----
