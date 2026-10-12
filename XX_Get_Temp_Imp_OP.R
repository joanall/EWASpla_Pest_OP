#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: April - October 2026                                                                              
#==============================================================================#

# [DESCRIPTION]: Perform temporary imputation on OP metabolites with <LOD 70 % 
# to be able to compute the molar sum. Temporary because Hana, who did the QC in BiSC
# will check the original QC code and update the db so they are available. 
# Participants in BiSC at 12w (N = 695) and at 32w (N= 755) 

# 5 OP metabolites 
# - DMP   --> fill-in, because detection frequency ≥ ~70% (i.e., <30% <LOD)
# - DMTP  --> Already imputed with fill-in and log transformed 
# - DEP   --> Already imputed with fill-in and log transformed 
# - DETP  --> LOD/√2, because detection frequency is low (~20%, i.e. 80% <LOD)
# Raw Values: nanograml / ml urine 
# Values for molar sum: µmoles/l 


# [INPUT]: .rdata file with chemcial ... XXX 

# [OUPUT]: Temporary files .rds with imputed OP metabolites. 2 files one for week 12w and
# one for week 32. 

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
# Clean enviornment
rm(list = ls())

# Load libraries
library(NADA)
library(msm)

# Define input path files 
in_file_op_12w <- "../db/chemicals/organophosphorous_urine/BISC_db_ops_m_12w_20240219.rdata"
in_file_op_32w <- "../db/chemicals/organophosphorous_urine/BISC_db_ops_m_32w_20240219.rdata"

# Define output path files 
out_file_op_12w <- "../db/chemicals/organophosphorous_urine/BISC_db_ops_m_12w_20240219_temp.rdata"
out_file_op_32w <- "../db/chemicals/organophosphorous_urine/BISC_db_ops_m_32w_20240219_temp.rdata"

# Set seed for  reproducibility
set.seed(666)  

# 3. Define LODs (from BiSC QC)
DMP_LOD_raw <- 1.20
DMP_LOQ_raw <- 4.00
DETP_LOD_raw <- 0.33
DETP_LOQ_raw <- 1.00

#==============================================================================#----
# OP metabolites 12W                                                           #----
#==============================================================================#----

# 1. Load OP Metabolites at 12 weeks (T2)
load(in_file_op_12w) # dataAll 695 x 37

# Impute DMP -------------------------------------------------------------------
# 2. Prepare data 
DMP_raw_num  <- as.numeric(as.character(dataAll$cc.DMP_m_12w_raw))
DMP_log  <- log(DMP_raw_num)
DMP_LOD <- log(DMP_LOD_raw)
DMP_LOQ <- log(DMP_LOQ_raw)

# 3. Flag data to be imputed
cens_DMP_log <- is.na(DMP_log) # NA = TRUE = censored (<LOD)
cens_DMP_loq_log <- !is.na(DMP_log) & DMP_log >= DMP_LOD & DMP_log <= DMP_LOQ  # between LOD and LOQ
obs_DMP_log <- ifelse(cens_DMP_log, DMP_LOD, DMP_log) # if TRUE (<LOD), put LOD-DEP value

# 4. Compute distribution parameters
cenros_DMP_log <- NADA::cenros(obs_DMP_log, cens_DMP_log, forwardT = NULL) #"forwardT = NULL": no log transform
mean_DMP_log <- NADA::mean(cenros_DMP_log)
sd_DMP_log   <- NADA::sd(cenros_DMP_log)

imputed_DMP_log <- DMP_log

# 5.1 Imputation if data is below LOD 
imputed_DMP_log[cens_DMP_log] <- msm::rtnorm(sum(cens_DMP_log), 
                                      mean = mean_DMP_log, 
                                      sd = sd_DMP_log,
                                      lower = -Inf, 
                                      upper = DMP_LOD)

# 5. 2  Imputation if is between LOD and LOQ
if (!is.null(DMP_LOQ)) {
  imputed_DMP_log[cens_DMP_loq_log] <- msm::rtnorm(
    n     = sum(cens_DMP_loq_log),
    mean  = mean_DMP_log,
    sd    = sd_DMP_log,
    lower = DMP_LOD,
    upper = DMP_LOQ
  )
}


# 6. Add as log2 to be consistent with previous imputations of other metabolites
dataAll$cc.DMP_m_12w_log2i <- imputed_DMP_log / log(2)

# Impute DETP ------------------------------------------------------------------
# 2. Prepare data 
DETP_raw_num  <- as.numeric(as.character(dataAll$cc.DETP_m_12w_raw))
DETP_raw_num[DETP_raw_num == 0] <- NA   # treat 0 as non-detect, not a real log-scale value - we have only 1 case 
DETP_log  <- log(DETP_raw_num)
DETP_LOD <- log(DETP_LOD_raw)
DETP_LOQ <- log(DETP_LOQ_raw)

# 3. Flag data to be imputed
cens_DETP_log <- is.na(DETP_log) # NA = TRUE = censored (<LOD)
cens_DETP_loq_log <- !is.na(DETP_log) & DETP_log >= DETP_LOD & DETP_log <= DETP_LOQ  # between LOD and LOQ
obs_DETP_log <- ifelse(cens_DETP_log, DETP_LOD, DETP_log) # if TRUE (<LOD), put LOD-DEP value

# 4. Compute distribution parameters
cenros_DETP_log <- NADA::cenros(obs_DETP_log, cens_DETP_log, forwardT = NULL) #"forwardT = NULL": no log transform
mean_DETP_log <- NADA::mean(cenros_DETP_log)
sd_DETP_log   <- NADA::sd(cenros_DETP_log)

imputed_DETP_log <- DETP_log

# 5.1 Imputation if data is below LOD 
imputed_DETP_log[cens_DETP_log] <- msm::rtnorm(sum(cens_DETP_log), 
                                             mean = mean_DETP_log, 
                                             sd = sd_DETP_log,
                                             lower = -Inf, 
                                             upper = DETP_LOD)

# 5. 2  Imputation if is between LOD and LOQ
if (!is.null(DETP_LOQ)) {
  imputed_DETP_log[cens_DETP_loq_log] <- msm::rtnorm(
    n     = sum(cens_DETP_loq_log),
    mean  = mean_DETP_log,
    sd    = sd_DETP_log,
    lower = DETP_LOD,
    upper = DETP_LOQ
  )
}


# 6. Add as log2 to be consistent with previous imputations of other metabolites
dataAll$cc.DETP_m_12w_log2i <- imputed_DETP_log / log(2)


# Save 12w DB ------------------------------------------------------------------
save(dataAll, file = out_file_op_12w)


#==============================================================================#----
# OP metabolites 32W                                                           #----
#==============================================================================#----

# 1. Load OP Metabolites at 32 weeks (T3)
load(in_file_op_32w) # dataAll 755 x 37

# Impute DMP -------------------------------------------------------------------
# 2. Prepare data 
DMP_raw_num  <- as.numeric(as.character(dataAll$cc.DMP_m_32w_raw))
DMP_log  <- log(DMP_raw_num)
DMP_LOD <- log(DMP_LOD_raw)
DMP_LOQ <- log(DMP_LOQ_raw)

# 3. Flag data to be imputed
cens_DMP_log <- is.na(DMP_log) # NA = TRUE = censored (<LOD)
cens_DMP_loq_log <- !is.na(DMP_log) & DMP_log >= DMP_LOD & DMP_log <= DMP_LOQ  # between LOD and LOQ
obs_DMP_log <- ifelse(cens_DMP_log, DMP_LOD, DMP_log) # if TRUE (<LOD), put LOD-DEP value

# 4. Compute distribution parameters
cenros_DMP_log <- NADA::cenros(obs_DMP_log, cens_DMP_log, forwardT = NULL) #"forwardT = NULL": no log transform
mean_DMP_log <- NADA::mean(cenros_DMP_log)
sd_DMP_log   <- NADA::sd(cenros_DMP_log)

imputed_DMP_log <- DMP_log

# 5.1 Imputation if data is below LOD 
imputed_DMP_log[cens_DMP_log] <- msm::rtnorm(sum(cens_DMP_log), 
                                             mean = mean_DMP_log, 
                                             sd = sd_DMP_log,
                                             lower = -Inf, 
                                             upper = DMP_LOD)

# 5. 2  Imputation if is between LOD and LOQ

if (!is.null(DMP_LOQ)) {
  imputed_DMP_log[cens_DMP_loq_log] <- msm::rtnorm(
    n     = sum(cens_DMP_loq_log),
    mean  = mean_DMP_log,
    sd    = sd_DMP_log,
    lower = DMP_LOD,
    upper = DMP_LOQ
  )
}


# 6. Add as log2 to be consistent with previous imputations of other metabolites
dataAll$cc.DMP_m_32w_log2i <- imputed_DMP_log / log(2)

# Impute DETP ------------------------------------------------------------------
# 2. Prepare data 
DETP_raw_num  <- as.numeric(as.character(dataAll$cc.DETP_m_32w_raw))
DETP_raw_num[DETP_raw_num == 0] <- NA   # treat 0 as non-detect, not a real log-scale value - we have only 1 case 
DETP_log  <- log(DETP_raw_num)
DETP_LOD <- log(DETP_LOD_raw)
DETP_LOQ <- log(DETP_LOQ_raw)

# 3. Flag data to be imputed
cens_DETP_log <- is.na(DETP_log) # NA = TRUE = censored (<LOD)
cens_DETP_loq_log <- !is.na(DETP_log) & DETP_log >= DETP_LOD & DETP_log <= DETP_LOQ  # between LOD and LOQ
obs_DETP_log <- ifelse(cens_DETP_log, DETP_LOD, DETP_log) # if TRUE (<LOD), put LOD-DEP value

# 4. Compute distribution parameters
cenros_DETP_log <- NADA::cenros(obs_DETP_log, cens_DETP_log, forwardT = NULL) #"forwardT = NULL": no log transform
mean_DETP_log <- NADA::mean(cenros_DETP_log)
sd_DETP_log   <- NADA::sd(cenros_DETP_log)

imputed_DETP_log <- DETP_log

# 5.1 Imputation if data is below LOD 
imputed_DETP_log[cens_DETP_log] <- msm::rtnorm(sum(cens_DETP_log), 
                                               mean = mean_DETP_log, 
                                               sd = sd_DETP_log,
                                               lower = -Inf, 
                                               upper = DETP_LOD)

# 5. 2  Imputation if is between LOD and LOQ
if (!is.null(DETP_LOQ)) {
  imputed_DETP_log[cens_DETP_loq_log] <- msm::rtnorm(
    n     = sum(cens_DETP_loq_log),
    mean  = mean_DETP_log,
    sd    = sd_DETP_log,
    lower = DETP_LOD,
    upper = DETP_LOQ
  )
}


# 6. Add as log2 to be consistent with previous imputations of other metabolites
dataAll$cc.DETP_m_32w_log2i <- imputed_DETP_log / log(2)


# Save 32w DB ------------------------------------------------------------------
save(dataAll, file = out_file_op_32w)

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----


