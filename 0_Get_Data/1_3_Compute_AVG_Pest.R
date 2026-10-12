#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: April - October 2026                                                                              
#==============================================================================#

# [DESCRIPTION]: Given data with chemicals in 2 time points for pesticides and op 
# metabolites, compute the averages between both time points. Also, compute
# the molar sum for the op metabolites. 

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

# Note: Geometric mean calculations require all values to be greater than 0. 
# If your time points contain zeros or negative numbers, the result will return NaN or 0. 

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

# Load libraries

# Define input & output path files 
input_file <- "db/all_pest_op.rds"
output_file <- "db/all_pest_op_avg.rds"

# Load data 
df <- readRDS(input_file) # 920 x 21 

#==============================================================================#----
# AVERAGES 12w - 32w                                                           #----
#==============================================================================#----

# 1. Compute Average Geometric Mean - TCPy, PBA - ln2 
df$TCPy_avg <- rowMeans(df[, c("TCPy_12w", "TCPy_32w")], na.rm = TRUE)
df$PBA_avg <- rowMeans(df[, c("PBA_12w", "PBA_32w")], na.rm = TRUE)

# 2. Compute Average Categorical Vars - 2,4-D
df$Diclor_avg[df$Diclor_12w == "detected" & df$Diclor_32w == "detected"] <- 2
df$Diclor_avg[(df$Diclor_12w == "detected" | df$Diclor_32w == "detected")] <- 1
df$Diclor_avg[df$Diclor_12w == "undetected" | df$Diclor_32w == "undetected"] <- 0

# 3. Compute average Specific Gravity 
df$SG_avg <- rowMeans(df[, c("SG_12w", "SG_32w")], na.rm = TRUE)


# 3. Compute Average Molar Sum OP metabolites
df$sum_DEP_avg <- rowMeans(df[, c("sum_DEP_12w", "sum_DEP_32w")], na.rm = TRUE)
df$sum_DMP_avg <- rowMeans(df[, c("sum_DMP_12w", "sum_DMP_32w")], na.rm = TRUE)
df$sum_DAP_avg <- rowMeans(df[, c("sum_DAP_12w", "sum_DAP_32w")], na.rm = TRUE)

# Save 
saveRDS(df, output_file)

#==============================================================================#----
# END OF CODE                                                                  #----
#==============================================================================#----



# get_avg_2tp <- function(df, tp1, tp2) { ----
#   
#   temp <- rep(NA_real_, nrow(df))                     
#   
#   x1 <- df[[tp1]]
#   x2 <- df[[tp2]]
#   
#   both  <- !is.na(x1) & !is.na(x2)
#   only1 <- !is.na(x1) &  is.na(x2)
#   only2 <-  is.na(x1) & !is.na(x2)
#   
#   temp[both] <- (x1[both] + x2[both]) / 2
#   temp[only1] <- x1[only1]
#   temp[only2] <- x2[only2]
#   
#   return(temp)
# }

