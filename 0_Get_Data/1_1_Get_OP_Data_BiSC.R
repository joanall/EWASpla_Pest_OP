#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: April - October 2026                                                                              
#==============================================================================#

# [DESCRIPTION]: Get OP Metabolites Exposure data from BiSC cohort from 12w (N = 695) 
# and 32w (N = 755) and built a unique dataframe. 

# [INPUT]: file.rdata with xxx xxxx 

# [OUTPUT]: Write output/s obtained after running code. 

# CHEMICAL COMPUNDS DATA: 8 compounds
# 2 pesticide biomarkers:
# - 2,4-Dichlorophenoxyacetic acid (2,4-D) — herbicide  
# - Pyrethorids - specific metabolite 3-Phenoxybenzoic acid (3-PBA) 
# 6 organophosphate pesticide (OP) metabolites
# - 1 specific biomarker: TCPy (metabolite of chlorpyrifos)
# - 5 non-specific OP metabolites: DMP, DMTP, DMDTP, DEP, DETP

# Values: Raw micromililiter/mol -> Because we will compute the molar sum later on.  
# In preliminary analyses we used: 
# - DMTP — continuous
# - DMDTP — dichotomized
# - DEP — continuous (and dichotomized)
# - DETP — dichotomized


#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

# Load libraries

# Define input path files 
file_op_12w <- "../db/chemicals/organophosphorous_urine/BISC_db_ops_m_12w_20240219_temp_molarsum.rdata"
file_op_32w <- "../db/chemicals/organophosphorous_urine/BISC_db_ops_m_32w_20240219_temp_molarsum.rdata"

# Define output path files 
output_file <- "db/all_op_12w_32w.rds"

#==============================================================================#----
# OP metabolites 12w (T2)                                                      #----
#==============================================================================#----

# 1. Load op metabolites at 12 weeks (T3)
load(file_op_12w) # dataAll 695 x 39

# 2. Subset compounds of interest 
op_12w <- dataAll[,c('ID_mother', 
                     "cc.DMP_m_12w_log2i" ,
                     "cc.DMTP_m_12w_log2i",
                     "cc.DEP_m_12w_log2i",
                     "cc.DETP_m_12w_log2i", 
                     'sum_DEP_12w', 
                     'sum_DMP_12w', 
                     'sum_DAP_12w')]

# 3. Rename for readability (new_name = old_name)
op_12w <- dplyr::rename(op_12w, 
                        DMP_12w = `cc.DMP_m_12w_log2i`, 
                        DMTP_12w = `cc.DMTP_m_12w_log2i`, 
                        DEP_12w = `cc.DEP_m_12w_log2i`,
                        DETP_12w = `cc.DETP_m_12w_log2i`)

head(op_12w)
# ID_mother    DMP_12w   DMTP_12w     DEP_12w  DETP_12w sum_DEP_12w sum_DMP_12w sum_DAP_12w 
# 1  10001611 -0.2281976 1.42223300  1.65535183 -4.346027 0.020730228  0.02563002  0.04636025
# 2  10002812 -1.4952266 0.05658353 -1.90064167 -2.771830 0.002598415  0.01013192  0.01273034 
# 3  10006511  1.5224593 1.98550043 -1.76337775 -2.230663 0.003163518  0.05065476  0.05381828 
# 4  10006911  0.1381090 1.73552218 -0.05889369 -5.380845 0.006370754  0.03216126  0.03853202 
# 5  10007411  2.7333543 3.11603199 -0.16812276 -3.726509 0.006219412  0.11376163  0.11998104 
# 6  10007911  2.1984942 2.27798475  0.22650853 -4.009708 0.007957289  0.07054021  0.07849750 

#==============================================================================#----
# OP metabolites 32w (T2)                                                      #----
#==============================================================================#----

# 1. Load op metabolites at 32 weeks (T3)
load(file_op_32w) # dataAll 755 x 39 

# 2. Subset compounds of interest 
op_32w <- dataAll[,c('ID_mother', 
                     "cc.DMP_m_32w_log2i" ,
                     "cc.DMTP_m_32w_log2i",
                     "cc.DEP_m_32w_log2i",
                     "cc.DETP_m_32w_log2i",
                     'sum_DEP_32w', 
                     'sum_DMP_32w', 
                     'sum_DAP_32w')]



# 3. Rename for readability (new_name = old_name)
op_32w <- dplyr::rename(op_32w, 
                        DMP_32w = `cc.DMP_m_32w_log2i`, 
                        DMTP_32w = `cc.DMTP_m_32w_log2i`, 
                        DEP_32w = `cc.DEP_m_32w_log2i`,
                        DETP_32w = `cc.DETP_m_32w_log2i`)
head(op_32w)
# ID_mother    DMP_32w  DMTP_32w    DEP_32w   DETP_32w sum_DEP_32w sum_DMP_32w sum_DAP_32w 
# 1  10000111  0.2182288 1.2986583  0.5753123 -1.1679693 0.012284356  0.02653823  0.03882259 
# 2  10000211  0.4456484 0.5655972  4.9818527 -4.2656892 0.205367152  0.02121838  0.22658554 
# 3  10000311 -0.1371868 0.6870607 -0.4941091 -2.9963757 0.005343805  0.01854218  0.02388599 
# 4  10000611  1.4483979 2.1342209  2.4828483  0.2868811 0.043444447  0.05253991  0.09598436 
# 5  10000711  2.2387869 3.3146965  0.8237494  1.9221978 0.033757894  0.10745672  0.14121461 
# 6  10000811 -0.3408193 0.8073549 -0.2688168 -2.2780521 0.006597701  0.01857766  0.02517536 

rm(dataAll)

#==============================================================================#----
# JOIN OP METABOLITES 12w + 32w                                                #----
#==============================================================================#----

# Merge pesticides data in whole single df, based on participant IDs
op_all <- merge(op_12w, op_32w, by = 'ID_mother', all = TRUE) # 920 


# Note: How many mothers that have op metabolites data in 12w (N = 695), have also data for 32w (N = 755)
length(intersect(op_12w$ID_mother, op_32w$ID_mother)) # 530 

# That means we have: 
# - 530 mothers with both time points 
# - 165 only at 12w
# - 225 only at 32w

saveRDS(op_all , output_file)


#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----





