#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: April - October 2026                                                                              
#==============================================================================#

# [DESCRIPTION]: Get "Pesticides" Exposure data from BiSC cohort from 12w (N = 695) 
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

# Values: 
# - 2,4-D — dichotomized
# - TCPy — continuous
# - 3-PBA — continuous and dichotomized

# [CONSIDERATIONS]: 
# - If possible better use the continuous variable. 
# - Compounds have been log transformed to follow normal distribution. 
# - No missing in the log transformed or dict vars, imputation done. Check QC. 

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

# Load libraries

# Define input path files 
file_pest_12w <- "../db/chemicals/pesticides_pah_urine/BISC_db_pest_m_12w_20240729.rdata"
file_pest_32w <- "../db/chemicals/pesticides_pah_urine/BISC_db_pest_m_32w_20240729.rdata"

# Define output path files 
output_file <- "db/all_pest_12w_32w.rds"

#==============================================================================#----
# PESTICIDES 12w (T2)                                                          #----
#==============================================================================#----

# 1. Load pesticides at 12 weeks (T2)
load(file_pest_12w) # dataAll 695 x 56

# 2. Subset compounds of interest 
# Note: If a compounds have both dichotomous and continuous, we pick the cont. 
# Note: PBA we choose cont. 
pest_12w <- dataAll[,c('ID_mother', 
                       "cc.2,4-D_m_12w_dic" ,
                       "cc.TCPy_m_12w_log2i",
                       "cc.3-PBA_m_12w_log2i", 
                       "ur_pool_grav_m_12w")]


# 3. Rename for readability (new_name = old_name)
pest_12w <- dplyr::rename(pest_12w, 
                          Diclor_12w = `cc.2,4-D_m_12w_dic`, 
                          TCPy_12w = `cc.TCPy_m_12w_log2i`, 
                          PBA_12w =  `cc.3-PBA_m_12w_log2i`, 
                          SG_12w = `ur_pool_grav_m_12w`) 



head(pest_12w)
# ID_mother Diclor_12w   TCPy_12w     PBA_12w SG_12w
# 1  10001611 undetected  0.8900525 -0.01525817     NA
# 2  10002812 undetected -1.5017996  0.25380834 1.0132
# 3  10006511   detected  0.6184253 -1.18736483 1.0194
# 4  10006911   detected -1.0474206 -1.24601212 1.0205
# 5  10007411 undetected  1.3518927  0.32914592     NA
# 6  10007911 undetected  0.7165077  0.26392610 1.0165

#==============================================================================#----
# PESTICIDES 32w (T2)                                                          #----
#==============================================================================#----

# 1. Load pesticides at 32 weeks (T3)
load(file_pest_32w) # dataAll 755 x 37

# 2. Subset compounds of interest 
# Note: If a compounds have both dichotomous and continous, we pick the cont. 
# Note: PBA we choose cont. 

pest_32w <- dataAll[ ,c('ID_mother',  
                       "cc.2,4-D_m_32w_dic" ,
                       "cc.TCPy_m_32w_log2i",
                       "cc.3-PBA_m_32w_log2i", 
                       "ur_pool_grav_m_32w")]

# 3. Rename for readability  (new_name = old_name)
pest_32w <- dplyr::rename(pest_32w, 
                          Diclor_32w = `cc.2,4-D_m_32w_dic`, 
                          TCPy_32w = `cc.TCPy_m_32w_log2i`, 
                          PBA_32w =  `cc.3-PBA_m_32w_log2i`, 
                          SG_32w = `ur_pool_grav_m_32w`)

head(pest_32w)
# ID_mother Diclor_32w   TCPy_32w     PBA_32w SG_32w
# 1  10000111 undetected  0.3724792 -4.17438903 1.0120
# 2  10000211   detected  1.0484010  0.98681908 1.0173
# 3  10000311 undetected -1.7176055 -3.47714779 1.0050
# 4  10000611   detected  0.9485731  0.08158072 1.0188
# 5  10000711   detected  2.0139919  0.28077865 1.0155
# 6  10000811   detected -0.6063034  0.35164410 1.0151

rm(dataAll)

#==============================================================================#----
# JOIN PESTICIDES 12w + 32w                                                    #----
#==============================================================================#----

# Merge pesticides data in whole single df, based on participant IDs
pest_all <- merge(pest_12w, pest_32w, by = 'ID_mother', all = TRUE) # 920 


# Note: How many mothers that have pesticides data in 12w (N=695), have also data for 32w (N = 755)
length(intersect(pest_12w$ID_mother, pest_32w$ID_mother)) # 530 

# That means we have: 
# - 530 mothers with both time points 
# - 165 only at 12w
# - 225 only at 32w

saveRDS(pest_all , output_file)

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----





