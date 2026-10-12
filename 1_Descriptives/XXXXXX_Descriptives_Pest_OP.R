#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: April - October 2026                                                                              
#==============================================================================#

# [DESCRIPTION]: Get Pesticides Exposure data from BiSC cohort from 12w (N = 695) 
# and 32w (N = 755) and built a unique dataframe. Additionally, compute the averages. 

# [INPUT]: file.rdata with xxx xxxx 


# [OUTPUT]: Write output/s obtained after running code. 



#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

# Load libraries

# Define input path files 
file_pest_12w <- "../db/chemicals/pesticides_pah_urine/BISC_db_pest_m_12w_20240729.rdata"
file_pest_32w <- "../db/chemicals/pesticides_pah_urine/BISC_db_pest_m_32w_20240729.rdata"

# Define output path files 
output_file <- ""

#==============================================================================#----
# OVERLAP SAMPLES X EXPOSURE 12 w                                            #----
#==============================================================================#----

#==============================================================================#----
# OVERLAP SAMPLES X EXPOSURE 32 w                                            #----
#==============================================================================#----



#==============================================================================#----
# OVERLAP SAMPLES X EXPOSURE AVG                                             #----
#==============================================================================#----



# Descro`toves `
# How many mothers that have pesticides data in 12w (N=695), have also data for 32w (530)
length(intersect(pest_12w$ID_mother, pest_32w$ID_mother)) # 530 

# That means we have: 
# - 530 mothers with both time points 
# - 165 only at 12w
# - 225 only at 32w
