#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: April - October 2026                                                                              
#==============================================================================#

# [DESCRIPTION]: EWAS analyses to explore the assoctiation between the exposure
# to Pesticides and OP during pregnancy and placental methylation 

# [INPUT]: Write input type file needed and formatting.  
# [OUTPUT]: Write output/s obtained after running code. 





#================================================================================
# [DESCRIPTION]: EWAS analyses to explore the assoctiation between the exposure
# to Pesticides and OP during pregnancy and placental methylation

# Note: The same code will be used for SEPAGES cohort. Scripts from part 0 will 
# be adapted to the 

# [INPUT BISC]:
Pest_OP + EWAS  -> ../db/Merged_Pest_OP/Pest_OP_12w_or_32w_Full_EWAS.rds
DNA meth data -> ../db/methyaltion_data/bVals_BISC_noout_nomask.RDS
General Covars -> "db/general_covars_data/covariables_data_02012025.rds"



#==============================================================================#
#                                                                              #
#==============================================================================#