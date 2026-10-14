#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: April - October 2026                                                                              
#==============================================================================#

# [DESCRIPTION]: Join all DBs needed to perfom EWAS analsysis on Pesticides - DNAm.
# Join exposure data (pest + op), covariates (general) + methylation covars
# [INPUT]: 2 .rds files, 
#  - ct has: basename as rownames (sample ID used in methylation files)
#  - dnam_cov has: ID_mother (mother ID)
#  - ids has: BOTH basename and ID_mother (mapping table)


# [OUTPUT]: 1.rds file, Partiipants and coavrs dna m ct and others 


#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

# Define input & output path files 
infile_chem <- "db/all_pest_op_avg.rds"
infile_cov <- "db/all_covars.rds"
infile_dnam <- "db/all_dnam_cov.rds"
outfile <- "db/final_df_incomplete_cases.rds"

# Load files 
chem <- readRDS(infile_chem) # 920  27 
cov <- readRDS(infile_cov) # 1080 15
dnam <- readRDS(infile_dnam) # 565 9

# Join all - Keep all participants 
final_df <- merge(chem, cov, by.x = 'ID_mother', by.y = "id_bisc", all = TRUE) # 1080 41
final_df <- merge(final_df, dnam, by.x = 'ID_mother', by.y = 'id_bisc', all = TRUE) # 1080 49 

# Save 
saveRDS(final_df, outfile)

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----

