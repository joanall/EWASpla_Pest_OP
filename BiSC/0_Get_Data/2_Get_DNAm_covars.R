#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: April - October 2026                                                                              
#==============================================================================#

# [DESCRIPTION]: Join cell-type proportions (ct) + DNAm covariates (dnam_cov) using ids dictionary

# [INPUT]: 2 .rds files, 
#  - ct has: basename as rownames (sample ID used in methylation files)
#  - dnam_cov has: ID_mother (mother ID)
#  - ids has: BOTH basename and ID_mother (mapping table)


# [OUTPUT]: 1.rds file, Participants and covars dna m ct and others 


#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

# Define input path files 
infile_ct <- "../db/methylation_data/celltype_BISC.RDS"
infile_dnam_cov <- "../db/methylation_data/covars_meth_BISC.RDS"
infile_ids <- "../db/methylation_data/ids_dictionary.RDS"

# Define output path files 
outfile <- "db/0_Get_Data/all_dnam_cov.rds"

#==============================================================================#----
# JOIN DBs                                                                     #----
#==============================================================================#----

# Load data
ct  <- readRDS(infile_ct) # 565 
dnam_cov  <- readRDS(infile_dnam_cov) # 565 
ids <- data.frame(readRDS(infile_ids)) # 565 

head(ids)
# id_mother            Basename id_child id_participant.c
# 1  12017711 206116800107_R07C01 12017731    BISC_12017731
# 2  12014111 206116800144_R02C01 12014131    BISC_12014131
# 3  12020311 206129780018_R04C01 12020331    BISC_12020331
# 4  12007411 206129780004_R06C01 12007431    BISC_12007431
# 5  12015511 206119350019_R05C01 12015531    BISC_12015531
# 6  10011811 206129780134_R04C01 10011831    BISC_10011831


# Add ID_mother to ct using basename
ct$Basename <- rownames(ct)
ct$id_mother <- ids$id_mother[match(ct$Basename, ids$Basename)]

# Merge ct + dnam_cov by Basename
dnam_cov <- dnam_cov[, c('id_mother', 'dna_contamination')]
df <- merge(ct, dnam_cov, by = 'id_mother') # 565 13

# Rename id_mother to "id_bisc"
names(df)[names(df) == "id_mother"] <- "id_bisc"

head(df) 
# id_bisc Trophoblasts      Stromal     Hofbauer Endothelial        nRBC Syncytiotrophoblast            Basename dna_contamination
# 1 10000311   0.13271524 8.791965e-21 0.000000e+00  0.09314076 0.037171405           0.7343509 206250280051_R01C01         -4.896945
# 2 10000611   0.22049676 0.000000e+00 2.704910e-19  0.10624197 0.013446614           0.6469980 206129780091_R05C01         -5.206538
# 3 10000811   0.05093993 0.000000e+00 6.002672e-20  0.04586278 0.004284509           0.8739253 206129780074_R07C01         -5.082517
# 4 10000911   0.14810667 4.042094e-20 0.000000e+00  0.08377731 0.004317056           0.7433552 206116790020_R04C01         -5.148032
# 5 10001011   0.00000000 0.000000e+00 4.078985e-19  0.05319040 0.016716840           0.9084622 206119350020_R02C01         -5.245660
# 6 10001211   0.01193852 1.008901e-03 8.712626e-19  0.03815969 0.023970944           0.9249219 206250280057_R01C01         -5.120396

# Save 
saveRDS(df, outfile)

#==============================================================================#----
# END OF CODE                                                                  #----
#==============================================================================#----