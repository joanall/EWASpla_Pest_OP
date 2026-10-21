#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: Write the list of cleaned EWAS result files (one path per
# line) and an index table describing each file, for the downstream scripts
# (results tables, plots, meta-analysis).
#
# Argument (optional): pattern for the run folders to include (regular
# expression; default: all), e.g. "sepages_2tp_t[23]_old"
#   Rscript script/3_Clean_Results/3_Get_Cleaned_Files.R "sepages_2tp_t[23]_old"
#
# [INPUT]:  results/3_Clean_Results/files_cleaned.txt
# [OUTPUT]: dmrs 

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

library(dmrff)







meth_data <- as.matrix(meth_data)

# Get DMRs 
# p.cutoff Unadjusted p-value cutoff for membership in a candidate DMR
# (Default: 0.05).

dmr <- dmrff (estimate = result_match$Effect,
              se = result_match$SE,
              p.value = result_match$pvalue,
              methylation = meth_data,
              chr = annotation$chr,
              pos = annotation$pos,
              p.cutoff = 0.05,
              maxgap = 500,
              verbose=T)

#write_rds(dmr.model, paste0("results/6_dmr.ewas.model.",ewaspcutoff,"cutoff.3.B.rds"))

saveRDS(dmr, file =  save_to[[which_file]])
