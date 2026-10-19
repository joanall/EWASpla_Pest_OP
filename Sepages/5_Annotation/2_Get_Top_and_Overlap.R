#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#

# [DESCRIPTION]: Manhattan plot for each EWAS result. Bonferroni-significant
# CpGs coloured by direction; dashed line = Bonferroni, dotted = suggestive
# (p < 1e-5).
#
# [INPUT]:  results/4_Annotation/main_t{2,3}/<exposure>_t{2,3}_annotated.rds
#           (all CpGs, ordered by p-value)

# [OUTPUT]: results/4_Analyze_Results/main_t{2,3}/<exposure>_t{2,3}_manhattan.png
#           results/4_Analyze_Results/manhattan_<exposure>_t2_t3.png (t2 above t3)

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
source("script/00_template_plot.R")
source("script/00_functions.R")


#==============================================================================#----
#   GET TOP (suggestive CpGs for each exposure)                                #----
#==============================================================================#----

# outputs:
# 1. rds and csv with suggestive cpgs, pvalue, effect
# 2. rds and csv with suggestive cpgs, pvalue, effect, gene that overlap between exposure
# 3. table 2 but in venn plot to see overlaprs or dotlpot 
# 4. for table 1 look at the top cpgs in the senseictivity analsysi so 
# we will obtain a table with 



# sav datafra,e of top cpgs for each exposure 


# lpok at cpg overlap o the top cpgs and at the gene overlap between exposures 


infile <- "results/4_Annotation/"

exposures <- c("24d", "3pba", "tcpy", "dmp", "dep", "dap")
