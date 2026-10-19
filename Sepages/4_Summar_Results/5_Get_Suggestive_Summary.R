#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#

# [DESCRIPTION]: Built a summary tables showing suggestive resutls 

# [INPUT]:  results/4_Annotation/main_t{2,3}/<exposure>_t{2,3}_annotated.rds



# Table 1-6 rds and csv with suggestive cpgs, pvalue, effect etc one table for each exposure
# Table and  plots 7-11. For each suggestive cpgs look at. the sesnetivities

# Table 7. rds and csv with suggestive cpgs, pvalue, effect and gene that overlap between exposures

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_template_plot.R")
source("script/00_functions.R")


#           results/4_Annotation/main_t{2,3}/<exposure>_t{2,3}_annotated.rds
