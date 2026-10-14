#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: Functions shared by the EWAS scripts. Load with:
# source("script/00_functions.R")

#==============================================================================#----
# EWAS FUNCTION                                                                #----
#==============================================================================#----
# Robust linear regression (MASS::rlm) for each CpG in a chunk.
# - cpg_chunk:  matrix, rows = CpGs, columns = samples (same order as metadata)
# - exposures:  character vector of exposure column names in metadata
# - covariates: character vector of covariate column names in metadata
# - metadata:   data.frame, rows = samples
# Returns a list with one data.frame per exposure:
#   Effect, SE, tvalue, pvalue (robust F-test, sfsmisc::f.robftest), CpG

rlm_single_cpg <- function(cpg_chunk, exposures, covariates, metadata) {
  
  df_exp_cov <- metadata[, c(exposures, covariates), drop = FALSE]
  f <- as.formula(paste("meth ~", paste(c(exposures, covariates), collapse = " + ")))
  
  # Pre-allocate one results data.frame per exposure
  n_cpg <- nrow(cpg_chunk)
  exposure_results <- lapply(exposures, function(e) {
    data.frame(Effect = rep(NA_real_, n_cpg), SE = NA_real_, tvalue = NA_real_,
               pvalue = NA_real_, CpG = rownames(cpg_chunk))
  })
  names(exposure_results) <- exposures
  
  for (i in seq_len(n_cpg)) {
    
    df_for_model <- df_exp_cov
    df_for_model$meth <- as.numeric(cpg_chunk[i, , drop = FALSE])
    
    # A failing CpG gets NA instead of stopping the chunk
    fit <- tryCatch(MASS::rlm(f, data = df_for_model, maxit = 200),
                    error = function(e) NULL)
    if (is.null(fit)) next
    
    coefs <- summary(fit)$coefficients
    for (e in exposures) {
      exposure_results[[e]][i, c("Effect", "SE", "tvalue")] <- coefs[e, 1:3]
      exposure_results[[e]]$pvalue[i] <- sfsmisc::f.robftest(fit, var = e)$p.value
    }
  }
  
  return(exposure_results)
}


