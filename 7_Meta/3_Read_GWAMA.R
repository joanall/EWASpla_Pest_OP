#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date:  October - November 2026
#==============================================================================#
# [DESCRIPTION]:  Reads the GWAMA outputs, keeps CpGs present in both cohorts, adds
# Bonferroni, and saves one .rds per meta-analysis (named list by exposure,
# same structure as the cohort clean results) plus a summary table.

# INPUT : results/meta/2_Run_GWAMA/<meta_id>/meta_<meta_id>_<exposure>.out
# OUTPUT: results/meta/3_Read_GWAMA/meta_<meta_id>.rds
#         results/meta/3_Read_GWAMA/meta_summary.csv
#==============================================================================#

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
dir_input <- "results/7_Meta/1_Input_GWAMA"
dir_run   <- "results/7_Meta/2_Run_GWAMA"
dir_read  <- "results/7_Meta/3_Read_GWAMA"
dir.create(dir_read, recursive = TRUE, showWarnings = FALSE)


files_out <- list.files(dir_run, pattern = "\\.out$", recursive = TRUE, full.names = TRUE)
files_out <- files_out[grepl("^meta_[^.]+\\.out$", basename(files_out))]   # drops .err.out, .gc.out
stopifnot(length(files_out) > 0)
cat(length(files_out), "GWAMA outputs found\n")                             # expected: 18
print(basename(files_out))

#==============================================================================#----
# READ GWAMA RESULTS                                                           #----
#==============================================================================#----
meta_list    <- list()
summary_rows <- list()

for (file_out in files_out) {
  meta_id  <- basename(dirname(file_out))
  exposure <- sub(paste0("^meta_", meta_id, "_(.*)\\.out$"), "\\1", basename(file_out))
  
  meta <- as.data.frame(data.table::fread(file_out))
  n_cpg_total <- nrow(meta)
  meta <- meta[meta$n_studies == 2, ]        # CpGs in both cohorts only
  n_cpg_both  <- nrow(meta)
  
  # Rename GWAMA columns (underscores, no special characters)
  names(meta)[names(meta) == "rs_number"]       <- "probeID"
  names(meta)[names(meta) == "p-value"]         <- "p_value"
  names(meta)[names(meta) == "_-log10_p-value"] <- "log10_p"
  names(meta)[names(meta) == "q_p-value"]       <- "q_p_value"
  names(meta)[names(meta) == "beta_95L"]        <- "beta_95_low"
  names(meta)[names(meta) == "beta_95U"]        <- "beta_95_up"
  meta <- meta[, !names(meta) %in% c("reference_allele", "other_allele", "eaf")]  # dummy alleles
  cols_needed  <- c("probeID", "beta", "se", "p_value", "q_p_value", "i2", "effects")
  cols_missing <- setdiff(cols_needed, names(meta))
  if (length(cols_missing) > 0) {
    stop(basename(file_out), " - missing: ", paste(cols_missing, collapse = ", "),
         "\nColumns in file: ", paste(names(meta), collapse = ", "))
  }  
  # Multiple testing within each meta-analysis x exposure
  meta$p_bonf <- p.adjust(meta$p_value, method = "bonferroni")
  meta <- meta[order(meta$p_value), ]
  rownames(meta) <- meta$probeID
  
  # Inline check: GWAMA beta vs manual inverse-variance weighted beta (top CpG)
  gwama_bisc    <- read.delim(file.path(dir_input, meta_id, paste0("bisc_gwamain_",    exposure, ".txt")))
  gwama_sepages <- read.delim(file.path(dir_input, meta_id, paste0("sepages_gwamain_", exposure, ".txt")))
  cpg_top    <- meta$probeID[1]
  in_bisc    <- gwama_bisc[gwama_bisc$MARKERNAME == cpg_top, ]
  in_sepages <- gwama_sepages[gwama_sepages$MARKERNAME == cpg_top, ]
  w_bisc     <- 1 / in_bisc$SE^2
  w_sepages  <- 1 / in_sepages$SE^2
  beta_ivw   <- (in_bisc$BETA * w_bisc + in_sepages$BETA * w_sepages) / (w_bisc + w_sepages)
  stopifnot(abs(beta_ivw - meta$beta[1]) < 1e-4)
  
  # Genomic inflation
  lambda <- median(qchisq(meta$p_value, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1)
  
  if (is.null(meta_list[[meta_id]])) meta_list[[meta_id]] <- list()
  meta_list[[meta_id]][[exposure]] <- meta
  
  summary_rows[[paste(meta_id, exposure)]] <- data.frame(
    meta_id      = meta_id,
    exposure     = exposure,
    n_cpg_total  = n_cpg_total,
    n_cpg_both   = n_cpg_both,
    lambda       = round(lambda, 3),
    n_p_1e5      = sum(meta$p_value < 1e-05),
    n_bonf_05    = sum(meta$p_bonf  < 0.05),
    top_cpg      = cpg_top,
    top_p        = signif(meta$p_value[1], 3),
    stringsAsFactors = FALSE)
}

#==============================================================================#----
# SAVE                                                                         #----
#==============================================================================#----
for (meta_id in names(meta_list)) {
  saveRDS(meta_list[[meta_id]], file.path(dir_read, paste0("meta_", meta_id, ".rds")))
}

meta_summary <- do.call(rbind, summary_rows)
rownames(meta_summary) <- NULL
write.csv(meta_summary, file.path(dir_read, "meta_summary.csv"), row.names = FALSE)
print(meta_summary)