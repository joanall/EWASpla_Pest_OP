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
# LIBRARIES.                                                                  #----
#==============================================================================#----

load("db/filter_DNAm/filter_EPIC.RData")
source("script/00_template_plot.R")

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

#==============================================================================#----
# FILTER UNRELIABLE PROBES.                                                    #----
#==============================================================================#----

# Function to remove unreliable CpGs 
exclude_CpGs <- function(cohort, cpgid, exclude, ethnic, artype='EPIC', filename = NULL, fileresume = NULL)
{
  # Test if filter data exists and gets it
  if( toupper(artype) == '450K' ) {
    #try(data("filter_450K") )
    filters <- filter_450K
  }else if(toupper(artype) == 'EPIC') {
    #try(data("filter_EPIC"))
    filters <- filter_EPIC
  }else{
    stop( paste0( "Unknown array type ", toupper(artype) ) )
  }
  
  
  # In order to minimize data size we merge only the selected ethnic column
  # Get first ethnic column position in filters
  firstEthnicPosition <-  min(grep(paste0("MASK_snp5_.*[^GMAF1p&^common]"), colnames(filters), perl = TRUE))-1
  
  # Select ethnic columns not related with our data from filters (ONLY IF ethnic != '' or NA)
  if(ethnic!='' && !is.na(ethnic)) {
    # Get first ethnic column position in filters
    fieldstodelete <- grep(paste0("MASK_snp5_.*[^",ethnic,"]"),  colnames(filters)[grep(paste0("MASK_snp5_.*[^GMAF1p&^common]"), colnames(filters), perl = TRUE)], perl = TRUE) + firstEthnicPosition
    
  } else {
    #.Remove from code.# fieldstodelete <- grep(paste0("MASK_snp5_.*[^EUR]"),  colnames(filters)[grep(paste0("MASK_snp5_.*[^GMAF1p&^common]"), colnames(filters), perl = TRUE)], perl = TRUE) + firstEthnicPosition
    fieldstodelete <- ""
  }
  
  fieldstomerge <- which(!seq(1:dim(filters)[2]) %in% fieldstodelete)
  
  # Merge cohort with CpGs filters
  cohort <- merge(cohort, filters[,fieldstomerge], by.x= cpgid, by.y = "probeID", all.x = TRUE )
  
  
  # CpGs id to remove
  if( is.null(exclude)) {
    excludeid <- NULL
  } else {
    excludeid <- cohort[eval(parse(text=getCritera(exclude, ethnic))), cpgid]
  }
  
  
  # Report descriptive exclussions to a descriptive file
  if(!is.null(fileresume)) {
    write(sprintf('\n# %s', strrep("-",16)), file = fileresume, append = TRUE)
    write(sprintf('# Remove "problematic"  CpGs : '), file = fileresume, append = TRUE)
    write(sprintf('# %s\n', strrep("-",16)), file = fileresume, append = TRUE)
    write(sprintf('# Criteria : \n\tArray type : %s \n\tEthnia : %s \n', toupper(artype), toupper(ethnic)), file = fileresume, append = TRUE)
    write(sprintf('# Mask : %s', exclude), file = fileresume, append = TRUE)
    write(sprintf('\n'), file = fileresume, append = TRUE)
    write(sprintf('# Total CpGs in data : %d', dim(cohort)[1]), file = fileresume, append = TRUE)
    write(sprintf('# Number of excluded CpGs: %d', length(excludeid)), file = fileresume, append = TRUE)
    write(sprintf('# Total CpGs after exclusion : %d\n', (dim(cohort)[1]) - length(excludeid)), file = fileresume, append = TRUE)
    write(sprintf('# Percent excluded CpGs: %f %%\n', ((length(excludeid)/dim(cohort)[1])*100 )), file = fileresume, append = TRUE)
    
    suppressWarnings(
      write.table( cohort[eval(parse(text=getCritera(exclude, ethnic))),],
                   filename, col.names = TRUE, row.names = FALSE, sep = '\t', append = FALSE, dec='.'))
  }
  
  # Report CpG excluded and reason to a file
  if(!is.null(filename)) {
    # write(sprintf('# Criteria : %s\n', toupper(artype)), file = filename)
    # write(sprintf('# Total CpGs in data : %d', dim(cohort)[1]), file = filename, append = TRUE)
    # write(sprintf('# Number of excluded CpGs: %d', length(excludeid)), file = filename, append = TRUE)
    # write(sprintf('# Total CpGs after exclusion : %d\n', (dim(cohort)[1]) - length(excludeid)), file = filename, append = TRUE)
    # write(sprintf('# Percent excluded CpGs: %f %%\n', ((length(excludeid)/dim(cohort)[1])*100 )), file = filename, append = TRUE)
    # Report exclusion reason
    suppressWarnings(
      write.table( cohort[eval(parse(text=getCritera(exclude, ethnic))),],
                   filename, col.names = TRUE, row.names = FALSE, sep = '\t', append = FALSE, dec='.'))
  }
  
  # Remove CpGs with exclusion parameters
  cohort <- cohort[ !cohort[,cpgid] %in% excludeid,]
  
  warning( length(excludeid), ' CpGs have been excluded')
  
  return(cohort)
  
}


getCritera <- function(exclude, ethnic) {
  
  criter = ''
  possible_crit <- c( 'MASK_sub25_copy', 'MASK_sub30_copy', 'MASK_sub35_copy', 'MASK_sub40_copy',
                      'MASK_mapping', 'MASK_extBase', 'MASK_typeINextBaseSwitch', 'MASK_snp5_common', 'MASK_snp5_GMAF1p',
                      'MASK_general', 'cpg_probes', 'noncpg_probes', 'control_probes', 'Unrel_450_EPIC_blood', 'MASK_rmsk15',
                      'Sex', 'Unrel_450_EPIC_pla_restrict', 'Unrel_450_EPIC_pla', 'MASK_snp5_ethnic')
  
  if(!is.null(exclude[1]) && !is.na(exclude[1]) && exclude[1] != '') {
    # Test if all parameters are allowed
    if(length(which(! exclude %in% possible_crit))>=1)
      stop(paste0('Parameter(s) : ',paste(exclude[which(! exclude %in% possible_crit)], sep = ','),' not valid. Possible values are ',possible_crit ))
    
    # Get all exclude variables
    ##.. Works with lists -> convert list to matrix ..## exclusion.crit <- do.call(cbind, lapply( ls(patt="exclude"), get) )
    
    # Create formula with exclusion criteria
    ##.. Works with lists -> gets only parameters with value 'Exclude'..## criter <- paste(paste("cohort$",attributes(exclusion.crit)$dimnames[[1]][which(exclusion.crit=='Exclude')], sep = ""),"TRUE | ",sep=" == ", collapse = '')
    criter <- paste(paste("cohort$",exclude, sep = ""),"TRUE | ",sep=" == ", collapse = '')
    
    # If  MASK_general = 'Exclude' --> Change to :  ("MASK.sub30.copy", "MASK.mapping", "MASK.extBase", "MASK.typeINextBaseSwitch" and "MASK.snp5.GMAF1p") == TRUE
    if( length(grep("cohort$MASK_general == TRUE",criter, fixed = TRUE))>0 )
      criter <- sub("cohort$MASK_general == TRUE", " cohort$MASK_sub30_copy == TRUE | cohort$MASK_mapping == TRUE | cohort$MASK_extBase == TRUE | cohort$MASK_typeINextBaseSwitch == TRUE | cohort$MASK_snp5_GMAF1p == TRUE " , criter, fixed = TRUE)
    
    # If Sex = 'Exclude' --> Change to : ( CpG_chrm %in% "chrX" | CpG_chrm %in% "chrY" )
    if( length(grep("cohort$Sex == TRUE",criter, fixed = TRUE))>0 )
      criter <- sub("cohort$Sex == TRUE", " cohort$CpG_chrm %in% 'chrX' | cohort$CpG_chrm %in% 'chrY' " , criter, fixed = TRUE)
    
    newCpGcond <- vector()
    if( length( grep("| cohort$cpg_probes == TRUE",criter, fixed = TRUE)) >0 )
      newCpGcond <- append(newCpGcond, c('cg'))
    if( length(grep("| cohort$noncpg_probes == TRUE",criter, fixed = TRUE))>0 )
      newCpGcond <- append(newCpGcond, c('ch'))
    if( length(grep("| cohort$control_probes == TRUE",criter, fixed = TRUE))>0 )
      newCpGcond <- append(newCpGcond, c('rs'))  
    
    if(!is.null(vector())) {
      # Remove invalid criteria
      criter <- sub("| cohort$cpg_probes == TRUE ", "", criter, fixed = TRUE)
      criter <- sub("| cohort$noncpg_probes == TRUE ", "", criter, fixed = TRUE)
      criter <- sub("| cohort$control_probes == TRUE ", "", criter, fixed = TRUE)
      criter <- sub("| cohort$MASK_snp5_ethnic == TRUE ", "", criter, fixed = TRUE)
      criter <- sub("cohort$cpg_probes == TRUE |", "", criter, fixed = TRUE)
      criter <- sub("cohort$noncpg_probes == TRUE |", "", criter, fixed = TRUE)
      criter <- sub("cohort$control_probes == TRUE |", "", criter, fixed = TRUE)
      criter <- sub("cohort$MASK_snp5_ethnic == TRUE |", "", criter, fixed = TRUE)
      
      # Add new criteria
      criter <- paste0(criter, paste0(" cohort$probeType %in% '",newCpGcond,"' | ", collapse = ''))
    }
    
    # Add ethnicity exlcusion criteria
    criter <- paste0(criter, paste0(" cohort$probeType %in% '",newCpGcond,"' | ", collapse= ''))
    if(is.na(ethnic) | ethnic == ''){
      criter <- paste("which(", substr(criter,0,nchar(criter)-3),")",sep="")
    }
  }
  
  if(!is.na(ethnic) & ethnic != ''){
    criter <- paste0(criter, " cohort$MASK_snp5_", ethnic," == TRUE | ")
    criter <- paste("which(", substr(criter, 0, nchar(criter)-2), ")", sep="")
  }
  criter
}

#==============================================================================#----
# INFLATION - LAMBDA - QQ PLOT                                                 #----
#==============================================================================#----

# Genomic inflation factor (lambda) with 95% CI
inflation <- function(pvector) {
  pvector <- pvector[!is.na(pvector)]
  chisq  <- qchisq(pvector, df = 1, lower.tail = FALSE)
  lambda <- median(chisq) / qchisq(0.5, df = 1)
  se_median <- 1.253 * sd(chisq) / sqrt(length(chisq))
  lower_ci  <- lambda - qnorm(0.975) * se_median / qchisq(0.5, df = 1)
  upper_ci  <- lambda + qnorm(0.975) * se_median / qchisq(0.5, df = 1)
  c(lambda = lambda, lower_ci = lower_ci, upper_ci = upper_ci)
}

# Expected vs observed -log10(p); thinned by rounding (same look, far fewer points)
qq_data <- function(pvector) {
  pvector <- sort(pvector[!is.na(pvector)])
  qqdata <- data.frame(e = -log10(ppoints(length(pvector))),
                       o = -log10(pvector))
  unique(round(qqdata, 3))
}

# QQ plot with lambda [95% CI] in the top-left corner
qq_plot <- function(pvector, title = NULL) {
  lambda <- inflation(pvector)
  lambda_txt <- sprintf("lambda = %.2f [%.2f, %.2f]", lambda[1], lambda[2], lambda[3])
  ggplot(qq_data(pvector), aes(x = e, y = o)) +
    geom_point(size = 0.6) +
    geom_abline(intercept = 0, slope = 1, colour = pal_bright[["red"]]) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
             label = lambda_txt, size = 3.5) +
    labs(title = title,
         x = expression(Expected~-log[10](p)),
         y = expression(Observed~-log[10](p))) +
    theme_ewas() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
}



#==============================================================================#----
# VOLCANO PLOT                                                                 #----
#==============================================================================#----
# Effect (x) vs -log10(p) (y); dashed line at the Bonferroni threshold
# (alpha / n unique CpGs). CpGs coloured by direction (hyper / hypo):
#   - p_cutoff = NULL (default): CpGs above the Bonferroni threshold
#   - p_cutoff = 1e-5:           CpGs with p < p_cutoff (e.g. suggestive),
#   - df:        EWAS results (one or more time points stacked)
#   - effect_col, p_col: column names
#   - annot:     TRUE to label coloured CpGs

volcano_plot <- function(df, effect_col = "Effect", p_col = "pvalue",
                         title = NULL, annot = FALSE,
                         alpha = 0.05, p_cutoff = NULL) {
  
  df$x <- df[[effect_col]]
  df$y <- -log10(df[[p_col]])
  
  # Bonferroni threshold (unique CpGs, so stacked time points don't double it)
  thr <- -log10(alpha / length(unique(df$CpG)))
  
  # Colouring threshold: Bonferroni (default) or cutoff_col
  thr_col <- if (is.null(p_cutoff)) thr else -log10(p_cutoff)
  
  # Direction of coloured CpGs
  sig <- !is.na(df$y) & df$y > thr_col
  df$direction <- "NS"
  df$direction[sig & df$x > 0] <- "Hypermethylated"
  df$direction[sig & df$x < 0] <- "Hypomethylated"
  df$direction <- factor(df$direction, levels = names(pal_direction))
  df <- df[order(df$direction == "NS", decreasing = TRUE), ]   # coloured on top
  
  p <- ggplot(df, aes(x = x, y = y, colour = direction)) +
    geom_point(size = 0.6, alpha = 0.6) +
    geom_hline(yintercept = thr, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = 0,   linetype = "dashed", colour = "grey70") +
    scale_colour_manual(values = pal_direction, drop = FALSE) +
    labs(title = title, x = "Effect (beta coefficient)",
         y = expression(-log[10](p)), colour = NULL) +
    theme_ewas() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  # Colouring threshold line (dotted), only when cutoff_col is used
  if (!is.null(p_cutoff)) {
    p <- p + geom_hline(yintercept = thr_col, linetype = "dotted", colour = "grey40")
  }
  
  if (annot) {
    p <- p + geom_text(data = df[df$direction != "NS", ], aes(label = CpG),
                       size = 2.5, hjust = -0.1, colour = "black",
                       check_overlap = TRUE, show.legend = FALSE)
  }
  return(p)
}

#==============================================================================#----
# MANHATTAN PLOT                                                               #----
#==============================================================================#----
# Genomic position (x) vs -log10(p) (y); chromosomes in alternating greys.
# Dashed line = Bonferroni (alpha / n unique CpGs), dotted line = suggestive.
# CpGs coloured by direction (hyper / hypo):
#   - colour_sugg = FALSE (default): CpGs above the Bonferroni line
#   - colour_sugg = TRUE:            CpGs above the suggestive line
# - df:         EWAS results with columns CpG, Effect, p-value (+ gene_col)
# - anno:       annotation, rownames = CpGs, columns chr and pos
# - suggestive: p-value for the suggestive line (default 1e-5)
# - annot:      TRUE to label coloured CpGs as "CpG (gene)"; gene omitted if
#               empty. Uses ggrepel if installed (no hidden labels).

manhattan_plot <- function(df, anno, effect_col = "Effect", p_col = "pvalue",
                           gene_col = "gene_name", title = NULL, annot = FALSE,
                           suggestive = 1e-5, alpha = 0.05, colour_sugg = FALSE) {
  
  # Position from the annotation (sex chromosomes already excluded -> NA dropped)
  idx <- match(df$CpG, rownames(anno))
  df$chr_num <- suppressWarnings(as.numeric(sub("chr", "", anno$chr[idx])))
  df$pos     <- as.numeric(anno$pos[idx])
  df$y       <- -log10(df[[p_col]])
  df <- df[!is.na(df$chr_num) & !is.na(df$pos) & !is.na(df$y), ]
  
  # Cumulative position along the genome
  chr_len    <- tapply(df$pos, df$chr_num, max)
  chr_offset <- c(0, cumsum(as.numeric(chr_len))[-length(chr_len)])
  names(chr_offset) <- names(chr_len)
  df$pos_cum <- df$pos + chr_offset[as.character(df$chr_num)]
  axis_df <- data.frame(chr    = names(chr_len),
                        center = tapply(df$pos_cum, df$chr_num, function(x) (min(x) + max(x)) / 2))
  
  # Thresholds (unique CpGs, so stacked time points don't double it)
  thr_bonf <- -log10(alpha / length(unique(df$CpG)))
  thr_sugg <- -log10(suggestive)
  thr_col  <- if (colour_sugg) thr_sugg else thr_bonf
  
  # Colour: alternating greys, CpGs above the colouring threshold by direction
  df$sig <- df$y > thr_col
  df$colour_group <- ifelse(df$chr_num %% 2 == 0, "even", "odd")
  df$colour_group[df$sig & df[[effect_col]] > 0] <- "Hypermethylated"
  df$colour_group[df$sig & df[[effect_col]] < 0] <- "Hypomethylated"
  df <- df[order(df$sig), ]                                    # coloured on top
  
  colour_values <- c(odd = "grey35", even = "grey70",
                     pal_direction[c("Hypermethylated", "Hypomethylated")])
  
  p <- ggplot(df, aes(x = pos_cum, y = y, colour = colour_group)) +
    geom_point(size = 0.5) +
    geom_hline(yintercept = thr_bonf, linetype = "dashed", colour = "grey20") +
    geom_hline(yintercept = thr_sugg, linetype = "dotted", colour = "grey50") +
    scale_colour_manual(values = colour_values,
                        breaks = c("Hypermethylated", "Hypomethylated")) +
    scale_x_continuous(breaks = axis_df$center, labels = axis_df$chr,
                       expand = c(0.01, 0)) +
    labs(title = title, x = "Chromosome", y = expression(-log[10](p)), colour = NULL) +
    theme_ewas() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7))
  
  if (annot && any(df$sig)) {
    df_lab <- df[df$sig, ]
    
    # Label: "CpG (gene)"; genes de-duplicated ("A;A;B" -> "A;B"), omitted if empty
    if (gene_col %in% colnames(df_lab)) {
      genes <- as.character(df_lab[[gene_col]])
      genes[is.na(genes)] <- ""
      genes <- sapply(strsplit(genes, ";"), function(g) paste(unique(g[g != ""]), collapse = ";"))
    } else {
      genes <- rep("", nrow(df_lab))
    }
    df_lab$label <- ifelse(genes == "", df_lab$CpG, paste0(df_lab$CpG, " (", genes, ")"))
    
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(data = df_lab, aes(label = label), size = 2.3,
                                        colour = "black", segment.colour = "grey60",
                                        segment.size = 0.2, min.segment.length = 0,
                                        max.overlaps = Inf, show.legend = FALSE)
    } else {
      p <- p + geom_text(data = df_lab, aes(label = label), size = 2.3,
                         vjust = -0.8, colour = "black",
                         check_overlap = TRUE, show.legend = FALSE)
    }
  }
  return(p)
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----

