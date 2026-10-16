#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#

# [DESCRIPTION]: Explore inflation of the EWAS results: genomic inflation
# factor (lambda) with 95% CI and QQ plots, for each exposure and time point.

# [INPUT]:  results/3_Clean_Results/main_t{2,3}/{24d,3pba,tcpy,dmp,dep,dap}_t{2,3}.rds
#           data.frame, rows = CpGs: Effect, SE, tvalue, pvalue, CpG, ...

# [OUTPUT]: results/3_Clean_Results/main_t{2,3}/<exposure>_t{2,3}_qqplot.png
#           results/3_Clean_Results/qqplots_all_t2.png, qqplots_all_t3.png
#           results/3_Clean_Results/lambda_main.csv / .png


#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
source("script/00_template_plot.R")
source("script/00_functions.R")

input_folder <- "results/3_Clean_Results/"

exposures  <- c("24d",   "3pba",  "tcpy", "dmp",       "dep",       "dap")
exp_labels <- c("2,4-D", "3-PBA", "TCPy", "DMP (sum)", "DEP (sum)", "DAP (sum)")

# Check: all input files exist
files_in <- c(paste0(input_folder, "main_t2/", exposures, "_t2.rds"),
              paste0(input_folder, "main_t3/", exposures, "_t3.rds"))
files_in[!file.exists(files_in)]   # should be character(0)

#==============================================================================#----
# LAMBDA & QQ PLOTS                                                            #----
#==============================================================================#----

lambda_tab <- NULL
qq_all     <- NULL


for (tp in c("t2", "t3")) {
  for (i in seq_along(exposures)) {
    
    file_in <- paste0(input_folder, "main_", tp, "/", exposures[i], "_", tp, ".rds")
    cat("Exposure:", exposures[i], "| Time:", tp, "\n")
    ewas_results <- readRDS(file_in)
    
    # Lambda and 95% CI
    lambda_ci <- inflation(ewas_results$pvalue)
    lambda_tab <- rbind(lambda_tab, data.frame(
      Time     = tp,
      Exposure = exp_labels[i],
      N_CpGs   = sum(!is.na(ewas_results$pvalue)),
      Lambda   = round(lambda_ci[["lambda"]],   3),
      Lower_CI = round(lambda_ci[["lower_ci"]], 3),
      Upper_CI = round(lambda_ci[["upper_ci"]], 3)
    ))
    
    # Plot 1. QQ plot of observed (y) vs expected (x) -log10 p-values for one
    # exposure and time point; lambda [95% CI] in the corner.
    p_qq <- qq_plot(ewas_results$pvalue, title = paste(exp_labels[i], "-", tp))
    save_plot(p_qq, paste0(input_folder, "main_", tp, "/", exposures[i], "_", tp, "_qqplot.png"),
              width = 12, height = 12)
    
    # Data for the combined figure
    qq_all <- rbind(qq_all, data.frame(Time = tp, Exposure = exp_labels[i],
                                       qq_data(ewas_results$pvalue)))
  }
}

lambda_tab

# Table. Genomic inflation factor (lambda) and 95% CI per exposure and time point.
save_table(lambda_tab, paste0(input_folder, "lambda_main"))

#==============================================================================#----
# COMBINED QQ PLOTS - ONE FIGURE PER TIME POINT                                #----
#==============================================================================#----

qq_all$Exposure     <- factor(qq_all$Exposure,     levels = exp_labels)
lambda_tab$Exposure <- factor(lambda_tab$Exposure, levels = exp_labels)
lambda_tab$label    <- sprintf("lambda = %.2f [%.2f, %.2f]",
                               lambda_tab$Lambda, lambda_tab$Lower_CI, lambda_tab$Upper_CI)

for (tp in c("t2", "t3")) {
  
  # Plot 2. QQ plots of the six exposures at one time point (one panel each);
  # lambda [95% CI] in each panel.
  p_all <- ggplot(qq_all[qq_all$Time == tp, ], aes(x = e, y = o)) +
    geom_point(size = 0.5) +
    geom_abline(intercept = 0, slope = 1, colour = pal_bright[["red"]]) +
    geom_text(data = lambda_tab[lambda_tab$Time == tp, ], aes(label = label),
              x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 3) +
    facet_wrap(~ Exposure, ncol = 3) +
    labs(title = paste("QQ plots - main models", tp),
         subtitle = "SEPAGES, robust linear regression",
         x = expression(Expected~-log[10](p)),
         y = expression(Observed~-log[10](p))) +
    theme_ewas() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  print(p_all)
  save_plot(p_all, paste0(input_folder, "qqplots_all_", tp, ".png"), width = 24, height = 16)
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----