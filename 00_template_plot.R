#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#

# [DESCRIPTION : Shared style for all plots and tables in the project, so every
# figure has the same layout, fonts, colours and size. Load it at the top of
# any plotting script with:
# source("script/00_template_plot.R")
# Then add theme_ewas() to each ggplot, use the palettes below for colours
# (no colours defined inside the scripts), save plots with save_plot() (one
# figure) or save_plots_pdf() (several pages), and tables with save_table()
# (.xlsx, never .csv: csv breaks with the Spanish/French decimal comma in Excel).
# Only files that other scripts read back (indexes, GWAMA input, percentile
# grid, top-hit lists) stay as plain text.

#==============================================================================#----
# PLOT                                                                         #----
#==============================================================================#----
# Example: how each plot in a script should look
#
# # Plot 1. Short description: what is on the x-axis, y-axis, colour and text.
# p1 <- ggplot(df, aes(x = Covariate, y = Component, fill = negLog10P)) +
#   geom_tile(colour = "white") +
#   fill_heat +
#   labs(title = "Title", subtitle = "Subtitle", x = "x-axis", y = "y-axis") +
#   theme_ewas()
# save_plot(p1, "results/<folder>/<name>.png")

# Example: how each table should be saved
#
# # Table 1. Short description: rows, columns and what the values are.
# save_table(tab1, "results/<folder>/<name>")   # writes .xlsx (+ .png if small)
# save_table(list(sheet_a = tab_a, sheet_b = tab_b), "results/<folder>/<name>")   # one .xlsx, two sheets


#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

library(ggplot2)
library(khroma)
library(colorspace)
library(gridExtra)
library(tableone)
library(writexl)     # .xlsx tables

# Default output sizes (cm) and resolution
plot_width  <- 16     # single-column figure ~ 8.5 cm, full width ~ 17 cm
plot_height <- 12
plot_dpi    <- 300

#==============================================================================#----
# THEME - PLOTS                                                                #----
#==============================================================================#----
# Base theme for all plots: white background, no grid, consistent text sizes
theme_ewas <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 2),
      plot.subtitle    = element_text(size = base_size - 1, colour = "grey30"),
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 1, colour = "black"),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      panel.grid       = element_blank(),
      panel.border     = element_rect(colour = "grey70", fill = NA, linewidth = 0.4),
      legend.position  = "right",
      legend.title     = element_text(size = base_size - 1),
      legend.text      = element_text(size = base_size - 2),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      plot.margin      = margin(10, 10, 10, 10)
    )
}

#==============================================================================#----
# THEME - TABLES                                                               #----
#==============================================================================#----
# Same font sizes as plots; bold header, light grey header background,
# alternating white / very light grey rows
theme_table <- function(base_size = 11) {
  ttheme_minimal(
    base_size = base_size - 1,
    core    = list(fg_params = list(hjust = 0, x = 0.05, col = "black"),
                   bg_params = list(fill = c("white", "grey97"), col = NA)),
    colhead = list(fg_params = list(fontface = "bold", hjust = 0, x = 0.05),
                   bg_params = list(fill = "grey90", col = NA)),
    rowhead = list(fg_params = list(fontface = "plain", hjust = 0, x = 0.05))
  )
}

#==============================================================================#----
# COLOURS                                                                      #----
#==============================================================================#----
# Palettes from khroma (Paul Tol); all available palettes and hex codes in
# khroma_palettes_hex.txt. Mothers' palettes are used for mp and m.
# Colour-blind safe. Use the named vectors below so the same category always
# has the same colour across figures.

# Paul Tol "bright" (qualitative)
pal_bright <- c(blue = "#4477AA", cyan = "#66CCEE", green = "#228833",
                yellow = "#CCBB44", red = "#EE6677", purple = "#AA3377",
                grey = "#BBBBBB")

# Cohorts
pal_cohort <- c(BiSC = "#4477AA", SEPAGES = "#EE6677")

# Time points (avg = average of t2 and t3)
pal_time <- c(t2 = "#66CCEE", t3 = "#AA3377", avg = "grey20")

# Source of an estimate in forest plots: meta-analysis + the two cohorts
pal_source <- c("Meta-analysis" = "black", pal_cohort)

# Models: main + sensitivity analyses s1-s6
pal_models <- c(main = "grey20",
                s1 = "#4477AA", s2 = "#66CCEE", s3 = "#228833",
                s4 = "#CCBB44", s5 = "#AA3377", s6 = "#EE6677", s7 = "#882255")

# Quantile cut-offs (mixture models)
pal_cutoff <- c("Tertile cut-off" = "#EE6677", "Quartile cut-off" = "#4477AA")

# Direction of effect (EWAS: hypo / hyper / not significant)
pal_direction <- c(Hypomethylated  = "#4477AA",
                   Hypermethylated = "#EE6677",
                   NS              = "#BBBBBB")

# Exposures (one colour each, fixed across all figures); names = labels used in plots
exp_labels   <- c("24d" = "2,4-D", "3pba" = "3-PBA", "tcpy" = "TCPy",
                  "dmp" = "ΣDMP",  "dep"  = "ΣDEP",  "dap"  = "ΣDAP", "mix" = "Mixture")
pal_exposure <- c("2,4-D" = "#4477AA", "3-PBA" = "#66CCEE",
                  "TCPy"  = "#228833", "ΣDMP"  = "#CCBB44",
                  "ΣDEP"  = "#EE6677", "ΣDAP"  = "#AA3377", "Mixture" = "grey20")

# Lighter versions (e.g. fill behind points / confidence bands)
pal_cohort_light <- lighten(pal_cohort, amount = 0.5)

# Ready-to-use scales
fill_heat      <- scale_fill_YlOrBr()                   # sequential (e.g. -log10 p)
fill_diverging <- scale_fill_sunset(midpoint = 0)       # diverging (e.g. correlations, betas)
colour_cohort    <- scale_colour_manual(values = pal_cohort)
fill_cohort      <- scale_fill_manual(values = pal_cohort)
colour_time      <- scale_colour_manual(values = pal_time)
colour_direction <- scale_colour_manual(values = pal_direction)
colour_exposure  <- scale_colour_manual(values = pal_exposure)
fill_exposure    <- scale_fill_manual(values = pal_exposure)

#==============================================================================#----
# SAVE                                                                         #----
#==============================================================================#----

# Plots: same size and resolution for every figure (override per plot if needed)
save_plot <- function(p, file, width = plot_width, height = plot_height, dpi = plot_dpi) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(file, plot = p, width = width, height = height, units = "cm",
         dpi = dpi, bg = "white")
}

# Several ggplots in one PDF, one plot per page (e.g. paginated panels)
save_plots_pdf <- function(plots, file, width = plot_width, height = plot_height) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  pdf(file, width = width / 2.54, height = height / 2.54)   # cm -> inches
  for (p in plots) print(p)
  dev.off()
  invisible(file)
}

# Tables: always .xlsx (opens in Excel without separator / decimal problems)
# + .png of small tables (for slides / quick checks).
# - tab:  data.frame or matrix (e.g. print(CreateTableOne(...))), or a named
#         list of them -> one sheet each
# - file: path without extension (an extension, if given, is dropped)
# - png:  NULL = only for a single table of <= 40 rows; TRUE / FALSE to force
# Row names are kept as a first column when they are not just 1, 2, 3...
save_table <- function(tab, file, base_size = 11, png = NULL) {
  file <- sub("\\.(xlsx|csv|png|txt)$", "", file)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  
  one_table <- is.data.frame(tab) || is.matrix(tab)
  sheets <- if (one_table) list(Sheet1 = tab) else tab
  stopifnot(is.list(sheets), !is.null(names(sheets)))
  
  sheets_xlsx <- lapply(sheets, function(x) {
    x_df <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
    if (!identical(rownames(x_df), as.character(seq_len(nrow(x_df))))) {
      x_df <- cbind(data.frame(" " = rownames(x_df), check.names = FALSE), x_df)
    }
    x_df
  })
  write_xlsx(sheets_xlsx, paste0(file, ".xlsx"))
  
  if (is.null(png)) png <- one_table && nrow(tab) <= 40
  if (png) {
    g <- tableGrob(sheets[[1]], theme = theme_table(base_size))
    png_height <- 0.6 + 0.6 * nrow(sheets[[1]])          # cm, grows with rows
    ggsave(paste0(file, ".png"), plot = g, width = plot_width,
           height = png_height, units = "cm", dpi = plot_dpi, bg = "white")
  }
  invisible(paste0(file, ".xlsx"))
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----
