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
# Then add theme_ewas() to each ggplot, use the palettes below for colours,
# save plots with save_plot() and tables with save_table().

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
# save_table(tab1, "results/<folder>/<name>")   # writes .csv and .png


#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----

library(ggplot2)
library(khroma)
library(colorspace)
library(gridExtra)
library(tableone)

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

# Time points
pal_time <- c(t2 = "#66CCEE", t3 = "#AA3377")

# Direction of effect (EWAS: hypo / hyper / not significant)
pal_direction <- c(Hypomethylated  = "#4477AA",
                   Hypermethylated = "#EE6677",
                   NS              = "#BBBBBB")

# Exposures (one colour each, fixed across all figures)
pal_exposure <- c("2,4-D"  = "#4477AA", "3-PBA"  = "#66CCEE",
                  "TCPy"   = "#228833", "DMP"    = "#CCBB44",
                  "DEP"    = "#EE6677", "DAP"    = "#AA3377")

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

# Tables: always .csv (for the paper) + .png (for slides / quick checks)
save_table <- function(tab, file, base_size = 11) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.csv(tab, paste0(file, ".csv"), row.names = TRUE)
  g <- tableGrob(tab, theme = theme_table(base_size))
  png_height <- 0.6 + 0.6 * nrow(tab)                   # cm, grows with rows
  ggsave(paste0(file, ".png"), plot = g, width = plot_width,
         height = png_height, units = "cm", dpi = plot_dpi, bg = "white")
}

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----
