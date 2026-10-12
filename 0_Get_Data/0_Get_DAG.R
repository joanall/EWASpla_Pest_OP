library(dagitty)
library(ggdag)
library(ggplot2)

# ── Define DAG ────────────────────────────────────────────────────────────────
dag <- dagitty('dag {

  # --- Nodes ---
  parental_age       [exposure, pos="0,0"]
  dna_methylation    [outcome,  pos="4,0"]

  # Mother-related confounders
  maternal_age       [pos="-3, 3"]
  education          [pos="-1, 3"]
  smoking            [pos=" 1, 3"]
  bmi                [pos=" 3, 3"]
  parity             [pos=" 5, 3"]

  # Child-related confounders
  child_sex          [pos="-2, 2"]
  child_ethnicity    [pos=" 0, 2"]
  gestational_age    [pos=" 2, 2"]
  preg_complications [pos=" 4, 2"]

  # Study design (BiSC only)
  study_center       [pos=" 6, 2"]

  # Technical covariates (adjusting outcome measurement, not causal)
  cell_types         [pos=" 2,-2"]
  dna_contamination  [pos=" 5,-2"]

  # --- Edges ---

  # Exposure -> Outcome
  parental_age -> dna_methylation

  # Mother-related -> Exposure & Outcome
  maternal_age       -> parental_age
  maternal_age       -> dna_methylation
  education          -> parental_age
  education          -> dna_methylation
  smoking            -> parental_age
  smoking            -> dna_methylation
  bmi                -> parental_age
  bmi                -> dna_methylation
  parity             -> parental_age
  parity             -> dna_methylation

  # Child-related -> Outcome only
  child_sex          -> dna_methylation
  child_ethnicity    -> dna_methylation
  gestational_age    -> dna_methylation
  preg_complications -> dna_methylation

  # Study center -> Outcome (BiSC only)
  study_center       -> dna_methylation

  # Technical covariates -> Outcome (measurement adjustment)
  cell_types         -> dna_methylation
  dna_contamination  -> dna_methylation
}')

# ── Tidy DAG for ggdag ────────────────────────────────────────────────────────
tdag <- tidy_dagitty(dag)

# ── Node labels (clean names for plot) ───────────────────────────────────────
label_map <- c(
  parental_age       = "Parental age",
  dna_methylation    = "DNA methylation",
  maternal_age       = "Maternal age",
  education          = "Education",
  smoking            = "Smoking",
  bmi                = "Maternal BMI",
  parity             = "Parity",
  child_sex          = "Child sex",
  child_ethnicity    = "Child ethnicity",
  gestational_age    = "Gestational age",
  preg_complications = "Preg. complications",
  study_center       = "Study center\n(BiSC only)",
  cell_types         = "Cell type\nproportions",
  dna_contamination  = "DNA\ncontamination"
)

# ── Node role categories (for colour) ────────────────────────────────────────
role_map <- c(
  parental_age       = "Exposure",
  dna_methylation    = "Outcome",
  maternal_age       = "Confounder (mother)",
  education          = "Confounder (mother)",
  smoking            = "Confounder (mother)",
  bmi                = "Confounder (mother)",
  parity             = "Confounder (mother)",
  child_sex          = "Confounder (child)",
  child_ethnicity    = "Confounder (child)",
  gestational_age    = "Confounder (child)",
  preg_complications = "Confounder (child)",
  study_center       = "Study design",
  cell_types         = "Technical covariate",
  dna_contamination  = "Technical covariate"
)

tdag$data$label <- label_map[tdag$data$name]
tdag$data$role  <- role_map[tdag$data$name]

# ── Colour palette ────────────────────────────────────────────────────────────
role_colours <- c(
  "Exposure"              = "#D85A30",
  "Outcome"               = "#1D9E75",
  "Confounder (mother)"   = "#7F77DD",
  "Confounder (child)"    = "#AFA9EC",
  "Study design"          = "#888780",
  "Technical covariate"   = "#B4B2A9"
)

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggdag(tdag, use_labels = "label", text = FALSE) +
  geom_dag_point(aes(colour = role), size = 18) +
  geom_dag_label_repel(
    aes(label = label),
    size          = 2.8,
    box.padding   = 0.3,
    label.padding = 0.2,
    label.size    = 0.2,
    fill          = "white",
    colour        = "grey20",
    segment.color = "grey60",
    max.overlaps  = 20
  ) +
  geom_dag_edges(edge_colour = "grey50", arrow_directed = grid::arrow(
    length = unit(6, "pt"), type = "closed"
  )) +
  scale_colour_manual(
    values = role_colours,
    name   = "Node role"
  ) +
  theme_dag() +
  theme(
    legend.position  = "bottom",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9),
    plot.title       = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 9,  hjust = 0.5, colour = "grey40"),
    plot.margin      = margin(10, 10, 10, 10)
  ) +
  labs(
    title    = "DAG — EWAS of parental age on placental DNA methylation",
    subtitle = "Technical covariates (cell types, DNA contamination) adjust the outcome measurement and are not causal nodes"
  )

# ── Save ──────────────────────────────────────────────────────────────────────
ggsave(
  filename = "dag_parental_age.pdf",
  plot     = p,
  width    = 14,
  height   = 9,
  device   = cairo_pdf
)

ggsave(
  filename = "dag_parental_age.png",
  plot     = p,
  width    = 14,
  height   = 9,
  dpi      = 300
)

message("DAG saved as dag_parental_age.pdf and dag_parental_age.png")

# ── Optional: check adjustment sets ──────────────────────────────────────────
cat("\nMinimal adjustment sets for estimating parental_age -> dna_methylation:\n")
print(adjustmentSets(dag, exposure = "parental_age", outcome = "dna_methylation"))