library(dagitty)

dag_pest <- dagitty('dag {
  pesticide [exposure]
  placental_dnam [outcome]

  ancestry -> diet
  ancestry -> placental_dnam
  maternal_age -> diet
  maternal_age -> parity
  maternal_age -> bmi
  maternal_age -> placental_dnam
  education -> diet
  education -> smoking
  education -> bmi
  education -> placental_dnam
  diet -> pesticide
  diet -> placental_dnam
  parity -> placental_dnam
  bmi -> placental_dnam
  smoking -> placental_dnam

  pesticide -> placental_dnam
  pesticide -> preg_compl
  pesticide -> gest_age
  pesticide -> cell_comp
  preg_compl -> gest_age
  preg_compl -> placental_dnam
  gest_age -> cell_comp
  gest_age -> placental_dnam
  cell_comp -> placental_dnam

  child_sex -> placental_dnam
  spec_gravity -> pesticide
  dna_contam -> placental_dnam
  hospital -> placental_dnam
}')

# minimal sufficient adjustment sets
adjustmentSets(dag_pest, exposure = "pesticide", outcome = "placental_dnam")

# what the DAG implies, to test against the data
impliedConditionalIndependencies(dag_pest)

# check the adjustment set actually used
isAdjustmentSet(dag_pest, c("ancestry", "maternal_age", "education", "diet"),
                exposure = "pesticide", outcome = "placental_dnam")