#!/bin/bash
#==============================================================================#
# Run all EWAS sequentially - SEPAGES
# 1) Main t2  2) Main t3  3) Sensitivity t2  4) Sensitivity t3
# Usage (from the project folder):
#   nohup bash script/2_Run_Models/Run_All_Models.sh > results/2_EWAS_Results/Run_All_Models.log 2>&1 &
# to check tail -f results/2_EWAS_Results/Run_All_Models.log.
# to go out of the tail ctrl + c
#==============================================================================#

models=(1 2 3 4 5 6 7 8 9 10 11 12)
failed=()

run_step() {   # $1 = label, rest = command
  label=$1; shift
  echo "[$(date '+%Y-%m-%d %H:%M')] Starting $label"
  if "$@"; then
    echo "[$(date '+%Y-%m-%d %H:%M')] $label finished successfully."
  else
    echo "[$(date '+%Y-%m-%d %H:%M')] $label FAILED - check its log."
    failed+=("$label")
  fi
}

# 1) Main t2
for m in "${models[@]}"; do
  run_step "Main t2 model $m" Rscript script/2_Run_Models/1_Run_Main_Model_T2.R "$m"
done

# 2) Main t3
for m in "${models[@]}"; do
  run_step "Main t3 model $m" Rscript script/2_Run_Models/1_Run_Main_Model_T3.R "$m"
done

# 3) Sensitivity t2 and 4) Sensitivity t3 (all 5 sensitivities x 6 exposures each)
run_step "Sensitivity t2" Rscript script/2_Run_Models/3_Run_Sens_Models_T2.R
run_step "Sensitivity t3" Rscript script/2_Run_Models/4_Run_Sens_Models_T3.R

if [ ${#failed[@]} -eq 0 ]; then
  echo "All steps completed successfully."
else
  echo "Failed steps: ${failed[*]}"
  exit 1
fi