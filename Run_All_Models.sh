#!/bin/bash
#==============================================================================#
# Run all EWAS sequentially - SEPAGES
# 1) Main t2  2) Main t3  3) Sensitivity t2  4) Sensitivity t3
# Usage (from the project folder):
#   nohup bash script/2_Run_Models/Run_All_Models.sh > results/2_EWAS_Results/Run_All_Models.log 2>&1 &
# To follow progress: tail -f results/2_EWAS_Results/Run_All_Models.log
# To exit tail: Ctrl + C (the runs keep going)
#==============================================================================#

models=(1 2 3 4 5 6)
time_points=(t2 t3)
run_date=$(date '+%Y%m%d')   # one date for the whole batch
failed=()

echo "Run date: $run_date"

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

# 1) and 2) Main models, t2 then t3
for tp in "${time_points[@]}"; do
  for m in "${models[@]}"; do
    run_step "Main $tp model $m" Rscript script/2_Run_Models/1_Run_Main_Model.R "$m" "$tp" "$run_date"
  done
done

# 3) and 4) Sensitivity models, t2 then t3 (5 sensitivities x 6 exposures each)
for tp in "${time_points[@]}"; do
  run_step "Sensitivity $tp" Rscript script/2_Run_Models/2_Run_Sens_Models.R all "$tp" "$run_date"
done

if [ ${#failed[@]} -eq 0 ]; then
  echo "All steps completed successfully."
else
  echo "Failed steps: ${failed[*]}"
  exit 1
fi