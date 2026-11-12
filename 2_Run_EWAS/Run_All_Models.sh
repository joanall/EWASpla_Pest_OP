#!/bin/bash
#==============================================================================#
# Run all EWAS sequentially for one cohort: every time point, all analyses
# (main + s1-s6), with one run date for the whole batch.
# Usage (from the project root, same command on both PCs):
#   nohup bash script/2_Run_EWAS/Run_All_Models.sh sepages 2tp > results/sepages/2_Run_EWAS/Run_All_Models.log 2>&1 &
#   nohup bash script/2_Run_EWAS/Run_All_Models.sh bisc tp     > results/bisc/2_Run_EWAS/Run_All_Models.log 2>&1 &
# Arguments: 1) cohort (sepages / bisc)  2) sample (2tp / tp; default 2tp)
#            3) time points (default "t2 t3 avg")
# To follow progress: tail -f results/<cohort>/2_Run_EWAS/Run_All_Models.log
# To exit tail: Ctrl + C (the runs keep going)
#==============================================================================#

cohort=${1:?"cohort missing: sepages or bisc"}
sample_set=${2:-2tp}
time_points=(${3:-t2 t3 avg})
run_date=$(date '+%Y%m%d')   # one date for the whole batch
failed=()

mkdir -p "results/${cohort}/2_Run_EWAS"
echo "Cohort: $cohort | Sample: $sample_set | Time points: ${time_points[*]} | Run date: $run_date"

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

for tp in "${time_points[@]}"; do
  run_step "$cohort $sample_set $tp (all analyses)" \
    Rscript script/2_Run_EWAS/Run_EWAS.R "$cohort" all "$tp" "$sample_set" "$run_date"
done

if [ ${#failed[@]} -eq 0 ]; then
  echo "All steps completed successfully."
else
  echo "Failed steps: ${failed[*]}"
  exit 1
fi
