#!/bin/bash
#==============================================================================#
# Run EWAS models sequentially - SEPAGES t2
# Usage (from the project folder):
#   bash script/4_Run_Main_Model_T2.sh
# Each model writes its own log in db/results/main_t2/logs/
#==============================================================================#

# Models: 1-6 with cell types, 7-12 without cell types
models=(1 2 3 4 5 6 7 8 9 10 11 12)

failed=()
for which_model in "${models[@]}"; do
  echo "[$(date '+%Y-%m-%d %H:%M')] Starting model $which_model (t2)"
  if Rscript script/4_Run_Main_Model_T2.R "$which_model"; then
    echo "[$(date '+%Y-%m-%d %H:%M')] Model $which_model finished successfully."
  else
    echo "[$(date '+%Y-%m-%d %H:%M')] Model $which_model FAILED - check its log."
    failed+=("$which_model")
  fi
done

if [ ${#failed[@]} -eq 0 ]; then
  echo "All t2 models completed successfully."
else
  echo "Failed models: ${failed[*]}"
  exit 1
fi

# Alternative: run a single model
# Rscript script/4_Run_Main_Model_T2.R 1