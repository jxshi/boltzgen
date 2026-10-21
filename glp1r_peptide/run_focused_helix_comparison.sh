#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
# Two sequential arms, each using both visible GPUs; 50 new designs in total.
for arm in focused focused_midhelix; do
  run="runs/round4_${arm}"
  if [ -e "$run" ]; then
    echo "Output already exists: $run. Choose a fresh output directory." >&2
    exit 1
  fi
  boltzgen check "glp1r_${arm}.yaml"
done
for arm in focused focused_midhelix; do
  boltzgen run "glp1r_${arm}.yaml" \
    --output "runs/round4_${arm}" \
    --protocol peptide-anything \
    --num_designs 25 --budget 25 --devices 2
done
