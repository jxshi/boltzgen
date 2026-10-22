#!/usr/bin/env bash
set -euo pipefail

# Run from the glp1r_peptide directory. Compares the established focused
# baseline with one sequence-only ablation that disallows Pro at positions 1-6.
baseline_yaml="glp1r_focused.yaml"
test_yaml="glp1r_focused_nopro_nterm6.yaml"
baseline_out="runs/round11_focused_control"
test_out="runs/round11_focused_nopro_nterm6"

for file in "$baseline_yaml" "$test_yaml"; do
  if [[ ! -f "$file" ]]; then
    printf 'Missing design specification: %s\n' "$file" >&2
    exit 1
  fi
done

for dir in "$baseline_out" "$test_out"; do
  if [[ -e "$dir" ]]; then
    printf 'Output path already exists; choose fresh paths before running: %s\n' "$dir" >&2
    exit 1
  fi
done

# Validate both configurations before spending GPU time.
boltzgen check "$baseline_yaml"
boltzgen check "$test_yaml"

# Keep all candidates for a fair structural comparison. Run arms sequentially
# to avoid GPU-memory competition; --devices 2 requests both available GPUs.
boltzgen run "$baseline_yaml" \
  --output "$baseline_out" \
  --protocol peptide-anything \
  --num_designs 50 \
  --budget 50 \
  --devices 2

boltzgen run "$test_yaml" \
  --output "$test_out" \
  --protocol peptide-anything \
  --num_designs 50 \
  --budget 50 \
  --devices 2
