#!/usr/bin/env bash
set -euo pipefail

# Run from this directory after activating the user's configured BoltzGen environment.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# A check creates a viewable CIF. Inspect that chain R and the highlighted residues
# correspond to Q234/W306/R380 before spending GPU time.
boltzgen check glp1r_peptide_28_32.yaml

# Small technical/geometry smoke test: 50 intermediates, 10 diverse finalists.
boltzgen run glp1r_peptide_28_32.yaml \
  --output runs/glp1r_28_32_smoke \
  --protocol peptide-anything \
  --num_designs 50 \
  --budget 10 \
  --devices 1

# After inspecting the smoke-test orientation, launch a larger central-length run.
# Uncomment and adjust --devices for the local server.
# boltzgen run glp1r_peptide_28_32.yaml \
#   --output runs/glp1r_28_32_prod \
#   --protocol peptide-anything \
#   --num_designs 10000 \
#   --budget 100 \
#   --devices 1

# Optional wider length exploration; do not run merely to increase candidate count.
# boltzgen run glp1r_peptide_25_35.yaml \
#   --output runs/glp1r_25_35_prod \
#   --protocol peptide-anything \
#   --num_designs 5000 \
#   --budget 50 \
#   --devices 1

# Resume an interrupted run by repeating its original command with --reuse.
# Example:
# boltzgen run glp1r_peptide_28_32.yaml \
#   --output runs/glp1r_28_32_prod \
#   --protocol peptide-anything \
#   --num_designs 10000 \
#   --budget 100 \
#   --devices 1 \
#   --reuse

