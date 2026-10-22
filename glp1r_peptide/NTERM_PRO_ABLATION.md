# GLP-1R N-terminal Pro sequence ablation

## What BoltzGen supports

The reviewed design-spec parser supports per-position `residue_constraints`
with allowed or disallowed amino-acid identities, and `secondary_structure`
annotations for designed residues. The source stores residue constraints as an
amino-acid mask. Neither the public design-spec keys nor the parser expose a
non-covalent distance, coordinate, or peptide-endpoint restraint between a
designed residue and a target residue. `structure_groups` controls structural
visibility and relative placement of entities; it is not an endpoint restraint.
The top-level `bond` constraint is covalent and must not be used to emulate a
non-covalent GLP-1R binding pose.

Sources: [BoltzGen design-spec documentation](https://github.com/HannesStark/boltzgen#how-to-make-a-design-specification-yaml), [schema parser](https://github.com/jxshi/boltzgen/blob/main/src/boltzgen/data/parse/schema.py), and [design-info validation](https://github.com/jxshi/boltzgen/blob/main/src/boltzgen/data/data.py).

The maintainer response to [issue #247](https://github.com/HannesStark/boltzgen/issues/247) confirms that `residue_constraints` are hard sequence constraints during inverse folding only. The pre-inverse-folding structures in `intermediate_designs/` will not show their effect. Use a fresh output directory, do not skip inverse folding, and inspect `intermediate_designs_inverse_folded/` or final sequences. The feature was introduced in v0.3.1; check the installed version on the server.

## Why test the Pro-free N terminus

Across the available rounds 6–10 structural review, the retrospective direction
screen marked 15 of 250 candidates as native-like. All 15 lacked Pro in the
first six sequence positions. In rounds 7–10, 12/26 Pro-free N termini passed
that geometry screen, versus 0/174 with at least one Pro in positions 1–6;
round 6 added 3/8 and 0/42, respectively. These are retrospective associations
in generated/model-filtered data, not independent experimental outcomes and
not evidence that Pro causes the reversed pose. The total counts are small and
the sequence may be correlated with other design features.

The new `glp1r_focused_nopro_nterm6.yaml` tests one narrow change: disallow Pro
at peptide positions 1–6 during inverse folding. It keeps the existing
28–32-residue length range, target chain, binding-site hints, and exclusion
regions. It does **not** force the N terminus toward H7, alter diffusion
backbone generation, or guarantee a binder or agonist. The tested arm should
therefore be interpreted as a sequence-level ablation that may or may not
change the refolded orientation.

Do not add a helix annotation based on the earlier mid-helix comparison: that
ablation did not show a reproducible improvement in the full-batch direction
screen.

## Controlled pilot

Copy the new YAML to the same `glp1r_peptide` directory that contains
`reference/6X18_R_canonical.cif`, then run:

```bash
bash run_focused_nopro_comparison.sh
```

The script checks both specifications before generation, uses fresh output
directories, runs 50 designs per arm with `budget=50`, and requests two devices.
It runs the arms sequentially to avoid competing for GPU memory. Save the
generated `config/`, `steps.yaml`, BoltzGen version, model artifacts, and full
outputs for both arms. Verify the constraint in the inverse-folded/final
sequences; the six N-terminal residues must contain no P in the experimental
arm. Do not use `--reuse` or skip inverse folding for this comparison.

Compare both complete batches with the same analysis code and criteria:
native-like terminal direction and local N-terminal geometry, built-in filter
status, interchain short contacts, key-site distances, and sequence diversity.
The direction screen is a structural prioritization heuristic, not a measured
binding or GLP-1R agonism rate. A better computational geometry fraction would
justify follow-up validation, not a claim of increased experimental hit rate.
