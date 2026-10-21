# Focused / focused-midhelix comparison

Run from the repository's `glp1r_peptide` directory:

```bash
bash run_focused_helix_comparison.sh
```

The control is the existing `glp1r_focused.yaml`. The experimental file is
`glp1r_focused_midhelix.yaml`. Both design linear 28–32 aa peptides using
Q234/W306/R380 target hints (label indices 206/278/352) and the same reference CIF.
Only the experimental arm adds `secondary_structure: {helix: 8..24}` to protein B.
Positions are 1-based and inclusive. Positions 1–7 and 25 through the sampled
chain end remain unspecified, including the N terminus. This interval is a
conservative exploratory choice, not an experimentally established optimum.

Each arm generates 25 designs, uses both visible GPUs, and exports a budget of
25. Run them sequentially. Keep all metrics, structures, generated pipeline
configs, and software/model versions. Budget is not a guarantee of passing
filters: inspect `pass_filters` explicitly. Existing output paths are rejected.

The secondary-structure field was checked against this repository's parser.
Local YAML/shell syntax checks do not replace server-side `boltzgen check`;
the script checks both specifications before starting generation. No GPU design
run was performed while preparing these files.

Compare the whole generated batches, not just the top ten. Report the fraction
with native-like direction/N-terminal placement, unusually short interchain
heavy-atom contacts, key-site coverage, and sequence diversity alongside ipTM/PAE.
The helix hint does not enforce N/C direction or guarantee improved activity.
The previous top-ten comparison (4/10 versus 2/10) is exploratory and cannot
establish an improvement in experimental hit rate. Neither file supplies a
native or pilot ligand sequence/backbone template.
