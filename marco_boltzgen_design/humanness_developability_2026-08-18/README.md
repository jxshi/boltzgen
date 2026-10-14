# MARCO Nanobody Humanness and Developability Triage

Date: 2026-08-18

This folder adds an early sequence-level humanness/developability check for the
50 MARCO nanobody candidates selected from the Boltz API validation and
counter-screening workflow.

## Purpose

The goal is to flag risks before wet-lab construction while preserving discovery
sensitivity. The recommendation is:

1. Construct and test the original 50 candidates first if the cost is acceptable.
2. Use this table to prioritize assay order and anticipate engineering risk.
3. Humanize only experimentally confirmed binders, then retest binding after
   humanization.

## Files

- `scripts/check_humanness_developability.py`: dependency-free Python triage script.
- `outputs/marco_top50_humanness_developability.csv`: full per-candidate metrics.
- `outputs/humanization_priority_list.csv`: compact humanization priority sheet.
- `outputs/marco_top50_nanobody_sequences.fasta`: FASTA for external tools.

## What Is Checked

The script computes:

- Length
- Hydrophobic fraction
- Aromatic fraction
- Net charge at pH 7.4
- Approximate pI
- N-glycosylation sequons
- Deamidation motifs
- DG isomerization motifs
- Methionine oxidation liability
- Cysteine count
- CDR3-proxy proline count
- Heuristic humanness proxy score
- Developability score and risk label
- Humanization priority

## Important Limitations

The `humanness_proxy_score` is a heuristic triage score, not a clinical
immunogenicity prediction. It should not be used as the final basis for lead
selection or regulatory decisions.

For confirmed binders, follow up with antibody-focused tools and databases such
as ANARCI/AbNatiV for numbering and structure-aware assessment, and BioPhi/OASis
or equivalent methods for humanness and potential immunogenicity assessment.

## Suggested Experimental Interpretation

- Keep high-confidence MARCO binders even if they need later humanization.
- Deprioritize candidates with both weak binding predictions and high sequence
  liabilities.
- For confirmed binders, design a small humanization panel:
  - original sequence
  - conservative humanized variant
  - more aggressive humanized variant
  - back-mutated rescue variant preserving binding-support residues

Humanization should be followed by repeat binding and specificity testing against
human MARCO SRCR, mouse Marco SRCR, MSR1/SCARA1 SRCR, SCARA5 SRCR, and CD163 SRCR3.

## Re-run

From this folder:

```bash
python3 scripts/check_humanness_developability.py
```

The default input is:

`../marco_top50_validation_counter_screen_2026-07-21/candidate_selectivity_summary.csv`
