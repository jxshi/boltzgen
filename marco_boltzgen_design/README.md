# MARCO Nanobody Design with BoltzGen

De novo VHH (nanobody) design against the SRCR domain of **human MARCO** (UniProt Q9UEW3) and **mouse Marco** (PDB 2OYA), powered by [BoltzGen](https://github.com/jxshi/boltzgen) and aligned with the **BoltzProt-1 Technical Report** (arXiv:2512.00000).

---

## Table of Contents

1. [Quick-Start](#1-quick-start)
2. [Pipeline Overview](#2-pipeline-overview)
3. [Step 1 — Clone & Setup](#step-1--clone--setup)
4. [Step 2 — Validate Specs](#step-2--validate-specs)
5. [Step 3 — Local Pilot](#step-3--local-pilot)
6. [Step 4 — HPC Production](#step-4--hpc-production)
7. [Step 5 — Collect & Merge Metrics](#step-5--collect--merge-metrics)
8. [Step 6 — Rank & Filter](#step-6--rank--filter)
9. [Step 7 — In-silico Specificity Screen](#step-7--in-silico-specificity-screen)
10. [Step 8 — CDR3 Novelty Check](#step-8--cdr3-novelty-check)
11. [Step 9 — AF2 Validation](#step-9--af2-validation)
12. [Interface Strategy Sets](#interface-strategy-sets)
13. [Key Scripts Reference](#key-scripts-reference)
14. [BoltzProt-1 Developability Flags](#boltzprot-1-developability-flags)
15. [SAbDab Novelty Cache](#sabdab-novelty-cache)
16. [Troubleshooting](#troubleshooting)

---

## 1. Quick-Start

```bash
# 1. Clone
git clone https://github.com/jxshi/boltzgen.git
cd boltzgen/marco_boltzgen_design

# 2. Validate
conda activate boltzgen
boltzgen check specs/crossreactive_marco_nanobody_setD_beta_pairing.yaml

# 3. Local pilot (50 designs)
NUM_DESIGNS=50 BUDGET=10 ./runs/run_nanobody_campaign.sh \
  specs/crossreactive_marco_nanobody_setD_beta_pairing.yaml runs/pilot

# 4. HPC production (60,000 designs)
NUM_DESIGNS=60000 BUDGET=150 sbatch scripts/run_hpc_campaign.sh \
  specs/crossreactive_marco_nanobody_setD_beta_pairing.yaml runs/setD_prod
```

---

## 2. Pipeline Overview

```
Stage 0       Stage 1       Stage 2       Stage 3       Stage 4       Stage 5       Stage 6
────────      ────────      ────────      ────────      ────────      ────────      ────────
Validate  →   HPC Design →  Collect &  →  Rank &    →  Specificity → Novelty   →  AF2
(specs)      (SLURM)        Merge         Filter       Screen        Check         Validation
```

| Stage | What | Where | Output |
|-------|------|-------|--------|
| **0 — Validate** | Spec check, target check | Local | PASS/FAIL |
| **1 — Design** | `boltzgen run` 5-step pipeline on SLURM | HPC (GPU) | CIF + NPZ in `runs/<name>/` |
| **2 — Collect** | Gather outputs from HPC, merge metrics | Local | `results/all_metrics.csv` |
| **3 — Rank** | `rank_designs.py` + developability filters | Local | `results/ranked_candidates.csv` |
| **4 — Specificity** | Re-predict against MARCO and a negative target panel | HPC + Local | `results/specificity_ranked.csv` |
| **5 — Novelty** | `novelty_check.py` — CDR3 edit distance ≥ 4 from SAbDab | Local | `results/novelty_checked.csv` |
| **6 — Validate** | `validate_designs.py` (AF2 backfold) | Local/HPC | `results/af_validation.csv` |

**Automatic post-processing:** Both `run_hpc_campaign.sh` and `run_nanobody_campaign.sh` automatically apply two hard-gate developability filters after generation — removing any design with an **N-glycosylation sequon** (NXS/T motif) and any design with a **proline in CDR3** (last 18% of sequence) — before ranking.

---

## Step 1 — Clone & Setup

```bash
git clone https://github.com/jxshi/boltzgen.git
cd boltzgen/marco_boltzgen_design
conda activate boltzgen
boltzgen --version
```

**First time on HPC only** — download BoltzGen models:

```bash
boltzgen run --force_download specs/mouse_marco_nanobody_setD_beta_pairing.yaml \
  --output /tmp/test_model_dl --num_designs 1 --budget 1
```

---

## Step 2 — Validate Specs

Always validate before submitting HPC jobs:

```bash
# Recommended: start with Set D (beta-pairing, highest priority)
boltzgen check specs/crossreactive_marco_nanobody_setD_beta_pairing.yaml

# Validate all priority specs
boltzgen check specs/mouse_marco_nanobody_setD_beta_pairing.yaml
boltzgen check specs/human_marco_nanobody_setD_beta_pairing.yaml
boltzgen check specs/crossreactive_marco_nanobody_setC_hybrid.yaml
boltzgen check specs/mouse_marco_nanobody_setA_so4_pocket.yaml
boltzgen check specs/human_marco_nanobody_setA_so4_pocket.yaml
boltzgen check specs/human_marco_nanobody_setB_patent_epitope.yaml
```

---

## Step 3 — Local Pilot

Run a small batch locally to verify the pipeline before submitting HPC jobs. All scripts automatically set `OPENBLAS_NUM_THREADS=1` to avoid RLIMIT_NPROC exhaustion.

> ⚠️ **Local Mac OOM risk:** Use `NUM_DESIGNS ≤ 100` and `DEVICES=1` locally. Production runs must go to HPC.

**Output:** `runs/<name>/final_ranked_designs/all_designs_metrics.csv`

---

### Priority 1 — Highest (cross-reactive binders targeting beta-sheet edge)

```bash
# ══ Set D: Beta-edge strand pairing (HIGHEST PRIORITY) ══

# Cross-reactive (recommended first run)
NUM_DESIGNS=50 BUDGET=10 ./runs/run_nanobody_campaign.sh \
  specs/crossreactive_marco_nanobody_setD_beta_pairing.yaml runs/pilot_setD_xr

# Human MARCO
NUM_DESIGNS=50 BUDGET=10 ./runs/run_nanobody_campaign.sh \
  specs/human_marco_nanobody_setD_beta_pairing.yaml runs/pilot_setD_human

# Mouse MARCO
NUM_DESIGNS=50 BUDGET=10 ./runs/run_nanobody_campaign.sh \
  specs/mouse_marco_nanobody_setD_beta_pairing.yaml runs/pilot_setD_mouse

# ══ Set C: Hybrid interface (maximum epitope breadth, cross-reactive only) ══
NUM_DESIGNS=50 BUDGET=10 ./runs/run_nanobody_campaign.sh \
  specs/crossreactive_marco_nanobody_setC_hybrid.yaml runs/pilot_setC_xr
```

---

### Priority 2 — High (species-specific interface targeting)

```bash
# ══ Set A: SO₄/pocket blocking ══
NUM_DESIGNS=50 BUDGET=10 ./runs/run_nanobody_campaign.sh \
  specs/human_marco_nanobody_setA_so4_pocket.yaml runs/pilot_setA_human

NUM_DESIGNS=50 BUDGET=10 ./runs/run_nanobody_campaign.sh \
  specs/mouse_marco_nanobody_setA_so4_pocket.yaml runs/pilot_setA_mouse

# ══ Set B: Patent antibody epitope (human MARCO only) ══
NUM_DESIGNS=50 BUDGET=10 ./runs/run_nanobody_campaign.sh \
  specs/human_marco_nanobody_setB_patent_epitope.yaml runs/pilot_setB
```

---

### Priority 3 — Medium (conserved hotspot targeting)

```bash
NUM_DESIGNS=50 BUDGET=10 ./runs/run_nanobody_campaign.sh \
  specs/crossreactive_marco_nanobody_hotspot.yaml runs/pilot_hotspot_xr

NUM_DESIGNS=50 BUDGET=10 ./runs/run_nanobody_campaign.sh \
  specs/human_marco_nanobody_hotspot.yaml runs/pilot_hotspot_human

NUM_DESIGNS=50 BUDGET=10 ./runs/run_nanobody_campaign.sh \
  specs/mouse_marco_nanobody_hotspot.yaml runs/pilot_hotspot_mouse
```

---

### Priority 4 — Exploratory (unconstrained surface, other binder types)

```bash
# ── Anywhere (nanobody, broadest exploration) ──
NUM_DESIGNS=50 BUDGET=10 ./runs/run_nanobody_campaign.sh \
  specs/human_marco_nanobody_anywhere.yaml runs/pilot_anywhere_human

NUM_DESIGNS=50 BUDGET=10 ./runs/run_nanobody_campaign.sh \
  specs/mouse_marco_nanobody_anywhere.yaml runs/pilot_anywhere_mouse

# ── Binder type variants (nanobody-anything protocol) ──
NUM_DESIGNS=50 BUDGET=10 ./runs/run_nanobody_campaign.sh \
  specs/human_marco_binder_anywhere.yaml runs/pilot_binder_human

NUM_DESIGNS=50 BUDGET=10 ./runs/run_nanobody_campaign.sh \
  specs/mouse_marco_binder_anywhere.yaml runs/pilot_binder_mouse

NUM_DESIGNS=50 BUDGET=10 ./runs/run_nanobody_campaign.sh \
  specs/human_marco_binder_hotspot.yaml runs/pilot_binder_hotspot_human

NUM_DESIGNS=50 BUDGET=10 ./runs/run_nanobody_campaign.sh \
  specs/mouse_marco_binder_hotspot.yaml runs/pilot_binder_hotspot_mouse

# ── Peptide (short linear motifs) ──
NUM_DESIGNS=50 BUDGET=10 ./runs/run_nanobody_campaign.sh \
  specs/human_marco_peptide_anywhere.yaml runs/pilot_peptide_human

NUM_DESIGNS=50 BUDGET=10 ./runs/run_nanobody_campaign.sh \
  specs/mouse_marco_peptide_anywhere.yaml runs/pilot_peptide_mouse

# ── Conserved surface (cross-reactive nanobody) ──
NUM_DESIGNS=50 BUDGET=10 ./runs/run_nanobody_campaign.sh \
  specs/crossreactive_conserved_surface.yaml runs/pilot_conserved_surface
```

---

## Step 4 — HPC Production

Choose your GPU cluster and submit accordingly. All scripts auto-set `OPENBLAS_NUM_THREADS=1` to prevent RLIMIT_NPROC exhaustion on shared nodes.

---

### Choose your script

|| GPU node | Use this script | Key settings |
|---|---|---|---|
| **4× A100 80GB** (recommended) | `scripts/run_a100hpc_campaign.sh` | `GPUS=4`, `BUDGET=200`, `DIFFUSION_BATCH_SIZE=16` |
| **2× RTX 5000 16GB** | `scripts/run_hpc_campaign.sh` | `GPUS=2`, `BUDGET=150`, `DIFFUSION_BATCH_SIZE=2` |

> ⚠️ **A100 script auto-sets `A100_MODE=1`** with optimal defaults. RTX script requires manual env var overrides for the same settings.

---

### A100 ×4 — Recommended production commands

```bash
# ── Priority 1: Set D beta-pairing (HIGHEST) ─────────────────────────────
# Cross-reactive (highest priority — targets both human + mouse MARCO)
NUM_DESIGNS=60000 sbatch scripts/run_a100hpc_campaign.sh \
  specs/crossreactive_marco_nanobody_setD_beta_pairing.yaml runs/setD_xr

# Human MARCO
NUM_DESIGNS=60000 sbatch scripts/run_a100hpc_campaign.sh \
  specs/human_marco_nanobody_setD_beta_pairing.yaml runs/setD_human

# Mouse MARCO
NUM_DESIGNS=60000 sbatch scripts/run_a100hpc_campaign.sh \
  specs/mouse_marco_nanobody_setD_beta_pairing.yaml runs/setD_mouse

# ── Priority 1: Set C hybrid (cross-reactive only) ──────────────────────
NUM_DESIGNS=60000 sbatch scripts/run_a100hpc_campaign.sh \
  specs/crossreactive_marco_nanobody_setC_hybrid.yaml runs/setC_xr

# ── Priority 2: Set A SO₄/pocket ─────────────────────────────────────────
NUM_DESIGNS=60000 sbatch scripts/run_a100hpc_campaign.sh \
  specs/human_marco_nanobody_setA_so4_pocket.yaml runs/setA_human

NUM_DESIGNS=60000 sbatch scripts/run_a100hpc_campaign.sh \
  specs/mouse_marco_nanobody_setA_so4_pocket.yaml runs/setA_mouse

# ── Priority 2: Set B patent epitope (human only) ───────────────────────
NUM_DESIGNS=60000 sbatch scripts/run_a100hpc_campaign.sh \
  specs/human_marco_nanobody_setB_patent_epitope.yaml runs/setB

# ── Priority 3: Hotspot ─────────────────────────────────────────────────
NUM_DESIGNS=60000 sbatch scripts/run_a100hpc_campaign.sh \
  specs/crossreactive_marco_nanobody_hotspot.yaml runs/hotspot_xr

NUM_DESIGNS=60000 sbatch scripts/run_a100hpc_campaign.sh \
  specs/human_marco_nanobody_hotspot.yaml runs/hotspot_human

NUM_DESIGNS=60000 sbatch scripts/run_a100hpc_campaign.sh \
  specs/mouse_marco_nanobody_hotspot.yaml runs/hotspot_mouse
```

---

### A100 ×4 — Speed mode (2–3× faster, for screening)

```bash
SPEED_MODE=1 NUM_DESIGNS=60000 sbatch scripts/run_a100hpc_campaign.sh \
  specs/crossreactive_marco_nanobody_setD_beta_pairing.yaml runs/setD_screen

SPEED_MODE=1 NUM_DESIGNS=60000 sbatch scripts/run_a100hpc_campaign.sh \
  specs/crossreactive_marco_nanobody_setC_hybrid.yaml runs/setC_screen
```

Speed mode sets: `recycling_steps=1`, `compile_pairformer=true`, `compile_structure=true`, `inverse_fold precision=bf16-mixed`, `diffusion_batch_size=32`.

---

### A100 ×4 — Run all priority specs in parallel

```bash
# Submit 7 jobs simultaneously, each on a dedicated 4×A100 node
for SPEC in crossreactive_marco_nanobody_setD_beta_pairing \
            human_marco_nanobody_setD_beta_pairing \
            mouse_marco_nanobody_setD_beta_pairing \
            crossreactive_marco_nanobody_setC_hybrid \
            human_marco_nanobody_setA_so4_pocket \
            mouse_marco_nanobody_setA_so4_pocket \
            human_marco_nanobody_setB_patent_epitope; do
  NUM_DESIGNS=60000 sbatch scripts/run_a100hpc_campaign.sh \
    specs/${SPEC}.yaml runs/$(basename $SPEC)_a100 &
done
wait
```

---

### RTX 5000 ×2 — Standard quality mode

```bash
# 60,000 designs (BoltzProt-1 production standard)
NUM_DESIGNS=60000 BUDGET=150 GPUS=2 sbatch scripts/run_hpc_campaign.sh \
  specs/crossreactive_marco_nanobody_setD_beta_pairing.yaml \
  runs/setD_rtx
```

### RTX 5000 ×2 — Speed mode

```bash
SPEED_MODE=1 NUM_DESIGNS=60000 BUDGET=150 GPUS=2 sbatch scripts/run_hpc_campaign.sh \
  specs/crossreactive_marco_nanobody_setD_beta_pairing.yaml \
  runs/setD_rtx_fast
```

---

### Key environment variables

| Variable | A100 default | RTX default | Meaning |
|----------|-------------|-------------|---------|
| `NUM_DESIGNS` | 60000 | 60000 | Designs per job |
| `BUDGET` | 200 | 150 | Inference steps per design |
| `GPUS` | 4 | 2 | GPU count |
| `DIFFUSION_BATCH_SIZE` | 16 (32 in speed) | 2 (8 in speed) | Batch size per step |
| `SPEED_MODE` | 0 | 0 | 1 = fast screening mode |
| `EXCLUDE_NGLYC` | 1 | 1 | Auto-filter N-glyc sequons |
| `FILTER_PROLINE` | 1 | 1 | Auto-filter proline-in-CDR3 |

---

### Monitor SLURM jobs

```bash
squeue -u $USER
tail -f logs/boltzgen_<spec>_<JOB_ID>.out
ls runs/<name>/intermediate_designs/*.cif 2>/dev/null | wc -l   # count done so far
```

**A100 ×4 timing estimates (BUDGET=200, diffusion_batch_size=16):**

| Mode | 60k designs | 100k designs |
|------|-------------|-------------|
| Quality (`SPEED_MODE=0`) | ~18–24 hours | ~30–40 hours |
| Speed (`SPEED_MODE=1`) | ~8–12 hours | ~14–20 hours |

---

## Step 5 — Collect & Merge Metrics

After all HPC jobs finish, copy results to local machine and merge:

```bash
# Copy from HPC to local
rsync -avz hpc:boltzgen/marco_boltzgen_design/runs/ runs/

# Merge all metrics CSVs
./scripts/collect_campaign.sh \
  --runs runs/setA_mouse runs/setA_human runs/setB runs/setC runs/setD \
  --out results/all_metrics.csv
```

---

## Step 6 — Rank & Filter

Rank by confidence, developability, and epitope coverage, with optional quality pre-filter gates:

```bash
# Standard ranking (no pre-filter gates):
python scripts/rank_designs.py \
  --metrics results/all_metrics.csv \
  --human-conserved "A:423,A:424,A:431,A:460,A:466,A:468,A:488,A:499" \
  --mouse-conserved "A:6,A:8,A:15,A:44,A:50,A:52,A:72,A:83" \
  --max-len 120 \
  --out results/ranked_candidates.csv

# With quality pre-filter gates (recommended for production):
python scripts/rank_designs.py \
  --metrics results/all_metrics.csv \
  --human-conserved "A:423,A:424,A:431,A:460,A:466,A:468,A:488,A:499" \
  --mouse-conserved "A:6,A:8,A:15,A:44,A:50,A:52,A:72,A:83" \
  --max-len 120 \
  --min_iptm 0.25 \
  --max_pae 15 \
  --max_gly_ala_frac 0.35 \
  --out results/ranked_candidates.csv

# Tune alpha for more cross-reactive bias (default: 0.5):
python scripts/rank_designs.py \
  --metrics results/all_metrics.csv \
  --alpha_crossreactivity 0.8 \
  --affinity_weight 0.3 \
  --min_iptm 0.25 \
  --out results/ranked_candidates.csv
```

**Quality pre-filter gates** (applied before scoring, to remove non-binders early):

| Flag | What it does | Recommended threshold |
|------|-------------|----------------------|
| `--min_iptm` | Remove designs below ipTM threshold | `0.25–0.30` (Boltz confirmed binders: ipTM > 0.5) |
| `--max_pae` | Remove designs above this PAE (Å) | `12–15` (close interface geometry) |
| `--max_gly_ala_frac` | Remove overly Gly/Ala-rich CDR3s | `0.30–0.35` (natural VHH: < 30%) |

What this does:
- Sorts by `final_score = base_confidence + alpha x crossreactivity_score + beta x interface_quality − developability_penalties`
- `interface_quality` is a min-max normalised composite of PLIP metrics (ipTM, H-bonds, salt bridges, buried SA) already present in BoltzGen output
- Applies 9 developability flag columns (Cys in CDR regions, length, charge, hydrophobicity, aromatic fraction, pI region, Pro in CDR3, N-glyc motifs, Gly/Ala-rich CDR3)
- Auto-detects contact columns in metrics CSV for cross-reactivity scoring

**Inspect top candidates:**

```python
import pandas as pd
df = pd.read_csv('results/ranked_candidates.csv')
print(df[['design_id','final_score','iptm','min_interaction_pae',
          'cdr3_gly_ala_frac','developability_flags','binder_sequence']].head(20))
```

---

## Step 7 — In-silico Specificity Screen

High MARCO confidence alone does **not** demonstrate specificity. Re-predict the
top 200–500 VHH sequences against the same positive and negative targets using
the same complex-prediction settings. A useful initial panel is human and mouse
MARCO as positives and MSR1/SCARA1, SCARA3, SCARA5, CD163, CD5L, MRC1, TREM2,
MERTK, and AXL as off-targets. Choose constructs containing the homologous,
solvent-exposed domain rather than comparing unrelated full-length constructs.

### Recommended funnel for the existing 8,000 + 3,700 designs

Do not immediately run all 11,700 sequences against every off-target. Use a
staged funnel so that the expensive counter-screen is spent on credible and
diverse sequences:

1. **Unify and deduplicate by exact amino-acid sequence.** Preserve a stable
   `design_id`, `design_source` (`boltzgen_8000` or `boltz_api_3700`), original
   rank, design target/species, epitope set, and paths to CIF/metadata. Do not
   merge records only by candidate name; identifiers can collide across runs.
2. **Apply source-neutral quality gates.** Recalculate sequence descriptors from
   the sequence for both sources, rather than trusting differently named source
   columns. As a balanced first pass use design ipTM ≥0.45, design interface PAE
   ≤12–15 Å, no N-glycosylation sequon, no extra CDR cysteine, acceptable length,
   charge and hydrophobic/aromatic patches. Retain failed rows with failure
   reasons for audit instead of deleting them.
3. **Cluster for diversity before counter-screening.** Cluster full VHHs at about
   90% identity and CDR3s at about 70–80% identity, then keep several members per
   cluster across source, epitope and human/mouse/cross-reactive strata. A useful
   first target is 1,000–2,000 quality-passing sequences, followed by 300–500
   diverse representatives. These are capacity targets, not biological cutoffs.
4. **Run a cheap homolog screen on the 1,000–2,000 pool.** Predict human MARCO,
   mouse MARCO, MSR1/SCARA1 and SCARA5 with identical settings and multiple
   seeds. Advance approximately 200–500 candidates with strong MARCO confidence,
   positive ΔipTM/ΔPAE, no obvious homolog hit, and broad cluster coverage.
5. **Run the full negative panel on the 200–500 shortlist.** Add CD163, CD5L,
   MRC1, TREM2, MERTK and AXL, use at least two independent model families where
   feasible, and keep per-seed rows rather than only the best seed. Use
   `rank_specificity.py --min-offtargets 9` for the production summary.
6. **Select a Pareto-balanced experimental panel.** Allocate the final 30–100
   constructs across high specificity margin, human/mouse cross-reactivity,
   epitope diversity, sequence-cluster diversity, both design sources and
   developability. Do not simply synthesize the top 30 values from one scalar
   score. Confirm by MARCO binding, homolog counter-binding, irrelevant-protein
   controls and MARCO-positive/negative cell assays.

The checked-in Boltz API package already contains a four-target counter-screen
for its top 50 (`MARCO`, `MSR1/SCARA1`, `SCARA5`, and `CD163`). Its
`validation_metrics.csv` can now be passed directly to `rank_specificity.py`:

```bash
python scripts/rank_specificity.py \
  --predictions boltz_api_results_package_2026-07-22/specificity_top50/tables/validation_metrics.csv \
  --min-offtargets 3 \
  --out results/api_top50_specificity.csv
```

The script recognizes both its canonical schema and the package aliases
`candidate_id`, `target_id`, `intended_target`, `off_target`, and
`interface_pae_mean`. The existing four-target result is useful evidence but is
not a complete nine-off-target panel. The 8,000 BoltzGen designs should first be
exported to the same long format after their positive/negative re-predictions;
their original design-time MARCO score is not an off-target specificity result.

Collect all predictions in a long-form CSV. Repeat a target with different
`model` values when using Boltz, AF3, or Chai-1; the script averages replicates
or model predictions for each design/target pair before calculating the margin.

```csv
design_id,target,target_role,model,iptm,interface_pae,buried_surface_area,hbonds,salt_bridges,shape_complementarity
vhh_001,human_MARCO,positive,boltz,0.72,7.1,840,7,2,0.71
vhh_001,mouse_MARCO,positive,boltz,0.68,8.0,790,6,2,0.68
vhh_001,MSR1,offtarget,boltz,0.21,18.2,310,2,0,0.42
```

Rank the resulting panel:

```bash
python scripts/rank_specificity.py \
  --predictions results/target_panel_predictions.csv \
  --positive-reducer min \
  --offtarget-iptm-threshold 0.30 \
  --min-offtargets 9 \
  --out results/specificity_ranked.csv
```

The fixed unit-scaled interface score is:

```text
0.35 × ipTM + 0.20 × inverted interface PAE + 0.15 × buried surface area
+ 0.10 × H-bonds + 0.10 × salt bridges + 0.10 × shape complementarity
```

`delta_specificity` is the conservative (lowest) positive MARCO score minus
the best off-target score. `offtarget_hit_count` is a simple promiscuity index.
For a cross-reactive campaign, keep the default `--positive-reducer min`; for a
species-selective campaign, label only the desired species as `positive` and
the other species as `offtarget`. Missing optional interface measurements are
excluded and the available weights are renormalized; inspect
`minimum_metric_coverage` rather than treating a partial score as equally
informative. `iptm` and `interface_pae` are always required. A candidate with no
off-target rows is marked as an incomplete panel and cannot pass. Do not interpret
the output as proof of specificity: prioritize a large margin and zero hits,
then confirm the shortlist experimentally with orthogonal binding and cell-panel
assays. Thresholds such as positive ipTM ≥0.5, interface PAE ≤10–12 Å, and
off-target ipTM ≤0.3 are triage heuristics that should be calibrated to controls.

### Are the defaults optimal?

No universal cutoff is optimal across predictors, target constructs, and model
versions. The defaults are deliberately moderate **starting points**, not fitted
MARCO decision boundaries:

| Parameter | Default | Recommended use |
|-----------|---------|-----------------|
| `--positive-reducer` | `min` | Keep for human/mouse cross-reactivity; use one declared positive for species selectivity. |
| `--min-positive-iptm` | `0.50` | Exploratory: 0.45; balanced: 0.50; strict: 0.55–0.60. |
| `--max-positive-pae` | `12 Å` | Tighten to 10 Å after confirming enough positive controls survive. |
| `--offtarget-iptm-threshold` | `0.30` | Lower to 0.25 for conservative polyspecificity triage; do not raise without controls. |
| `--min-specificity-margin` | `0.15` | Use 0.10 for broad discovery and 0.20 for a strict shortlist. This is a composite-score margin, not ΔipTM. |
| `--min-offtargets` | `1` | Minimal smoke-test default; set to the actual panel size (9 for the panel above) in production. |

Tune cutoffs on known MARCO binders, non-binders, irrelevant VHHs, and deliberately
challenging homologs. Apply the identical target construct, MSA/template policy,
seed count, and model ensemble to positives and negatives. Prefer selecting a
Pareto set over repeatedly changing weights to promote favored candidates. With
experimental labels, choose thresholds using held-out precision/recall or ROC
analysis; without labels, report sensitivity analyses at exploratory, balanced,
and strict settings instead of claiming that a single default is validated.

---

## Step 8 — CDR3 Novelty Check

Flags any design whose CDR3 is **edit-distance < 4** from a known SAbDab antibody. Per BoltzProt-1: *"every recovered design has a minimum CDR3 edit distance of at least four to its closest SAbDab match."*

The pre-built cache (4,466 unique CDR3s from 32k SAbDab PDBs) is committed to the repo at `.sabdab_reference.json` — no rebuild needed on most machines.

```bash
# First time: build the cache from local SAbDab zip
# (only needed if .sabdab_reference.json is missing or stale)
python scripts/novelty_check.py \
  --build_cache \
  --sabdab_zip ~/Downloads/all_structures.zip

# Check designs against the reference
python scripts/novelty_check.py \
  --designs results/ranked_candidates.csv \
  --out results/novelty_checked.csv
```

**Output columns:**

| Column | Meaning |
|--------|---------|
| `cdr3_edit_distance` | Min edit distance of CDR3 to any SAbDab entry |
| `cdrs_edit_distance` | Min edit distance of CDR1+CDR2+CDR3 combined |
| `novelty_flag` | `low_novelty` if edit distance < 4 |

**Filter to novel designs only:**

```python
import pandas as pd
df = pd.read_csv('results/novelty_checked.csv')
novel = df[df['novelty_flag'] != 'low_novelty']
novel.to_csv('results/novel_candidates.csv', index=False)
print(f"Novel candidates: {len(novel)} / {len(df)}")
```

---

## Step 9 — AF2 Validation

Validates that designed binder sequences back-fold correctly to the predicted structures:

```bash
python scripts/validate_designs.py \
  --complexes results/candidate_cifs \
  --metrics results/novel_candidates.csv \
  --top_n 50 \
  --method colabfold \
  --out results/af_validation.csv
```

**Thresholds:** CA RMSD < 2.5 Å **AND** mean PAE < 5.0 Å

**Merge AF2 results:**

```python
import pandas as pd
ranked = pd.read_csv('results/novel_candidates.csv')
af2 = pd.read_csv('results/af_validation.csv')
merged = ranked.merge(af2[['design_id','af2_rmsd','af2_pae','af2_plddt','flag_ok']], on='design_id', how='left')
passing = merged[merged['flag_ok'] == True]
print(f"AF2-passing designs: {len(passing)} / {len(merged)}")
print(passing[['design_id','plddt','final_score','af2_rmsd','af2_pae']].head(20))
```

**Designs passing AF2 are ready for experimental characterization.**

---

## Interface Strategy Sets

Four strategy groups cover distinct SRCR surfaces. **Set D** and **Set C** are the highest priority.

| Set | Strategy | Specs | Priority |
|-----|----------|-------|----------|
| **D** | Beta-edge strand targeting — polar beta-sheet face | `*_setD_beta_pairing.yaml` | 🔴 Highest |
| **C** | Hybrid interface — Sets A + B union | `*_setC_hybrid.yaml` | 🔴 Highest |
| **A** | SO₄/pocket blocking — ligand-binding crevice | `*_setA_so4_pocket.yaml` | 🟡 High |
| **B** | Patent antibody epitope — human-only | `*_setB_patent_epitope.yaml` | 🟡 High |
| Hotspot | ARG-rich basic patch (conserved) | `*_hotspot.yaml` | 🟢 Medium |
| Anywhere | Unconstrained surface exploration | `*_anywhere.yaml` | 🔵 Exploratory |

**Set D — Beta-Pairing (highest priority):**

| Species | Residues | Notes |
|---------|----------|-------|
| Mouse | 7,8,10,12,14,15,17,20,21,22,50,52,54,98,101,102 | 2OYA label_seq |
| Human | 423,424,426,428,430,431,433,436,437,438,466,468,470,514,517,518 | Q9UEW3 label_seq |

**Set C — Hybrid (cross-reactive, maximum breadth):**

| Species | Residues | Notes |
|---------|----------|-------|
| Mouse | 12,14,21,50,56,58,78,89 | 2OYA label_seq |
| Human | 429,431,438,450,452,467,472,473,475,487,495,499,505,506,507,509,511 | Q9UEW3 label_seq |

**Set A — SO₄/Pocket Blocking:**

| Species | Residues | Notes |
|---------|----------|-------|
| Mouse | 12,14,21,50,56,58,78,89 | 2OYA label_seq |
| Human | 429,431,438,467,473,475,495,506 | Q9UEW3 label_seq |

**Set B — Patent Epitope (human-only):**

| Species | Residues | Notes |
|---------|----------|-------|
| Human | 450,452,472,473,487,499,505,507,509,511 | Q9UEW3 label_seq |

> ⚠️ **Numbering:** Spec files use **mmCIF `label_seq`** (not Q9UEW3 full-sequence positions). Mouse uses 2OYA `label_seq` directly.

---

## Key Scripts Reference

| Script | What it does |
|--------|-------------|
| `scripts/run_hpc_campaign.sh` | SLURM submission — full pipeline, N-glyc + proline pre-filters, then ranking |
| `runs/run_nanobody_campaign.sh` | Local campaign runner (Mac/HPC login node) |
| `scripts/filter_developability.py` | Unified N-glyc sequon filter, proline-in-CDR3 filter, and Gly/Ala-rich CDR3 filter (last 18% heuristic); use `--filter_gly_ala --gly_ala_threshold 0.35` |
| `scripts/rank_designs.py` | Rank by `base_conf + alpha x crossreactivity + beta x interface_quality − penalties`; output `interface_quality` and `cdr3_gly_ala_frac` columns; has quality pre-filter gates `--min_iptm`, `--max_pae`, `--max_gly_ala_frac` |
| `scripts/novelty_check.py` | Check CDR3 edit distance against SAbDab reference |
| `scripts/validate_designs.py` | AF2 backfold validation |
| `scripts/collect_campaign.sh` | Merge metrics from multiple runs |

---

## BoltzProt-1 Developability Flags

The `rank_designs.py` script applies 10 sequence-based developability flags aligned with the BoltzProt-1 six-assay panel:

| Flag | Threshold | Risk |
|------|-----------|------|
| `has_cys` | Any Cys in CDR1/2/3 regions | Disulfide scrambling / oxidative aggregation (framework C's at pos ~22 and ~95 are structural and excluded) |
| `nglyc_motif` | N[^P][ST] pattern present | Glycan heterogeneity during expression |
| `too_long` | Binder length > 120 aa | High-risk for expression |
| `excess_positive_charge` | Net charge > +8 | Self-association / AC-SINS risk |
| `hydrophobic_patch_flag` | Hydrophobic fraction > 0.42 | HIC retention / aggregation |
| `aromatic_high` | Aromatic fraction > 0.14 | Polyspecificity / BVP ELISA risk |
| `pi_acidic` | Net charge < −5 | HIC retention / acidic pI risk |
| `pi_basic` | Net charge > +5 | HIC retention / basic pI risk |
| `proline_cdr3` | Pro in last 18% of sequence (~CDR3 region) | Thermal stability / Tm disruption |
| `gly_ala_rich_cdr3` | Gly+Ala fraction > 35% in CDR3 | Synthetic loop / non-native paratope (diffusion-model failure mode) |

**Penalty weights in final score:** `nglyc_motif`, `proline_cdr3`, and `gly_ala_rich_cdr3` incur a **2× penalty** (stronger weight); all others incur 1×.

**Developability tiers (from BoltzProt-1):**

| Tier | Score range | Interpretation |
|------|-------------|----------------|
| Tier-1 | 0 penalties | Best — proceed to experimental validation |
| Tier-2 | 1–2 penalties | Acceptable with characterization |
| Screening-Hit | 3–4 penalties | Needs developability screening assays |
| Problematic | 5–6 penalties | High risk — consider redesign |

---

## SAbDab Novelty Cache

The file `.sabdab_reference.json` contains **4,466 unique CDR3 sequences** (length 6–22 aa) extracted from 32,000+ IMGT-renumbered PDBs in the SAbDab archive.

- **Location:** `.sabdab_reference.json` (committed to repo, 277 KB)
- **Rebuild only if:** The SAbDab version changes or the cache is deleted
- **Rebuild command:** `python scripts/novelty_check.py --build_cache --sabdab_zip ~/Downloads/all_structures.zip`

CDR3 is extracted from IMGT positions 105–117 (inclusive, 1-indexed → Python slice `[104:117]`). Verified against 6xul (15-aa CDR3) and 7tlz (22-aa camelid VHH CDR3).

**Novelty threshold:** min edit distance ≥ 4 (BoltzProt-1 standard). This was the minimum distance observed across all BoltzProt-1 confirmed binders.

---

## Troubleshooting

**"ERROR: spec not found"**
→ Run from `boltzgen/marco_boltzgen_design/` directory.

**OOM / CUDA out of memory on HPC**
→ Reduce `NUM_DESIGNS` or set `SPEED_MODE=1` (lowers `diffusion_batch_size` footprint).

**SLURM job killed before completing**
→ Re-run with `--reuse` flag (already set in `run_hpc_campaign.sh`). Already-generated designs are skipped.

**"conda: command not found"**
→ `eval "$(conda shell.bash hook)" && conda activate boltzgen`

**`all_designs_metrics.csv` not found after job finishes**
→ The job was likely killed early. Check `logs/<spec>_<JOB_ID>.out` for the last completed design and re-run with `--reuse`.

**Novelty check very slow**
→ On first run, `novelty_check.py` parses `~/Downloads/all_structures.zip`. This takes ~30–60 seconds. Subsequent runs use `.sabdab_reference.json` cache (~1 second).

---

## Complete End-to-End Workflow

```bash
# ═══════════════════════════════════════════════════════════════
# STAGE 0: Validate
# ═══════════════════════════════════════════════════════════════
conda activate boltzgen
boltzgen check specs/crossreactive_marco_nanobody_setD_beta_pairing.yaml

# ═══════════════════════════════════════════════════════════════
# STAGE 1: HPC production (60k designs per spec)
# ═══════════════════════════════════════════════════════════════
NUM_DESIGNS=60000 BUDGET=150 sbatch scripts/run_hpc_campaign.sh \
  specs/crossreactive_marco_nanobody_setD_beta_pairing.yaml runs/setD_prod

# ═══════════════════════════════════════════════════════════════
# STAGE 2: Copy & merge (run locally after HPC jobs finish)
# ═══════════════════════════════════════════════════════════════
rsync -avz hpc:boltzgen/marco_boltzgen_design/runs/setD_prod/ runs/setD_prod/
./scripts/collect_campaign.sh --runs runs/setD_prod --out results/all_metrics.csv

# ═══════════════════════════════════════════════════════════════
# STAGE 3: Rank & filter
# ═══════════════════════════════════════════════════════════════
python scripts/rank_designs.py \
  --metrics results/all_metrics.csv \
  --human-conserved "A:423,A:424,A:431,A:460,A:466,A:468,A:488,A:499" \
  --mouse-conserved "A:6,A:8,A:15,A:44,A:50,A:52,A:72,A:83" \
  --max-len 120 \
  --out results/ranked_candidates.csv

# ═══════════════════════════════════════════════════════════════
# STAGE 4: CDR3 novelty check
# ═══════════════════════════════════════════════════════════════
python scripts/novelty_check.py \
  --designs results/ranked_candidates.csv \
  --out results/novelty_checked.csv

# Filter to novel candidates only
python3 -c "
import pandas as pd
df = pd.read_csv('results/novelty_checked.csv')
novel = df[df['novelty_flag'] != 'low_novelty']
novel.to_csv('results/novel_candidates.csv', index=False)
print(f'Novel candidates: {len(novel)} / {len(df)}')
"

# ═══════════════════════════════════════════════════════════════
# STAGE 5: AF2 backfold validation
# ═══════════════════════════════════════════════════════════════
python scripts/validate_designs.py \
  --complexes results/candidate_cifs \
  --metrics results/novel_candidates.csv \
  --top_n 50 \
  --method colabfold \
  --out results/af_validation.csv
```
