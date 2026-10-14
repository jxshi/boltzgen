#!/usr/bin/env python3
"""Early humanness/developability triage for MARCO nanobody candidates.

This script is intentionally dependency-free. It provides reproducible first-pass
sequence checks for the 50 Boltz API-designed MARCO nanobody candidates. The
humanness score is a heuristic proxy only; clinical developability should be
followed up with antibody-specific tools such as ANARCI/AbNatiV and BioPhi/OASis.
"""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path


AA = set("ACDEFGHIKLMNPQRSTVWY")
HYDROPHOBIC = set("AILMFWVY")
AROMATIC = set("FWY")
POSITIVE = set("KRH")
NEGATIVE = set("DE")

PKA = {
    "Cterm": 3.55,
    "Nterm": 7.50,
    "C": 8.50,
    "D": 3.90,
    "E": 4.07,
    "H": 6.04,
    "K": 10.54,
    "R": 12.48,
    "Y": 10.46,
}

HUMAN_LIKE_MOTIFS = [
    "VESGGG",
    "TESGGG",
    "TVSGGG",
    "VTLT",
    "SVTL",
    "WG.GT",
    "QVTV",
    "TVV",
    "TVS",
]


def fraction(seq: str, residues: set[str]) -> float:
    return sum(1 for aa in seq if aa in residues) / len(seq)


def count_regex(seq: str, pattern: str) -> int:
    return len(re.findall(pattern, seq))


def net_charge(seq: str, ph: float = 7.4) -> float:
    charge = 1 / (1 + 10 ** (ph - PKA["Nterm"]))
    charge -= 1 / (1 + 10 ** (PKA["Cterm"] - ph))
    for aa in seq:
        if aa in "KRH":
            charge += 1 / (1 + 10 ** (ph - PKA[aa]))
        elif aa in "CDEY":
            charge -= 1 / (1 + 10 ** (PKA[aa] - ph))
    return charge


def estimate_pi(seq: str) -> float:
    lo, hi = 2.0, 12.5
    for _ in range(50):
        mid = (lo + hi) / 2
        if net_charge(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def motif_score(seq: str) -> int:
    score = 0
    for motif in HUMAN_LIKE_MOTIFS:
        if re.search(motif, seq):
            score += 1
    return score


def humanness_proxy(seq: str) -> tuple[float, str]:
    """Return heuristic humanness proxy and evidence note.

    This is not an immunogenicity predictor. It rewards broad antibody-like
    framework motifs and penalizes unusual sequence liabilities.
    """
    motifs = motif_score(seq)
    length = len(seq)
    hyd = fraction(seq, HYDROPHOBIC)
    cys = seq.count("C")
    nonstandard = sum(1 for aa in seq if aa not in AA)

    score = 55 + motifs * 4
    if 110 <= length <= 140:
        score += 10
    else:
        score -= 8
    if hyd <= 0.42:
        score += 8
    elif hyd <= 0.46:
        score += 2
    else:
        score -= 8
    if cys == 0:
        score += 5
    elif cys in (2, 4):
        score += 2
    else:
        score -= 8
    score -= nonstandard * 10
    score = max(0, min(100, score))

    if score >= 78:
        label = "higher proxy humanness"
    elif score >= 65:
        label = "moderate proxy humanness"
    else:
        label = "lower proxy humanness"
    return score, label


def developability(seq: str) -> tuple[float, list[str], dict[str, float]]:
    length = len(seq)
    hyd = fraction(seq, HYDROPHOBIC)
    aromatic = fraction(seq, AROMATIC)
    charge = net_charge(seq, 7.4)
    pi = estimate_pi(seq)
    nglyc = count_regex(seq, r"N[^P][ST]")
    deamidation = count_regex(seq, r"N[GST]")
    dg_isomerization = count_regex(seq, r"DG")
    methionine = seq.count("M")
    cysteine = seq.count("C")
    cdr3_proxy = seq[math.floor(length * 0.82) :]
    cdr3_proline = cdr3_proxy.count("P")

    flags: list[str] = []
    score = 100.0

    if not (110 <= length <= 140):
        flags.append("length_outside_typical_vhh_range")
        score -= 10
    if hyd > 0.46:
        flags.append("high_hydrophobic_fraction")
        score -= 18
    elif hyd > 0.42:
        flags.append("moderate_hydrophobic_fraction")
        score -= 8
    if aromatic > 0.16:
        flags.append("high_aromatic_fraction")
        score -= 8
    if abs(charge) > 12:
        flags.append("extreme_net_charge_pH7_4")
        score -= 8
    if pi < 5.0 or pi > 9.5:
        flags.append("extreme_predicted_pI")
        score -= 8
    if nglyc:
        flags.append("n_glycosylation_sequon")
        score -= 25
    if deamidation:
        flags.append("deamidation_motif")
        score -= 5 * deamidation
    if dg_isomerization:
        flags.append("dg_isomerization_motif")
        score -= 4 * dg_isomerization
    if methionine:
        flags.append("methionine_oxidation_liability")
        score -= 3 * methionine
    if cysteine not in (0, 2, 4):
        flags.append("odd_cysteine_count")
        score -= 12
    if cdr3_proline:
        flags.append("cdr3_proxy_proline")
        score -= 4 * cdr3_proline

    metrics = {
        "length": length,
        "hydrophobic_fraction": hyd,
        "aromatic_fraction": aromatic,
        "net_charge_pH7_4": charge,
        "predicted_pI": pi,
        "n_glyc_sequon_count": nglyc,
        "deamidation_motif_count": deamidation,
        "dg_isomerization_motif_count": dg_isomerization,
        "methionine_count": methionine,
        "cysteine_count": cysteine,
        "cdr3_proxy_proline_count": cdr3_proline,
    }
    return max(0, min(100, score)), flags, metrics


def risk_label(score: float) -> str:
    if score >= 85:
        return "low"
    if score >= 70:
        return "moderate"
    if score >= 55:
        return "elevated"
    return "high"


def humanization_priority(humanness_score: float, dev_score: float, validation_tier: str) -> str:
    if "selective" in validation_tier and (humanness_score < 70 or dev_score < 75):
        return "humanize_after_binding_confirmed_high_priority"
    if "selective" in validation_tier:
        return "humanize_after_binding_confirmed_standard"
    if humanness_score < 65 or dev_score < 65:
        return "do_not_humanize_until_binding_confirmed"
    return "optional_if_selected_as_backup"


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    default_input = root.parents[0] / "marco_top50_validation_counter_screen_2026-07-21" / "candidate_selectivity_summary.csv"
    input_csv = default_input if default_input.exists() else root / "inputs" / "candidate_selectivity_summary.csv"
    output_dir = root / "outputs"
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = list(csv.DictReader(input_csv.open()))
    out_rows = []
    priority_rows = []

    for idx, row in enumerate(rows, start=1):
        seq = row["binder_sequence"].strip().upper()
        dev_score, flags, metrics = developability(seq)
        human_score, human_label = humanness_proxy(seq)
        overall = round(0.55 * dev_score + 0.45 * human_score, 2)
        priority = humanization_priority(human_score, dev_score, row["validation_tier"])
        out = {
            "panel_rank": idx,
            "candidate_id": row["candidate_id"],
            "validation_tier": row["validation_tier"],
            "design_source": row["design_source"],
            "variant": row["variant"],
            "epitope_set": row["epitope_set"],
            "binder_sequence": seq,
            "marco_complex_iplddt": row["marco_complex_iplddt"],
            "marco_iptm": row["marco_iptm"],
            "iplddt_margin_vs_best_offtarget": row["iplddt_margin_vs_best_offtarget"],
            "iptm_margin_vs_best_offtarget": row["iptm_margin_vs_best_offtarget"],
            "developability_score": round(dev_score, 2),
            "developability_risk": risk_label(dev_score),
            "humanness_proxy_score": round(human_score, 2),
            "humanness_proxy_label": human_label,
            "overall_sequence_triage_score": overall,
            "humanization_priority": priority,
            "liability_flags": ";".join(flags) if flags else "none",
        }
        for key, value in metrics.items():
            out[key] = round(value, 4) if isinstance(value, float) else value
        out_rows.append(out)
        priority_rows.append(
            {
                "priority_rank": idx,
                "candidate_id": row["candidate_id"],
                "validation_tier": row["validation_tier"],
                "humanization_priority": priority,
                "overall_sequence_triage_score": overall,
                "developability_score": round(dev_score, 2),
                "humanness_proxy_score": round(human_score, 2),
                "recommended_next_step": (
                    "Keep original sequence for first binding screen; humanize only confirmed binders."
                    if idx <= 50
                    else "Backup only."
                ),
                "reason": (
                    f"{human_label}; developability risk {risk_label(dev_score)}; "
                    f"flags: {';'.join(flags) if flags else 'none'}."
                ),
            }
        )

    output_csv = output_dir / "marco_top50_humanness_developability.csv"
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        writer.writerows(out_rows)

    priority_csv = output_dir / "humanization_priority_list.csv"
    with priority_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(priority_rows[0].keys()))
        writer.writeheader()
        writer.writerows(priority_rows)

    fasta = output_dir / "marco_top50_nanobody_sequences.fasta"
    with fasta.open("w") as handle:
        for out in out_rows:
            handle.write(f">{out['candidate_id']}|rank={out['panel_rank']}|{out['validation_tier']}|{out['humanization_priority']}\n")
            seq = out["binder_sequence"]
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")

    print(f"Wrote {output_csv}")
    print(f"Wrote {priority_csv}")
    print(f"Wrote {fasta}")
    print(f"Rows: {len(out_rows)}")


if __name__ == "__main__":
    main()
