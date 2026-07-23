#!/usr/bin/env python3
"""Rank VHHs by their predicted MARCO-versus-off-target specificity margin.

The input is a long-form CSV with one row per design/target/model prediction.
Scores are first calculated per prediction, averaged across models for each
design/target pair, and then reduced to a conservative positive-target score
and the strongest predicted off-target score.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = {"design_id", "target", "target_role", "iptm", "interface_pae"}
COLUMN_ALIASES = {
    "candidate_id": "design_id",
    "target_id": "target",
    "interface_pae_mean": "interface_pae",
    "delta_sasa_refolded": "buried_surface_area",
    "plip_hbonds_refolded": "hbonds",
    "plip_saltbridges_refolded": "salt_bridges",
}
ROLE_ALIASES = {
    "positive": "positive",
    "intended_target": "positive",
    "target": "positive",
    "offtarget": "offtarget",
    "off_target": "offtarget",
    "negative": "offtarget",
}
OPTIONAL_METRICS = {
    "buried_surface_area": (0.15, 0.0, 1200.0),
    "hbonds": (0.10, 0.0, 10.0),
    "salt_bridges": (0.10, 0.0, 5.0),
    "shape_complementarity": (0.10, 0.0, 1.0),
}


@dataclass(frozen=True)
class TriageThresholds:
    """Starting-point gates; campaigns should calibrate these on controls."""

    min_positive_iptm: float = 0.50
    max_positive_pae: float = 12.0
    max_offtarget_iptm: float = 0.30
    min_specificity_margin: float = 0.15
    min_offtargets: int = 1


def normalize(series: pd.Series, low: float, high: float) -> pd.Series:
    """Linearly normalize a metric and clamp it to the unit interval."""
    return ((pd.to_numeric(series, errors="coerce").fillna(low) - low) / (high - low)).clip(0, 1)


def score_predictions(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate and add a unit-scaled interface score to prediction rows."""
    scored = frame.copy()
    for source, destination in COLUMN_ALIASES.items():
        if destination not in scored and source in scored:
            scored[destination] = scored[source]

    missing = REQUIRED_COLUMNS - set(scored.columns)
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")

    raw_roles = scored["target_role"].astype("string").str.strip().str.lower()
    invalid_roles = set(raw_roles.dropna()) - set(ROLE_ALIASES)
    if invalid_roles:
        raise ValueError(f"Invalid target_role values: {', '.join(sorted(invalid_roles))}")
    scored["target_role"] = raw_roles.map(ROLE_ALIASES)

    for column in OPTIONAL_METRICS:
        if column not in scored:
            scored[column] = pd.NA

    # These fixed, interpretable ranges make scores comparable between batches.
    # PAE is inverted because a smaller interface uncertainty is preferable.
    components = pd.DataFrame(index=scored.index)
    components["iptm"] = 0.35 * normalize(scored["iptm"], 0.0, 1.0)
    components["interface_pae"] = 0.20 * (1.0 - normalize(scored["interface_pae"], 0.0, 30.0))
    available_weight = pd.Series(0.55, index=scored.index)
    for column, (weight, low, high) in OPTIONAL_METRICS.items():
        present = pd.to_numeric(scored[column], errors="coerce").notna()
        components[column] = weight * normalize(scored[column], low, high).where(present, 0.0)
        available_weight += weight * present.astype(float)

    # Do not punish a row merely because an optional program/metric was absent.
    # Reweight the observed components to 1.0 and report their coverage instead.
    scored["metric_coverage"] = available_weight
    scored["interface_score"] = components.sum(axis=1) / available_weight
    return scored


def summarize_specificity(
    predictions: pd.DataFrame,
    *,
    positive_reducer: str = "min",
    offtarget_iptm_threshold: float = 0.30,
    thresholds: TriageThresholds | None = None,
) -> pd.DataFrame:
    """Return one specificity summary row per design.

    ``positive_reducer='min'`` requires every positive target (for example,
    human and mouse MARCO) to score well. Use ``max`` for a species-selective
    campaign where binding either declared positive is sufficient.
    """
    if positive_reducer not in {"min", "max", "mean"}:
        raise ValueError("positive_reducer must be one of: min, max, mean")

    thresholds = thresholds or TriageThresholds(max_offtarget_iptm=offtarget_iptm_threshold)
    if not 0 <= thresholds.min_positive_iptm <= 1:
        raise ValueError("min_positive_iptm must be between 0 and 1")
    if not 0 <= thresholds.max_offtarget_iptm <= 1:
        raise ValueError("max_offtarget_iptm must be between 0 and 1")
    if thresholds.max_positive_pae < 0 or thresholds.min_specificity_margin < 0:
        raise ValueError("PAE and specificity margin thresholds must be non-negative")
    if thresholds.min_offtargets < 1:
        raise ValueError("min_offtargets must be at least 1")
    scored = score_predictions(predictions)
    per_target = (
        scored.groupby(["design_id", "target", "target_role"], as_index=False)
        .agg(
            target_score=("interface_score", "mean"),
            target_iptm=("iptm", "mean"),
            target_pae=("interface_pae", "mean"),
            metric_coverage=("metric_coverage", "mean"),
            model_count=("interface_score", "size"),
        )
    )
    rows: list[dict[str, object]] = []
    for design_id, group in per_target.groupby("design_id", sort=False):
        positives = group[group["target_role"] == "positive"]
        negatives = group[group["target_role"] == "offtarget"]
        if positives.empty:
            raise ValueError(f"Design {design_id!r} has no positive target prediction")

        positive_score = float(getattr(positives["target_score"], positive_reducer)())
        if negatives.empty:
            best_offtarget_score = float("nan")
            best_offtarget = ""
            hit_count = 0
        else:
            best_index = negatives["target_score"].idxmax()
            best_offtarget_score = float(negatives.loc[best_index, "target_score"])
            best_offtarget = str(negatives.loc[best_index, "target"])
            hit_count = int((negatives["target_iptm"] >= thresholds.max_offtarget_iptm).sum())

        panel_complete = len(negatives) >= thresholds.min_offtargets
        positive_quality_pass = bool(
            (positives["target_iptm"] >= thresholds.min_positive_iptm).all()
            and (positives["target_pae"] <= thresholds.max_positive_pae).all()
        )
        margin = positive_score - best_offtarget_score
        specificity_pass = bool(
            panel_complete
            and positive_quality_pass
            and hit_count == 0
            and margin >= thresholds.min_specificity_margin
        )

        rows.append(
            {
                "design_id": design_id,
                "marco_score": positive_score,
                "best_offtarget_score": best_offtarget_score,
                "best_offtarget": best_offtarget,
                "delta_specificity": margin,
                "offtarget_hit_count": hit_count,
                "positive_quality_pass": positive_quality_pass,
                "panel_complete": panel_complete,
                "specificity_pass": specificity_pass,
                "positive_target_count": int(len(positives)),
                "offtarget_target_count": int(len(negatives)),
                "minimum_model_count": int(group["model_count"].min()),
                "minimum_metric_coverage": float(group["metric_coverage"].min()),
            }
        )

    return pd.DataFrame(rows).sort_values(
        ["specificity_pass", "panel_complete", "delta_specificity", "marco_score", "offtarget_hit_count"],
        ascending=[False, False, False, False, True],
    ).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", required=True, help="Long-form positive/off-target prediction CSV")
    parser.add_argument("--out", default="results/specificity_ranked.csv", help="Output summary CSV")
    parser.add_argument(
        "--positive-reducer",
        choices=("min", "max", "mean"),
        default="min",
        help="How to combine positive targets; min is conservative for human/mouse cross-reactivity",
    )
    parser.add_argument("--min-positive-iptm", type=float, default=0.50)
    parser.add_argument("--max-positive-pae", type=float, default=12.0)
    parser.add_argument("--min-specificity-margin", type=float, default=0.15)
    parser.add_argument(
        "--min-offtargets",
        type=int,
        default=1,
        help="Minimum evaluated off-target count; set to the full panel size in production",
    )
    parser.add_argument(
        "--offtarget-iptm-threshold",
        type=float,
        default=0.30,
        help="Mean target ipTM above which an off-target counts as a promiscuity hit",
    )
    args = parser.parse_args()

    thresholds = TriageThresholds(
        min_positive_iptm=args.min_positive_iptm,
        max_positive_pae=args.max_positive_pae,
        max_offtarget_iptm=args.offtarget_iptm_threshold,
        min_specificity_margin=args.min_specificity_margin,
        min_offtargets=args.min_offtargets,
    )
    ranked = summarize_specificity(
        pd.read_csv(args.predictions),
        positive_reducer=args.positive_reducer,
        offtarget_iptm_threshold=args.offtarget_iptm_threshold,
        thresholds=thresholds,
    )
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    ranked.to_csv(output, index=False)
    print(f"Wrote {output} ({len(ranked)} designs)")


if __name__ == "__main__":
    main()
