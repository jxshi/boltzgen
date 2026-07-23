import importlib.util
from pathlib import Path

import pandas as pd
import pytest


SCRIPT = Path(__file__).parents[1] / "marco_boltzgen_design" / "scripts" / "rank_specificity.py"
SPEC = importlib.util.spec_from_file_location("rank_specificity", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def prediction(design, target, role, iptm, pae=6.0):
    return {
        "design_id": design,
        "target": target,
        "target_role": role,
        "iptm": iptm,
        "interface_pae": pae,
        "buried_surface_area": 800,
        "hbonds": 6,
        "salt_bridges": 2,
        "shape_complementarity": 0.7,
    }


def test_specific_binder_ranks_above_promiscuous_binder():
    rows = [
        prediction("specific", "human_MARCO", "positive", 0.75),
        prediction("specific", "mouse_MARCO", "positive", 0.70),
        prediction("specific", "MSR1", "offtarget", 0.15, 20),
        prediction("promiscuous", "human_MARCO", "positive", 0.80),
        prediction("promiscuous", "mouse_MARCO", "positive", 0.76),
        prediction("promiscuous", "MSR1", "offtarget", 0.72),
    ]

    result = MODULE.summarize_specificity(pd.DataFrame(rows))

    assert result.iloc[0]["design_id"] == "specific"
    assert result.iloc[0]["offtarget_hit_count"] == 0
    assert result.iloc[1]["offtarget_hit_count"] == 1
    assert result.iloc[0]["delta_specificity"] > result.iloc[1]["delta_specificity"]


def test_model_predictions_are_averaged_before_margin():
    rows = [
        prediction("vhh1", "human_MARCO", "positive", 0.8),
        prediction("vhh1", "human_MARCO", "positive", 0.6),
        prediction("vhh1", "MSR1", "offtarget", 0.2),
    ]

    result = MODULE.summarize_specificity(pd.DataFrame(rows))

    assert result.iloc[0]["minimum_model_count"] == 1
    assert result.iloc[0]["positive_target_count"] == 1


def test_missing_required_columns_are_reported():
    with pytest.raises(ValueError, match="interface_pae"):
        MODULE.score_predictions(pd.DataFrame({"design_id": ["x"]}))


def test_missing_offtarget_does_not_look_specific():
    result = MODULE.summarize_specificity(
        pd.DataFrame([prediction("vhh1", "human_MARCO", "positive", 0.8)])
    )

    assert not result.iloc[0]["panel_complete"]
    assert not result.iloc[0]["specificity_pass"]
    assert pd.isna(result.iloc[0]["delta_specificity"])


def test_optional_metrics_are_reweighted_instead_of_scored_as_zero():
    complete = pd.DataFrame([prediction("vhh1", "human_MARCO", "positive", 0.7)])
    minimal = complete.drop(
        columns=["buried_surface_area", "hbonds", "salt_bridges", "shape_complementarity"]
    )

    complete_score = MODULE.score_predictions(complete).iloc[0]
    minimal_score = MODULE.score_predictions(minimal).iloc[0]

    assert complete_score["metric_coverage"] == pytest.approx(1.0)
    assert minimal_score["metric_coverage"] == pytest.approx(0.55)
    assert minimal_score["interface_score"] > 0


def test_boltz_api_column_and_role_aliases_are_accepted():
    api_rows = pd.DataFrame(
        {
            "candidate_id": ["api_vhh", "api_vhh"],
            "target_id": ["MARCO_SRCR", "MSR1_SRCR"],
            "target_role": ["intended_target", "off_target"],
            "iptm": [0.8, 0.2],
            "interface_pae_mean": [5.0, 20.0],
        }
    )

    result = MODULE.summarize_specificity(api_rows)

    assert result.iloc[0]["design_id"] == "api_vhh"
    assert result.iloc[0]["best_offtarget"] == "MSR1_SRCR"
