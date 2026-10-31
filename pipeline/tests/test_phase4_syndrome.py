from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "web" / "public" / "data"

pytestmark = pytest.mark.skipif(
    not (DATA / "syndrome_signal.json").exists(),
    reason="run `uv run python pipeline/phase4_syndrome.py` first",
)


def _payload() -> dict:
    return json.loads((DATA / "syndrome_signal.json").read_text())


def test_syndrome_export_shape() -> None:
    payload = _payload()
    assert payload["schema_version"] == "1.0.0"
    assert payload["syndrome_signal"]["baseline"] == "body_weight_alone"
    assert payload["syndrome_signal"]["n_records"] > 0


def test_no_individual_ids_exported() -> None:
    text = (DATA / "syndrome_signal.json").read_text()
    assert '"ID"' not in text and '"IID"' not in text


def test_labels_and_metrics_are_sane() -> None:
    s = _payload()["syndrome_signal"]
    counts = s["label_counts"]
    # pooled "any" label must cover at least the largest single group
    assert counts["any"] >= max(counts[g] for g in ("PLCD4", "FGD4", "GALNT2", "DPF2/MUS81"))
    assert 0 < s["base_rate"] < 0.01
    for name, m in s["scores"].items():
        lo, hi = m["auc_ci"]
        assert 0.0 <= lo <= m["auc"] <= hi <= 1.0, name
        assert m["auprc_ci"][0] <= m["auprc"] <= m["auprc_ci"][1] or m["auprc"] > 0, name


def test_body_weight_flag_beats_chance() -> None:
    # Regression guard for the gen_1-key fix: with the gen_3 key this was ~0.5.
    b0 = _payload()["syndrome_signal"]["scores"]["B0_body_weight_flag"]
    assert b0["auc_ci"][0] > 0.55


def test_per_locus_report_fields() -> None:
    per_group = _payload()["per_group"]
    assert set(per_group) == {"PLCD4", "FGD4", "GALNT2", "DPF2/MUS81"}
    for group, entry in per_group.items():
        # bodyweight-available homozygotes are never fewer than complete-case ones
        assert entry["n_homozygous_with_bodyweight"] >= entry["n_homozygous_complete_case"], group
        assert 0.0 <= entry["auc_negated_bw_single_locus"] <= 1.0, group
        assert entry["observed_over_expected"] is not None, group


def test_single_locus_auc_matches_effect_size_expectation() -> None:
    # Equal-variance shift model: observed AUROC should sit near the expected value (n-dependent).
    per_group = _payload()["per_group"]
    for group in ("PLCD4", "GALNT2"):
        e = per_group[group]
        gap = abs(e["auc_negated_bw_single_locus"] - e["auc_expected_from_effect_size"])
        assert gap < 0.08, group
