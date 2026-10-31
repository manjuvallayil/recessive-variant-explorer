from __future__ import annotations

import json
from pathlib import Path

import pytest

from pipeline.phase3_inbreeding import PAPER_INBREEDING_DEPRESSION_KG_PER_1PCT

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "web" / "public" / "data"

pytestmark = pytest.mark.skipif(
    not (DATA / "inbreeding.json").exists(),
    reason="run `uv run python pipeline/phase3_inbreeding.py` first",
)


def _payload() -> dict:
    return json.loads((DATA / "inbreeding.json").read_text())


def test_inbreeding_export_shape() -> None:
    payload = _payload()
    assert payload["schema_version"] == "1.0.0"
    assert payload["inbreeding"]["n_animals"] > 0
    assert payload["inbreeding"]["slope_kg_per_unit_F"] != 0.0
    assert (
        payload["paper_comparison"]["paper_kg_per_1pct_inbreeding"]
        == PAPER_INBREEDING_DEPRESSION_KG_PER_1PCT
    )


def test_per_1pct_is_slope_divided_by_100() -> None:
    # Regression guard: an earlier version multiplied by 100 (gave -7402 kg per 1%).
    inb = _payload()["inbreeding"]
    assert inb["slope_kg_per_1pct_F"] == pytest.approx(inb["slope_kg_per_unit_F"] / 100.0)
    lo, hi = inb["slope_kg_per_1pct_F_ci"]
    assert lo < inb["slope_kg_per_1pct_F"] < hi
    assert -2.0 < inb["slope_kg_per_1pct_F"] < 0.0


def test_recessive_attenuation_is_small_and_het_model_is_flagged() -> None:
    na = _payload()["nonadditive_model"]
    assert na["attenuation_fraction_homozygote_only"] < na["attenuation_fraction"]
    assert abs(na["attenuation_fraction_homozygote_only"]) < 0.10
    assert any("NOT evidence" in note for note in _payload()["notes"])
