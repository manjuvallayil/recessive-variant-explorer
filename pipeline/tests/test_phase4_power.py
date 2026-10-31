from __future__ import annotations

import json
from pathlib import Path

import pytest

from pipeline.phase4_power import PAPER_CLAIM

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "web" / "public" / "data"

pytestmark = pytest.mark.skipif(
    not (DATA / "power_reference.json").exists(),
    reason="run `uv run python pipeline/phase4_power.py` first",
)


def test_power_export_shape() -> None:
    payload = json.loads((DATA / "power_reference.json").read_text())
    assert payload["schema_version"] == "1.0.0"
    assert payload["power"]["sample_size"] == PAPER_CLAIM["n_animals"]
    assert "source_ext_fig1" in payload["validation"]


def test_power_validation_passes_against_ext_fig1() -> None:
    payload = json.loads((DATA / "power_reference.json").read_text())
    validation = payload["validation"]
    assert validation["passes"] is True
    assert validation["max_abs_diff_vs_source"] <= validation["tolerance_abs"]
    assert set(validation["compared_models"]) == {"Cow - Recessive", "Cow - Genotype Class"}
