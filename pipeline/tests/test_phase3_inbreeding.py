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


def test_inbreeding_export_shape() -> None:
    payload = json.loads((DATA / "inbreeding.json").read_text())
    assert payload["schema_version"] == "1.0.0"
    assert payload["inbreeding"]["n_animals"] > 0
    assert payload["inbreeding"]["inbreeding_depression_slope_kg_per_unit_f"] != 0.0
    assert (
        payload["paper_comparison"]["paper_kg_per_1pct_inbreeding"]
        == PAPER_INBREEDING_DEPRESSION_KG_PER_1PCT
    )
    assert "paper_nonadditive_fraction" in payload["paper_comparison"]
