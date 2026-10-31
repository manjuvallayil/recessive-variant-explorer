from __future__ import annotations

import json
from pathlib import Path

import pytest

from pipeline.phase3_farm import FARM_GEN7, FARM_GEN8

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "web" / "public" / "data"

pytestmark = pytest.mark.skipif(
    not (DATA / "farm_profiles.json").exists(),
    reason="run `uv run python pipeline/phase3_farm.py` first",
)


def test_farm_profiles_export_shape() -> None:
    payload = json.loads((DATA / "farm_profiles.json").read_text())
    assert payload["schema_version"] == "1.0.0"
    profiles = payload["farm_profiles"]
    assert FARM_GEN8 in profiles
    assert FARM_GEN7 in profiles
    for entry in profiles.values():
        assert isinstance(entry, dict)
        assert "small_n_caveat" in payload
