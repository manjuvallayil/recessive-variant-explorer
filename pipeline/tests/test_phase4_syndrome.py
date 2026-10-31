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


def test_syndrome_export_shape() -> None:
    payload = json.loads((DATA / "syndrome_signal.json").read_text())
    assert payload["schema_version"] == "1.0.0"
    assert payload["syndrome_signal"]["baseline"] == "body_weight_alone"
    assert payload["syndrome_signal"]["n_records"] > 0
