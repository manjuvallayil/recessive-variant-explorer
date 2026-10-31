from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "web" / "public" / "data"

pytestmark = pytest.mark.skipif(
    not (DATA / "detection_summary.json").exists(),
    reason="run `uv run python pipeline/phase4_detection.py` first",
)


def test_detection_export_shape() -> None:
    payload = json.loads((DATA / "detection_summary.json").read_text())
    assert payload["schema_version"] == "1.0.0"
    assert payload["detection"]["n_total"] > 0
    assert payload["detection"]["n_detected"] >= 0
