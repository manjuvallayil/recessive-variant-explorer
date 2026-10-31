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


def test_qtl_truth_table_shape() -> None:
    truth = json.loads((DATA / "detection_summary.json").read_text())["qtl_truth"]
    rows = truth["rows"]
    assert len(rows) == 30
    assert all(0.0 < r["maf"] <= 0.06 for r in rows)  # truth MAF up to ~5.1%
    assert all(sum(r["genotype_counts"].values()) > 70_000 for r in rows)
    assert sum(r["n_qtl"] for r in truth["by_maf_bin"]) == 30
    assert sum(r["n_qtl"] for r in truth["by_effect_bin"]) == 30


def test_dominance_detects_more_truth_qtl_than_additive() -> None:
    # Paper claim: the additive model misses recessive QTL that the dominance test finds.
    truth = json.loads((DATA / "detection_summary.json").read_text())["qtl_truth"]
    n_dom = sum(r["n_detected_dominance"] for r in truth["by_maf_bin"])
    n_add = sum(r["n_detected_additive"] for r in truth["by_maf_bin"])
    assert n_dom > n_add
