from __future__ import annotations

import json
from pathlib import Path

import pytest

from pipeline.phase3_rna import EXPRESSION_GENES, RNA_PROXY_GENE, SPLICING_GENES

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "web" / "public" / "data"

pytestmark = pytest.mark.skipif(
    not (DATA / "rna_phase3.json").exists(),
    reason="run `uv run python pipeline/phase3_rna.py` first",
)


def test_rna_phase3_export_shape() -> None:
    payload = json.loads((DATA / "rna_phase3.json").read_text())
    assert payload["schema_version"] == "1.0.0"
    assert set(payload["rna"]["splice"]) == set(SPLICING_GENES)
    assert set(payload["rna"]["expression"]) == set(EXPRESSION_GENES)
    assert payload["proxy"]["MUS81"]["proxy_gene"] == RNA_PROXY_GENE
    for section in ("splice", "expression"):
        for _gene, stats in payload["rna"][section].items():
            assert stats["n_animals"] > 0
            assert stats["mean_value"] > 0.0
