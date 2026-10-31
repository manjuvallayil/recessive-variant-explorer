from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "web" / "public" / "data"
FASTA = ROOT / "pipeline" / "data" / "refseq_proteins.fasta"

pytestmark = pytest.mark.skipif(
    not (DATA / "plm_scores.json").exists(),
    reason="run `uv run python pipeline/phase4_plm.py` first",
)


def _payload() -> dict:
    return json.loads((DATA / "plm_scores.json").read_text())


def test_three_candidates_arsi_dropped() -> None:
    p = _payload()
    assert [g["gene"] for g in p["genes"]] == ["PLCD4", "MTRF1", "DPF2"]
    assert "ARSI" in p["dropped"]


def test_matrix_shape_and_wild_type_rows() -> None:
    for g in _payload()["genes"]:
        assert len(g["matrix"]) == g["length"], g["gene"]
        assert len(g["matrix_columns"]) == 22
        pos, wt = g["position"], g["wt"]
        assert g["matrix"][pos - 1][1] == wt, g["gene"]


def test_percentiles_in_range_and_no_ids() -> None:
    for g in _payload()["genes"]:
        assert 0.0 <= g["damage_percentile"] <= 100.0, g["gene"]
    assert '"ID"' not in (DATA / "plm_scores.json").read_text()


def test_fasta_residues_match_wt() -> None:
    seqs: dict[str, str] = {}
    name = None
    for line in FASTA.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        else:
            seqs[name] += line.strip()
    for g in _payload()["genes"]:
        assert seqs[g["accession"]][g["position"] - 1] == g["wt"], g["gene"]
