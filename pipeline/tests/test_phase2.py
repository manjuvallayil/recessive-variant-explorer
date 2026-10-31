from __future__ import annotations

import json
import math
from pathlib import Path

import pytest
from scipy import stats

from pipeline.phase2 import neg_log10_p_from_z

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "web" / "public" / "data"

PHASE2_FILES = [
    "manhattan_dom.json",
    "manhattan_add.json",
    "qtl_contrast.json",
    "chip_vs_sequence.json",
    "imputation_by_maf.json",
    "regional/index.json",
    "regional/chr2_22.json",
    "regional/plcd4.json",
    "regional/fgd4.json",
    "regional/arsi.json",
    "regional/mtrf1.json",
    "regional/plag1.json",
    "regional/myh1.json",
    "regional/galnt2.json",
    "regional/dpf2.json",
    "regional/mus81.json",
]

pytestmark = pytest.mark.skipif(
    not (DATA / "qtl_contrast.json").exists(),
    reason="run `uv run python pipeline/phase2.py` first",
)


def load(rel_path: str) -> dict:
    return json.loads((DATA / rel_path).read_text())


def test_log_ndtr_extreme_z_no_underflow() -> None:
    # Moderate Z matches standard normal survival function
    z_mod = 5.45131
    expected = -math.log10(2.0 * float(stats.norm.sf(z_mod)))
    assert abs(neg_log10_p_from_z(z_mod) - expected) < 1e-6

    # Extreme Z (PLAG1 P ~ 4.97e-442, Z ~ 44.991) where float64 P underflows to 0.0
    nlp_extreme = neg_log10_p_from_z(44.991)
    assert 440.0 < nlp_extreme < 443.0


def test_all_phase2_exports_exist_and_fit_size_budget() -> None:
    manifest = load("manifest.json")
    by_file = {entry["file"]: entry for entry in manifest["files"]}
    for rel_path in PHASE2_FILES:
        assert (DATA / rel_path).exists(), f"missing {rel_path}"
        assert rel_path in by_file, f"{rel_path} missing from manifest.json"
        payload = load(rel_path)
        assert payload["schema_version"] == "1.0.0"
        assert "generated_at" in payload
        assert "code_version" in payload
        assert "source_files" in payload
        assert by_file[rel_path]["gzip_bytes"] < 500_000, f"{rel_path} exceeds 500 KB gz"
    assert manifest["total_bytes"] < 10_000_000


def test_manhattan_thinning_and_unclipped_peaks() -> None:
    dom = load("manhattan_dom.json")
    add = load("manhattan_add.json")

    for payload in (dom, add):
        total_pts = 0
        for trait in ("bodyweight", "stature"):
            t = payload["traits"][trait]
            meta = t["meta"]
            cols = t["columns"]
            n = len(cols["chr"])
            total_pts += n
            assert n == meta["n_thinned_points"]
            assert n <= 60_000
            assert len(cols["bp_umd"]) == n
            assert len(cols["bp_ars"]) == n
            assert len(cols["neg_log10_p_x100"]) == n
            assert 0.001 < meta["raw_unmapped_ars_fraction"] < 0.01
        assert total_pts <= 60_000

    # PLAG1 standard-additive peak on BTA14 exceeds int16 max (32767) and is unclipped
    max_add_bw = max(add["traits"]["bodyweight"]["columns"]["neg_log10_p_x100"])
    assert max_add_bw > 32_000
    max_dom_bw = max(dom["traits"]["bodyweight"]["columns"]["neg_log10_p_x100"])
    assert max_dom_bw > 5_500


def test_qtl_contrast_7_of_8_invisible_to_additive() -> None:
    qtl = load("qtl_contrast.json")
    s = qtl["summary"]
    assert s["n_primary_non_additive_qtl"] == 8
    assert s["detected_by_dominance_gws"] == 8
    assert s["candidate_detected_by_standard_additive_gws"] == 1
    assert s["detected_in_paper_cojo_standard_additive"] == 1
    assert s["invisible_to_standard_additive_at_candidate_and_cojo"] == 7

    by_name = {r["qtl"]: r for r in qtl["rows"]}
    assert by_name["PLAG1"]["bodyweight"]["candidate_detected_by_additive_gws"] is True
    for rec in ("Chr2:22", "PLCD4", "FGD4", "ARSI", "MTRF1", "GALNT2", "DPF2/MUS81"):
        bw = by_name[rec]["bodyweight"]
        assert bw["detected_by_dominance_gws"] is True
        assert bw["candidate_detected_by_additive_gws"] is False
        assert bw["candidate_add_neg_log10_p"] < 3.0


def test_regional_plots_and_extfig3_crosscheck() -> None:
    idx = load("regional/index.json")
    assert len(idx["loci"]) == 10
    for item in idx["loci"]:
        reg = load(item["file"])
        cand = reg["candidate_variant"]
        c_idx = cand["index_in_columns"]
        assert c_idx is not None
        assert reg["columns"]["r2_x1000"][c_idx] == 1000
        assert reg["columns"]["bp_umd"][c_idx] == cand["bp_umd"]
        if reg["gene"] is not None:
            symbols = {g["symbol"] for g in reg["genes"]}
            assert reg["gene"] in symbols, f"{reg['gene']} not in {item['file']} gene track"

    myh1 = load("regional/myh1.json")
    ext3 = myh1["extfig3_comparison"]
    assert ext3["n_matched"] == 7862
    assert ext3["r2_correlation"] > 0.99


def test_chip_blind_spot() -> None:
    chip = load("chip_vs_sequence.json")
    assert chip["n_chip_distinct_positions"] == 45_865
    assert chip["fraction_of_chip_in_sequence"] == 1.0
    assert 0.002 < chip["fraction_of_sequence_on_chip"] < 0.004
    assert chip["summary"]["major_recessive_syndromes_significant_on_50k_chip"] == ["DPF2/MUS81"]
    assert set(chip["summary"]["primary_8_loci_significant_on_50k_chip"]) == {
        "Chr2:22",
        "PLAG1",
        "DPF2/MUS81",
    }


def test_imputation_by_maf() -> None:
    imp = load("imputation_by_maf.json")
    assert imp["summary"]["all_bw_recessive_candidates_gt_0_95"] is True
    assert imp["summary"]["min_ar2_bw_recessive_candidates"] >= 0.979
    bins = imp["maf_bins"]
    assert len(bins) == 9
    assert bins[0]["mean_ar2"] < bins[-1]["mean_ar2"]
