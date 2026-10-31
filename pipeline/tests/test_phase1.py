from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "web" / "public" / "data"

PHASE1_FILES = [
    "loci.json",
    "genotype_class_stats.json",
    "effects.json",
    "validation_milk.json",
    "pleiotropy.json",
    "frequencies.json",
    "depletion.json",
    "carrier_cross.json",
    "national_burden.json",
    "sire_amplification.json",
    "prospective_calves.json",
    "manifest.json",
]

pytestmark = pytest.mark.skipif(
    not (DATA / "effects.json").exists(),
    reason="run `uv run python pipeline/phase1.py` first",
)


def load(name: str) -> dict:
    return json.loads((DATA / name).read_text())


def test_all_phase1_exports_exist_and_have_metadata() -> None:
    manifest = load("manifest.json")
    manifest_files = {f["file"] for f in manifest["files"]}
    for fname in PHASE1_FILES:
        path = DATA / fname
        assert path.exists(), f"missing {fname}"
        payload = load(fname)
        assert payload.get("schema_version") == "1.0.0"
        assert "generated_at" in payload
        assert "code_version" in payload
        if fname != "manifest.json":
            assert fname in manifest_files
            assert "source_files" in payload


def test_no_animal_or_sire_ids_in_exports() -> None:
    forbidden_keys = {"ID", "IID", "SIRE", "DAM", "SIRE_ID", "MGS_ID", "KEY", "anml_key", "animal"}
    for fname in PHASE1_FILES:
        text = (DATA / fname).read_text()
        for key in forbidden_keys:
            assert not re.search(rf'"{key}"\s*:', text), f"forbidden key {key} in {fname}"


def test_small_n_suppression_rule() -> None:
    for fname in ("genotype_class_stats.json", "validation_milk.json"):
        payload = load(fname)
        for row in payload["rows"]:
            for k, cstat in row["class_stats"].items():
                assert "n" in cstat
                if cstat["n"] < 5:
                    assert cstat["mean"] is None, f"un-suppressed mean in {fname} class {k}"
                    assert cstat["sd"] is None
                    assert cstat["se"] is None
                    assert all(v is None for v in cstat["quantiles"].values())


def test_table1_reproduction_estimator_b_vs_paper() -> None:
    effects = load("effects.json")
    summary = effects["table1_comparison_summary"]
    assert summary["closest_by_hom_alt_mae"] == "B"
    assert summary["B"]["mae_bodyweight_hom_alt"] < 1.5

    by_key = {(r["trait"], r["locus"], r["estimator"]): r for r in effects["rows"]}
    targets = effects["table1_targets"]

    # Large bodyweight recessive effects: same sign and within 15% of Table 1
    for locus in ("PLCD4", "FGD4", "MTRF1", "GALNT2", "DPF2", "MUS81"):
        row_b = by_key[("bodyweight", locus, "B")]
        _, paper_hom = targets["bodyweight"][locus]
        assert row_b["hom_alt"] * paper_hom > 0
        rel_err = abs(row_b["hom_alt"] - paper_hom) / abs(paper_hom)
        assert rel_err < 0.10, f"{locus} bodyweight rel_err={rel_err:.3f}"

    # Moderate bodyweight recessive loci (Chr2:22, ARSI): within 15%
    for locus in ("Chr2:22", "ARSI"):
        row_b = by_key[("bodyweight", locus, "B")]
        _, paper_hom = targets["bodyweight"][locus]
        assert row_b["hom_alt"] * paper_hom > 0
        rel_err = abs(row_b["hom_alt"] - paper_hom) / abs(paper_hom)
        assert rel_err < 0.15

    # Stature (scaled 5 cm per TOP score unit): major stature loci within 10%
    for locus in ("PLCD4", "GALNT2", "FGD4", "DPF2", "MUS81", "ARSI", "MTRF1"):
        row_b = by_key[("stature", locus, "B")]
        _, paper_hom = targets["stature"][locus]
        assert row_b["hom_alt"] * paper_hom > 0
        rel_err = abs(row_b["hom_alt"] - paper_hom) / abs(paper_hom)
        assert rel_err < 0.10, f"{locus} stature rel_err={rel_err:.3f}"

    # BCS: PLCD4, FGD4, MTRF1, MYH1 hom_alt same sign and within 0.03 units
    for locus in ("PLCD4", "FGD4", "MTRF1", "MYH1"):
        row_b = by_key[("bcs", locus, "B")]
        _, paper_hom = targets["bcs"][locus]
        assert row_b["hom_alt"] * paper_hom > 0
        assert abs(row_b["hom_alt"] - paper_hom) < 0.03


def test_unweighted_vs_rinverse_sensitivity_reported() -> None:
    effects = load("effects.json")
    estimators = {r["estimator"] for r in effects["rows"]}
    assert {"A", "A_rinverse", "B"} <= estimators


def test_extfig2c_aggregate_comparison_matches_loso() -> None:
    effects = load("effects.json")
    ext2c = {r["locus"]: r for r in effects["extfig2c_comparison"]["rows"]}
    for locus in ("Chr2:22", "PLCD4", "ARSI", "MTRF1", "GALNT2", "DPF2", "chr14:25015640"):
        assert abs(ext2c[locus]["delta_vs_B_primary_hom_alt"]) < 1e-2
        assert abs(ext2c[locus]["delta_vs_B_primary_het"]) < 1e-2
    # FGD4 in ExtFig2c uses window 6; verify window_sensitivity_rows matches it
    win_sens = {
        (r["trait"], r["locus"], r["window"]): r for r in effects["window_sensitivity_rows"]
    }
    fgd4_w6 = win_sens[("bodyweight", "FGD4", 6)]
    assert abs(ext2c["FGD4"]["hom_alt"] - fgd4_w6["hom_alt"]) < 1e-2


def test_frequencies_and_depletion() -> None:
    freq = load("frequencies.json")
    for locus, payload in freq["loci"].items():
        hf = payload["breeds"]["HF"]
        j = payload["breeds"]["J"]
        assert abs(hf["alt_freq"] - hf["paper_hf"]) < 0.005, f"{locus} HF freq mismatch"
        assert abs(j["alt_freq"] - j["paper_j"]) < 0.011, f"{locus} J freq mismatch"
    assert freq["cumulative"]["HF"]["sum_of_carrier_frequencies_6_unique_loci"] > 0.40
    assert freq["cumulative"]["J"]["sum_of_carrier_frequencies_6_unique_loci"] > 0.40

    dep = load("depletion.json")
    by_locus = {r["locus"]: r for r in dep["rows"]}
    for locus in ("PLCD4", "FGD4", "MTRF1", "GALNT2"):
        assert by_locus[locus]["observed_over_expected"] < 0.75
        assert by_locus[locus]["lrt_p_one_sided"] < 1e-3


def test_national_burden_and_prospective_calves() -> None:
    nb = load("national_burden.json")
    assert nb["n_animals"] == 2_799_022
    bw_rate = nb["any_bw_significant_mutation"]["rate_per_animal"]
    bw_scaled = nb["any_bw_significant_mutation"]["scaled_per_year"]
    assert 0.009 < bw_rate < 0.012
    assert 8500 < bw_scaled < 10500

    calves = load("prospective_calves.json")
    assert calves["n_calves"] == 568
    assert calves["effects"]["PLCD4"]["effect_vs_control"] < -40
    assert calves["effects"]["PLCD4"]["p"] < 1e-10
    assert calves["effects"]["FGD4"]["effect_vs_control"] < -5
    assert calves["effects"]["FGD4"]["p"] < 0.05
