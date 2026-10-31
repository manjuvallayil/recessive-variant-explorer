import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SUMMARY = ROOT / "work" / "phase0_summary.json"

pytestmark = pytest.mark.skipif(
    not SUMMARY.exists(), reason="run `uv run python pipeline/phase0_verify.py` first"
)


@pytest.fixture(scope="module")
def checks() -> dict:
    return json.loads(SUMMARY.read_text())["checks"]


def test_all_checks_pass(checks: dict) -> None:
    failed = [
        name
        for name, payload in checks.items()
        if name != "environment" and not payload.get("pass")
    ]
    assert failed == []


def test_v1_iid_join_reproduces_recessive_effect(checks: dict) -> None:
    v1 = checks["V1"]
    assert v1["iid_join_adjusted"]["hom_alt_minus_ref"] < -100
    assert abs(v1["row_order_join_raw"]["hom_alt_minus_ref"]) < 20


def test_v2_allele_coding_matches_paper(checks: dict) -> None:
    v2 = checks["V2"]
    assert v2["max_abs_delta_all"] < 0.001
    assert v2["max_abs_delta_breed"] < 0.02
    assert all(row["a1_freq"] > 0 for row in v2["lead_variants"])


def test_v3_id_spaces_and_gen3_row_order(checks: dict) -> None:
    v3 = checks["V3"]
    assert v3["other_5_rows"] == 2197
    assert v3["other_5_union_n"] == 4276
    assert all(row["row_identity"] == 1.0 for row in v3["gen3_row_order"])
    assert all(not row["fam_order_equal"] for row in v3["gen3_row_order"])
    by_pair = {(row["a"], row["b"]): row for row in v3["comparisons"]}
    assert by_pair[("gen_1", "phen_1.1")]["only_a"] == 82
    assert by_pair[("gen_3", "phen_4.1")]["only_a"] == 368


def test_v4_weighting_choice_is_documented(checks: dict) -> None:
    v4 = checks["V4"]
    assert v4["chosen_for_pipeline"] == "unweighted"
    assert all(error < 10 for error in v4["mean_abs_error_hom_alt"].values())


def test_v5_inbreeding_units_and_slope(checks: dict) -> None:
    v5 = checks["V5"]
    assert 0 <= v5["F_min"] <= v5["F_max"] <= 5
    assert -1.0 < v5["slope_yd_per_1pct"] < -0.5
    assert v5["corr_value_yd"] > v5["corr_value_adj"]


def test_v6_carrier_cross_column_order(checks: dict) -> None:
    v6 = checks["V6"]
    assert v6["other_5_rows"] == 2197
    assert v6["best_permutation"] in {"(0, 1, 2)", "(0, 2, 1)"}
    assert v6["best_rate"] < 0.01


def test_v7_national_burden_keys(checks: dict) -> None:
    v7 = checks["V7"]
    assert v7["gen_10_shape"] == [5550, 7]
    assert v7["sire_in_key_frac"] == 1.0
    assert v7["mgs_in_key_frac"] == 1.0
    assert v7["birth_year_min"] == 2010
    assert v7["birth_year_max"] == 2019


def test_v8_loso_trait_mapping(checks: dict) -> None:
    v8 = checks["V8"]
    mapping = {row["other_2_trait"]: row["matches_phen_1"] for row in v8["trait_mapping"]}
    assert mapping == {1: 1, 2: 3, 3: 2}
    assert len(v8["windows"]) == 13
    assert all(row["ncols"] == 501 for row in v8["windows"])


def test_v9_summary_stat_coverage(checks: dict) -> None:
    v9 = checks["V9"]
    assert v9["other_1_1"]["rows"] == 16_128_394
    assert v9["coverage_frac"]["other_9"] == 1.0
    assert v9["coverage_frac"]["other_10"] > 0.99
    assert v9["coverage_frac"]["other_15"] < 0.01
    assert v9["other_15_rows_in_other_1_frac"] == 1.0
    assert v9["other_1_chroms"] == list(range(1, 30))


def test_v10_rna_and_chip_id_links(checks: dict) -> None:
    v10 = checks["V10"]
    assert v10["gen6_iids_equal_phen7_key"]
    assert v10["gen6_iids_equal_phen7_2_id"]
    assert v10["gen_7_ids_equal_phen_6"]
    assert v10["gen_8_ids_equal_phen_5"]
    assert v10["gen_6_1_bim_rows"] == 17_180_153
    nearest = v10["candidate_nearest"]
    for locus in ["PLCD4", "FGD4", "ARSI", "DPF2", "GALNT2", "MYH1", "Chr2:22"]:
        assert nearest[locus]["distance"] == 0
    assert nearest["MTRF1"]["distance"] <= 150
    assert nearest["MUS81"]["distance"] <= 150
    ld = v10["mus81_proxy_ld"]
    assert ld["n"] == 350
    assert ld["distance_bp"] == 103
    assert 0 <= ld["r2"] <= 1


def test_v11_source_data_inventory(checks: dict) -> None:
    v11 = checks["V11"]
    assert v11["exists"]
    assert v11["n_csv"] == 10
    names = {row["file"] for row in v11["files"]}
    assert "source_data.ExtFig9b.csv" in names
    assert "source_data.ExtFig2c.csv" in names
