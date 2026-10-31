"""Phase 5 inputs: parameters for the client-side mating-risk simulator.

Reads only validated exports (frequencies, depletion, effects, national_burden) plus one aggregate
from dryad/other/other_11.txt (mean breed-weighted allele frequency over animals, used as the
national calibration input). Writes web/public/data/mating_inputs.json. No IDs are exported.

Loci (ARSI excluded: no breed frequency available; see docs/ASSUMPTIONS.md)
- PLCD4, FGD4, MTRF1, GALNT2, DPF2/MUS81 (body-weight significant set)
- MYH1 (BCS-suggestive, body-weight effect not significant; off by default in the simulator)
- DPF2/MUS81 is one locus: DPF2 frequencies, effects and depletion ratio are used as the proxy.
"""

from __future__ import annotations

import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web" / "public" / "data"
OTHER11 = ROOT / "dryad" / "other" / "other_11.txt"
FEMALES_PER_YEAR = 940_100  # from national_burden.json (paper value)

LOCI = [
    # id, label, frequency/effect/depletion source key, in_bw_set
    ("PLCD4", "PLCD4", "PLCD4", True),
    ("FGD4", "FGD4", "FGD4", True),
    ("MTRF1", "MTRF1", "MTRF1", True),
    ("GALNT2", "GALNT2", "GALNT2", True),
    ("DPF2_MUS81", "DPF2/MUS81", "DPF2", True),
    ("MYH1", "MYH1", "MYH1", False),
]


def code_version() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


def load(name: str) -> dict:
    return json.loads((OUT / name).read_text())


_OTHER11_CACHE: dict[str, pd.DataFrame] = {}


def _other11() -> pd.DataFrame:
    if "df" not in _OTHER11_CACHE:
        _OTHER11_CACHE["df"] = pd.read_csv(OTHER11, sep=r"\s+", header=0)
    return _OTHER11_CACHE["df"]


def national_mean_q(q_hf: float, q_j: float, q_all: float) -> dict:
    """Mean breed-weighted allele frequency over animals (other_11 dam breed composition).

    Animal breed composition is proxied by the dam's (dam_HF_16, dam_J_16) in sixteenths; the
    remainder uses the all-breed frequency. Only the aggregate is returned.
    """
    df = _other11()
    f_hf = df["dam_HF_16"].to_numpy(dtype=float) / 16.0
    f_j = df["dam_J_16"].to_numpy(dtype=float) / 16.0
    f_hf = np.clip(f_hf, 0, 1)
    f_j = np.clip(f_j, 0, 1)
    rem = np.clip(1.0 - f_hf - f_j, 0, 1)
    q_i = f_hf * q_hf + f_j * q_j + rem * q_all
    return {
        "n_animals": int(len(df)),
        "mean_breed_weighted_q": float(q_i.mean()),
        "mean_fraction_HF": float(f_hf.mean()),
        "mean_fraction_J": float(f_j.mean()),
    }


def compute() -> dict:
    freq = load("frequencies.json")["loci"]
    dep_rows = load("depletion.json")["rows"]
    effects = load("effects.json")["rows"]
    nb = load("national_burden.json")

    def freq_of(src: str, breed: str) -> float:
        return float(freq[src]["breeds"][breed]["alt_freq"])

    def depletion_of(src: str, breed: str) -> float | None:
        for r in dep_rows:
            if r["gene"] == src and r["breed"] == breed:
                return float(r["observed_over_expected"])
        return None

    def effect_of(src: str) -> dict:
        for r in effects:
            if r["trait"] == "bodyweight" and r["estimator"] == "B" and r["locus"] == src:
                return {"hom_kg": float(r["hom_alt"]), "hom_se_kg": float(r["se_hom_alt"]),
                        "het_kg": float(r["het"]), "het_se_kg": float(r["se_het"]),
                        "p_hom": float(r["p_hom_alt"])}
        raise KeyError(src)

    loci_out = []
    for lid, label, src, in_bw in LOCI:
        q_hf = freq_of(src, "HF")
        q_j = freq_of(src, "J")
        q_all_nb = float(nb["mutations"][src]["freq_all"])
        nat = national_mean_q(q_hf, q_j, q_all_nb)
        loci_out.append(
            {
                "id": lid,
                "label": label,
                "in_body_weight_set": in_bw,
                "q_HF": q_hf,
                "q_J": q_j,
                "depletion_ratio_HF": depletion_of(src, "HF"),
                "depletion_ratio_J": depletion_of(src, "J"),
                # Breed where the allele is common (HF: PLCD4/FGD4/GALNT2/DPF2; J: MTRF1/MYH1).
                "depletion_ratio_primary": depletion_of(src, "HF" if q_hf >= q_j else "J"),
                "depletion_primary_breed": "HF" if q_hf >= q_j else "J",
                "effect": effect_of(src),
                "national_mean_q": nat["mean_breed_weighted_q"],
            }
        )

    national = {
        "females_born_per_year": FEMALES_PER_YEAR,
        "n_animals_other_11": nat["n_animals"],
        "mean_fraction_HF": nat["mean_fraction_HF"],
        "mean_fraction_J": nat["mean_fraction_J"],
        "paper_homozygous_rate": nb["paper_homozygous_rate"],
        "paper_scaled_homozygotes": nb["paper_scaled_homozygotes"],
        "export_bw_significant_rate": nb["any_bw_significant_mutation"]["rate_per_animal"],
        "export_bw_significant_scaled": nb["any_bw_significant_mutation"]["scaled_per_year"],
        "export_bw_significant_mutations": nb["any_bw_significant_mutation"]["mutations_included"],
        "export_six_locus_rate": nb["any_mutation"]["rate_per_animal"],
        "export_six_locus_scaled": nb["any_mutation"]["scaled_per_year"],
        "export_per_mutation_scaled": {
            k: float(v["scaled_per_year"]) for k, v in nb["mutations"].items()
        },
    }

    return {
        "schema_version": "1.0.0",
        "generated_at": datetime.now(UTC).isoformat(),
        "code_version": code_version(),
        "source_files": {
            "frequencies": "web/public/data/frequencies.json",
            "depletion": "web/public/data/depletion.json",
            "effects": "web/public/data/effects.json",
            "national_burden": "web/public/data/national_burden.json",
            "other_11_aggregate": "dryad/other/other_11.txt (mean breed-weighted q only)",
        },
        "loci": loci_out,
        "national": national,
        "notes": [
            "Transmission: non-carrier 0, carrier 0.5, affected 1, unknown = breed-weighted q.",
            "Depletion ratios (observed/expected homozygotes, purebreds) are applied to "
            "homozygous calf probabilities when the toggle is on.",
            "Effects: LOSO-adjusted estimator B (paper-like). Het effect used only for DPF2/MUS81.",
            "DPF2/MUS81 treated as one locus using DPF2 frequency, effect and depletion ratio. "
            "DPF2-MUS81 LD is weak (r2 = 0.018 in 350 animals), so independence is an "
            "approximation.",
            "Calibration uses mean breed-weighted q over animals (dam breed composition as proxy); "
            "it ignores pedigree structure, which the national export uses.",
        ],
    }


def main() -> None:
    payload = compute()
    (OUT / "mating_inputs.json").write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
