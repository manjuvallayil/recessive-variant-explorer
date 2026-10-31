"""Phase 3c: inbreeding and inbreeding depression (trimmed)."""

from __future__ import annotations

import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.api import OLS

ROOT = Path(__file__).resolve().parents[1]
DRYAD = ROOT / "dryad"
PHENS = DRYAD / "phens"
OTHER = DRYAD / "other"
OUT = ROOT / "web" / "public" / "data"

PAPER_INBREEDING_DEPRESSION_KG_PER_1PCT = -0.77
PAPER_NONADDITIVE_FRACTION = 0.235


def code_version() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


def compute_inbreeding() -> dict:
    other7 = pd.read_csv(OTHER / "other_7.txt", sep=r"\s+")
    other7.columns = ["ID", "F", "value"]
    phen1 = pd.read_csv(PHENS / "phen_1.1.txt", sep=r"\s+")
    merged = phen1.merge(other7, on="ID", how="inner")
    merged = merged.dropna(subset=["yd_trait", "F"])
    x = merged["F"].to_numpy(dtype=float)
    y = merged["yd_trait"].to_numpy(dtype=float)
    model = OLS(y, np.column_stack([np.ones(len(x)), x])).fit()
    slope_kg_per_unit_f = float(model.params[1])
    r2_f = float(model.rsquared)

    return {
        "schema_version": "1.0.0",
        "generated_at": datetime.now(UTC).isoformat(),
        "code_version": code_version(),
        "source_files": {
            "other_7": "dryad/other/other_7.txt",
            "phen_1.1": "dryad/phens/phen_1.1.txt",
        },
        "inbreeding": {
            "n_animals": int(len(merged)),
            "inbreeding_coefficient_mean": float(merged["F"].mean()),
            "inbreeding_depression_slope_kg_per_unit_f": slope_kg_per_unit_f,
            "inbreeding_depression_slope_kg_per_1pct_inbreeding": slope_kg_per_unit_f * 100,
            "r_squared_weight_on_f": r2_f,
        },
        "paper_comparison": {
            "paper_kg_per_1pct_inbreeding": PAPER_INBREEDING_DEPRESSION_KG_PER_1PCT,
            "paper_nonadditive_fraction": PAPER_NONADDITIVE_FRACTION,
        },
        "notes": [
            "Phase 3c trimmed: inbreeding-depression slope from OLS of bodyweight on inbreeding.",
            "Joint non-additive attenuation not computed here; paper 23.5%.",
        ],
    }


def main() -> None:
    payload = compute_inbreeding()
    (OUT / "inbreeding.json").write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
