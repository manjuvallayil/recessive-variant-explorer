"""Phase 4d: analytic power calculator (noncentral chi-square) for recessive/class tests."""

from __future__ import annotations

import csv
import json
import math
import subprocess
from datetime import UTC, datetime
from pathlib import Path

from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web" / "public" / "data"
EXT_FIG1 = ROOT / "dryad" / "source_data" / "source_data.ExtFig1.csv"

PAPER_CLAIM = {
    "maf": 0.025,
    "variance_explained": 0.001,
    "heritability": 0.25,
    "n_animals": 50_000,
    "power_recessive_class": 0.90,
}
VARIANCE_EXPLAINED = 0.001


def code_version() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


def normal_cdf(z: float) -> float:
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def power_for(
    variance_explained: float,
    sample_size: int,
    df: int = 1,
    alpha: float = 5e-8,
) -> float:
    lam = sample_size * variance_explained / (1.0 - variance_explained)
    crit = float(chi2.ppf(1.0 - alpha, df))
    root_lam = math.sqrt(lam)
    root_x = math.sqrt(crit)
    sd = math.sqrt(df)
    return float(
        1.0 - normal_cdf((root_x - root_lam) / sd) + normal_cdf(-((root_x + root_lam) / sd))
    )


def compute_validation() -> dict:
    rows = []
    max_abs_diff = 0.0
    worst = None
    with open(EXT_FIG1, newline="") as fh:
        for row in csv.DictReader(fh):
            n = int(round(float(row["N"])))
            model = row["model"]
            p_csv = float(row["Power"])
            df = 2 if "Genotype Class" in model else 1
            p_an = float(power_for(VARIANCE_EXPLAINED, n, df=df))
            diff = abs(p_an - p_csv)
            if diff > max_abs_diff:
                max_abs_diff = diff
                worst = {"model": model, "n": n, "analytic": p_an, "source": p_csv, "diff": diff}
            rows.append(
                {
                    "model": model,
                    "n": n,
                    "source_power": p_csv,
                    "analytic_power": p_an,
                }
            )
    return {
        "source_ext_fig1": {
            "variance_explained": VARIANCE_EXPLAINED,
            "n": 50_000,
            "cow_recessive_power": 0.947,
            "cow_genotype_class_power": 0.911,
            "analytic_recessive_power": float(power_for(VARIANCE_EXPLAINED, 50_000, df=1)),
            "analytic_genotype_class_power": float(power_for(VARIANCE_EXPLAINED, 50_000, df=2)),
        },
        "max_abs_diff_vs_source": max_abs_diff,
        "worst_row": worst,
        "paper_claim": PAPER_CLAIM,
        "notes": [
            "Phase 4d trimmed: analytic power via Z-transformation of noncentral chi-square.",
        ],
    }


def compute_power() -> dict:
    n = float(PAPER_CLAIM["n_animals"])
    return {
        "schema_version": "1.0.0",
        "generated_at": datetime.now(UTC).isoformat(),
        "code_version": code_version(),
        "source_files": {"ext_fig1": "dryad/source_data/source_data.ExtFig1.csv"},
        "power": {
            "test": "recessive_class",
            "variance_explained": VARIANCE_EXPLAINED,
            "sample_size": int(n),
            "lambda": n * VARIANCE_EXPLAINED / (1.0 - VARIANCE_EXPLAINED),
            "gws_alpha": 5e-08,
            "critical_chi2_df1": float(chi2.ppf(1.0 - 5e-8, 1)),
            "power": float(power_for(VARIANCE_EXPLAINED, n, df=1)),
        },
        "paper_claim": PAPER_CLAIM,
        "validation": compute_validation(),
        "notes": [
            "Phase 4d trimmed: analytic power via Z-transformation of noncentral chi-square.",
            "Paper claim: MAF 2.5%, variance 0.1%, h2=0.25, n=50,000 -> ~90% power.",
        ],
    }


def main() -> None:
    payload = compute_power()
    (OUT / "power_reference.json").write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
