"""Phase 4d: analytic power calculator (noncentral chi-square) for recessive/class tests."""

from __future__ import annotations

import csv
import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web" / "public" / "data"
EXT_FIG1 = ROOT / "dryad" / "source_data" / "source_data.ExtFig1.csv"

PAPER_CLAIM = {
    "maf": 0.025,
    "variance_explained": 0.10,
    "heritability": 0.25,
    "n_animals": 50_000,
    "power_recessive_class": 0.90,
}


def code_version() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


def noncentral_chi2_power(
    effect_size_sd: float, sample_size: int, df: int = 1, alpha: float = 5e-8
) -> float:
    """Power for a one-shot test with noncentrality lambda = N*delta^2/sigma^2."""
    lam = sample_size * effect_size_sd**2
    crit = float(chi2.ppf(1.0 - alpha, df))
    return float(chi2.sf(crit, df, lam))


def compute_validation() -> dict:
    n = float(PAPER_CLAIM["n_animals"])
    h2 = float(PAPER_CLAIM["heritability"])
    q = float(PAPER_CLAIM["maf"])
    delta_sd = float(np.sqrt(h2 / (2.0 * q * (1.0 - q))))
    analytic_recessive = float(noncentral_chi2_power(delta_sd, n, df=1))
    analytic_genotype_class = float(noncentral_chi2_power(delta_sd, n, df=2))
    cow_recessive = None
    cow_gc = None
    with open(EXT_FIG1, newline="") as fh:
        for row in csv.DictReader(fh):
            if float(row["N"]) != n:
                continue
            if row["model"] == "Cow - Recessive":
                cow_recessive = float(row["Power"])
            if row["model"] == "Cow - Genotype Class":
                cow_gc = float(row["Power"])
    return {
        "source_ext_fig1": {
            "n": n,
            "cow_recessive_power": cow_recessive,
            "cow_genotype_class_power": cow_gc,
            "analytic_recessive_power": analytic_recessive,
            "analytic_genotype_class_power": analytic_genotype_class,
        },
        "paper_claim": PAPER_CLAIM,
        "notes": [
            "Phase 4d trimmed: analytic power vs source_data.ExtFig1.csv and paper claim.",
        ],
    }


def compute_power() -> dict:
    n = float(PAPER_CLAIM["n_animals"])
    h2 = float(PAPER_CLAIM["heritability"])
    q = float(PAPER_CLAIM["maf"])
    delta_sd = float(np.sqrt(h2 / (2.0 * q * (1.0 - q))))
    power_recessive = float(noncentral_chi2_power(delta_sd, n, df=1))
    return {
        "schema_version": "1.0.0",
        "generated_at": datetime.now(UTC).isoformat(),
        "code_version": code_version(),
        "source_files": {"ext_fig1": "dryad/source_data/source_data.ExtFig1.csv"},
        "power": {
            "test": "recessive_class",
            "effect_size_sd": delta_sd,
            "sample_size": int(n),
            "noncentrality": n * delta_sd**2,
            "gws_alpha": 5e-08,
            "power": power_recessive,
        },
        "paper_claim": PAPER_CLAIM,
        "validation": compute_validation(),
        "notes": [
            "Phase 4d trimmed: analytic noncentral chi-square power for recessive/class test.",
            "Paper claim: MAF 2.5%, variance 0.1%, h2=0.25, n=50,000 -> ~90% power.",
        ],
    }


def main() -> None:
    payload = compute_power()
    (OUT / "power_reference.json").write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
