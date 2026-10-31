"""Phase 4d: analytic power calculator (noncentral chi-square) for recessive/class tests."""

from __future__ import annotations

import csv
import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

from scipy.stats import chi2, ncx2

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


def noncentrality(variance_explained: float, sample_size: int) -> float:
    return sample_size * variance_explained / (1.0 - variance_explained)


def power_for(
    variance_explained: float,
    sample_size: int,
    df: int = 1,
    alpha: float = 5e-8,
) -> float:
    """Exact power: P(noncentral chi2(df, lambda) > central chi2 critical value at alpha)."""
    lam = noncentrality(variance_explained, sample_size)
    crit = float(chi2.ppf(1.0 - alpha, df))
    return float(ncx2.sf(crit, df, lam))


# Only the cow recessive and genotype-class curves are the closed-form test this module
# implements (1-df and 2-df noncentral chi2 at a single variance explained). The additive
# curves use a different coding (variance captured by the additive term is not 0.1%), and
# the bull/sire curves include progeny-based designs. They are reported, not validated.
COMPARABLE_MODELS = {"Cow - Recessive": 1, "Cow - Genotype Class": 2}


def compute_validation() -> dict:
    rows = []
    max_abs_diff = 0.0
    worst = None
    not_compared: dict[str, float] = {}
    with open(EXT_FIG1, newline="") as fh:
        for row in csv.DictReader(fh):
            n = int(round(float(row["N"])))
            model = row["model"]
            p_csv = float(row["Power"])
            if model not in COMPARABLE_MODELS:
                not_compared[model] = max(not_compared.get(model, 0.0), abs(p_csv))
                continue
            df = COMPARABLE_MODELS[model]
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
        "compared_models": sorted(COMPARABLE_MODELS),
        "not_compared_models": sorted(not_compared),
        "tolerance_abs": 0.005,
        "passes": max_abs_diff <= 0.005,
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
            "Analytic power: exact noncentral chi-square (scipy.stats.ncx2), df = 1 or 2.",
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
            "Analytic power: exact noncentral chi-square (scipy.stats.ncx2), df = 1 or 2.",
            "Paper claim: MAF 2.5%, variance 0.1%, h2=0.25, n=50,000 -> ~90% power.",
        ],
    }


def main() -> None:
    payload = compute_power()
    (OUT / "power_reference.json").write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
