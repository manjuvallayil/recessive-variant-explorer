"""Phase 3a: RNA-level follow-up for the recessive variant exploration.

Splice junction means from ``phen_7.1`` (FGD4, GALNT2) and expression means from
``phen_7.2`` (FGD4, GALNT2, MUS81), with DPF2 as the labelled proxy gene for MUS81.
"""

from __future__ import annotations

import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DRYAD = ROOT / "dryad"
PHENS = DRYAD / "phens"
OUT = ROOT / "web" / "public" / "data"

SPLICING_GENES = ["FGD4", "GALNT2"]
EXPRESSION_GENES = ["FGD4", "GALNT2", "MUS81"]
RNA_PROXY_GENE = "DPF2"


def code_version() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


def _mean_by_key(df: pd.DataFrame, value_cols: list[str]) -> dict[str, float]:
    means: dict[str, float] = {}
    for _, row in df.iterrows():
        key = str(row.iloc[0])
        vals = [float(row[c]) for c in value_cols if pd.notna(row[c])]
        means[key] = float(sum(vals) / len(vals)) if vals else 0.0
    return means


def _summary(means: dict[str, float]) -> dict:
    vals = list(means.values())
    return {
        "n_animals": len(vals),
        "mean_value": float(sum(vals) / len(vals)) if vals else 0.0,
    }


def compute_rna() -> dict:
    phen7_1 = pd.read_csv(PHENS / "phen_7.1.txt", sep=r"\s+")
    phen7_2 = pd.read_csv(PHENS / "phen_7.2.txt", sep=r"\s+")

    splice_means = {
        gene: _mean_by_key(phen7_1, [c for c in phen7_1.columns if c.startswith(f"{gene}_")])
        for gene in SPLICING_GENES
    }
    expr_means = {
        gene: _mean_by_key(phen7_2, [c for c in phen7_2.columns if c.strip() == gene])
        for gene in EXPRESSION_GENES
    }

    return {
        "schema_version": "1.0.0",
        "generated_at": datetime.now(UTC).isoformat(),
        "code_version": code_version(),
        "source_files": {
            "phen_7.1": "dryad/phens/phen_7.1.txt",
            "phen_7.2": "dryad/phens/phen_7.2.txt",
            "proxy_gene_for_mus81": "dryad/gens/gen_1.9.bim",
        },
        "rna": {
            "splice": {gene: _summary(splice_means[gene]) for gene in SPLICING_GENES},
            "expression": {gene: _summary(expr_means[gene]) for gene in EXPRESSION_GENES},
        },
        "proxy": {
            "MUS81": {
                "proxy_gene": RNA_PROXY_GENE,
                "note": "MUS81 exact variant absent from gen_6.1; DPF2 used as labelled proxy.",
            }
        },
        "notes": [
            "Phase 3a trimmed: splice/expression means for the recessive RNA targets.",
            "MUS81 expression reported via DPF2 labelled proxy.",
        ],
    }


def main() -> None:
    payload = compute_rna()
    (OUT / "rna_phase3.json").write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
