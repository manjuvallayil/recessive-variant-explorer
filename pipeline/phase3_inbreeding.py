"""Phase 3c: inbreeding depression and the non-additive attenuation (OLS approximation).

Model
- A: bodyweight yd ~ 1 + F. Slope b_A is the inbreeding-depression slope.
- B: bodyweight yd ~ 1 + F + (het, hom dummies for the 8 non-additive loci).
- Attenuation = (b_A - b_B) / b_A: the share of the F slope absorbed by the non-additive loci.

This is an OLS approximation of the paper's BayesC0 (see docs/ASSUMPTIONS.md). It is not the
same estimand as the paper's variance fraction (23.5%), so the two are reported side by side.

Keys and units
- phen_1.1 and other_7 share the discovery ID space; genotypes use the gen_1 windows, which
  are the phen_1 key (docs/ASSUMPTIONS.md item 5 / item 32).
- F units are assumed to be fractions (the -0.77 kg per 1% comparison depends on this).
  Sensitivity: F > 1 rows (179 animals, max 3.1) are excluded in a second fit.
- Genomic cross-check of the scale was inconclusive: it used one 4,447-variant window only.
"""

from __future__ import annotations

import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
from bed_reader import open_bed

ROOT = Path(__file__).resolve().parents[1]
DRYAD = ROOT / "dryad"
PHENS = DRYAD / "phens"
OTHER = DRYAD / "other"
GENS = DRYAD / "gens"
OUT = ROOT / "web" / "public" / "data"

PAPER_INBREEDING_DEPRESSION_KG_PER_1PCT = -0.77
PAPER_NONADDITIVE_FRACTION = 0.235
BOOTSTRAP_B = 200
SEED = 0

# (label, gen_1 window, chrom, [bp...]) ; DPF2 and MUS81 are one locus (hom for either)
LOCI = [
    ("Chr2:22", "gen_1.1", 2, [22836659]),
    ("PLCD4", "gen_1.2", 2, [107313998]),
    ("FGD4", "gen_1.3", 5, [77632752]),
    ("ARSI", "gen_1.4", 7, [63649071]),
    ("MTRF1", "gen_1.5", 12, [11463981]),
    ("MYH1", "gen_1.7", 19, [30114249]),
    ("GALNT2", "gen_1.8", 28, [1312334]),
    ("DPF2/MUS81", "gen_1.9", 29, [44213160, 44645469]),
]


def code_version() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


def _variant_genotypes(stem: str, chrom: int, bp: int) -> pd.Series:
    bim = pd.read_csv(GENS / f"{stem}.bim", sep=r"\s+", header=None)
    hit = bim.index[(bim[0] == chrom) & (bim[3] == bp)]
    if len(hit) == 0:
        raise ValueError(f"variant {chrom}:{bp} not in {stem}")
    g = open_bed(GENS / f"{stem}.bed").read((slice(None), [int(hit[0])]), dtype="int8")[:, 0]
    iids = pd.read_csv(GENS / f"{stem}.fam", sep=r"\s+", header=None)[1].astype(int).to_numpy()
    keep = g >= 0
    return pd.Series(g[keep].astype(int), index=iids[keep])


def locus_class(stem: str, chrom: int, bps: list[int]) -> pd.Series:
    """Genotype class per IID: 2 if homozygous for the mutant at any listed variant, else 1 if
    heterozygous at any, else 0."""
    cls: pd.Series | None = None
    for bp in bps:
        g = _variant_genotypes(stem, chrom, bp)
        cls = g if cls is None else pd.concat([cls, g], axis=1).max(axis=1)
    assert cls is not None
    return cls.astype(int)


def _ols(y: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, float]:
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    resid = y - x @ beta
    r2 = 1.0 - float(resid @ resid) / float(((y - y.mean()) ** 2).sum())
    return beta, r2


def _design(frame: pd.DataFrame, loci_cols: list[str], hom_only: bool = False) -> np.ndarray:
    cols = [np.ones(len(frame)), frame["F"].to_numpy(dtype=float)]
    for c in loci_cols:
        if not hom_only:
            cols.append(frame[f"{c}_het"].to_numpy(dtype=float))
        cols.append(frame[f"{c}_hom"].to_numpy(dtype=float))
    return np.column_stack(cols)


def compute_inbreeding() -> dict:
    other7 = pd.read_csv(OTHER / "other_7.txt", sep=r"\s+", header=None, names=["ID", "F", "value"])
    phen1 = pd.read_csv(PHENS / "phen_1.1.txt", sep=r"\s+")
    frame = phen1.merge(other7[["ID", "F"]], on="ID", how="inner").dropna(subset=["yd_trait", "F"])
    frame = frame.set_index("ID")

    kept_loci: list[str] = []
    for label, stem, chrom, bps in LOCI:
        cls = locus_class(stem, chrom, bps)
        col = label.replace("/", "_")
        frame[f"{col}_het"] = (frame.index.map(cls) == 1).astype(float)
        frame[f"{col}_hom"] = (frame.index.map(cls) == 2).astype(float)
        frame = frame[frame.index.map(cls).notna()]
        if frame[f"{col}_hom"].sum() >= 1:
            kept_loci.append(col)
    frame = frame.reset_index()

    y = frame["yd_trait"].to_numpy(dtype=float)
    xa = _design(frame, [])
    xb = _design(frame, kept_loci)
    xh = _design(frame, kept_loci, hom_only=True)
    beta_a, r2_a = _ols(y, xa)
    beta_b, r2_b = _ols(y, xb)
    beta_h, r2_h = _ols(y, xh)
    att_hom_only = 1.0 - float(beta_h[1]) / float(beta_a[1])
    b_a = float(beta_a[1])
    b_b = float(beta_b[1])
    attenuation = 1.0 - b_b / b_a

    rng = np.random.default_rng(SEED)
    boot_att: list[float] = []
    boot_slope_a: list[float] = []
    n = len(y)
    for _ in range(BOOTSTRAP_B):
        idx = rng.integers(0, n, n)
        ba, _ = _ols(y[idx], xa[idx])
        bb, _ = _ols(y[idx], xb[idx])
        boot_slope_a.append(float(ba[1]))
        boot_att.append(1.0 - float(bb[1]) / float(ba[1]))

    # Sensitivity: drop F > 1 (not valid under a fraction definition)
    ok = frame["F"].to_numpy(dtype=float) <= 1.0
    beta_a_ok, _ = _ols(y[ok], xa[ok])
    beta_b_ok, _ = _ols(y[ok], xb[ok])
    att_ok = 1.0 - float(beta_b_ok[1]) / float(beta_a_ok[1])

    return {
        "schema_version": "1.0.0",
        "generated_at": datetime.now(UTC).isoformat(),
        "code_version": code_version(),
        "source_files": {
            "other_7": "dryad/other/other_7.txt",
            "phen_1_1": "dryad/phens/phen_1.1.txt",
            "gen_1_windows": "dryad/gens/gen_1.1 .. gen_1.9 (8 non-additive loci)",
        },
        "inbreeding": {
            "n_animals": int(n),
            "F_mean": float(frame["F"].mean()),
            "F_median": float(frame["F"].median()),
            "F_n_above_1": int((~ok).sum()),
            "slope_kg_per_unit_F": b_a,
            "slope_kg_per_1pct_F": b_a / 100.0,
            "slope_kg_per_1pct_F_ci": [
                float(np.percentile(boot_slope_a, 2.5)) / 100.0,
                float(np.percentile(boot_slope_a, 97.5)) / 100.0,
            ],
            "r2_F_only": r2_a,
        },
        "nonadditive_model": {
            "loci_in_model": kept_loci,
            "slope_kg_per_unit_F_with_loci": b_b,
            "slope_kg_per_1pct_F_with_loci": b_b / 100.0,
            "attenuation_fraction": attenuation,
            "attenuation_fraction_ci": [
                float(np.percentile(boot_att, 2.5)),
                float(np.percentile(boot_att, 97.5)),
            ],
            "attenuation_fraction_excluding_F_gt_1": att_ok,
            "r2_loci_and_F": r2_b,
            "r2_gain_from_loci": r2_b - r2_a,
            "attenuation_fraction_homozygote_only": att_hom_only,
            "slope_kg_per_1pct_F_homozygote_only": float(beta_h[1]) / 100.0,
            "r2_homozygote_only_and_F": r2_h,
            "method": "OLS approximation (not BayesC0); attenuation = 1 - b_with_loci / b_F_only",
        },
        "paper_comparison": {
            "paper_kg_per_1pct_inbreeding": PAPER_INBREEDING_DEPRESSION_KG_PER_1PCT,
            "paper_nonadditive_fraction": PAPER_NONADDITIVE_FRACTION,
            "paper_nonadditive_fraction_is_variance": True,
        },
        "notes": [
            "Phenotype: raw yield deviation (yd_trait), not LOSO-adjusted.",
            "F assumed to be a fraction; the -0.77 kg per 1% comparison depends on this.",
            "Attenuation is a slope comparison, not the paper's variance fraction; "
            "the two numbers are not directly comparable.",
            "Genomic scale cross-check inconclusive (single 4,447-variant window).",
            "Loci: 8 non-additive loci; DPF2 and MUS81 combined (hom for either).",
            "PRIMARY recessive estimate = homozygote-only attenuation (about 2%).",
            "The het+hom attenuation (about 24%) is NOT evidence of a recessive effect: het "
            "dummies on raw yd absorb population structure (PLCD4 het is +24 kg raw vs -1 kg "
            "LOSO-adjusted). Do not present it as reproducing the paper's 23.5%.",
            "Paper 23.5% is a variance fraction from BayesC0; not reproduced here.",
        ],
    }


def main() -> None:
    payload = compute_inbreeding()
    (OUT / "inbreeding.json").write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
