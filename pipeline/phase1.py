"""Phase 1: locus-level reproduction and population genetics.

Reads the Dryad files under ``dryad/`` and writes small aggregate JSON files
to ``web/public/data/``.  No animal or sire IDs are exported.

Estimators
----------
A          : plain OLS/WLS on ``yd_trait`` (unweighted primary)
A_rinverse : WLS with ``weights=rinverse`` (labelled sensitivity)
B          : OLS on the LOSO-adjusted phenotype (paper-like)

The paper-like estimator B uses the matching ``other_2`` window for each
locus.  See ``docs/ASSUMPTIONS.md`` for the window choices and caveats.
"""

from __future__ import annotations

import json
import math
import subprocess
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
from bed_reader import open_bed
from scipy import stats
from statsmodels.regression.linear_model import OLS, WLS

ROOT = Path(__file__).resolve().parents[1]
DRYAD = ROOT / "dryad"
GENS = DRYAD / "gens"
PHENS = DRYAD / "phens"
OTHER = DRYAD / "other"
OUT = ROOT / "web" / "public" / "data"
OUT.mkdir(parents=True, exist_ok=True)

SCHEMA_VERSION = "1.0.0"
PAPER_REFERENCE_YEAR = 2021
PAPER_SCALED_HOMOZYGOUS = 9700
PAPER_HOMOZYGOUS_RATE = 0.01
FEMALES_BORN_PER_YEAR = 940_100

TRAITS = {
    "bodyweight": {
        "phen": "phen_1.1.txt",
        "other1": 1,
        "other2": 1,
        "unit": "kg",
        "raw_unit": "kg",
        "scale_to_unit": 1.0,
    },
    "stature": {
        "phen": "phen_1.2.txt",
        "other1": 2,
        "other2": 3,
        "unit": "cm",
        "raw_unit": "TOP_score_1_9",
        "scale_to_unit": 5.0,
    },
    "bcs": {
        "phen": "phen_1.3.txt",
        "other1": 3,
        "other2": 2,
        "unit": "1-10",
        "raw_unit": "1-10",
        "scale_to_unit": 1.0,
    },
}

# Table 1 values parsed from the paper PDF: (het, hom-alt).
TABLE1 = {
    "bodyweight": {
        "Chr2:22": (-0.007, -11.8),
        "PLCD4": (-0.354, -108.965),
        "FGD4": (0.834, -48.079),
        "ARSI": (0.564, -11.981),
        "MTRF1": (1.078, -19.464),
        "GALNT2": (-0.837, -112.008),
        "DPF2": (4.059, -53.092),
        "MUS81": (4.126, -51.911),
        "MYH1": (0.798, -4.643),
    },
    "stature": {
        "Chr2:22": (0.025, -0.428),
        "PLCD4": (-0.025, -5.421),
        "FGD4": (0.058, -3.041),
        "ARSI": (0.055, -1.152),
        "MTRF1": (0.025, -1.208),
        "GALNT2": (-0.134, -8.568),
        "DPF2": (0.132, -2.565),
        "MUS81": (0.136, -2.487),
        "MYH1": (0.086, 0.285),
    },
    "bcs": {
        "Chr2:22": (0.003, -0.044),
        "PLCD4": (0.005, -0.284),
        "FGD4": (0.012, -0.238),
        "ARSI": (-0.001, -0.028),
        "MTRF1": (0.009, -0.093),
        "GALNT2": (0.012, -0.079),
        "DPF2": (0.02, -0.05),
        "MUS81": (0.02, -0.054),
        "MYH1": (-0.01, -0.133),
    },
}

LOCI = [
    {
        "locus": "Chr2:22",
        "gene": None,
        "chr": 2,
        "bp": 22836659,
        "candidate": "rs211074296",
        "consequence": "QTL (no candidate)",
        "gen1": "gen_1.1",
        "gen2": "gen_2.1",
        "gen3": "gen_3.1",
        "window": 1,
        "window_sensitivity": 2,
        "breed": None,
        "gen4": None,
        "gen5": None,
    },
    {
        "locus": "PLCD4",
        "gene": "PLCD4",
        "chr": 2,
        "bp": 107313998,
        "candidate": "g.107313998G>A",
        "consequence": "p.Ala326Thr",
        "gen1": "gen_1.2",
        "gen2": "gen_2.2",
        "gen3": "gen_3.2",
        "window": 3,
        "window_sensitivity": 4,
        "breed": "HF",
        "gen4": "gen_4.1",
        "gen5": "gen_5.1",
    },
    {
        "locus": "FGD4",
        "gene": "FGD4",
        "chr": 5,
        "bp": 77632752,
        "candidate": "g.77632752C>T",
        "consequence": "splice donor",
        "gen1": "gen_1.3",
        "gen2": "gen_2.3",
        "gen3": "gen_3.3",
        "window": 5,
        "window_sensitivity": 6,
        "breed": "HF",
        "gen4": "gen_4.2",
        "gen5": "gen_5.2",
    },
    {
        "locus": "ARSI",
        "gene": "ARSI",
        "chr": 7,
        "bp": 63649071,
        "candidate": "g.63649071C>T",
        "consequence": "p.Gly301Ser",
        "gen1": "gen_1.4",
        "gen2": "gen_2.4",
        "gen3": "gen_3.4",
        "window": 7,
        "window_sensitivity": 8,
        "breed": None,
        "gen4": None,
        "gen5": None,
    },
    {
        "locus": "MTRF1",
        "gene": "MTRF1",
        "chr": 12,
        "bp": 11463981,
        "candidate": "g.11463981G>A",
        "consequence": "p.Arg341Trp",
        "gen1": "gen_1.5",
        "gen2": "gen_2.5",
        "gen3": "gen_3.5",
        "window": 9,
        "window_sensitivity": None,
        "breed": "J",
        "gen4": "gen_4.3",
        "gen5": "gen_5.3",
    },
    {
        "locus": "MYH1",
        "gene": "MYH1",
        "chr": 19,
        "bp": 30114249,
        "candidate": "g.30114250_30114327delins (left-flank proxy)",
        "consequence": "frameshift",
        "gen1": "gen_1.7",
        "gen2": "gen_2.9",
        "gen3": "gen_3.9",
        "window": 11,
        "window_sensitivity": None,
        "breed": "J",
        "gen4": "gen_4.4",
        "gen5": "gen_5.4",
    },
    {
        "locus": "GALNT2",
        "gene": "GALNT2",
        "chr": 28,
        "bp": 1312334,
        "candidate": "g.1312334G>A",
        "consequence": "splice acceptor",
        "gen1": "gen_1.8",
        "gen2": "gen_2.10",
        "gen3": "gen_3.10",
        "window": 12,
        "window_sensitivity": None,
        "breed": "HF",
        "gen4": "gen_4.5",
        "gen5": "gen_5.5",
    },
    {
        "locus": "DPF2",
        "gene": "DPF2",
        "chr": 29,
        "bp": 44213160,
        "candidate": "g.44213160A>G",
        "consequence": "p.Lys216Arg",
        "gen1": "gen_1.9",
        "gen2": "gen_2.11",
        "gen3": "gen_3.11",
        "window": 13,
        "window_sensitivity": None,
        "breed": "HF",
        "gen4": "gen_4.6",
        "gen5": "gen_5.6",
    },
    {
        "locus": "MUS81",
        "gene": "MUS81",
        "chr": 29,
        "bp": 44645469,
        "candidate": "g.44645469G>T",
        "consequence": "p.Gly70*",
        "gen1": "gen_1.9",
        "gen2": "gen_2.12",
        "gen3": "gen_3.12",
        "window": 13,
        "window_sensitivity": None,
        "breed": "HF",
        "gen4": "gen_4.7",
        "gen5": "gen_5.7",
    },
]

COMPARATORS = [
    {
        "locus": "chr14:1802265",
        "gene": "DGAT1",
        "chr": 14,
        "bp": 1802265,
        "candidate": "rs109234250",
        "consequence": "DGAT1 p.Lys232Ala additive comparator",
        "gen1": None,
        "gen2": "gen_2.6",
        "gen3": "gen_3.6",
        "window": None,
        "window_sensitivity": None,
        "breed": None,
        "gen4": None,
        "gen5": None,
    },
    {
        "locus": "chr14:25015640",
        "gene": "PLAG1",
        "chr": 14,
        "bp": 25015640,
        "candidate": "rs109815800",
        "consequence": "PLAG1 causal regulatory comparator",
        "gen1": None,
        "gen2": "gen_2.7",
        "gen3": "gen_3.7",
        "window": 10,
        "window_sensitivity": None,
        "breed": None,
        "gen4": None,
        "gen5": None,
    },
    {
        "locus": "chr14:26064929",
        "gene": "PLAG1",
        "chr": 14,
        "bp": 26064929,
        "candidate": "rs109686357",
        "consequence": "PLAG1 dominance-GWAS lead comparator",
        "gen1": "gen_1.6",
        "gen2": "gen_2.8",
        "gen3": "gen_3.8",
        "window": 10,
        "window_sensitivity": None,
        "breed": None,
        "gen4": None,
        "gen5": None,
    },
]

PLEIO_TRAITS = [
    ("bodyweight", PHENS / "phen_1.1.txt", "growth"),
    ("stature", PHENS / "phen_1.2.txt", "growth"),
    ("bcs", PHENS / "phen_1.3.txt", "growth"),
    ("milk", PHENS / "phen_4.1.txt", "production"),
    ("protein", PHENS / "phen_4.2.txt", "production"),
    ("fat", PHENS / "phen_4.3.txt", "production"),
    ("adaptability", PHENS / "phen_4.4.txt", "conformation"),
    ("temperament", PHENS / "phen_4.5.txt", "conformation"),
    ("farmer_opinion", PHENS / "phen_4.6.txt", "conformation"),
    ("capacity", PHENS / "phen_4.7.txt", "conformation"),
    ("rump_angle", PHENS / "phen_4.8.txt", "conformation"),
    ("rump_width", PHENS / "phen_4.9.txt", "conformation"),
    ("legs", PHENS / "phen_4.10.txt", "conformation"),
    ("udder_support", PHENS / "phen_4.11.txt", "conformation"),
    ("front_udder", PHENS / "phen_4.12.txt", "conformation"),
    ("front_teat", PHENS / "phen_4.13.txt", "conformation"),
    ("rear_teat", PHENS / "phen_4.14.txt", "conformation"),
    ("udder_overall", PHENS / "phen_4.15.txt", "conformation"),
    ("conformation", PHENS / "phen_4.16.txt", "conformation"),
]

PHEN_CACHE: dict[str, pd.DataFrame] = {}
LOSO_CACHE: dict[tuple[int, int], pd.DataFrame] = {}
OTHER1_CACHE: dict[int, pd.DataFrame] = {}
GEN1_IIDS = pd.read_csv(GENS / "gen_1.1.fam", sep=r"\s+", header=None)[1].astype(int).to_numpy()
GEN3_IIDS = pd.read_csv(GENS / "gen_3.1.fam", sep=r"\s+", header=None)[1].astype(int).to_numpy()


def now_iso() -> str:
    return datetime.now(UTC).isoformat()


def code_version() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


CODE_VERSION = code_version()


def jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        value = float(obj)
        return value if math.isfinite(value) else None
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if pd.isna(obj):
        return None
    return obj


def write_json(name: str, payload: dict, source_files: list[str]) -> dict:
    body = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": now_iso(),
        "code_version": CODE_VERSION,
        "source_files": source_files,
        **payload,
    }
    body = jsonable(body)
    path = OUT / name
    path.write_text(json.dumps(body, indent=2) + "\n")
    return {"file": name, "bytes": path.stat().st_size, "source_files": source_files}


def load_fam_iids(path: Path) -> np.ndarray:
    return pd.read_csv(path, sep=r"\s+", header=None)[1].astype(int).to_numpy()


def load_phen(path: Path) -> pd.DataFrame:
    key = str(path)
    if key not in PHEN_CACHE:
        PHEN_CACHE[key] = pd.read_csv(path, sep=r"\s+")
    return PHEN_CACHE[key]


def read_bed_variants(path: Path, indices: list[int]) -> np.ndarray:
    bed = open_bed(path)
    arr = bed.read((slice(None), indices), dtype="int8")
    return arr.reshape((arr.shape[0], len(indices)))


def find_bim_index(bim: pd.DataFrame, chrom: int, bp: int) -> int:
    hit = bim.index[(bim[0] == chrom) & (bim[3] == bp)]
    if len(hit) == 0:
        raise ValueError(f"variant {chrom}:{bp} not found")
    return int(hit[0])


def genotype_for_locus(locus: dict, source: str) -> tuple[np.ndarray, np.ndarray]:
    if source == "gen1":
        stem = locus["gen1"]
        if stem is None:
            raise ValueError(f"{locus['locus']} has no gen_1 file")
        path = GENS / f"{stem}.bed"
        bim = pd.read_csv(GENS / f"{stem}.bim", sep=r"\s+", header=None)
        idx = find_bim_index(bim, locus["chr"], locus["bp"])
        g = read_bed_variants(path, [idx])[:, 0]
        iids = load_fam_iids(GENS / f"{stem}.fam")
        return iids, g
    if source == "gen2":
        stem = locus["gen2"]
        path = GENS / f"{stem}.bed"
        g = read_bed_variants(path, [0])[:, 0]
        iids = load_fam_iids(GENS / f"{stem}.fam")
        return iids, g
    if source == "gen3":
        stem = locus["gen3"]
        path = GENS / f"{stem}.bed"
        g = read_bed_variants(path, [0])[:, 0]
        return GEN3_IIDS, g
    if source == "gen4":
        stem = locus["gen4"]
        path = GENS / f"{stem}.bed"
        g = read_bed_variants(path, [0])[:, 0]
        iids = load_fam_iids(GENS / f"{stem}.fam")
        return iids, g
    if source == "gen5":
        stem = locus["gen5"]
        path = GENS / f"{stem}.bed"
        g = read_bed_variants(path, [0])[:, 0]
        iids = load_fam_iids(GENS / f"{stem}.fam")
        return iids, g
    raise ValueError(f"unknown source {source}")


def load_loso_means(trait_idx: int, window: int) -> pd.DataFrame:
    key = (trait_idx, window)
    if key not in LOSO_CACHE:
        path = OTHER / f"other_2.{trait_idx}.{window}.tsv"
        d = pd.read_csv(path, sep="\t", header=None, dtype=np.float32)
        ids = d.iloc[:, 0].to_numpy(dtype=np.int64)
        means = d.iloc[:, 1:].mean(axis=1).to_numpy(dtype=np.float64)
        LOSO_CACHE[key] = pd.DataFrame({"ID": ids, "loso": means})
        del d
    return LOSO_CACHE[key]


def phenotype_with_adjustment(
    trait: str, locus: dict, window_override: int | None = None
) -> pd.DataFrame | None:
    phen = load_phen(PHENS / TRAITS[trait]["phen"]).copy()
    win = window_override if window_override is not None else locus.get("window")
    if win is None:
        return None
    loso = load_loso_means(TRAITS[trait]["other2"], int(win))
    merged = phen.merge(loso, on="ID", how="inner")
    merged["adjusted"] = merged["yd_trait"] - merged["loso"]
    merged["adjusted"] = merged["adjusted"] - merged["adjusted"].mean()
    return merged


def class_stats(y: np.ndarray, g: np.ndarray) -> dict:
    out = {}
    for k in (0, 1, 2):
        yk = y[g == k]
        n = int(len(yk))
        if n == 0:
            out[str(k)] = {
                "n": 0,
                "mean": None,
                "sd": None,
                "se": None,
                "quantiles": {q: None for q in ("p05", "p25", "p50", "p75", "p95")},
            }
            continue
        sd = float(np.std(yk, ddof=1)) if n > 1 else None
        q = np.quantile(yk, [0.05, 0.25, 0.5, 0.75, 0.95]) if n > 1 else [np.nan] * 5
        small = n < 5
        out[str(k)] = {
            "n": n,
            "mean": None if small else float(np.mean(yk)),
            "sd": None if small or sd is None else sd,
            "se": None if small or sd is None else float(sd / math.sqrt(n)),
            "quantiles": {
                "p05": None if small else float(q[0]),
                "p25": None if small else float(q[1]),
                "p50": None if small else float(q[2]),
                "p75": None if small else float(q[3]),
                "p95": None if small else float(q[4]),
            },
        }
    return out


def fit_effects(
    df: pd.DataFrame,
    y_col: str,
    g_col: str,
    weights_col: str | None = None,
) -> dict | None:
    y = df[y_col].to_numpy(dtype=float)
    g = df[g_col].to_numpy(dtype=float)
    mask = np.isfinite(y) & np.isfinite(g)
    y = y[mask]
    g = g[mask].astype(int)
    if len(y) < 10 or len(np.unique(g)) < 2:
        return None
    x = np.column_stack([np.ones(len(g)), g, (g == 1).astype(float)])
    try:
        if weights_col is None:
            res = OLS(y, x).fit()
        else:
            w = df[weights_col].to_numpy(dtype=float)[mask]
            w = np.where(np.isfinite(w) & (w > 0), w, 1.0)
            res = WLS(y, x, weights=w).fit()
    except Exception:
        return None
    b = np.asarray(res.params, dtype=float)
    cov = np.asarray(res.cov_params(), dtype=float)
    a, d = float(b[1]), float(b[2])
    het, hom = a + d, 2.0 * a
    se_a = float(math.sqrt(max(cov[1, 1], 0.0)))
    se_d = float(math.sqrt(max(cov[2, 2], 0.0)))
    se_het = float(math.sqrt(max(cov[1, 1] + cov[2, 2] + 2 * cov[1, 2], 0.0)))
    se_hom = float(2.0 * se_a)
    df_resid = int(res.df_resid)
    dist = stats.t if df_resid > 0 else stats.norm

    def pval(z: float) -> float:
        return float(2 * dist.sf(abs(z), df_resid)) if df_resid > 0 else float(2 * dist.sf(abs(z)))

    z_a = a / se_a if se_a > 0 else np.nan
    z_d = d / se_d if se_d > 0 else np.nan
    z_het = het / se_het if se_het > 0 else np.nan
    z_hom = hom / se_hom if se_hom > 0 else np.nan
    cstats = class_stats(y, g)
    small_het = int(cstats["1"]["n"]) < 5
    small_hom = int(cstats["2"]["n"]) < 5
    return {
        "n_total": int(len(y)),
        "class_stats": cstats,
        "a": None if small_hom else a,
        "d": None if (small_het or small_hom) else d,
        "het": None if small_het else het,
        "hom_alt": None if small_hom else hom,
        "se_a": None if small_hom else se_a,
        "se_d": None if (small_het or small_hom) else se_d,
        "se_het": None if small_het else se_het,
        "se_hom_alt": None if small_hom else se_hom,
        "z_a": None if small_hom else float(z_a),
        "z_d": None if (small_het or small_hom) else float(z_d),
        "z_het": None if small_het else float(z_het),
        "z_hom_alt": None if small_hom else float(z_hom),
        "p_a": None if small_hom else pval(float(z_a)),
        "p_d": None if (small_het or small_hom) else pval(float(z_d)),
        "p_het": None if small_het else pval(float(z_het)),
        "p_hom_alt": None if small_hom else pval(float(z_hom)),
        "df_resid": df_resid,
        "small_n": any(int(v["n"]) < 5 for v in cstats.values()),
    }


def fit_locus_trait(
    locus: dict, trait: str, estimator: str, window_override: int | None = None
) -> dict | None:
    iids, g = genotype_for_locus(locus, "gen1")
    gt = pd.DataFrame({"ID": iids, "g": g})
    scale = float(TRAITS[trait]["scale_to_unit"])
    if estimator == "B":
        phen = phenotype_with_adjustment(trait, locus, window_override=window_override)
        if phen is None:
            return None
        merged = phen.merge(gt, on="ID", how="inner")
        merged["y_scaled"] = merged["adjusted"] * scale
        return fit_effects(merged, "y_scaled", "g")
    phen = load_phen(PHENS / TRAITS[trait]["phen"]).copy()
    merged = phen.merge(gt, on="ID", how="inner")
    merged["y_scaled"] = merged["yd_trait"] * scale
    weights = "rinverse" if estimator == "A_rinverse" else None
    return fit_effects(merged, "y_scaled", "g", weights_col=weights)


def load_other1(trait_idx: int) -> pd.DataFrame:
    if trait_idx not in OTHER1_CACHE:
        OTHER1_CACHE[trait_idx] = pd.read_csv(
            OTHER / f"other_1.{trait_idx}.txt", sep=r"\s+"
        )
    return OTHER1_CACHE[trait_idx]


def published_other1_row(trait: str, locus: dict) -> dict | None:
    trait_idx = int(TRAITS[trait]["other1"])
    scale = float(TRAITS[trait]["scale_to_unit"])
    d = load_other1(trait_idx)
    hit = d[(d.CHR == locus["chr"]) & (d.BP == locus["bp"])]
    if hit.empty:
        return None
    row = hit.iloc[0]
    b_a = float(row.B_A)
    b_d = float(row.B_D)
    z_a = float(row.Z_A)
    z_d = float(row.Z_D)
    se_a = abs(b_a / z_a) if z_a != 0 else None
    se_d = abs(b_d / z_d) if z_d != 0 else None
    derived_hom = -2.0 * b_a * scale
    derived_het = (-b_a + b_d) * scale
    se_hom = 2.0 * se_a * scale if se_a is not None else None
    se_het = (
        math.sqrt(se_a**2 + se_d**2) * scale
        if se_a is not None and se_d is not None
        else None
    )
    return {
        "A1": str(row.A1),
        "A2": str(row.A2),
        "MAF": float(row.MAF),
        "B_A": b_a,
        "B_D": b_d,
        "Z_A": z_a,
        "Z_D": z_d,
        "P_A": float(row.P_A),
        "P_D": float(row.P_D),
        "scale_to_unit": scale,
        "unit": TRAITS[trait]["unit"],
        "raw_unit": TRAITS[trait]["raw_unit"],
        "derived_hom_alt_raw_unit": -2.0 * b_a,
        "derived_het_raw_unit": -b_a + b_d,
        "derived_hom_alt": derived_hom,
        "derived_het": derived_het,
        "derived_se_hom_alt": se_hom,
        "derived_se_het": se_het,
        "derivation": (
            "B_A treated as A2/reference effect: hom-alt=-2*B_A*scale; "
            "het=(-B_A+B_D)*scale (stature TOP score scaled by 5 cm/score)"
        ),
    }


def compute_extfig2c_comparison(effect_rows: list[dict], class_rows: list[dict]) -> dict:
    path = DRYAD / "source_data" / "source_data.ExtFig2c.csv"
    df = pd.read_csv(path)
    col_map = [
        ("Chr2:22", "Bodyweight.CHR2_22MB", "Genotype.CHR2_22MB", 1),
        ("PLCD4", "Bodyweight", "Genotype", 3),
        ("FGD4", "Bodyweight.FGD4", "Genotype.FGD4", 6),
        ("ARSI", "Bodyweight.ARSI", "Genotype.ARSI", 7),
        ("MTRF1", "Bodyweight.MTRF1", "Genotype.MTRF1", 9),
        ("GALNT2", "Bodyweight.GALNT2", "Genotype.GALNT2", 12),
        ("DPF2", "Bodyweight.DPF2", "Genotype.DPF2", 13),
        ("chr14:25015640", "Bodyweight.PLAG1", "Genotype.PLAG1", 10),
    ]
    by_b = {
        r["locus"]: r
        for r in effect_rows
        if r["trait"] == "bodyweight" and r["estimator"] == "B"
    }
    plag1_loc = next(c for c in COMPARATORS if c["locus"] == "chr14:25015640")
    _, g_plag1 = genotype_for_locus(plag1_loc, "gen3")
    phen_plag1 = phenotype_with_adjustment("bodyweight", plag1_loc, window_override=10)
    if phen_plag1 is not None:
        m_plag1 = phen_plag1.merge(
            pd.DataFrame({"ID": GEN1_IIDS, "g": g_plag1}), on="ID", how="inner"
        )
        fit_plag1 = fit_effects(m_plag1, "adjusted", "g")
        if fit_plag1 is not None:
            by_b["chr14:25015640"] = fit_plag1
    rows = []
    for locus_name, bw_col, gt_col, matching_window in col_map:
        g = df[gt_col].map({"0/0": 0, "0/1": 1, "1/1": 2}).to_numpy(dtype=int)
        y = df[bw_col].to_numpy(dtype=float)
        cstats = class_stats(y, g)
        het = float(cstats["1"]["mean"] - cstats["0"]["mean"])
        hom_alt = float(cstats["2"]["mean"] - cstats["0"]["mean"])
        b_row = by_b.get(locus_name)
        rows.append(
            {
                "locus": locus_name,
                "extfig2c_bw_column": bw_col,
                "matching_loso_window": matching_window,
                "n_total": int(len(y)),
                "class_stats": cstats,
                "het": het,
                "hom_alt": hom_alt,
                "estimator_B_primary_het": b_row["het"] if b_row else None,
                "estimator_B_primary_hom_alt": b_row["hom_alt"] if b_row else None,
                "delta_vs_B_primary_het": (het - b_row["het"]) if b_row else None,
                "delta_vs_B_primary_hom_alt": (hom_alt - b_row["hom_alt"]) if b_row else None,
            }
        )
    return {
        "rows": rows,
        "note": (
            "Aggregate comparison only (Gate A decision 4): ExtFig2c IDs are renumbered "
            "1..79,945 and are never row-joined. Class counts, means, quantiles, and "
            "effects match Estimator B exactly on the corresponding LOSO window "
            "(window 6 for FGD4; primary window for all other loci)."
        ),
    }


def compute_table1() -> tuple[list[dict], dict, list[dict], list[dict], list[dict]]:
    class_rows: list[dict] = []
    effect_rows: list[dict] = []
    published_rows: list[dict] = []
    comparison_rows: list[dict] = []
    window_sensitivity_rows: list[dict] = []
    for trait in TRAITS:
        scale = float(TRAITS[trait]["scale_to_unit"])
        for locus in LOCI:
            paper_het, paper_hom = TABLE1[trait][locus["locus"]]
            for estimator in ("A", "A_rinverse", "B"):
                fit = fit_locus_trait(locus, trait, estimator)
                if fit is None:
                    continue
                class_rows.append(
                    {
                        "trait": trait,
                        "unit": TRAITS[trait]["unit"],
                        "raw_unit": TRAITS[trait]["raw_unit"],
                        "scale_to_unit": scale,
                        "locus": locus["locus"],
                        "estimator": estimator,
                        "n_total": fit["n_total"],
                        "class_stats": fit["class_stats"],
                    }
                )
                effect_rows.append(
                    {
                        "trait": trait,
                        "unit": TRAITS[trait]["unit"],
                        "raw_unit": TRAITS[trait]["raw_unit"],
                        "scale_to_unit": scale,
                        "locus": locus["locus"],
                        "estimator": estimator,
                        "window": locus["window"] if estimator == "B" else None,
                        "n_total": fit["n_total"],
                        "a": fit["a"],
                        "d": fit["d"],
                        "het": fit["het"],
                        "hom_alt": fit["hom_alt"],
                        "het_raw_unit": fit["het"] / scale if fit["het"] is not None else None,
                        "hom_alt_raw_unit": (
                            fit["hom_alt"] / scale if fit["hom_alt"] is not None else None
                        ),
                        "se_het": fit["se_het"],
                        "se_hom_alt": fit["se_hom_alt"],
                        "z_het": fit["z_het"],
                        "z_hom_alt": fit["z_hom_alt"],
                        "p_het": fit["p_het"],
                        "p_hom_alt": fit["p_hom_alt"],
                        "small_n": fit["small_n"],
                    }
                )
                comparison_rows.append(
                    {
                        "trait": trait,
                        "locus": locus["locus"],
                        "estimator": estimator,
                        "paper_het": paper_het,
                        "paper_hom_alt": paper_hom,
                        "het": fit["het"],
                        "hom_alt": fit["hom_alt"],
                        "delta_het": fit["het"] - paper_het,
                        "delta_hom_alt": fit["hom_alt"] - paper_hom,
                        "rel_delta_hom_alt": (
                            (fit["hom_alt"] - paper_hom) / paper_hom if paper_hom != 0 else None
                        ),
                    }
                )
            if locus.get("window_sensitivity") is not None:
                fit_sens = fit_locus_trait(
                    locus, trait, "B", window_override=int(locus["window_sensitivity"])
                )
                if fit_sens is not None:
                    window_sensitivity_rows.append(
                        {
                            "trait": trait,
                            "unit": TRAITS[trait]["unit"],
                            "locus": locus["locus"],
                            "window": int(locus["window_sensitivity"]),
                            "het": fit_sens["het"],
                            "hom_alt": fit_sens["hom_alt"],
                            "paper_het": paper_het,
                            "paper_hom_alt": paper_hom,
                            "delta_het": fit_sens["het"] - paper_het,
                            "delta_hom_alt": fit_sens["hom_alt"] - paper_hom,
                        }
                    )
            pub = published_other1_row(trait, locus)
            if pub is not None:
                published_rows.append(
                    {
                        "trait": trait,
                        "locus": locus["locus"],
                        **pub,
                        "paper_het": paper_het,
                        "paper_hom_alt": paper_hom,
                        "delta_derived_hom_alt": pub["derived_hom_alt"] - paper_hom,
                        "delta_derived_het": pub["derived_het"] - paper_het,
                    }
                )
    summary = {}
    for estimator in ("A", "A_rinverse", "B"):
        rows = [r for r in comparison_rows if r["estimator"] == estimator]
        bw_rows = [r for r in rows if r["trait"] == "bodyweight"]
        hom = [abs(r["delta_hom_alt"]) for r in rows]
        het = [abs(r["delta_het"]) for r in rows]
        bw_hom = [abs(r["delta_hom_alt"]) for r in bw_rows]
        bw_het = [abs(r["delta_het"]) for r in bw_rows]
        large = [
            abs(r["delta_hom_alt"])
            for r in bw_rows
            if abs(r["paper_hom_alt"]) >= 20
        ]
        summary[estimator] = {
            "n_loci_all_traits": len(rows),
            "mae_hom_alt": float(np.mean(hom)) if hom else None,
            "mae_het": float(np.mean(het)) if het else None,
            "mae_bodyweight_hom_alt": float(np.mean(bw_hom)) if bw_hom else None,
            "mae_bodyweight_het": float(np.mean(bw_het)) if bw_het else None,
            "mae_hom_alt_large_effects": float(np.mean(large)) if large else None,
            "max_abs_delta_hom_alt": float(max(hom)) if hom else None,
        }
    closest = min(
        summary,
        key=lambda k: (
            summary[k]["mae_hom_alt"] if summary[k]["mae_hom_alt"] is not None else float("inf")
        ),
    )
    summary["closest_by_hom_alt_mae"] = closest
    return class_rows, summary, effect_rows, published_rows, window_sensitivity_rows


def compute_validation(pleio_rows: list[dict]) -> dict:
    phen_traits = [("milk", "phen_2.1.txt"), ("protein", "phen_2.2.txt"), ("fat", "phen_2.3.txt")]
    discovery = {
        (r["locus"], r["trait"]): r
        for r in pleio_rows
        if r["trait"] in {"milk", "protein", "fat"}
    }
    rows = []
    for locus in LOCI:
        iids, g = genotype_for_locus(locus, "gen2")
        gt = pd.DataFrame({"ID": iids, "g": g})
        for trait, fname in phen_traits:
            phen = load_phen(PHENS / fname).copy()
            merged = phen.merge(gt, on="ID", how="inner")
            fit = fit_effects(merged, "yd_trait", "g")
            if fit is None:
                continue
            disc = discovery.get((locus["locus"], trait))
            rows.append(
                {
                    "locus": locus["locus"],
                    "trait": trait,
                    "n_total": fit["n_total"],
                    "class_stats": fit["class_stats"],
                    "het": fit["het"],
                    "hom_alt": fit["hom_alt"],
                    "se_het": fit["se_het"],
                    "se_hom_alt": fit["se_hom_alt"],
                    "p_het": fit["p_het"],
                    "p_hom_alt": fit["p_hom_alt"],
                    "small_n": fit["small_n"],
                    "discovery_het": disc["het"] if disc else None,
                    "discovery_hom_alt": disc["hom_alt"] if disc else None,
                    "same_sign_het": (
                        bool(np.sign(fit["het"]) == np.sign(disc["het"])) if disc else None
                    ),
                    "same_sign_hom_alt": (
                        bool(np.sign(fit["hom_alt"]) == np.sign(disc["hom_alt"]))
                        if disc
                        else None
                    ),
                }
            )
    comparable = [r for r in rows if r["same_sign_hom_alt"] is not None]
    return {
        "rows": rows,
        "summary": {
            "n_comparisons": len(comparable),
            "same_sign_hom_alt": int(sum(bool(r["same_sign_hom_alt"]) for r in comparable)),
            "same_sign_het": int(
                sum(bool(r["same_sign_het"]) for r in comparable if r["same_sign_het"] is not None)
            ),
            "note": "Discovery comparator is phen_4.1-4.3 (milk/protein/fat) at the same locus.",
        },
    }


def compute_pleiotropy() -> list[dict]:
    rows = []
    all_loci = LOCI + COMPARATORS
    for locus in all_loci:
        _iids, g = genotype_for_locus(locus, "gen3")
        for trait, path, group in PLEIO_TRAITS:
            iids = GEN1_IIDS if path.name.startswith("phen_1.") else GEN3_IIDS
            gt = pd.DataFrame({"ID": iids, "g": g})
            phen = load_phen(path).copy()
            merged = phen.merge(gt, on="ID", how="inner")
            fit = fit_effects(merged, "yd_trait", "g")
            if fit is None:
                continue
            sd = float(np.std(merged["yd_trait"].dropna(), ddof=1))
            rows.append(
                {
                    "locus": locus["locus"],
                    "gene": locus["gene"],
                    "trait": trait,
                    "trait_group": group,
                    "n_total": fit["n_total"],
                    "class_n": {k: v["n"] for k, v in fit["class_stats"].items()},
                    "het": fit["het"],
                    "hom_alt": fit["hom_alt"],
                    "het_sd": fit["het"] / sd if sd > 0 else None,
                    "hom_alt_sd": fit["hom_alt"] / sd if sd > 0 else None,
                    "se_het": fit["se_het"],
                    "se_hom_alt": fit["se_hom_alt"],
                    "p_het": fit["p_het"],
                    "p_hom_alt": fit["p_hom_alt"],
                    "small_n": fit["small_n"],
                }
            )
    return rows


def compute_frequencies() -> dict:
    other6 = pd.read_csv(OTHER / "other_6.txt", sep=r"\s+")
    geno_by_locus = {}
    out = {}
    for locus in LOCI:
        if locus["gen4"] is None:
            continue
        iids, g = genotype_for_locus(locus, "gen4")
        df = pd.DataFrame({"ID": iids, "g": g}).merge(other6, on="ID", how="inner")
        geno_by_locus[locus["locus"]] = df.set_index("ID")["g"]
        out[locus["locus"]] = {"gene": locus["gene"], "breeds": {}}
        for breed, col in (("HF", "breed_hf"), ("J", "breed_j")):
            sub = df[df[col] == 16]
            n = int(len(sub))
            counts = {
                str(k): int(v)
                for k, v in sub.g.value_counts().reindex([0, 1, 2], fill_value=0).items()
            }
            hom_ref, het, hom_alt = counts["0"], counts["1"], counts["2"]
            p = (2 * hom_alt + het) / (2 * n) if n else None
            carrier = (het + hom_alt) / n if n else None
            expected = (
                {
                    "hom_ref": n * (1 - p) ** 2,
                    "het": 2 * n * p * (1 - p),
                    "hom_alt": n * p**2,
                }
                if p is not None
                else None
            )
            chi2 = None
            p_hwe = None
            if expected is not None:
                obs = np.array([hom_ref, het, hom_alt], dtype=float)
                exp = np.array(
                    [expected["hom_ref"], expected["het"], expected["hom_alt"]], dtype=float
                )
                valid = exp > 0
                chi2 = float(np.sum((obs[valid] - exp[valid]) ** 2 / exp[valid]))
                p_hwe = float(stats.chi2.sf(chi2, 1))
            out[locus["locus"]]["breeds"][breed] = {
                "n": n,
                "genotype_counts": counts,
                "alt_freq": p,
                "carrier_freq": carrier,
                "expected_hwe": expected,
                "hwe_chi2": chi2,
                "hwe_p": p_hwe,
                "paper_hf": None,
                "paper_j": None,
            }
    paper = {
        "PLCD4": (0.039, 0.001),
        "FGD4": (0.039, 0.001),
        "MTRF1": (0.002, 0.114),
        "MYH1": (0.006, 0.109),
        "GALNT2": (0.055, 0.001),
        "DPF2": (0.066, 0.001),
        "MUS81": (0.066, 0.001),
    }
    for locus, (ph, pj) in paper.items():
        if locus in out:
            out[locus]["breeds"]["HF"]["paper_hf"] = ph
            out[locus]["breeds"]["J"]["paper_j"] = pj
    cumulative = {}
    unique_loci_6 = ["PLCD4", "FGD4", "MTRF1", "MYH1", "GALNT2", "DPF2"]
    extra_table1_aaf = {"HF": [0.124, 0.124], "J": [0.002, 0.002]}  # Chr2:22, ARSI
    for breed, col in (("HF", "breed_hf"), ("J", "breed_j")):
        ids = None
        any_carrier = None
        any_carrier_6 = None
        for loc_name, series in geno_by_locus.items():
            if ids is None:
                ids = series.index.to_numpy()
                any_carrier = np.zeros(len(ids), dtype=bool)
                any_carrier_6 = np.zeros(len(ids), dtype=bool)
            vals = series.reindex(ids).to_numpy() > 0
            any_carrier |= vals
            if loc_name in unique_loci_6:
                any_carrier_6 |= vals
        if ids is None:
            continue
        breed_ids = set(other6.loc[other6[col] == 16, "ID"].astype(int))
        keep = np.array([int(i) in breed_ids for i in ids])
        emp_any_7 = float(any_carrier[keep].mean()) if keep.any() else None
        emp_any_6 = float(any_carrier_6[keep].mean()) if keep.any() else None
        sum_carrier_6 = float(
            sum(out[m]["breeds"][breed]["carrier_freq"] for m in unique_loci_6)
        )
        p_non_extra = 1.0
        for q_extra in extra_table1_aaf[breed]:
            c_extra = 1.0 - (1.0 - q_extra) ** 2
            p_non_extra *= 1.0 - c_extra
        expected_any_8 = (
            float(1.0 - (1.0 - emp_any_6) * p_non_extra) if emp_any_6 is not None else None
        )
        cumulative[breed] = {
            "n": int(keep.sum()),
            "any_carrier_fraction": emp_any_7,
            "any_carrier_fraction_6_unique_loci": emp_any_6,
            "sum_of_carrier_frequencies_6_unique_loci": sum_carrier_6,
            "expected_any_carrier_8_table1_loci": expected_any_8,
            "loci": list(geno_by_locus),
            "note": (
                "gen_4 includes the 7 candidate mutations (DPF2 and MUS81 are in strong LD). "
                "Including Chr2:22 and ARSI at Table 1 AAFs yields >40% cumulative carrier "
                "rate in both HF and Jersey."
            ),
        }
    return {"loci": out, "cumulative": cumulative}


def compute_depletion(freq: dict) -> dict:
    rows = []
    for locus in LOCI:
        if locus["gen4"] is None:
            continue
        breed = locus["breed"]
        if breed is None:
            continue
        data = freq["loci"][locus["locus"]]["breeds"][breed]
        n = data["n"]
        p = data["alt_freq"]
        k = data["genotype_counts"]["2"]
        expected = n * p**2 if p is not None else None
        if expected is None or expected <= 0 or n == 0:
            continue
        binom = stats.binomtest(k, n, p**2, alternative="less")
        # Binomial log-likelihood ratio, one-sided in the depletion direction.
        phat = k / n
        ll0 = stats.binom.logpmf(k, n, p**2)
        ll1 = stats.binom.logpmf(k, n, phat) if 0 < phat < 1 else (0.0 if k == 0 else -np.inf)
        lrt = max(0.0, 2 * (ll1 - ll0))
        lrt_p_two = float(stats.chi2.sf(lrt, 1))
        lrt_p_one = lrt_p_two / 2 if k < expected else 1 - lrt_p_two / 2
        rows.append(
            {
                "locus": locus["locus"],
                "gene": locus["gene"],
                "breed": breed,
                "n": n,
                "alt_freq": p,
                "observed_hom_alt": k,
                "expected_hom_alt": expected,
                "observed_over_expected": k / expected,
                "binom_p_one_sided": float(binom.pvalue),
                "lrt_stat": float(lrt),
                "lrt_p_one_sided": float(lrt_p_one),
            }
        )
    return {
        "rows": rows,
        "note": "One-sided tests for mutant-homozygote depletion vs HWE expectation.",
    }


def compute_carrier_cross() -> dict:
    o5 = pd.read_csv(OTHER / "other_5.txt", sep=r"\s+", header=None, names=["off", "sire", "dam"])
    rows = []
    for locus in LOCI:
        if locus["gen5"] is None:
            continue
        iids, g = genotype_for_locus(locus, "gen5")
        mp = pd.Series(g, index=iids)
        off = o5["off"].map(mp)
        sire = o5["sire"].map(mp)
        dam = o5["dam"].map(mp)
        mask = (sire == 1) & (dam == 1) & off.notna()
        counts = (
            off[mask].astype(int).value_counts().reindex([0, 1, 2], fill_value=0).astype(int)
        )
        n = int(counts.sum())
        expected = [n / 4, n / 2, n / 4]
        chi2 = None
        p = None
        if n > 0:
            chi2, p = stats.chisquare(counts.to_numpy(dtype=float), f_exp=expected)
            chi2, p = float(chi2), float(p)
        rows.append(
            {
                "locus": locus["locus"],
                "gene": locus["gene"],
                "n_trios_both_het": n,
                "offspring_counts": {str(k): int(v) for k, v in counts.items()},
                "expected_1_2_1": [float(v) for v in expected],
                "chi2": chi2,
                "p": p,
            }
        )
    return {
        "rows": rows,
        "selection": (
            "Trios with sire genotype=1 and dam genotype=1 for the same variant; "
            "offspring genotype non-missing."
        ),
        "note": "Small n; chi-square approximate. Counts are aggregated only.",
    }


def compute_national_burden(freq: dict) -> dict:
    gen10 = pd.read_csv(GENS / "gen_10.1.txt", sep=r"\s+").set_index("KEY")
    o11 = pd.read_csv(OTHER / "other_11.txt", sep=r"\s+")
    mutations = ["PLCD4", "FGD4", "MTRF1", "MYH1", "GALNT2", "DPF2"]
    bw_significant_mutations = {"PLCD4", "FGD4", "MTRF1", "GALNT2", "DPF2"}
    per_mutation = {}
    p_any = np.zeros(len(o11), dtype=float)
    p_any_bw = np.zeros(len(o11), dtype=float)
    for mut in mutations:
        if mut not in gen10.columns:
            continue
        f_hf = freq["loci"][mut]["breeds"]["HF"]["alt_freq"]
        f_j = freq["loci"][mut]["breeds"]["J"]["alt_freq"]
        all_counts = freq["loci"][mut]["breeds"]["HF"]["genotype_counts"]
        all_counts_j = freq["loci"][mut]["breeds"]["J"]["genotype_counts"]
        n_hf = freq["loci"][mut]["breeds"]["HF"]["n"]
        n_j = freq["loci"][mut]["breeds"]["J"]["n"]
        n_all = n_hf + n_j
        f_all = (
            (2 * all_counts["2"] + all_counts["1"] + 2 * all_counts_j["2"] + all_counts_j["1"])
            / (2 * n_all)
            if n_all
            else None
        )
        sire_g = o11["SIRE_ID"].map(gen10[mut]).astype(float)
        mgs_g = o11["MGS_ID"].map(gen10[mut]).astype(float)
        hf = o11["mgd_HF_16"].astype(float) / 16
        j = o11["mgd_J_16"].astype(float) / 16
        remainder = (1 - hf - j).clip(lower=0)
        q_mgd = hf * f_hf + j * f_j + remainder * (f_all if f_all is not None else 0.0)
        p_sire = sire_g / 2
        p_dam = (mgs_g / 2 + q_mgd) / 2
        p_homo = (p_sire * p_dam).clip(lower=0, upper=1)
        p_any = 1 - (1 - p_any) * (1 - p_homo)
        if mut in bw_significant_mutations:
            p_any_bw = 1 - (1 - p_any_bw) * (1 - p_homo)
        by_year = (
            pd.DataFrame({"birth_year": o11["birth_year"], "p_homo": p_homo})
            .groupby("birth_year")["p_homo"]
            .agg(["sum", "mean", "size"])
            .reset_index()
        )
        per_mutation[mut] = {
            "freq_HF": f_hf,
            "freq_J": f_j,
            "freq_all": f_all,
            "expected_total": float(p_homo.sum()),
            "rate_per_animal": float(p_homo.mean()),
            "scaled_per_year": float(p_homo.mean() * FEMALES_BORN_PER_YEAR),
            "by_birth_year": [
                {
                    "birth_year": int(r.birth_year),
                    "expected_homozygotes": float(r.sum),
                    "rate": float(r.mean),
                    "n_animals": int(r.size),
                }
                for r in by_year.itertuples()
            ],
        }
    any_by_year = (
        pd.DataFrame({"birth_year": o11["birth_year"], "p_any": p_any})
        .groupby("birth_year")["p_any"]
        .agg(["sum", "mean", "size"])
        .reset_index()
    )
    any_summary = {
        "expected_total": float(p_any.sum()),
        "rate_per_animal": float(p_any.mean()),
        "scaled_per_year": float(p_any.mean() * FEMALES_BORN_PER_YEAR),
        "by_birth_year": [
            {
                "birth_year": int(r.birth_year),
                "expected_homozygotes": float(r.sum),
                "rate": float(r.mean),
                "n_animals": int(r.size),
            }
            for r in any_by_year.itertuples()
        ],
    }
    any_bw_by_year = (
        pd.DataFrame({"birth_year": o11["birth_year"], "p_any_bw": p_any_bw})
        .groupby("birth_year")["p_any_bw"]
        .agg(["sum", "mean", "size"])
        .reset_index()
    )
    any_bw_summary = {
        "mutations_included": sorted(bw_significant_mutations),
        "expected_total": float(p_any_bw.sum()),
        "rate_per_animal": float(p_any_bw.mean()),
        "scaled_per_year": float(p_any_bw.mean() * FEMALES_BORN_PER_YEAR),
        "by_birth_year": [
            {
                "birth_year": int(r.birth_year),
                "expected_homozygotes": float(r.sum),
                "rate": float(r.mean),
                "n_animals": int(r.size),
            }
            for r in any_bw_by_year.itertuples()
        ],
    }
    return {
        "n_animals": int(len(o11)),
        "mutations": per_mutation,
        "any_mutation": any_summary,
        "any_bw_significant_mutation": any_bw_summary,
        "paper_homozygous_rate": PAPER_HOMOZYGOUS_RATE,
        "paper_scaled_homozygotes": PAPER_SCALED_HOMOZYGOUS,
        "females_born_per_year": FEMALES_BORN_PER_YEAR,
        "mendelian_derivation": (
            "P(calf AA) = P(sire transmits A) * P(dam transmits A), where "
            "P(sire transmits A) = g_sire / 2 for sire genotype g_sire in {0,1,2}, and "
            "P(dam transmits A) = 0.5 * (g_MGS / 2) + 0.5 * q_MGD because the dam inherits "
            "one allele from the genotyped maternal grandsire (MGS) and one allele from the "
            "ungenotyped maternal granddam (MGD) with breed-weighted allele frequency q_MGD."
        ),
        "q_mgd_method": (
            "q_MGD = (mgd_HF_16/16)*f_HF + (mgd_J_16/16)*f_J + remainder*f_all; "
            "DPF2 column represents DPF2/MUS81."
        ),
        "note": (
            "Independence across mutations assumed for the any-mutation rate; "
            "DPF2/MUS81 treated as one locus. Excluding MYH1 (which is not genome-wide "
            "significant for body weight and is introduced later in the paper) gives "
            "rate=1.015% (~9,543 females/year), matching the paper's ~1% (~9,700) figure; "
            "including MYH1 across all 6 gen_10 columns gives 1.315% (~12,367 females/year)."
        ),
    }


def compute_sire_amplification() -> dict:
    gen10 = pd.read_csv(GENS / "gen_10.1.txt", sep=r"\s+").set_index("KEY")
    o11 = pd.read_csv(OTHER / "other_11.txt", sep=r"\s+")
    progeny = o11.groupby(["birth_year", "SIRE_ID"]).size().reset_index(name="n_progeny")
    total_by_year = progeny.groupby("birth_year")["n_progeny"].sum()
    per_mutation = {}
    for mut in ["PLCD4", "FGD4", "MTRF1", "MYH1", "GALNT2", "DPF2"]:
        if mut not in gen10.columns:
            continue
        carrier_keys = set(gen10.index[gen10[mut] >= 1].astype(int))
        carrier = progeny[progeny["SIRE_ID"].isin(carrier_keys)]
        carrier_by_year = carrier.groupby("birth_year")["n_progeny"].sum()
        share = {
            str(int(y)): float(carrier_by_year.get(y, 0) / total_by_year.get(y, np.nan))
            for y in total_by_year.index
        }
        top = (
            carrier.groupby("SIRE_ID")["n_progeny"]
            .sum()
            .sort_values(ascending=False)
            .head(10)
            .reset_index()
        )
        top_rows = [
            {
                "rank": i + 1,
                "n_progeny": int(r.n_progeny),
                "genotype": int(gen10.loc[int(r.SIRE_ID), mut]),
            }
            for i, r in enumerate(top.itertuples())
        ]
        per_mutation[mut] = {
            "n_carrier_sires": len(carrier_keys),
            "total_progeny": int(carrier["n_progeny"].sum()),
            "carrier_share_of_all_progeny": float(
                carrier["n_progeny"].sum() / progeny["n_progeny"].sum()
            ),
            "carrier_share_by_birth_year": share,
            "top_carrier_sires_ranked": top_rows,
        }
    return {
        "per_mutation": per_mutation,
        "note": "Ranks only; no sire IDs exported. Progeny counts from other_11 by SIRE_ID.",
    }


def compute_prospective_calves() -> dict:
    d = pd.read_csv(PHENS / "phen_3.1.txt", sep=r"\s+")
    d = d.dropna(subset=["weight", "age", "HERD", "breed_hf", "breed_j", "VARIANT"]).copy()
    d["HERD"] = d["HERD"].astype(str)
    variant_dummies = pd.get_dummies(d["VARIANT"], prefix="V", dtype=float)
    ref = "V_CONTROL"
    variant_cols = [c for c in variant_dummies.columns if c != ref]
    herd_dummies = pd.get_dummies(d["HERD"], prefix="H", drop_first=True, dtype=float)
    x = pd.concat(
        [
            pd.Series(1.0, index=d.index, name="const"),
            variant_dummies[variant_cols],
            d[["age", "breed_hf", "breed_j"]].astype(float),
            herd_dummies,
        ],
        axis=1,
    )
    res = OLS(d["weight"].astype(float), x).fit()
    ci = res.conf_int()
    effects = {}
    for col in variant_cols:
        variant = col.removeprefix("V_")
        effects[variant] = {
            "n": int((d["VARIANT"] == variant).sum()),
            "mean_weight_unadjusted": float(d.loc[d["VARIANT"] == variant, "weight"].mean()),
            "effect_vs_control": float(res.params[col]),
            "se": float(res.bse[col]),
            "ci_low": float(ci.loc[col, 0]),
            "ci_high": float(ci.loc[col, 1]),
            "p": float(res.pvalues[col]),
        }
    return {
        "n_calves": int(len(d)),
        "n_herds": int(d["HERD"].nunique()),
        "control_n": int((d["VARIANT"] == "CONTROL").sum()),
        "control_mean_weight": float(d.loc[d["VARIANT"] == "CONTROL", "weight"].mean()),
        "effects": effects,
        "model": "weight ~ VARIANT + age + HERD + breed_hf + breed_j; CONTROL reference",
        "note": "Small groups (PLCD4 n=9, FGD4 n=22, DPF2/MUS81 n=11); CIs are wide.",
    }


def compute_loci_metadata() -> dict:
    loci = []
    for locus in LOCI + COMPARATORS:
        iids, g = genotype_for_locus(locus, "gen1" if locus["gen1"] else "gen3")
        counts = {
            str(k): int(v)
            for k, v in pd.Series(g).value_counts().reindex([0, 1, 2], fill_value=0).items()
        }
        a1_freq = (2 * counts["2"] + counts["1"]) / (2 * len(g))
        entry = {
            "locus": locus["locus"],
            "gene": locus["gene"],
            "chr": locus["chr"],
            "bp": locus["bp"],
            "candidate": locus["candidate"],
            "consequence": locus["consequence"],
            "gen1": locus["gen1"],
            "gen2": locus["gen2"],
            "gen3": locus["gen3"],
            "window": locus["window"],
            "window_sensitivity": locus["window_sensitivity"],
            "breed": locus["breed"],
            "genotype_counts_all": counts,
            "a1_freq_all": a1_freq,
        }
        if locus["gen3"] is not None and locus["gen1"] is not None:
            _, g3 = genotype_for_locus(locus, "gen3")
            entry["gen3_matches_gen1"] = bool(np.array_equal(g, g3))
        loci.append(entry)
    return {"loci": loci, "table1_targets": TABLE1}


def main() -> None:
    manifest = []
    loci_meta = compute_loci_metadata()
    manifest.append(write_json("loci.json", loci_meta, ["dryad/gens/*.bim", "dryad/gens/*.fam"]))

    (
        class_rows,
        table1_summary,
        effect_rows,
        published_rows,
        window_sensitivity_rows,
    ) = compute_table1()
    extfig2c = compute_extfig2c_comparison(effect_rows, class_rows)
    manifest.append(
        write_json(
            "genotype_class_stats.json",
            {"rows": class_rows, "extfig2c_comparison": extfig2c},
            [
                "dryad/phens/phen_1.*.txt",
                "dryad/other/other_2.*.tsv",
                "dryad/gens/gen_1.*",
                "dryad/source_data/source_data.ExtFig2c.csv",
            ],
        )
    )
    manifest.append(
        write_json(
            "effects.json",
            {
                "rows": effect_rows,
                "window_sensitivity_rows": window_sensitivity_rows,
                "table1_targets": TABLE1,
                "table1_comparison_summary": table1_summary,
                "published_other_1": published_rows,
                "extfig2c_comparison": extfig2c,
            },
            [
                "dryad/other/other_1.*.txt",
                "dryad/phens/phen_1.*.txt",
                "dryad/gens/gen_1.*",
                "dryad/source_data/source_data.ExtFig2c.csv",
            ],
        )
    )

    pleio_rows = compute_pleiotropy()
    manifest.append(
        write_json(
            "pleiotropy.json",
            {"rows": pleio_rows, "n_rows": len(pleio_rows)},
            ["dryad/phens/phen_1.*.txt", "dryad/phens/phen_4.*.txt", "dryad/gens/gen_3.*"],
        )
    )

    validation = compute_validation(pleio_rows)
    manifest.append(
        write_json(
            "validation_milk.json",
            validation,
            ["dryad/phens/phen_2.*.txt", "dryad/phens/phen_4.1-3.txt", "dryad/gens/gen_2.*"],
        )
    )

    freq = compute_frequencies()
    manifest.append(
        write_json(
            "frequencies.json",
            freq,
            ["dryad/gens/gen_4.*", "dryad/other/other_6.txt"],
        )
    )
    depletion = compute_depletion(freq)
    manifest.append(
        write_json(
            "depletion.json",
            depletion,
            ["dryad/gens/gen_4.*", "dryad/other/other_6.txt"],
        )
    )
    carrier = compute_carrier_cross()
    manifest.append(
        write_json(
            "carrier_cross.json",
            carrier,
            ["dryad/gens/gen_5.*", "dryad/other/other_5.txt"],
        )
    )
    burden = compute_national_burden(freq)
    manifest.append(
        write_json(
            "national_burden.json",
            burden,
            ["dryad/gens/gen_10.1.txt", "dryad/other/other_11.txt", "dryad/gens/gen_4.*"],
        )
    )
    amplification = compute_sire_amplification()
    manifest.append(
        write_json(
            "sire_amplification.json",
            amplification,
            ["dryad/gens/gen_10.1.txt", "dryad/other/other_11.txt"],
        )
    )
    calves = compute_prospective_calves()
    manifest.append(
        write_json(
            "prospective_calves.json",
            calves,
            ["dryad/phens/phen_3.1.txt"],
        )
    )

    manifest_payload = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": now_iso(),
        "code_version": CODE_VERSION,
        "files": manifest,
        "total_bytes": int(sum(m["bytes"] for m in manifest)),
    }
    (OUT / "manifest.json").write_text(json.dumps(jsonable(manifest_payload), indent=2) + "\n")
    print(f"wrote {len(manifest)} exports to {OUT} ({manifest_payload['total_bytes']} bytes)")


if __name__ == "__main__":
    main()
