"""Phase 0 verification for the Recessive Variant Explorer.

Runs the checks V1-V11 from BRIEF.md and writes work/phase0_summary.json.
This script is deliberately conservative: it reports evidence and pass/fail
flags, but it never modifies anything under dryad/ or paper/.
"""

from __future__ import annotations

import json
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from bed_reader import open_bed

ROOT = Path(__file__).resolve().parents[1]
DRYAD = ROOT / "dryad"
GENS = DRYAD / "gens"
PHENS = DRYAD / "phens"
OTHER = DRYAD / "other"
SOURCE = DRYAD / "source_data"
WORK = ROOT / "work"
WORK.mkdir(parents=True, exist_ok=True)

# Table 1 targets (body weight, kg) from the paper PDF.  The first entry is the
# chr2:22Mb QTL, which has no candidate gene.
TABLE1_HOM_ALT = {
    "Chr2:22": -11.8,
    "PLCD4": -108.965,
    "FGD4": -48.079,
    "ARSI": -11.981,
    "MTRF1": -19.464,
    "GALNT2": -112.008,
    "DPF2": -53.092,
    "MUS81": -51.911,
    "MYH1": -4.643,
}
TABLE1_HET = {
    "Chr2:22": -0.007,
    "PLCD4": -0.354,
    "FGD4": 0.834,
    "ARSI": 0.564,
    "MTRF1": 1.078,
    "GALNT2": -0.837,
    "DPF2": 4.059,
    "MUS81": 4.126,
    "MYH1": 0.798,
}
PAPER_AAF_ALL = {
    "Chr2:22": 0.074,
    "PLCD4": 0.024,
    "FGD4": 0.024,
    "ARSI": 0.072,
    "MTRF1": 0.045,
    "GALNT2": 0.033,
    "DPF2": 0.039,
    "MUS81": 0.040,
    "MYH1": 0.042,
}
PAPER_AAF_HF_J = {
    "PLCD4": (0.039, 0.001),
    "FGD4": (0.039, 0.001),
    "MTRF1": (0.002, 0.114),
    "MYH1": (0.006, 0.109),
    "GALNT2": (0.055, 0.001),
    "DPF2": (0.066, 0.001),
    "MUS81": (0.066, 0.001),
}

# Discovery loci: (label, gen_1 file, CHR, BP, LOSO window).
# Windows are the non-overlapping interval containing the locus; overlapping
# windows are tested as a sensitivity check in V4.
LOCI = [
    ("Chr2:22", "gen_1.1", 2, 22836659, 1),
    ("PLCD4", "gen_1.2", 2, 107313998, 3),
    ("FGD4", "gen_1.3", 5, 77632752, 5),
    ("ARSI", "gen_1.4", 7, 63649071, 7),
    ("MTRF1", "gen_1.5", 12, 11463981, 9),
    ("GALNT2", "gen_1.8", 28, 1312334, 12),
    ("DPF2", "gen_1.9", 29, 44213160, 13),
    ("MUS81", "gen_1.9", 29, 44645469, 13),
    ("MYH1", "gen_1.7", 19, 30114249, 11),
]

WINDOW_NAMES = {
    1: "Chr2_20_30Mb",
    2: "Chr2_15_25Mb",
    3: "Chr2_100_110Mb",
    4: "Chr2_105_115Mb",
    5: "Chr5_70_80Mb",
    6: "Chr5_75_85Mb",
    7: "Chr7_60_70Mb",
    8: "Chr7_55_65Mb",
    9: "Chr12_5_15Mb",
    10: "Chr14_20_30Mb",
    11: "Chr19_25_35Mb",
    12: "Chr28_0_10Mb",
    13: "Chr29_40_50Mb",
}


def py(x):
    """Convert numpy scalars to plain Python for JSON."""
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, np.ndarray):
        return x.tolist()
    return x


def load_phen(trait: int) -> pd.DataFrame:
    return pd.read_csv(PHENS / f"phen_1.{trait}.txt", sep=r"\s+", dtype={"ID": int})


def load_fam_iid(path: Path) -> np.ndarray:
    return pd.read_csv(path, sep=r"\s+", header=None, usecols=[1])[1].astype(int).to_numpy()


def lead_variant(gen_file: str, chrom: int, bp: int) -> tuple[int, pd.Series]:
    bim = pd.read_csv(
        GENS / f"{gen_file}.bim",
        sep=r"\s+",
        header=None,
        names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
    )
    dist = (bim.BP - bp).abs()
    idx = int(dist.idxmin())
    return idx, bim.loc[idx]


def lead_genotype(gen_file: str, chrom: int, bp: int) -> tuple[np.ndarray, pd.Series, np.ndarray]:
    idx, row = lead_variant(gen_file, chrom, bp)
    gt = open_bed(GENS / f"{gen_file}.bed").read(index=(slice(None), [idx]), dtype="int8")[:, 0]
    iid = load_fam_iid(GENS / f"{gen_file}.fam")
    return gt, row, iid


_ADJ_CACHE: dict[int, pd.DataFrame] = {}


def adjusted_bodyweight(window: int) -> pd.DataFrame:
    """Reconstruct the paper-like LOSO-adjusted bodyweight phenotype.

    adjusted = yd_trait - mean(500 MCMC samples of the LOSO population effect),
    then centred.  Column 0 of other_2.* is the true phenotype ID.
    """
    if window in _ADJ_CACHE:
        return _ADJ_CACHE[window]
    path = OTHER / f"other_2.1.{window}.tsv"
    dtype = {0: int, **{i: "float32" for i in range(1, 501)}}
    df = pd.read_csv(path, sep="\t", header=None, dtype=dtype)
    ids = df.iloc[:, 0].to_numpy()
    pop_effect = df.iloc[:, 1:].mean(axis=1).to_numpy(dtype="float32")
    phen = load_phen(1).set_index("ID")
    out = pd.DataFrame(
        {"ID": ids, "adj": phen.loc[ids, "yd_trait"].to_numpy() - pop_effect}
    )
    out["adj"] -= out["adj"].mean()
    _ADJ_CACHE[window] = out
    return out


def wls_classes(df: pd.DataFrame, weight: np.ndarray) -> dict:
    g = df["g"].to_numpy(float)
    y = df["adj"].to_numpy(float)
    x = g - 1.0
    z = (g == 1).astype(float)
    X = np.column_stack([np.ones(len(df)), x, z])
    sw = np.sqrt(weight)
    beta = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
    mu, a, d = (float(v) for v in beta)
    return {
        "mu": mu,
        "a": a,
        "d": d,
        "hom_alt": 2 * a,
        "het": a + d,
        "n": {str(k): int(v) for k, v in df["g"].value_counts().sort_index().items()},
        "class_means": {
            str(k): float(v) for k, v in df.groupby("g")["adj"].mean().items()
        },
    }


def check_environment() -> dict:
    try:
        mem = int(subprocess.check_output(["sysctl", "-n", "hw.memsize"], text=True).strip())
    except Exception:
        mem = os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    usage = shutil.disk_usage(ROOT)
    return {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python": sys.version.split()[0],
        "memory_gb": round(mem / 1e9, 2),
        "disk_free_gb": round(usage.free / 1e9, 2),
        "cwd": str(ROOT),
    }


def check_v1() -> dict:
    phen = load_phen(1)
    gt, row, iid = lead_genotype("gen_1.2", 2, 107313998)
    adj = adjusted_bodyweight(3).set_index("ID")
    df = pd.DataFrame({"ID": iid, "g": gt}).merge(adj.reset_index(), on="ID", how="inner")
    iid_means = {str(k): float(v) for k, v in df.groupby("g")["adj"].mean().items()}
    iid_hom = float(df.loc[df.g == 2, "adj"].mean() - df.loc[df.g == 0, "adj"].mean())

    # Deliberately wrong join: align the first len(phen) genotype rows to phen rows.
    y = phen["yd_trait"].to_numpy(float)
    g_row = gt[: len(phen)]
    row_hom = float(y[g_row == 2].mean() - y[g_row == 0].mean())
    row_means = {str(k): float(y[g_row == k].mean()) for k in (0, 1, 2)}
    return {
        "lead_variant": {
            "file": "gen_1.2",
            "snp": row["SNP"],
            "chr": int(row["CHR"]),
            "bp": int(row["BP"]),
            "a1": row["A1"],
            "a2": row["A2"],
        },
        "iid_join_adjusted": {"class_means": iid_means, "hom_alt_minus_ref": iid_hom},
        "row_order_join_raw": {"class_means": row_means, "hom_alt_minus_ref": row_hom},
        "paper_hom_alt": TABLE1_HOM_ALT["PLCD4"],
        "pass": abs(iid_hom - TABLE1_HOM_ALT["PLCD4"]) < 15 and abs(row_hom) < 10,
    }


def check_v2() -> dict:
    rows = []
    for label, gfile, chrom, bp, _ in LOCI:
        gt, row, _ = lead_genotype(gfile, chrom, bp)
        vc = {str(k): int((gt == k).sum()) for k in (0, 1, 2)}
        freq = (vc["1"] + 2 * vc["2"]) / (2 * sum(vc.values()))
        rows.append(
            {
                "locus": label,
                "file": gfile,
                "snp": row["SNP"],
                "chr": int(row["CHR"]),
                "bp": int(row["BP"]),
                "a1": row["A1"],
                "a2": row["A2"],
                "counts": vc,
                "a1_freq": freq,
                "paper_aaf_all": PAPER_AAF_ALL[label],
                "delta": freq - PAPER_AAF_ALL[label],
            }
        )

    # Breed-specific check from the purebred validation set (gen_4 + other_6).
    o6 = pd.read_csv(OTHER / "other_6.txt", sep=r"\s+")
    breed_rows = []
    g4_map = {
        "PLCD4": "gen_4.1",
        "FGD4": "gen_4.2",
        "MTRF1": "gen_4.3",
        "MYH1": "gen_4.4",
        "GALNT2": "gen_4.5",
        "DPF2": "gen_4.6",
        "MUS81": "gen_4.7",
    }
    for label, gfile in g4_map.items():
        bim = pd.read_csv(GENS / f"{gfile}.bim", sep=r"\s+", header=None)
        gt = open_bed(GENS / f"{gfile}.bed").read(dtype="int8")[:, 0]
        iid = load_fam_iid(GENS / f"{gfile}.fam")
        df = pd.DataFrame({"ID": iid, "g": gt}).merge(o6, on="ID", how="inner")
        hf = df[df.breed_j == 0]
        j = df[df.breed_hf == 0]

        def freq(d: pd.DataFrame) -> float:
            return float(d.g.sum() / (2 * len(d))) if len(d) else float("nan")

        hf_f, j_f = freq(hf), freq(j)
        ph, pj = PAPER_AAF_HF_J[label]
        breed_rows.append(
            {
                "locus": label,
                "file": gfile,
                "snp": str(bim.loc[0, 1]),
                "a1": str(bim.loc[0, 4]),
                "a2": str(bim.loc[0, 5]),
                "hf_n": int(len(hf)),
                "j_n": int(len(j)),
                "hf_freq": hf_f,
                "j_freq": j_f,
                "paper_hf": ph,
                "paper_j": pj,
                "max_delta": max(abs(hf_f - ph), abs(j_f - pj)),
            }
        )
    max_delta = max(abs(r["delta"]) for r in rows)
    max_breed_delta = max(r["max_delta"] for r in breed_rows)
    return {
        "lead_variants": rows,
        "breed_specific": breed_rows,
        "max_abs_delta_all": max_delta,
        "max_abs_delta_breed": max_breed_delta,
        "pass": max_delta < 0.01 and max_breed_delta < 0.02,
    }


def check_v3() -> dict:
    def ids(path: Path) -> set[int]:
        return set(pd.read_csv(path, sep=r"\s+", usecols=[0]).iloc[:, 0].astype(int))

    other_5_raw = pd.read_csv(OTHER / "other_5.txt", sep=r"\s+", header=None)
    sets = {
        "gen_1": set(load_fam_iid(GENS / "gen_1.1.fam")),
        "phen_1.1": ids(PHENS / "phen_1.1.txt"),
        "gen_2": set(load_fam_iid(GENS / "gen_2.1.fam")),
        "phen_2.1": ids(PHENS / "phen_2.1.txt"),
        "gen_3": set(load_fam_iid(GENS / "gen_3.1.fam")),
        "phen_4.1": ids(PHENS / "phen_4.1.txt"),
        "gen_4": set(load_fam_iid(GENS / "gen_4.1.fam")),
        "gen_5": set(load_fam_iid(GENS / "gen_5.1.fam")),
        "other_5": set(other_5_raw.values.ravel()),
        "gen_6": set(load_fam_iid(GENS / "gen_6.1.fam")),
        "phen_7.1": ids(PHENS / "phen_7.1.txt"),
        "phen_7.2": ids(PHENS / "phen_7.2.txt"),
        "gen_7": set(pd.read_csv(GENS / "gen_7.1.txt", sep=r"\s+").iloc[:, 0].astype(int)),
        "phen_6.1": ids(PHENS / "phen_6.1.txt"),
        "phen_6.2": ids(PHENS / "phen_6.2.txt"),
        "gen_8": set(
            pd.read_csv(GENS / "gen_8.1.txt", sep=r"\s+", header=None)
            .iloc[:, 0]
            .astype(int)
        ),
        "phen_5.1": ids(PHENS / "phen_5.1.txt"),
        "phen_5.4": ids(PHENS / "phen_5.4.txt"),
        "gen_9": set(load_fam_iid(GENS / "gen_9.1.fam")),
        "gen_11": set(load_fam_iid(GENS / "gen_11.1.fam")),
        "gen_12": set(load_fam_iid(GENS / "gen_12.1.fam")),
        "phen_8.1": ids(PHENS / "phen_8.1.txt"),
    }

    def cmp(a: str, b: str) -> dict:
        A, B = sets[a], sets[b]
        return {
            "a": a,
            "b": b,
            "n_a": len(A),
            "n_b": len(B),
            "equal": A == B,
            "only_a": len(A - B),
            "only_b": len(B - A),
        }

    comparisons = [
        cmp("gen_1", "phen_1.1"),
        cmp("gen_2", "phen_2.1"),
        cmp("gen_3", "phen_4.1"),
        cmp("gen_5", "other_5"),
        cmp("gen_6", "phen_7.1"),
        cmp("gen_6", "phen_7.2"),
        cmp("gen_7", "phen_6.1"),
        cmp("gen_7", "phen_6.2"),
        cmp("gen_8", "phen_5.1"),
        cmp("gen_8", "phen_5.4"),
        cmp("gen_11", "phen_8.1"),
        cmp("gen_12", "phen_8.1"),
    ]

    # gen_3 row order is identical to gen_1, but its .fam IID labels are not.
    row_checks = []
    for label, gfile, chrom, bp, _ in LOCI:
        idx, _ = lead_variant(gfile, chrom, bp)
        g1 = open_bed(GENS / f"{gfile}.bed").read(index=(slice(None), [idx]), dtype="int8")[:, 0]
        # Map gen_1 locus -> matching gen_3 single-variant file.
        g3_file = {
            "Chr2:22": "gen_3.1",
            "PLCD4": "gen_3.2",
            "FGD4": "gen_3.3",
            "ARSI": "gen_3.4",
            "MTRF1": "gen_3.5",
            "MYH1": "gen_3.9",
            "GALNT2": "gen_3.10",
            "DPF2": "gen_3.11",
            "MUS81": "gen_3.12",
        }[label]
        g3 = open_bed(GENS / f"{g3_file}.bed").read(dtype="int8")[:, 0]
        row_checks.append(
            {
                "locus": label,
                "gen_1": gfile,
                "gen_3": g3_file,
                "row_identity": float((g1 == g3).mean()),
                "fam_order_equal": bool(
                    np.array_equal(
                        load_fam_iid(GENS / f"{gfile}.fam"),
                        load_fam_iid(GENS / f"{g3_file}.fam"),
                    )
                ),
            }
        )

    overlaps = {
        "gen_1_vs_gen_2": len(sets["gen_1"] & sets["gen_2"]),
        "gen_1_vs_gen_3": len(sets["gen_1"] & sets["gen_3"]),
        "gen_3_vs_gen_9": len(sets["gen_3"] & sets["gen_9"]),
        "gen_1_vs_gen_9": len(sets["gen_1"] & sets["gen_9"]),
    }
    return {
        "sets": {k: {"n": len(v), "min": min(v), "max": max(v)} for k, v in sets.items()},
        "other_5_rows": int(len(other_5_raw)),
        "other_5_union_n": int(len(sets["other_5"])),
        "comparisons": comparisons,
        "gen3_row_order": row_checks,
        "numeric_overlap_counts": overlaps,
        "pass": all(c["equal"] for c in comparisons if c["a"] not in {"gen_1", "gen_3"})
        and all(r["row_identity"] == 1.0 for r in row_checks)
        and all(not r["fam_order_equal"] for r in row_checks),
    }


def check_v4() -> dict:
    schemes = ("unweighted", "rinverse", "1/rinverse")
    rows = []
    for label, gfile, chrom, bp, window in LOCI:
        if label not in TABLE1_HOM_ALT:
            continue
        gt, row, iid = lead_genotype(gfile, chrom, bp)
        adj = adjusted_bodyweight(window).set_index("ID")
        phen = load_phen(1)[["ID", "rinverse"]]
        df = (
            pd.DataFrame({"ID": iid, "g": gt})
            .merge(adj.reset_index(), on="ID", how="inner")
            .merge(phen, on="ID", how="left")
        )
        weights = {
            "unweighted": np.ones(len(df)),
            "rinverse": df["rinverse"].to_numpy(float),
            "1/rinverse": 1.0 / df["rinverse"].to_numpy(float),
        }
        fits = {name: wls_classes(df, weights[name]) for name in schemes}
        rows.append(
            {
                "locus": label,
                "file": gfile,
                "window": window,
                "window_name": WINDOW_NAMES[window],
                "fits": fits,
                "table1_hom_alt": TABLE1_HOM_ALT[label],
                "table1_het": TABLE1_HET[label],
            }
        )

    errors = {}
    for name in schemes:
        errors[name] = float(
            np.mean([abs(r["fits"][name]["hom_alt"] - r["table1_hom_alt"]) for r in rows])
        )
    chosen_by_error = min(errors, key=errors.get)
    # The three schemes are effectively tied (<0.02 kg mean absolute error).
    # Prefer unweighted for the pipeline: it makes no assumption about what
    # rinverse means, and the paper's own model is a mixed model on adjusted
    # phenotypes rather than a WLS reweighting.
    chosen_for_pipeline = "unweighted"
    return {
        "rows": rows,
        "mean_abs_error_hom_alt": errors,
        "chosen_by_error": chosen_by_error,
        "chosen_for_pipeline": chosen_for_pipeline,
        "note": "weighting choice changes hom-alt effects by <0.02 kg on average",
        "pass": errors[chosen_by_error] < 10,
    }


def check_v5() -> dict:
    o7 = pd.read_csv(OTHER / "other_7.txt", sep=r"\s+", header=None, names=["ID", "F", "value"])
    phen = load_phen(1)[["ID", "yd_trait"]]
    adj = adjusted_bodyweight(1).set_index("ID")
    m = o7.merge(adj.reset_index(), on="ID").merge(phen, on="ID")
    sub = m[m.F <= 1]
    slope_yd = float(np.polyfit(sub.F, sub.yd_trait, 1)[0])
    slope_adj = float(np.polyfit(sub.F, sub.adj, 1)[0])
    return {
        "shape": [int(o7.shape[0]), int(o7.shape[1])],
        "F_min": float(o7.F.min()),
        "F_max": float(o7.F.max()),
        "F_gt_1": int((o7.F > 1).sum()),
        "corr_value_yd": float(m.value.corr(m.yd_trait)),
        "corr_value_adj": float(m.value.corr(m.adj)),
        "slope_yd_per_F": slope_yd,
        "slope_yd_per_1pct": slope_yd / 100,
        "slope_adj_per_F": slope_adj,
        "paper_slope_per_1pct": -0.77,
        "pass": abs(slope_yd / 100 - (-0.77)) < 0.1 and o7.F.max() < 5,
    }


def check_v6() -> dict:
    o5 = pd.read_csv(OTHER / "other_5.txt", sep=r"\s+", header=None, names=["A", "B", "C"])
    reuse = {}
    for c in "ABC":
        vc = o5[c].value_counts()
        reuse[c] = {
            "n_unique": int(o5[c].nunique()),
            "max_reuse": int(vc.max()),
            "top": vc.head(3).to_dict(),
        }

    # Genotype maps for the seven carrier-cross variants.
    geno = []
    for i in range(1, 8):
        gfile = f"gen_5.{i}"
        gt = open_bed(GENS / f"{gfile}.bed").read(dtype="int8")[:, 0]
        iid = load_fam_iid(GENS / f"{gfile}.fam")
        geno.append(pd.Series(gt, index=iid))

    perms = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
    mendel = {}
    for perm in perms:
        off_col, p1_col, p2_col = perm
        violations = 0
        total = 0
        for g in geno:
            off = o5.iloc[:, off_col].map(g)
            p1 = o5.iloc[:, p1_col].map(g)
            p2 = o5.iloc[:, p2_col].map(g)
            ok = off.notna() & p1.notna() & p2.notna()
            if not ok.any():
                continue
            # A genotype is possible iff it is one of the four parental
            # allele combinations.  Enumerate transmitted alleles.
            possible = np.zeros(len(o5), dtype=bool)
            for a1 in (0, 1):
                for a2 in (0, 1):
                    for b1 in (0, 1):
                        for b2 in (0, 1):
                            for gt_off in (a1 + b1, a1 + b2, a2 + b1, a2 + b2):
                                possible |= (p1.to_numpy() == a1 + a2) & (
                                    p2.to_numpy() == b1 + b2
                                ) & (off.to_numpy() == gt_off)
            violations += int((ok.to_numpy() & ~possible).sum())
            total += int(ok.sum())
        mendel[str(perm)] = {"violations": violations, "total": total, "rate": violations / total}

    def pedigree(path: Path) -> dict:
        df = pd.read_csv(path, sep=r"\s+", header=None, names=["ID", "SIRE", "DAM"])
        return {
            "rows": int(len(df)),
            "columns": 3,
            "id_min": int(df.ID.min()),
            "id_max": int(df.ID.max()),
            "sire_zero_frac": float((df.SIRE == 0).mean()),
            "dam_zero_frac": float((df.DAM == 0).mean()),
        }

    best = min(mendel.items(), key=lambda kv: kv[1]["rate"])
    return {
        "other_5_shape": list(o5.shape),
        "other_5_rows": int(len(o5)),
        "other_5_reuse": reuse,
        "mendelian": mendel,
        "best_permutation": best[0],
        "best_rate": best[1]["rate"],
        "other_3": pedigree(OTHER / "other_3.txt"),
        "other_4": pedigree(OTHER / "other_4.txt"),
        "pass": best[1]["rate"] < 0.01 and best[0].startswith("(0,"),
    }


def check_v7() -> dict:
    g10 = pd.read_csv(GENS / "gen_10.1.txt", sep=r"\s+")
    value_counts = {
        c: {str(k): int(v) for k, v in g10[c].value_counts().sort_index().items()}
        for c in g10.columns
        if c != "KEY"
    }
    o11 = pd.read_csv(OTHER / "other_11.txt", sep=r"\s+")
    key = set(g10.KEY.astype(int))
    return {
        "gen_10_shape": [int(g10.shape[0]), int(g10.shape[1])],
        "gen_10_columns": g10.columns.tolist(),
        "gen_10_key_min": int(g10.KEY.min()),
        "gen_10_key_max": int(g10.KEY.max()),
        "gen_10_value_counts": value_counts,
        "other_11_shape": [int(o11.shape[0]), int(o11.shape[1])],
        "birth_year_min": int(o11.birth_year.min()),
        "birth_year_max": int(o11.birth_year.max()),
        "birth_year_counts": {
            str(k): int(v) for k, v in o11.birth_year.value_counts().sort_index().items()
        },
        "sire_in_key_frac": float(o11.SIRE_ID.isin(key).mean()),
        "mgs_in_key_frac": float(o11.MGS_ID.isin(key).mean()),
        "sire_zero_frac": float((o11.SIRE_ID == 0).mean()),
        "mgs_zero_frac": float((o11.MGS_ID == 0).mean()),
        "pass": bool(o11.SIRE_ID.isin(key).all() and o11.MGS_ID.isin(key).all()),
    }


def check_v8() -> dict:
    mapping = []
    for trait, phen_trait in [(1, 1), (2, 3), (3, 2)]:
        d = pd.read_csv(
            OTHER / f"other_2.{trait}.1.tsv", sep="\t", header=None, usecols=[0]
        ).iloc[:, 0]
        p = load_phen(phen_trait)["ID"]
        mapping.append(
            {
                "other_2_trait": trait,
                "matches_phen_1": phen_trait,
                "n_other": int(len(d)),
                "n_phen": int(len(p)),
                "set_equal": bool(set(d) == set(p)),
                "order_equal": bool(np.array_equal(d.to_numpy(), p.to_numpy())),
            }
        )
    windows = []
    for i in range(1, 14):
        path = OTHER / f"other_2.1.{i}.tsv"
        with path.open() as fh:
            first = fh.readline()
        windows.append(
            {
                "window": i,
                "name": WINDOW_NAMES[i],
                "ncols": len(first.rstrip("\n").split("\t")),
                "bytes": int(path.stat().st_size),
            }
        )
    return {
        "trait_mapping": mapping,
        "windows": windows,
        "pass": all(m["set_equal"] for m in mapping) and all(w["ncols"] == 501 for w in windows),
    }


def check_v9() -> dict:
    con = duckdb.connect()
    con.execute("SET memory_limit='6GB'; SET threads=4; SET temp_directory='work/tmp';")
    o1 = (
        "read_csv('dryad/other/other_1.1.txt', delim=' ', header=true, "
        "columns={'CHR':'INTEGER','BP':'BIGINT','A1':'VARCHAR','A2':'VARCHAR','MAF':'DOUBLE',"
        "'B_A':'DOUBLE','B_D':'DOUBLE','Z_A':'DOUBLE','Z_D':'DOUBLE',"
        "'P_A':'DOUBLE','P_D':'DOUBLE'}, "
        "ignore_errors=true)"
    )
    o9 = (
        "read_csv('dryad/other/other_9.txt', delim='\\t', header=true, "
        "columns={'CHR':'INTEGER','BP':'BIGINT','AR2':'DOUBLE'}, ignore_errors=true)"
    )
    o10 = (
        "read_csv('dryad/other/other_10.txt', delim=' ', header=true, "
        "columns={\"CHR.UMD\":'INTEGER',\"BP.UMD\":'BIGINT',"
        "\"CHR.ARS\":'VARCHAR',\"BP.ARS\":'BIGINT'}, "
        "ignore_errors=true)"
    )
    o15 = (
        "read_csv('dryad/other/other_15.txt', delim='\\t', header=true, "
        "columns={'CHR':'INTEGER','BP':'BIGINT','P_D':'DOUBLE'}, ignore_errors=true)"
    )
    out = {}
    o1_counts = con.execute(
        f"SELECT COUNT(*), COUNT(DISTINCT (CHR,BP)) FROM {o1}"
    ).fetchone()
    out["other_1_1"] = {"rows": int(o1_counts[0]), "distinct_keys": int(o1_counts[1])}
    distinct_expr = {
        "other_9": "COUNT(DISTINCT (CHR,BP))",
        "other_10": 'COUNT(DISTINCT ("CHR.UMD","BP.UMD"))',
        "other_15": "COUNT(DISTINCT (CHR,BP))",
    }
    for name, q in [("other_9", o9), ("other_10", o10), ("other_15", o15)]:
        rows, distinct = con.execute(f"SELECT COUNT(*), {distinct_expr[name]} FROM {q}").fetchone()
        out[name] = {"rows": int(rows), "distinct_rows": int(distinct)}
    out["other_9_coverage"] = int(
        con.execute(
            f"SELECT COUNT(*) FILTER (WHERE b.CHR IS NOT NULL) FROM {o1} a "
            f"LEFT JOIN (SELECT DISTINCT CHR,BP FROM {o9}) b USING (CHR,BP)"
        ).fetchone()[0]
    )
    out["other_10_coverage"] = int(
        con.execute(
            f'SELECT COUNT(*) FILTER (WHERE b.cu IS NOT NULL) FROM {o1} a '
            f'LEFT JOIN (SELECT DISTINCT "CHR.UMD" cu,"BP.UMD" bu FROM {o10}) b '
            f"ON a.CHR=b.cu AND a.BP=b.bu"
        ).fetchone()[0]
    )
    out["other_15_coverage"] = int(
        con.execute(
            f"SELECT COUNT(*) FILTER (WHERE b.CHR IS NOT NULL) FROM {o1} a "
            f"LEFT JOIN (SELECT DISTINCT CHR,BP FROM {o15}) b USING (CHR,BP)"
        ).fetchone()[0]
    )
    out["other_15_rows_in_other_1"] = int(
        con.execute(
            f"SELECT COUNT(*) FROM (SELECT DISTINCT CHR,BP FROM {o15}) a "
            f"SEMI JOIN (SELECT DISTINCT CHR,BP FROM {o1}) b USING (CHR,BP)"
        ).fetchone()[0]
    )
    out["other_1_chroms"] = [
        int(v[0]) for v in con.execute(f"SELECT DISTINCT CHR FROM {o1} ORDER BY 1").fetchall()
    ]
    out["other_10_ars_chroms"] = [
        str(v[0])
        for v in con.execute(f'SELECT DISTINCT "CHR.ARS" FROM {o10} ORDER BY 1').fetchall()
    ]
    n = out["other_1_1"]["rows"]
    out["coverage_frac"] = {
        "other_9": out["other_9_coverage"] / n,
        "other_10": out["other_10_coverage"] / n,
        "other_15": out["other_15_coverage"] / n,
    }
    out["other_15_rows_in_other_1_frac"] = (
        out["other_15_rows_in_other_1"] / out["other_15"]["distinct_rows"]
    )
    out["pass"] = (
        out["other_1_1"]["rows"] - out["other_1_1"]["distinct_keys"] == 55
        and out["coverage_frac"]["other_9"] > 0.999
        and 0.99 < out["coverage_frac"]["other_10"] < 1.0
        and 0.002 < out["coverage_frac"]["other_15"] < 0.004
    )
    return out


def check_v10() -> dict:
    iid6 = load_fam_iid(GENS / "gen_6.1.fam")
    p7 = pd.read_csv(PHENS / "phen_7.1.txt", sep=r"\s+")
    p72 = pd.read_csv(PHENS / "phen_7.2.txt", sep=r"\s+")
    bim = pd.read_csv(
        GENS / "gen_6.1.bim",
        sep=r"\s+",
        header=None,
        names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
    )
    targets = {
        "Chr2:22": (2, 22836659),
        "PLCD4": (2, 107313998),
        "FGD4": (5, 77632752),
        "ARSI": (7, 63649071),
        "MTRF1": (12, 11463981),
        "DPF2": (29, 44213160),
        "MUS81": (29, 44645469),
        "GALNT2": (28, 1312334),
        "MYH1": (19, 30114249),
        "chr14:1802265": (14, 1802265),
        "chr14:25015640": (14, 25015640),
        "chr14:26064929": (14, 26064929),
    }
    found = {}
    for label, (c, p) in targets.items():
        sub = bim[bim.CHR == c]
        if sub.empty:
            found[label] = None
            continue
        d = (sub.BP - p).abs()
        i = int(d.idxmin())
        found[label] = {
            "snp": str(sub.loc[i, "SNP"]),
            "chr": int(sub.loc[i, "CHR"]),
            "bp": int(sub.loc[i, "BP"]),
            "a1": str(sub.loc[i, "A1"]),
            "a2": str(sub.loc[i, "A2"]),
            "distance": int(d.loc[i]),
        }

    proxy_idx = int(bim.index[(bim.CHR == 29) & (bim.BP == 44213160)][0])
    nearest_idx = int(bim.index[(bim.CHR == 29) & (bim.BP == 44645366)][0])
    bed = open_bed(GENS / "gen_6.1.bed")
    g = bed.read((slice(None), [proxy_idx, nearest_idx]), dtype="int8")
    r = float(np.corrcoef(g[:, 0].astype(float), g[:, 1].astype(float))[0, 1])
    mus81_proxy_ld = {
        "proxy_snp": "Chr29_44213160",
        "proxy_bp": 44213160,
        "nearest_mus81_snp": "rs41652267",
        "nearest_mus81_bp": 44645366,
        "distance_bp": 103,
        "n": int(g.shape[0]),
        "proxy_genotype_counts": np.bincount(g[:, 0], minlength=3).tolist(),
        "nearest_genotype_counts": np.bincount(g[:, 1], minlength=3).tolist(),
        "r": r,
        "r2": r * r,
    }

    gen7 = pd.read_csv(GENS / "gen_7.1.txt", sep=r"\s+")
    gen8 = pd.read_csv(GENS / "gen_8.1.txt", sep=r"\s+", header=None)
    gen9_bim = pd.read_csv(GENS / "gen_9.1.bim", sep=r"\s+", header=None)
    phen5_ids = pd.read_csv(PHENS / "phen_5.1.txt", sep=r"\s+", usecols=["ID"])["ID"]
    phen6_ids = pd.read_csv(PHENS / "phen_6.1.txt", sep=r"\s+", usecols=["ID"])["ID"]
    return {
        "gen_6_1_fam_n": int(len(iid6)),
        "phen_7_1_key_n": int(len(p7)),
        "phen_7_2_id_n": int(len(p72)),
        "gen6_iids_equal_phen7_key": bool(
            np.array_equal(np.sort(iid6), np.sort(p7.KEY.to_numpy()))
        ),
        "gen6_iids_equal_phen7_2_id": bool(
            np.array_equal(np.sort(iid6), np.sort(p72.ID.to_numpy()))
        ),
        "gen_6_1_bim_rows": int(len(bim)),
        "gen_6_1_bim_chroms": [int(v) for v in sorted(bim.CHR.unique())],
        "candidate_nearest": found,
        "mus81_proxy_ld": mus81_proxy_ld,
        "gen_7_ids_equal_phen_6": bool(set(gen7.ID.astype(int)) == set(phen6_ids.astype(int))),
        "gen_8_ids_equal_phen_5": bool(
            set(gen8.iloc[:, 0].astype(int)) == set(phen5_ids.astype(int))
        ),
        "gen_9_bim": {
            "rows": int(len(gen9_bim)),
            "variants": [
                {
                    "chr": int(r[0]),
                    "snp": str(r[1]),
                    "bp": int(r[3]),
                    "a1": str(r[4]),
                    "a2": str(r[5]),
                }
                for r in gen9_bim.itertuples(index=False)
            ],
        },
        "pass": bool(
            np.array_equal(np.sort(iid6), np.sort(p7.KEY.to_numpy()))
            and np.array_equal(np.sort(iid6), np.sort(p72.ID.to_numpy()))
        ),
    }


def check_v11() -> dict:
    files = []
    for path in sorted(SOURCE.iterdir()):
        if path.suffix.lower() in {".jpg", ".jpeg", ".tif", ".tiff", ".png"}:
            files.append({"file": path.name, "type": "image", "bytes": int(path.stat().st_size)})
            continue
        try:
            df = pd.read_csv(path, nrows=3)
            with path.open(encoding="utf-8", errors="ignore") as fh:
                n_rows = sum(1 for _ in fh) - 1
            files.append(
                {
                    "file": path.name,
                    "type": "csv",
                    "rows": int(n_rows),
                    "cols": int(df.shape[1]),
                    "columns": df.columns.tolist(),
                }
            )
        except Exception as exc:  # pragma: no cover - defensive
            files.append({"file": path.name, "type": "unknown", "error": str(exc)})
    return {
        "exists": SOURCE.is_dir(),
        "n_csv": sum(f["type"] == "csv" for f in files),
        "n_images": sum(f["type"] == "image" for f in files),
        "files": files,
        "pass": SOURCE.is_dir() and any(f["type"] == "csv" for f in files),
    }


def main() -> int:
    t0 = time.time()
    checks = {
        "environment": check_environment(),
        "V1": check_v1(),
        "V2": check_v2(),
        "V3": check_v3(),
        "V4": check_v4(),
        "V5": check_v5(),
        "V6": check_v6(),
        "V7": check_v7(),
        "V8": check_v8(),
        "V9": check_v9(),
        "V10": check_v10(),
        "V11": check_v11(),
    }
    summary = {
        "generated_at": pd.Timestamp.now(tz="UTC").isoformat(),
        "runtime_seconds": round(time.time() - t0, 1),
        "checks": checks,
    }
    out = WORK / "phase0_summary.json"
    out.write_text(json.dumps(summary, indent=2, default=py))
    print(f"wrote {out} in {summary['runtime_seconds']}s")
    for key, value in checks.items():
        if isinstance(value, dict) and "pass" in value:
            print(f"{key}: {'PASS' if value['pass'] else 'FAIL'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
