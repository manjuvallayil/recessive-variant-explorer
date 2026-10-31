"""Phase 4b: syndrome-signal detection with AUROC/AUPRC metrics.

Key facts (see docs/ASSUMPTIONS.md item 5):
- phen_1.* (body weight, stature, BCS) is keyed by gen_1.*.fam IIDs. Labels for phen_1
  outcomes must therefore come from the gen_1 windows. Using gen_3.*.fam (the phen_4 key)
  for phen_1 gives no signal at all.
- Labels are homozygous-mutant status. DPF2 and MUS81 are one locus (hom for either).
- Pre-specified direction: mutant homozygotes are lighter (Phase 1 effects are negative),
  so the body-weight flag is -bw. The direction comes from Phase 1 effect signs, not from
  label-tuned thresholds.
- Labels with fewer than N_PER_LOCUS_MIN positives are not evaluated on their own. They
  are included in the pooled "any" label.
"""

from __future__ import annotations

import json
import subprocess
from datetime import UTC, datetime
from math import erf, sqrt
from pathlib import Path

import numpy as np
import pandas as pd
from bed_reader import open_bed
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
DRYAD = ROOT / "dryad"
PHENS = DRYAD / "phens"
GENS = DRYAD / "gens"
OUT = ROOT / "web" / "public" / "data"
N_PER_LOCUS_MIN = 30
BOOTSTRAP_B = 1000
CV_FOLDS = 5
CV_SEED = 0

# (variant name, gen_1 window, chrom, bp, label group)
VARIANTS = [
    ("PLCD4", "gen_1.2", 2, 107313998, "PLCD4"),
    ("FGD4", "gen_1.3", 5, 77632752, "FGD4"),
    ("GALNT2", "gen_1.8", 28, 1312334, "GALNT2"),
    ("DPF2", "gen_1.9", 29, 44213160, "DPF2/MUS81"),
    ("MUS81", "gen_1.9", 29, 44645469, "DPF2/MUS81"),
]
GROUPS = ["PLCD4", "FGD4", "GALNT2", "DPF2/MUS81"]


def code_version() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


def variant_homozygous(stem: str, chrom: int, bp: int) -> pd.Series:
    """Index: IID (gen_1.fam key). Value 1 if homozygous for A1 (mutant), else 0."""
    bim = pd.read_csv(GENS / f"{stem}.bim", sep=r"\s+", header=None)
    hit = bim.index[(bim[0] == chrom) & (bim[3] == bp)]
    if len(hit) == 0:
        raise ValueError(f"variant {chrom}:{bp} not found in {stem}")
    bed = open_bed(GENS / f"{stem}.bed")
    g = bed.read((slice(None), [int(hit[0])]), dtype="int8")[:, 0]
    iids = pd.read_csv(GENS / f"{stem}.fam", sep=r"\s+", header=None)[1].astype(int).to_numpy()
    keep = g >= 0
    return pd.Series((g[keep] == 2).astype(int), index=iids[keep])


def variant_counts(stem: str, chrom: int, bp: int) -> dict[str, int]:
    """Genotype counts (0/1/2, missing) for one variant in the gen_1 window."""
    bim = pd.read_csv(GENS / f"{stem}.bim", sep=r"\s+", header=None)
    hit = bim.index[(bim[0] == chrom) & (bim[3] == bp)]
    g = open_bed(GENS / f"{stem}.bed").read((slice(None), [int(hit[0])]), dtype="int8")[:, 0]
    return {
        "n0": int((g == 0).sum()),
        "n1": int((g == 1).sum()),
        "n2": int((g == 2).sum()),
        "missing": int((g < 0).sum()),
    }


def _discovery() -> tuple[pd.DataFrame, dict[str, np.ndarray], dict[str, int]]:
    bw = pd.read_csv(PHENS / "phen_1.1.txt", sep=r"\s+").rename(columns={"yd_trait": "bw"})
    sta = pd.read_csv(PHENS / "phen_1.2.txt", sep=r"\s+").rename(columns={"yd_trait": "stature"})
    bcs = pd.read_csv(PHENS / "phen_1.3.txt", sep=r"\s+").rename(columns={"yd_trait": "bcs"})
    frame = (
        bw[["ID", "bw"]]
        .merge(sta[["ID", "stature"]], on="ID", how="inner")
        .merge(bcs[["ID", "bcs"]], on="ID", how="inner")
        .dropna()
        .reset_index(drop=True)
    )
    hom: dict[str, pd.Series] = {}
    for name, stem, chrom, bp, _group in VARIANTS:
        hom[name] = variant_homozygous(stem, chrom, bp)
    labels: dict[str, np.ndarray] = {}
    ids = frame["ID"].astype(int)
    for group in GROUPS:
        members = [name for name, _stem, _chrom, _bp, grp in VARIANTS if grp == group]
        lab = np.zeros(len(frame), dtype=int)
        for name in members:
            lab = lab | ids.map(hom[name]).fillna(0).astype(int).to_numpy()
        labels[group] = lab
    pooled = np.zeros(len(frame), dtype=int)
    for group in GROUPS:
        pooled = pooled | labels[group]
    counts = {g: int(labels[g].sum()) for g in GROUPS}
    counts["any"] = int(pooled.sum())
    frame["label"] = pooled
    return frame, labels, counts


def _boot_ci(y: np.ndarray, score: np.ndarray, b: int = BOOTSTRAP_B) -> dict:
    rng = np.random.default_rng(0)
    aucs: list[float] = []
    auprcs: list[float] = []
    n = len(y)
    for _ in range(b):
        idx = rng.integers(0, n, n)
        yi = y[idx]
        if yi.min() == yi.max():
            continue
        aucs.append(float(roc_auc_score(yi, score[idx])))
        auprcs.append(float(average_precision_score(yi, score[idx])))
    a = np.array(aucs)
    p = np.array(auprcs)
    return {
        "auc": float(roc_auc_score(y, score)),
        "auc_ci": [float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))],
        "auprc": float(average_precision_score(y, score)),
        "auprc_ci": [float(np.percentile(p, 2.5)), float(np.percentile(p, 97.5))],
    }


def _top_precision(y: np.ndarray, score: np.ndarray, frac: float) -> dict:
    k = max(1, int(np.ceil(frac * len(y))))
    top = np.argsort(-score, kind="stable")[:k]
    prec = float(np.mean(y[top]))
    base = float(np.mean(y))
    return {"k": k, "precision": prec, "enrichment": prec / base if base > 0 else None}


def _mahalanobis(X: np.ndarray) -> np.ndarray:
    """Unsupervised: distance from the sample mean in the sample covariance metric."""
    mu = X.mean(axis=0)
    cov = np.cov(X, rowvar=False)
    inv = np.linalg.inv(cov)
    d = X - mu
    return np.sqrt(np.einsum("ij,jk,ik->i", d, inv, d))


def _oof_logistic(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    """5-fold stratified out-of-fold predicted probabilities (scaler fit inside each fold)."""
    oof = np.zeros(len(y), dtype=float)
    skf = StratifiedKFold(n_splits=CV_FOLDS, shuffle=True, random_state=CV_SEED)
    for tr, te in skf.split(X, y):
        model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000))
        model.fit(X[tr], y[tr])
        oof[te] = model.predict_proba(X[te])[:, 1]
    return oof


def compute_syndrome() -> dict:
    frame, labels, counts = _discovery()
    y = frame["label"].to_numpy(dtype=int)
    if y.min() == y.max():
        raise ValueError("syndrome labels must have at least two classes")
    bw = frame["bw"].to_numpy(dtype=float)
    X = frame[["bw", "stature", "bcs"]].to_numpy(dtype=float)

    scores = {
        "B0_body_weight_flag": -bw,
        "B1_mahalanobis": _mahalanobis(X),
        "S_logistic_cv": _oof_logistic(X, y),
    }
    metrics = {name: _boot_ci(y, s) for name, s in scores.items()}
    top = {
        name: {"top1": _top_precision(y, s, 0.01), "top5": _top_precision(y, s, 0.05)}
        for name, s in scores.items()
    }

    # Per-locus statistics use every genotyped record with body weight (no stature/BCS needed),
    # so that records lost to missing traits do not dilute the homozygote counts.
    bw_all = pd.read_csv(PHENS / "phen_1.1.txt", sep=r"\s+")[["ID", "yd_trait"]]
    bw_all = bw_all.rename(columns={"yd_trait": "bw"}).dropna().reset_index(drop=True)
    bw_vals = bw_all["bw"].to_numpy(dtype=float)
    per_group = {}
    any_hom_all = np.zeros(len(bw_all), dtype=bool)
    for group in GROUPS:
        members = [v for v in VARIANTS if v[4] == group]
        hom_all = np.zeros(len(bw_all), dtype=bool)
        hwe_expected = 0.0
        counts_out = {}
        for name, stem, chrom, bp, _grp in members:
            hom_map = variant_homozygous(stem, chrom, bp)
            hom_all |= bw_all["ID"].map(hom_map).fillna(0).astype(bool).to_numpy()
            c = variant_counts(stem, chrom, bp)
            n_tot = c["n0"] + c["n1"] + c["n2"]
            q = (c["n1"] + 2 * c["n2"]) / (2 * n_tot)
            hwe_expected += q * q * n_tot
            counts_out[name] = {**c, "alt_freq": round(q, 5)}
        any_hom_all |= hom_all
        n_hom = int(hom_all.sum())
        d_hom = float(bw_vals[hom_all].mean()) if n_hom else float("nan")
        d_ref = float(bw_vals[~hom_all].mean())
        sd_ref = float(bw_vals[~hom_all].std(ddof=1))
        auc_single = float(roc_auc_score(hom_all.astype(int), -bw_vals))
        rng = np.random.default_rng(CV_SEED)
        boot = []
        yb = hom_all.astype(int)
        for _ in range(BOOTSTRAP_B):
            idx = rng.integers(0, len(yb), len(yb))
            if yb[idx].min() == yb[idx].max():
                continue
            boot.append(float(roc_auc_score(yb[idx], -bw_vals[idx])))
        auc_ci = [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))]
        shift = (d_ref - d_hom) / (sd_ref * sqrt(2.0))
        auc_expected = 0.5 * (1.0 + erf(shift / sqrt(2.0)))
        complete_case_pos = int(labels[group].sum())
        entry = {
            "n_homozygous_with_bodyweight": n_hom,
            "n_homozygous_complete_case": complete_case_pos,
            "n_expected_hwe_homozygous": round(hwe_expected, 1),
            "observed_over_expected": round(n_hom / hwe_expected, 3) if hwe_expected else None,
            "mean_bw_homozygous": round(d_hom, 2),
            "mean_bw_reference": round(d_ref, 2),
            "auc_negated_bw_single_locus": round(auc_single, 4),
            "auc_negated_bw_single_locus_ci": [round(auc_ci[0], 4), round(auc_ci[1], 4)],
            "auc_expected_from_effect_size": round(auc_expected, 4),
            "variants": counts_out,
            "evaluated_in_pooled_multivariate": complete_case_pos >= N_PER_LOCUS_MIN,
        }
        per_group[group] = entry
    per_group_meta = {
        "pooled_any_homozygous_with_bodyweight": int(any_hom_all.sum()),
        "pooled_any_homozygous_complete_case": counts["any"],
        "n_bodyweight_records": int(len(bw_all)),
    }

    return {
        "schema_version": "1.0.0",
        "generated_at": datetime.now(UTC).isoformat(),
        "code_version": code_version(),
        "source_files": {
            "phen_1_1": "dryad/phens/phen_1.1.txt",
            "phen_1_2": "dryad/phens/phen_1.2.txt",
            "phen_1_3": "dryad/phens/phen_1.3.txt",
            "gen_1_windows": "dryad/gens/gen_1.2, gen_1.3, gen_1.8, gen_1.9 (.bed/.bim/.fam)",
        },
        "syndrome_signal": {
            "baseline": "body_weight_alone",
            "baseline_direction": "lower body weight flagged (from Phase 1 effect signs)",
            "n_records": int(len(frame)),
            "n_positive_labels": int(y.sum()),
            "base_rate": float(np.mean(y)),
            "label_counts": counts,
            "scores": metrics,
            "top_precision": top,
        },
        "per_group": per_group,
        "per_group_meta": per_group_meta,
        "notes": [
            "Labels: homozygous-mutant status, gen_1 key (matches phen_1 IIDs).",
            "DPF2 and MUS81 are one locus (homozygous for either).",
            "Per-locus statistics use all genotyped records with body weight (no stature/BCS).",
            "Pooled and multivariate metrics use complete cases (bw, stature, BCS).",
            "Expected AUROC from effect size assumes equal-variance normal classes (shift model).",
            "Logistic model: 5-fold stratified out-of-fold predictions. Sire grouping not "
            "possible (pedigree parents are outside the discovery ID space): leakage risk "
            "from related animals across folds.",
            "Mahalanobis is unsupervised (no labels used).",
            "Not yet done: unsupervised isolation forest / PCA / autoencoder, LightGBM, nested CV.",
        ],
    }


def main() -> None:
    payload = compute_syndrome()
    (OUT / "syndrome_signal.json").write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
