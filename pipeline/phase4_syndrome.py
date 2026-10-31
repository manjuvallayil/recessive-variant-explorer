"""Phase 4b: syndrome-signal detection with AUROC/AUPRC metrics (trimmed)."""

from __future__ import annotations

import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
from bed_reader import open_bed
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
DRYAD = ROOT / "dryad"
PHENS = DRYAD / "phens"
GENS = DRYAD / "gens"
OUT = ROOT / "web" / "public" / "data"
N_PER_LOCUS_MIN = 30
BOOTSTRAP_B = 1000

LOCI = [
    ("PLCD4", "gen_3.2", 2, 107313998),
    ("FGD4", "gen_3.3", 5, 77632752),
    ("GALNT2", "gen_3.10", 28, 1312334),
    ("DPF2/MUS81", "gen_3.11", 29, 44213160),
]


def code_version() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


def load_fam_iids(path: Path) -> np.ndarray:
    return pd.read_csv(path, sep=r"\s+", header=None)[1].astype(int).to_numpy()


def read_bed_variants(path: Path, indices: list[int]) -> np.ndarray:
    bed = open_bed(path)
    arr = bed.read((slice(None), indices), dtype="int8")
    return arr.reshape((arr.shape[0], len(indices)))


def find_bim_index(bim: pd.DataFrame, chrom: int, bp: int) -> int:
    hit = bim.index[(bim[0] == chrom) & (bim[3] == bp)]
    if len(hit) == 0:
        raise ValueError(f"variant {chrom}:{bp} not found")
    return int(hit[0])


def locus_genotypes(locus: str, stem: str, chrom: int, bp: int) -> dict[int, int]:
    bim = pd.read_csv(GENS / f"{stem}.bim", sep=r"\s+", header=None)
    idx = find_bim_index(bim, chrom, bp)
    g = read_bed_variants(GENS / f"{stem}.bed", [idx])[:, 0]
    iids = load_fam_iids(GENS / f"{stem}.fam")
    return {int(iid): int(geno) for iid, geno in zip(iids, g, strict=True)}


def _discovery() -> tuple[np.ndarray, pd.DataFrame]:
    bw_df = pd.read_csv(PHENS / "phen_1.1.txt", sep=r"\s+").rename(columns={"yd_trait": "bw"})
    sta_df = pd.read_csv(PHENS / "phen_1.2.txt", sep=r"\s+").rename(columns={"yd_trait": "stature"})
    bcs_df = pd.read_csv(PHENS / "phen_1.3.txt", sep=r"\s+").rename(columns={"yd_trait": "bcs"})
    frame = (
        bw_df[["ID", "bw"]]
        .merge(sta_df[["ID", "stature"]], on="ID", how="inner")
        .merge(bcs_df[["ID", "bcs"]], on="ID", how="inner")
    )
    frame = frame.dropna().reset_index(drop=True)
    per_locus: dict[str, np.ndarray] = {}
    pooled = np.zeros(len(frame), dtype=int)
    for locus, stem, chrom, bp in LOCI:
        gm = locus_genotypes(locus, stem, chrom, bp)
        lab = np.zeros(len(frame), dtype=int)
        for pos, (_, row) in enumerate(frame.iterrows()):
            iid = int(row["ID"])
            g = gm.get(iid)
            if g is not None and g == 2:
                lab[pos] = 1
        per_locus[locus] = lab
        print(f"{locus}: homozygous-mutant positives = {int(lab.sum())}")
        if int(lab.sum()) >= N_PER_LOCUS_MIN:
            pooled = pooled | lab.astype(bool)
    frame["label"] = pooled.astype(float)
    return frame[["ID", "bw", "stature", "bcs", "label"]].to_numpy(dtype=float), frame, per_locus


def _boot_auc(y: np.ndarray, pred: np.ndarray, b: int = BOOTSTRAP_B) -> tuple[float, float]:
    rng = np.random.default_rng(0)
    vals = []
    for _ in range(b):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            vals.append(float("nan"))
            continue
        vals.append(float(roc_auc_score(y[idx], pred[idx])))
    vals = np.array(vals)
    vals = vals[~np.isnan(vals)]
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def compute_syndrome() -> dict:
    data, frame, per_locus = _discovery()
    bw = data[:, 1]
    stature = data[:, 2]
    bcs = data[:, 3]
    y = data[:, 4]
    if len(np.unique(y)) < 2:
        raise ValueError("syndrome labels must have at least two classes")
    base_auc = float(roc_auc_score(y, bw)) if len(np.unique(y)) > 1 else float("nan")
    base_auprc = float(average_precision_score(y, bw)) if len(np.unique(y)) > 1 else float("nan")
    full_X = np.column_stack([np.ones(len(frame)), bw, stature, bcs])
    clf = LogisticRegression(max_iter=1000).fit(full_X, y)
    full_pred = clf.predict_proba(full_X)[:, 1]
    full_auc = float(roc_auc_score(y, full_pred))
    full_auprc = float(average_precision_score(y, full_pred))

    lo, hi = _boot_auc(y, bw)
    base_auc_ci = [float(lo), float(hi)]
    lo2, hi2 = _boot_auc(y, full_pred)
    full_auc_ci = [float(lo2), float(hi2)]

    order = np.argsort(-full_pred)
    n = len(y)
    top1 = order[: max(1, int(np.ceil(0.01 * n)))]
    top5 = order[: max(1, int(np.ceil(0.05 * n)))]
    prec1 = float(np.mean(y[top1]))
    prec5 = float(np.mean(y[top5]))
    base_rate = float(np.mean(y))

    per_locus_out = {}
    for locus, lab in per_locus.items():
        per_locus_out[locus] = {
            "n_positive": int(lab.sum()),
            "n_used": N_PER_LOCUS_MIN,
            "homozygous_mutant_auc": (
                float(roc_auc_score(y, lab.astype(float))) if lab.sum() >= N_PER_LOCUS_MIN else None
            ),
            "homozygous_mutant_auprc": (
                float(average_precision_score(y, lab.astype(float)))
                if lab.sum() >= N_PER_LOCUS_MIN
                else None
            ),
        }

    return {
        "schema_version": "1.0.0",
        "generated_at": datetime.now(UTC).isoformat(),
        "code_version": code_version(),
        "source_files": {
            "phen_1_1": "dryad/phens/phen_1.1.txt",
            "phen_1_2": "dryad/phens/phen_1.2.txt",
            "phen_1_3": "dryad/phens/phen_1.3.txt",
            "gen_3": "dryad/gens/gen_3.*.bed",
        },
        "syndrome_signal": {
            "baseline": "body_weight_alone",
            "baseline_auc": base_auc,
            "baseline_auprc": base_auprc,
            "baseline_auc_ci": base_auc_ci,
            "baseline_auprc_ci": None,
            "full_auc": full_auc,
            "full_auprc": full_auprc,
            "full_auc_ci": full_auc_ci,
            "full_auprc_ci": None,
            "base_rate": base_rate,
            "top1_precision": float(prec1),
            "top5_precision": float(prec5),
            "n_positive_labels": int(y.sum()),
            "n_records": int(len(frame)),
        },
        "per_locus": per_locus_out,
        "notes": [
            "Phase 4b trimmed: B0 body weight alone as flag; label from homozygous-mutant status.",
            "Full multi-trait logistic on bw/stature/BCS; paper-style negative results pending.",
        ],
    }


def main() -> None:
    payload = compute_syndrome()
    (OUT / "syndrome_signal.json").write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
