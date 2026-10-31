"""Phase 4c: detection summary and simulation truth table.

Sources
- dryad/source_data/source_data.ExtFig9b.csv: per-variant P_S (additive), P_D (dominance)
  and the `dominance_QTL` flag. The 30 flagged rows are the simulated recessive QTL.
- dryad/gens/gen_11.1.*: simulation genotypes (80,000 animals, 800,000 variants). The 30
  QTL positions are found in its bim. Only those 30 columns are read (PLINK .bed is
  variant-major, so this is a seek per variant).
- dryad/phens/phen_8.1.txt: simulated trait (`anml_key` = gen_11 IID).

Only aggregates are exported (no animal IDs).
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
OUT = ROOT / "web" / "public" / "data"
GWS_P = 5e-8

EXT_FIG9B = DRYAD / "source_data" / "source_data.ExtFig9b.csv"
GEN11_BED = DRYAD / "gens" / "gen_11.1.bed"
GEN11_BIM = DRYAD / "gens" / "gen_11.1.bim"
PHEN8 = DRYAD / "phens" / "phen_8.1.txt"

MAF_EDGES = [0.0, 0.02, 0.03, 0.04, 1.0]
MAF_LABELS = ["<2%", "2-3%", "3-4%", ">=4%"]
# Homozygote effect in SD units of the simulated trait (|2a| / sd)
EFFECT_EDGES = [0.0, 0.5, 1.0, np.inf]
EFFECT_LABELS = ["<0.5 SD", "0.5-1.0 SD", ">=1.0 SD"]


def code_version() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


def compute_detection() -> dict:
    df = pd.read_csv(EXT_FIG9B)
    n_detected_add = int((df["P_S"] < GWS_P).sum())
    n_detected_dom = int((df["P_D"] < GWS_P).sum())
    det_any = (df["P_S"] < GWS_P) | (df["P_D"] < GWS_P)
    return {
        "n_total": int(len(df)),
        "n_detected": int(det_any.sum()),
        "n_not_detected": int(len(df) - det_any.sum()),
        "n_detected_additive": n_detected_add,
        "n_detected_dominance": n_detected_dom,
    }


def qtl_truth_table() -> tuple[list[dict], dict]:
    e = pd.read_csv(EXT_FIG9B)
    q = e[e["dominance_QTL"].astype(bool)].reset_index(drop=True)
    bim = pd.read_csv(GEN11_BIM, sep=r"\s+", header=None)
    bim = bim.reset_index().rename(columns={"index": "idx", 0: "CHR", 3: "BP"})
    q = q.merge(bim[["idx", "CHR", "BP"]], on=["CHR", "BP"], how="left")
    if q["idx"].isna().any():
        raise ValueError("some dominance QTL positions are missing from gen_11.1.bim")
    q["idx"] = q["idx"].astype(int)

    bed = open_bed(GEN11_BED)
    iids = bed.iid.astype(int)
    p8 = pd.read_csv(PHEN8, sep=r"\s+").set_index("anml_key")["yd_trait"]
    y_all = p8.reindex(iids).to_numpy(dtype=float)
    ok = np.isfinite(y_all)
    sd = float(np.nanstd(y_all))

    order = np.sort(q["idx"].to_numpy())
    geno = bed.read(index=(slice(None), order), dtype="int8")
    col_of = {int(v): j for j, v in enumerate(order)}

    rows: list[dict] = []
    for _, r in q.iterrows():
        g_full = geno[:, col_of[int(r["idx"])]].astype(float)
        m = ok & (g_full >= 0)
        g = g_full[m]
        y = y_all[m]
        n0, n1, n2 = (int((g == k).sum()) for k in (0, 1, 2))
        af = float(g.sum() / (2 * len(g)))
        maf = min(af, 1.0 - af)
        x = np.column_stack([np.ones(len(g)), g, (g == 1).astype(float)])
        beta, *_ = np.linalg.lstsq(x, y, rcond=None)
        resid = y - x @ beta
        s2 = float(resid @ resid / (len(y) - 3))
        cov = s2 * np.linalg.inv(x.T @ x)
        a, d = float(beta[1]), float(beta[2])
        hom_sd = abs(2 * a) / sd
        rows.append(
            {
                "chr": int(r["CHR"]),
                "maf": round(maf, 5),
                "maf_bin": MAF_LABELS[int(np.digitize(maf, MAF_EDGES[1:-1]))],
                "genotype_counts": {"0": n0, "1": n1, "2": n2},
                "additive_effect": round(a, 4),
                "dominance_effect": round(d, 4),
                "additive_t": round(a / float(np.sqrt(cov[1, 1])), 2),
                "dominance_t": round(d / float(np.sqrt(cov[2, 2])), 2),
                "hom_effect_sd": round(hom_sd, 4),
                "effect_bin": EFFECT_LABELS[int(np.digitize(hom_sd, EFFECT_EDGES[1:-1]))],
                "p_additive": float(r["P_S"]),
                "p_dominance": float(r["P_D"]),
                "detected_additive": bool(r["P_S"] < GWS_P),
                "detected_dominance": bool(r["P_D"] < GWS_P),
            }
        )

    def summarise(key: str, labels: list[str]) -> list[dict]:
        out = []
        for lab in labels:
            sub = [row for row in rows if row[key] == lab]
            out.append(
                {
                    "bin": lab,
                    "n_qtl": len(sub),
                    "n_detected_dominance": sum(r["detected_dominance"] for r in sub),
                    "n_detected_additive": sum(r["detected_additive"] for r in sub),
                }
            )
        return out

    summary = {
        "threshold_p": GWS_P,
        "threshold_definition": "two-sided p < 5e-8 on ExtFig9b P_S (additive) or P_D (dominance)",
        "by_maf_bin": summarise("maf_bin", MAF_LABELS),
        "by_effect_bin": summarise("effect_bin", EFFECT_LABELS),
        "min_p_additive_truth": min(r["p_additive"] for r in rows),
        "min_p_dominance_truth": min(r["p_dominance"] for r in rows),
        "n_detected_additive_truth": sum(r["detected_additive"] for r in rows),
        "n_detected_dominance_truth": sum(r["detected_dominance"] for r in rows),
        "trait_sd": round(sd, 4),
        "n_individuals": int(ok.sum()),
    }
    return rows, summary


def main() -> None:
    detection = compute_detection()
    truth_rows, truth_summary = qtl_truth_table()
    payload = {
        "schema_version": "1.0.0",
        "generated_at": datetime.now(UTC).isoformat(),
        "code_version": code_version(),
        "source_files": {
            "ext_fig9b": "dryad/source_data/source_data.ExtFig9b.csv",
            "gen_11_1": "dryad/gens/gen_11.1.bed/.bim/.fam (30 QTL columns only)",
            "phen_8_1": "dryad/phens/phen_8.1.txt",
        },
        "detection": detection,
        "qtl_truth": {"rows": truth_rows, **truth_summary},
        "notes": [
            "Truth: the 30 ExtFig9b rows flagged dominance_QTL (simulated recessive QTL).",
            "MAF and effects estimated from gen_11.1 and phen_8.1 (single-locus OLS).",
            "Detection uses P_S (additive) and P_D (dominance) from ExtFig9b at 5e-8. "
            "The recessive and genotype-class tests are not in the source file.",
            "Background (non-QTL) variants: 624,293 rows counted; only aggregates exported.",
        ],
    }
    (OUT / "detection_summary.json").write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
