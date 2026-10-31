"""Phase 4c: detection summary from source_data.ExtFig9b.csv (trimmed)."""

from __future__ import annotations

import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DRYAD = ROOT / "dryad"
OUT = ROOT / "web" / "public" / "data"
GWS_P = 5e-8

EXT_FIG9B = DRYAD / "source_data" / "source_data.ExtFig9b.csv"


def code_version() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


def compute_detection() -> dict:
    df = pd.read_csv(EXT_FIG9B)
    rows = []
    n_detected_add = 0
    n_detected_dom = 0
    n_detected = 0
    for _, r in df.iterrows():
        p_s = float(r["P_S"])
        p_d = float(r["P_D"])
        dom_qtl = bool(r["dominance_QTL"])
        det_add = p_s < GWS_P
        det_dom = p_d < GWS_P
        det = det_add or det_dom
        n_detected_add += int(det_add)
        n_detected_dom += int(det_dom)
        n_detected += int(det)
        rows.append(
            {
                "chr": int(r["CHR"]),
                "bp": int(r["BP"]),
                "p_s": p_s,
                "p_d": p_d,
                "dominance_qtl": dom_qtl,
                "detected_additive": bool(det_add),
                "detected_dominance": bool(det_dom),
                "detected": bool(det),
            }
        )
    return {
        "schema_version": "1.0.0",
        "generated_at": datetime.now(UTC).isoformat(),
        "code_version": code_version(),
        "source_files": {"ext_fig9b": "dryad/source_data/source_data.ExtFig9b.csv"},
        "detection": {
            "n_total": int(len(df)),
            "n_detected": int(n_detected),
            "n_not_detected": int(len(df) - n_detected),
            "n_detected_additive": int(n_detected_add),
            "n_detected_dominance": int(n_detected_dom),
            "n_rows_exported": int(len(rows)),
            "rows": rows[:200],
        },
        "notes": [
            "Phase 4c trimmed: primary source = source_data.ExtFig9b.csv.",
            "MAF bin x effect size vs detected not computable (no MAF columns).",
            "Rescan of gen_11.1 skipped by decision.",
        ],
    }


def main() -> None:
    payload = compute_detection()
    (OUT / "detection_summary.json").write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
