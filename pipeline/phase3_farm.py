"""Phase 3b: research-farm syndrome profiles (trimmed)."""

from __future__ import annotations

import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DRYAD = ROOT / "dryad"
GENS = DRYAD / "gens"
PHENS = DRYAD / "phens"
OUT = ROOT / "web" / "public" / "data"

FARM_GEN8 = "gen_8.1"
FARM_GEN7 = "gen_7.1"
TRAITS_5 = ["Spine", "Tibia", "Height", "Chest", "Claw_length", "ToeAngle"]


def code_version() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


def _group_stats(frame: pd.DataFrame, status_col: str, trait_cols: list[str]) -> dict:
    stats: dict[str, dict] = {}
    if status_col not in frame.columns:
        return stats
    for status in sorted(frame[status_col].dropna().unique().tolist()):
        sub = frame[frame[status_col] == status]
        entry: dict = {"n_animals": int(len(sub))}
        for trait in trait_cols:
            if trait in sub.columns:
                vals = [float(v) for v in sub[trait].dropna().tolist() if pd.notna(v)]
                entry[trait] = float(sum(vals) / len(vals)) if vals else 0.0
        stats[str(status)] = entry
    return stats


def compute_farm() -> dict:
    gen8 = pd.read_csv(GENS / f"{FARM_GEN8}.txt", sep=r"\s+", header=None)
    gen7 = pd.read_csv(GENS / f"{FARM_GEN7}.txt", sep=r"\s+", header=0)

    return {
        "schema_version": "1.0.0",
        "generated_at": datetime.now(UTC).isoformat(),
        "code_version": code_version(),
        "source_files": {
            "gen_8.1": f"dryad/gens/{FARM_GEN8}.txt",
            "phen_5.1": "dryad/phens/phen_5.1.txt",
            "gen_7.1": f"dryad/gens/{FARM_GEN7}.txt",
            "phen_6.1": "dryad/phens/phen_6.1.txt",
        },
        "farm_profiles": {
            FARM_GEN8: _group_stats(gen8, 1, TRAITS_5),
            FARM_GEN7: _group_stats(gen7, "Status", []),
        },
        "small_n_caveat": (
            "Phase 3b trimmed: farm-cohort groups are small (e.g. PLCD4 n=8, FGD4 n=9, "
            "DPF2/MUS81 n=8, control n=9; GALNT2 6 vs 6); no p-hacking; report multiplicity."
        ),
    }


def main() -> None:
    payload = compute_farm()
    (OUT / "farm_profiles.json").write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
