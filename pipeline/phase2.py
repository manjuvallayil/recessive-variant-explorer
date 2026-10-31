"""Phase 2: genome-wide additive vs dominance GWAS, resolution, regional LD, and imputation.

Reads the zstd Parquet files under ``work/parquet/`` (creating them once from
``dryad/other/`` via DuckDB if missing) and the PLINK windows ``dryad/gens/gen_1.*``,
and writes aggregated JSON files to ``web/public/data/``.
"""

from __future__ import annotations

import gzip
import json
import math
import subprocess
import time
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from bed_reader import open_bed
from scipy.special import log_ndtr

ROOT = Path(__file__).resolve().parents[1]
DRYAD = ROOT / "dryad"
GENS = DRYAD / "gens"
OTHER = DRYAD / "other"
SOURCE = DRYAD / "source_data"
WORK = ROOT / "work"
PARQUET = WORK / "parquet"
TMP = WORK / "tmp"
OUT = ROOT / "web" / "public" / "data"
REGIONAL_OUT = OUT / "regional"

PARQUET.mkdir(parents=True, exist_ok=True)
TMP.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)
REGIONAL_OUT.mkdir(parents=True, exist_ok=True)

SCHEMA_VERSION = "1.0.0"
GWS_P = 5e-8
GWS_NEG_LOG10_P = -math.log10(GWS_P)  # 7.301029995663981

PARQUET_SPECS = [
    (
        "other_1.1",
        "other_1.1.txt",
        " ",
        (
            "{'CHR':'INTEGER','BP':'BIGINT','A1':'VARCHAR','A2':'VARCHAR','MAF':'DOUBLE',"
            "'B_A':'DOUBLE','B_D':'DOUBLE','Z_A':'DOUBLE','Z_D':'DOUBLE',"
            "'P_A':'DOUBLE','P_D':'DOUBLE'}"
        ),
    ),
    (
        "other_1.2",
        "other_1.2.txt",
        " ",
        (
            "{'CHR':'INTEGER','BP':'BIGINT','A1':'VARCHAR','A2':'VARCHAR','MAF':'DOUBLE',"
            "'B_A':'DOUBLE','B_D':'DOUBLE','Z_A':'DOUBLE','Z_D':'DOUBLE',"
            "'P_A':'DOUBLE','P_D':'DOUBLE'}"
        ),
    ),
    (
        "other_1.3",
        "other_1.3.txt",
        " ",
        (
            "{'CHR':'INTEGER','BP':'BIGINT','A1':'VARCHAR','A2':'VARCHAR','MAF':'DOUBLE',"
            "'B_A':'DOUBLE','B_D':'DOUBLE','Z_A':'DOUBLE','Z_D':'DOUBLE',"
            "'P_A':'DOUBLE','P_D':'DOUBLE'}"
        ),
    ),
    (
        "other_8.1",
        "other_8.1.txt",
        " ",
        (
            "{'CHR':'INTEGER','BP':'BIGINT','A1':'VARCHAR','A2':'VARCHAR','MAF':'DOUBLE',"
            "'B_A':'DOUBLE','Z_A':'DOUBLE','P_A':'VARCHAR'}"
        ),
    ),
    (
        "other_8.2",
        "other_8.2.txt",
        " ",
        (
            "{'CHR':'INTEGER','BP':'BIGINT','A1':'VARCHAR','A2':'VARCHAR','MAF':'DOUBLE',"
            "'B_A':'DOUBLE','Z_A':'DOUBLE','P_A':'VARCHAR'}"
        ),
    ),
    (
        "other_8.3",
        "other_8.3.txt",
        " ",
        (
            "{'CHR':'INTEGER','BP':'BIGINT','A1':'VARCHAR','A2':'VARCHAR','MAF':'DOUBLE',"
            "'B_A':'DOUBLE','Z_A':'DOUBLE','P_A':'VARCHAR'}"
        ),
    ),
    (
        "other_9",
        "other_9.txt",
        "\t",
        "{'CHR':'INTEGER','BP':'BIGINT','AR2':'DOUBLE'}",
    ),
    (
        "other_10",
        "other_10.txt",
        " ",
        "{'CHR.UMD':'INTEGER','BP.UMD':'BIGINT','CHR.ARS':'VARCHAR','BP.ARS':'BIGINT'}",
    ),
    (
        "other_15",
        "other_15.txt",
        "\t",
        "{'CHR':'INTEGER','BP':'BIGINT','P_D':'DOUBLE'}",
    ),
]

# The 8 non-additive body-weight QTL from Table 1 / Fig. 1 / ExtFig 2, plus MYH1 (BCS suggestive).
NON_ADDITIVE_QTLS = [
    {
        "qtl": "Chr2:22",
        "gene": None,
        "chr": 2,
        "candidate_bp": 22836659,
        "window_start": 21836659,
        "window_end": 23836659,
        "is_primary_8": True,
        "is_major_recessive_syndrome": False,
        "category": "modest_recessive_qtl",
        "gen1": "gen_1.1",
        "slug": "chr2_22",
    },
    {
        "qtl": "PLCD4",
        "gene": "PLCD4",
        "chr": 2,
        "candidate_bp": 107313998,
        "window_start": 106313998,
        "window_end": 108313998,
        "is_primary_8": True,
        "is_major_recessive_syndrome": True,
        "category": "novel_recessive_syndrome",
        "gen1": "gen_1.2",
        "slug": "plcd4",
    },
    {
        "qtl": "FGD4",
        "gene": "FGD4",
        "chr": 5,
        "candidate_bp": 77632752,
        "window_start": 76632752,
        "window_end": 78632752,
        "is_primary_8": True,
        "is_major_recessive_syndrome": True,
        "category": "novel_recessive_syndrome",
        "gen1": "gen_1.3",
        "slug": "fgd4",
    },
    {
        "qtl": "ARSI",
        "gene": "ARSI",
        "chr": 7,
        "candidate_bp": 63649071,
        "window_start": 62649071,
        "window_end": 64649071,
        "is_primary_8": True,
        "is_major_recessive_syndrome": False,
        "category": "known_stature_recessive_qtl",
        "gen1": "gen_1.4",
        "slug": "arsi",
    },
    {
        "qtl": "MTRF1",
        "gene": "MTRF1",
        "chr": 12,
        "candidate_bp": 11463981,
        "window_start": 10463981,
        "window_end": 12463981,
        "is_primary_8": True,
        "is_major_recessive_syndrome": True,
        "category": "novel_recessive_syndrome",
        "gen1": "gen_1.5",
        "slug": "mtrf1",
    },
    {
        "qtl": "PLAG1",
        "gene": "PLAG1",
        "chr": 14,
        "candidate_bp": 26064929,
        "causal_bp": 25015640,
        "window_start": 24015640,
        "window_end": 27064929,
        "is_primary_8": True,
        "is_major_recessive_syndrome": False,
        "category": "partial_dominance_additive_locus",
        "gen1": "gen_1.6",
        "slug": "plag1",
    },
    {
        "qtl": "GALNT2",
        "gene": "GALNT2",
        "chr": 28,
        "candidate_bp": 1312334,
        "window_start": 312334,
        "window_end": 2312334,
        "is_primary_8": True,
        "is_major_recessive_syndrome": True,
        "category": "novel_recessive_syndrome",
        "gen1": "gen_1.8",
        "slug": "galnt2",
    },
    {
        "qtl": "DPF2/MUS81",
        "gene": "DPF2/MUS81",
        "chr": 29,
        "candidate_bp": 44213160,
        "secondary_candidate_bp": 44645469,
        "window_start": 43213160,
        "window_end": 45645469,
        "is_primary_8": True,
        "is_major_recessive_syndrome": True,
        "category": "novel_recessive_syndrome",
        "gen1": "gen_1.9",
        "slug": "dpf2",
    },
    {
        "qtl": "MYH1",
        "gene": "MYH1",
        "chr": 19,
        "candidate_bp": 30114249,
        "window_start": 29114249,
        "window_end": 31114249,
        "is_primary_8": False,
        "is_major_recessive_syndrome": True,
        "category": "bcs_suggestive_recessive_syndrome",
        "gen1": "gen_1.7",
        "slug": "myh1",
    },
]

REGIONAL_LOCI = [
    ("chr2_22", "Chr2:22", None, 2, 22836659, "gen_1.1", "rs211074296"),
    ("plcd4", "PLCD4", "PLCD4", 2, 107313998, "gen_1.2", "g.107313998G>A"),
    ("fgd4", "FGD4", "FGD4", 5, 77632752, "gen_1.3", "g.77632752C>T"),
    ("arsi", "ARSI", "ARSI", 7, 63649071, "gen_1.4", "g.63649071C>T"),
    ("mtrf1", "MTRF1", "MTRF1", 12, 11463981, "gen_1.5", "g.11463981G>A"),
    (
        "plag1",
        "chr14 QTL (PLAG1 regulatory domain; gene ~263 kb outside the genotyped window)",
        None,
        14,
        26064929,
        "gen_1.6",
        "rs109686357 (dominance GWAS lead in PLAG1 regulatory domain)",
    ),
    ("myh1", "MYH1", "MYH1", 19, 30114249, "gen_1.7", "g.30114249C>CCCACCTGGG (proxy)"),
    ("galnt2", "GALNT2", "GALNT2", 28, 1312334, "gen_1.8", "g.1312334G>A"),
    ("dpf2", "DPF2", "DPF2", 29, 44213160, "gen_1.9", "g.44213160A>G"),
    ("mus81", "MUS81", "MUS81", 29, 44645469, "gen_1.9", "g.44645469G>T"),
]

# Candidate mutations, comparators, and window-best lead variants that must never be dropped
# by spatial bin thinning in manhattan_dom.json or manhattan_add.json.
MUST_KEEP_VARIANTS = [
    # 9 Table 1 candidate / lead positions
    (2, 22836659),
    (2, 107313998),
    (5, 77632752),
    (7, 63649071),
    (12, 11463981),
    (19, 30114249),
    (28, 1312334),
    (29, 44213160),
    (29, 44645469),
    # 3 BTA14 comparators (DGAT1, PLAG1 causal, PLAG1 dominance lead)
    (14, 1802265),
    (14, 25015640),
    (14, 26064929),
    # Window-best dominance and additive leads across bodyweight and stature
    (2, 21999730),
    (2, 23118213),
    (2, 108131543),
    (5, 77241202),
    (5, 106266836),
    (5, 106250171),
    (5, 106250892),
    (5, 106264863),
    (5, 106289128),
    (7, 62891843),
    (7, 62718020),
    (7, 64626068),
    (12, 11415931),
    (12, 11159910),
    (14, 24974811),
    (14, 24973953),
    (14, 24991209),
    (14, 25307116),
    (14, 25474773),
    (19, 30805217),
    (19, 30592432),
    (19, 29149312),
    (28, 680910),
    (28, 334652),
    (29, 44143856),
    (29, 45056385),
]


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


def neg_log10_p_from_z(z: np.ndarray | float) -> np.ndarray | float:
    """Compute exact two-sided -log10(P) from Z using scipy.special.log_ndtr.

    Never underflows or clips, even for |Z| > 40 (e.g. PLAG1 P ~ 5e-442).
    """
    arr = np.asarray(z, dtype=np.float64)
    out = -(math.log(2.0) + log_ndtr(-np.abs(arr))) / math.log(10.0)
    out = np.where(np.isfinite(out) & (out >= 0.0), out, 0.0)
    if np.ndim(z) == 0:
        return float(out)
    return out


def jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        v = float(obj)
        return v if math.isfinite(v) else None
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if pd.isna(obj):
        return None
    return obj


def write_json(
    rel_path: str, payload: dict, source_files: list[str], compact: bool = False
) -> dict:
    body = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": now_iso(),
        "code_version": CODE_VERSION,
        "source_files": source_files,
        **payload,
    }
    body = jsonable(body)
    path = OUT / rel_path
    path.parent.mkdir(parents=True, exist_ok=True)
    text = (
        json.dumps(body, separators=(",", ":")) + "\n"
        if compact
        else json.dumps(body, indent=2) + "\n"
    )
    raw_bytes = text.encode("utf-8")
    path.write_bytes(raw_bytes)
    gz_bytes = len(gzip.compress(raw_bytes, compresslevel=6))
    return {
        "file": rel_path,
        "bytes": len(raw_bytes),
        "gzip_bytes": gz_bytes,
        "source_files": source_files,
    }


def connect_duckdb() -> duckdb.DuckDBPyConnection:
    con = duckdb.connect()
    con.execute(
        f"SET memory_limit='8GB'; SET threads=6; SET temp_directory='{TMP.as_posix()}';"
    )
    return con


def ensure_parquet_files(con: duckdb.DuckDBPyConnection) -> None:
    # Rebuild other_8.1.parquet if it was created before '<1e-308' string handling (506 rows short)
    o81_pq = PARQUET / "other_8.1.parquet"
    if o81_pq.exists():
        n_o81 = con.execute(f"SELECT COUNT(*) FROM '{o81_pq.as_posix()}'").fetchone()[0]
        if n_o81 < 16_128_364:
            o81_pq.unlink()
            cov_cache = WORK / "ars_coverage_cache.json"
            if cov_cache.exists():
                cov_cache.unlink()

    for stem, fname, delim, cols in PARQUET_SPECS:
        out_path = PARQUET / f"{stem}.parquet"
        if out_path.exists():
            continue
        src = OTHER / fname
        if stem.startswith("other_8."):
            select_expr = (
                "SELECT CHR, BP, A1, A2, MAF, B_A, Z_A, "
                "TRY_CAST(REPLACE(P_A, '<', '') AS DOUBLE) AS P_A"
            )
        else:
            select_expr = "SELECT *"
        con.execute(
            f"COPY ({select_expr} FROM read_csv({src.as_posix()!r}, delim={delim!r}, "
            f"header=true, columns={cols}, ignore_errors=true)) "
            f"TO {out_path.as_posix()!r} (FORMAT PARQUET, COMPRESSION ZSTD)"
        )


def parquet_path(stem: str) -> str:
    return (PARQUET / f"{stem}.parquet").as_posix()


def read_bed_chunked_r2(bed_path: Path, lead_idx: int, chunk_size: int = 2000) -> np.ndarray:
    """Compute Pearson R^2 between lead_idx and all variants in bed_path in chunks."""
    bed = open_bed(bed_path)
    n_variants = bed.sid_count
    lead = bed.read((slice(None), [lead_idx]), dtype="int8")[:, 0].astype(np.float64)
    lead_valid = lead >= 0
    lead_mean = lead[lead_valid].mean()
    lead_c = np.where(lead_valid, lead - lead_mean, 0.0)
    lead_ss = float(np.sum(lead_c**2))

    r2 = np.zeros(n_variants, dtype=np.float64)
    for start in range(0, n_variants, chunk_size):
        end = min(start + chunk_size, n_variants)
        chunk = bed.read((slice(None), slice(start, end)), dtype="int8").astype(np.float64)
        valid = (chunk >= 0) & lead_valid[:, None]
        counts = np.maximum(valid.sum(axis=0), 1)
        means = np.where(valid, chunk, 0.0).sum(axis=0) / counts
        centered = np.where(valid, chunk - means[None, :], 0.0)
        ss = np.sum(centered**2, axis=0)
        num = np.sum(centered * lead_c[:, None], axis=0)
        denom = np.sqrt(lead_ss * ss)
        r = np.divide(num, denom, out=np.zeros_like(num), where=denom > 0)
        r2[start:end] = np.clip(r**2, 0.0, 1.0)
    r2[lead_idx] = 1.0
    return r2


def thin_manhattan_df(df: pd.DataFrame, max_points: int = 24_000) -> tuple[pd.DataFrame, dict]:
    """Thin genome-wide summary stats to <= max_points per trait (<= 48k per model file).

    Rule:
    1. Compute exact -log10(P) from Z via scipy.special.log_ndtr (no clipping).
    2. Keep all strong hits (-log10(P) >= 5.0), capped at top 3 per 20-kb bin (CHR, BP // 20_000).
    3. For background (-log10(P) < 5.0), first take the minimum-P (max |Z|) variant per 25-kb bin
       (CHR, BP // 25_000). Keep all 25-kb bin winners with -log10(P) >= 2.5, and for the null
       baseline (-log10(P) < 2.5) keep the maximum-|Z| winner per 125-kb group of 25-kb winners.
    """
    n_raw = int(len(df))
    z = df["Z"].to_numpy(dtype=np.float64)
    nlp = neg_log10_p_from_z(z)
    df = df.assign(nlp=nlp, abs_z=np.abs(z))

    must_keep_set = set(MUST_KEEP_VARIANTS)
    df["must_keep"] = [
        (int(c), int(b)) in must_keep_set
        for c, b in zip(df["CHR"].to_numpy(), df["BP"].to_numpy(), strict=False)
    ]
    must_keep_df = df[df["must_keep"]].copy()

    strong = df[df["nlp"] >= 5.0].copy()
    strong["bin20k"] = strong["BP"] // 20_000
    strong_kept = (
        strong.sort_values(["CHR", "bin20k", "abs_z"], ascending=[True, True, False])
        .groupby(["CHR", "bin20k"], as_index=False)
        .head(3)
        .drop(columns=["bin20k"])
    )

    bg = df[df["nlp"] < 5.0].copy()
    bg["bin25k"] = bg["BP"] // 25_000
    bg25 = (
        bg.sort_values("abs_z", ascending=False)
        .drop_duplicates(subset=["CHR", "bin25k"], keep="first")
        .drop(columns=["bin25k"])
    )
    bg_mid = bg25[bg25["nlp"] >= 2.5]
    bg_low = bg25[bg25["nlp"] < 2.5].copy()
    bg_low["bin125k"] = bg_low["BP"] // 125_000
    bg_low_kept = (
        bg_low.sort_values("abs_z", ascending=False)
        .drop_duplicates(subset=["CHR", "bin125k"], keep="first")
        .drop(columns=["bin125k"])
    )

    combined = (
        pd.concat([must_keep_df, strong_kept, bg_mid, bg_low_kept], ignore_index=True)
        .sort_values(["must_keep", "abs_z"], ascending=[False, False])
        .drop_duplicates(subset=["CHR", "BP"], keep="first")
    )
    if len(combined) > max_points:
        combined = combined.sort_values(["must_keep", "nlp"], ascending=[False, False]).head(
            max_points
        )
    combined = combined.sort_values(["CHR", "BP"]).reset_index(drop=True)

    meta = {
        "n_raw_variants": n_raw,
        "n_thinned_points": int(len(combined)),
        "n_strong_hits_kept": int(len(strong_kept)),
        "n_moderate_25kb_kept": int(len(bg_mid)),
        "n_baseline_kept": int(len(bg_low_kept)),
        "max_points_cap": max_points,
        "max_neg_log10_p": float(combined["nlp"].max()) if len(combined) else 0.0,
        "n_gws_raw": int((nlp >= GWS_NEG_LOG10_P).sum()),
        "rule": (
            "Exact -log10(P) from Z via scipy.special.log_ndtr (unclipped); "
            "all hits with -log10(P)>=5.0 capped at top 3 per 20-kb bin; "
            "background minimum-P per 25-kb bin for 2.5<=-log10(P)<5.0 and "
            "per 125-kb group of 25-kb winners for -log10(P)<2.5; <=24k points/trait "
            "(<=48k points/model, <500 KB gzipped)."
        ),
    }
    return combined, meta


def load_thinned_candidates_from_duckdb(
    con: duckdb.DuckDBPyConnection, pq_file: str, z_col: str, extra_z_col: str | None = None
) -> tuple[pd.DataFrame, dict]:
    """Use DuckDB hash aggregation for 25-kb bins and window functions only on strong hits."""
    total_rows, n_gws_raw = con.execute(
        f"""
        SELECT
            COUNT(*),
            COUNT(*) FILTER (WHERE ABS({z_col}) >= 5.45131)
        FROM '{pq_file}'
        """
    ).fetchone()

    extra_inner = f", a.{extra_z_col} AS Z_EXTRA" if extra_z_col else ""
    extra_outer = ", Z_EXTRA" if extra_z_col else ""
    extra_argmax = f", ARG_MAX(a.{extra_z_col}, ABS(a.{z_col})) AS Z_EXTRA" if extra_z_col else ""

    must_keep_values = ", ".join(f"({c}, {b})" for c, b in MUST_KEEP_VARIANTS)

    query = f"""
    WITH must_keep_keys(CHR, BP) AS (
        VALUES {must_keep_values}
    ),
    must_keep AS (
        SELECT a.CHR, a.BP, a.{z_col} AS Z {extra_inner}
        FROM '{pq_file}' a
        SEMI JOIN must_keep_keys k USING (CHR, BP)
    ),
    strong AS (
        SELECT CHR, BP, {z_col} AS Z {extra_outer}
        FROM (
            SELECT a.CHR, a.BP, a.{z_col} {extra_inner},
                   ROW_NUMBER() OVER (
                       PARTITION BY a.CHR, a.BP // 20000
                       ORDER BY ABS(a.{z_col}) DESC
                   ) AS rn
            FROM '{pq_file}' a
            WHERE ABS(a.{z_col}) >= 4.417
        )
        WHERE rn <= 3
    ),
    bg25 AS (
        SELECT
            a.CHR,
            ARG_MAX(a.BP, ABS(a.{z_col})) AS BP,
            ARG_MAX(a.{z_col}, ABS(a.{z_col})) AS Z
            {extra_argmax}
        FROM '{pq_file}' a
        GROUP BY a.CHR, a.BP // 25000
    ),
    unioned AS (
        SELECT * FROM must_keep
        UNION ALL
        SELECT * FROM strong
        UNION ALL
        SELECT * FROM bg25
    )
    SELECT CHR, BP, FIRST(Z) AS Z
           {', FIRST(Z_EXTRA) AS Z_EXTRA' if extra_z_col else ''}
    FROM unioned
    GROUP BY CHR, BP
    """
    df = con.execute(query).df()
    thinned, meta = thin_manhattan_df(df, max_points=24_000)

    o10 = parquet_path("other_10")
    con.register("thinned_pts", thinned[["CHR", "BP"]])
    mapped = con.execute(
        f"""
        SELECT t.CHR, t.BP,
               FIRST(m."CHR.ARS") AS CHR_ARS,
               FIRST(m."BP.ARS") AS BP_ARS
        FROM thinned_pts t
        LEFT JOIN '{o10}' m ON t.CHR = m."CHR.UMD" AND t.BP = m."BP.UMD"
        GROUP BY t.CHR, t.BP
        """
    ).df()
    con.unregister("thinned_pts")
    thinned = thinned.merge(mapped, on=["CHR", "BP"], how="left").sort_values(
        ["CHR", "BP"]
    ).reset_index(drop=True)

    cov_cache_path = WORK / "ars_coverage_cache.json"
    cov_cache = json.loads(cov_cache_path.read_text()) if cov_cache_path.exists() else {}
    cache_key = Path(pq_file).name
    if cache_key not in cov_cache:
        unmapped_raw = int(
            con.execute(
                f"""
                SELECT COUNT(*)
                FROM '{pq_file}' a
                ANTI JOIN '{o10}' m ON a.CHR = m."CHR.UMD" AND a.BP = m."BP.UMD"
                """
            ).fetchone()[0]
        )
        cov_cache[cache_key] = unmapped_raw
        cov_cache_path.write_text(json.dumps(cov_cache, indent=2))
    else:
        unmapped_raw = int(cov_cache[cache_key])

    meta["n_raw_variants"] = int(total_rows)
    meta["n_gws_raw"] = int(n_gws_raw)
    meta["raw_unmapped_ars_count"] = int(unmapped_raw)
    meta["raw_unmapped_ars_fraction"] = float(unmapped_raw / total_rows) if total_rows else 0.0

    ars_chr_num = pd.to_numeric(thinned["CHR_ARS"], errors="coerce").fillna(0).astype(int)
    ars_bp_num = thinned["BP_ARS"].fillna(0).astype(int)
    unmapped_thinned = int((ars_bp_num == 0).sum())
    meta["thinned_unmapped_ars_count"] = unmapped_thinned
    meta["thinned_unmapped_ars_fraction"] = (
        float(unmapped_thinned / len(thinned)) if len(thinned) else 0.0
    )

    columns = {
        "chr": thinned["CHR"].astype(int).tolist(),
        "bp_umd": thinned["BP"].astype(int).tolist(),
        "chr_ars": ars_chr_num.tolist(),
        "bp_ars": ars_bp_num.tolist(),
        "neg_log10_p_x100": np.rint(thinned["nlp"].to_numpy() * 100.0).astype(int).tolist(),
    }
    if extra_z_col:
        nlp_extra = neg_log10_p_from_z(thinned["Z_EXTRA"].to_numpy(dtype=np.float64))
        meta["max_neg_log10_pa"] = float(np.max(nlp_extra)) if len(nlp_extra) else 0.0
    return columns, meta


def compute_manhattan_exports(con: duckdb.DuckDBPyConnection) -> tuple[dict, dict]:
    dom_traits = {}
    add_traits = {}
    for trait_name, idx in (("bodyweight", 1), ("stature", 2)):
        cols_dom, meta_dom = load_thinned_candidates_from_duckdb(
            con, parquet_path(f"other_1.{idx}"), z_col="Z_D", extra_z_col="Z_A"
        )
        dom_traits[trait_name] = {"meta": meta_dom, "columns": cols_dom}

        cols_add, meta_add = load_thinned_candidates_from_duckdb(
            con, parquet_path(f"other_8.{idx}"), z_col="Z_A"
        )
        add_traits[trait_name] = {"meta": meta_add, "columns": cols_add}

    common_note = (
        "Columnar arrays sorted by (chr, bp_umd). neg_log10_p_x100 = round(-log10(P)*100) "
        "computed from Z via scipy.special.log_ndtr without clipping (stored as unclipped "
        "integers because PLAG1 reaches -log10(P)=441.3 -> 44130 > int16 max 32767). "
        "bp_ars=0 indicates unmapped in other_10."
    )
    dom_payload = {
        "model": "dominance",
        "gws_threshold_neg_log10_p": GWS_NEG_LOG10_P,
        "quantisation": "int_x100_unclipped",
        "note": common_note,
        "traits": dom_traits,
    }
    add_payload = {
        "model": "standard_additive",
        "gws_threshold_neg_log10_p": GWS_NEG_LOG10_P,
        "quantisation": "int_x100_unclipped",
        "note": common_note,
        "traits": add_traits,
    }
    return dom_payload, add_payload


def compute_qtl_contrast(con: duckdb.DuckDBPyConnection) -> dict:
    o1_bw = parquet_path("other_1.1")
    o8_bw = parquet_path("other_8.1")
    o1_st = parquet_path("other_1.2")
    o8_st = parquet_path("other_8.2")
    o1_bcs = parquet_path("other_1.3")
    o8_bcs = parquet_path("other_8.3")

    ext2a = pd.read_csv(SOURCE / "source_data.ExtFig2a.csv")
    std_cojo = ext2a[ext2a["Model"] == "Standard"].copy()

    # Precompute LD R^2 in gen_1.9 around DPF2 so we can annotate the chr29:45.06Mb SNP.
    bim_19 = pd.read_csv(
        GENS / "gen_1.9.bim",
        sep=r"\s+",
        header=None,
        names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
    )
    dpf2_idx = int(bim_19.index[bim_19["BP"] == 44213160][0])
    r2_19 = read_bed_chunked_r2(GENS / "gen_1.9.bed", dpf2_idx)
    r2_map_19 = dict(zip(bim_19["BP"].astype(int), r2_19, strict=False))

    rows = []
    for q in NON_ADDITIVE_QTLS:
        chrom = q["chr"]
        lo, hi = q["window_start"], q["window_end"]
        cand_bp = q["candidate_bp"]

        dom_best_d = con.execute(
            f"SELECT BP, A1, A2, MAF, B_A, B_D, Z_A, Z_D, P_A, P_D "
            f"FROM '{o1_bw}' WHERE CHR={chrom} AND BP BETWEEN {lo} AND {hi} "
            f"ORDER BY ABS(Z_D) DESC LIMIT 1"
        ).fetchone()
        dom_best_a = con.execute(
            f"SELECT BP, Z_A, Z_D FROM '{o1_bw}' WHERE CHR={chrom} AND BP BETWEEN {lo} AND {hi} "
            f"ORDER BY ABS(Z_A) DESC LIMIT 1"
        ).fetchone()
        add_best = con.execute(
            f"SELECT BP, A1, A2, MAF, B_A, Z_A, P_A "
            f"FROM '{o8_bw}' WHERE CHR={chrom} AND BP BETWEEN {lo} AND {hi} "
            f"ORDER BY ABS(Z_A) DESC LIMIT 1"
        ).fetchone()

        cand_dom = con.execute(
            f"SELECT MAF, B_A, B_D, Z_A, Z_D, P_A, P_D "
            f"FROM '{o1_bw}' WHERE CHR={chrom} AND BP={cand_bp} LIMIT 1"
        ).fetchone()
        cand_add = con.execute(
            f"SELECT MAF, B_A, Z_A, P_A "
            f"FROM '{o8_bw}' WHERE CHR={chrom} AND BP={cand_bp} LIMIT 1"
        ).fetchone()

        st_dom = con.execute(
            f"SELECT MAX(ABS(Z_D)), MAX(ABS(Z_A)) FROM '{o1_st}' "
            f"WHERE CHR={chrom} AND BP BETWEEN {lo} AND {hi}"
        ).fetchone()
        st_add = con.execute(
            f"SELECT MAX(ABS(Z_A)) FROM '{o8_st}' "
            f"WHERE CHR={chrom} AND BP BETWEEN {lo} AND {hi}"
        ).fetchone()
        bcs_dom = con.execute(
            f"SELECT MAX(ABS(Z_D)), MAX(ABS(Z_A)) FROM '{o1_bcs}' "
            f"WHERE CHR={chrom} AND BP BETWEEN {lo} AND {hi}"
        ).fetchone()
        bcs_add = con.execute(
            f"SELECT MAX(ABS(Z_A)) FROM '{o8_bcs}' "
            f"WHERE CHR={chrom} AND BP BETWEEN {lo} AND {hi}"
        ).fetchone()

        cojo_in_win = std_cojo[
            (std_cojo["CHR"] == chrom) & (std_cojo["BP"].between(lo, hi))
        ]

        best_dom_pd = float(neg_log10_p_from_z(dom_best_d[7]))
        best_dom_pa = float(neg_log10_p_from_z(dom_best_a[1]))
        best_add_pa = float(neg_log10_p_from_z(add_best[5]))
        cand_dom_pd = float(neg_log10_p_from_z(cand_dom[4])) if cand_dom else None
        cand_dom_pa = float(neg_log10_p_from_z(cand_dom[3])) if cand_dom else None
        cand_add_pa = float(neg_log10_p_from_z(cand_add[2])) if cand_add else None

        add_lead_bp = int(add_best[0])
        add_lead_r2_to_cand = r2_map_19.get(add_lead_bp) if chrom == 29 else None

        rows.append(
            {
                "qtl": q["qtl"],
                "gene": q["gene"],
                "chr": chrom,
                "candidate_bp": cand_bp,
                "window_start": lo,
                "window_end": hi,
                "is_primary_8": q["is_primary_8"],
                "category": q["category"],
                "bodyweight": {
                    "window_best_dom_d_neg_log10_p": best_dom_pd,
                    "window_best_dom_d_bp": int(dom_best_d[0]),
                    "window_best_dom_d_z": float(dom_best_d[7]),
                    "window_best_dom_a_neg_log10_p": best_dom_pa,
                    "window_best_dom_a_bp": int(dom_best_a[0]),
                    "window_best_add_neg_log10_p": best_add_pa,
                    "window_best_add_bp": add_lead_bp,
                    "window_best_add_z": float(add_best[5]),
                    "window_best_add_maf": float(add_best[3]),
                    "window_best_add_r2_to_candidate": add_lead_r2_to_cand,
                    "candidate_dom_d_neg_log10_p": cand_dom_pd,
                    "candidate_dom_a_neg_log10_p": cand_dom_pa,
                    "candidate_add_neg_log10_p": cand_add_pa,
                    "detected_by_dominance_gws": best_dom_pd >= GWS_NEG_LOG10_P,
                    "candidate_detected_by_additive_gws": (
                        cand_add_pa >= GWS_NEG_LOG10_P if cand_add_pa is not None else False
                    ),
                    "detected_in_cojo_standard_additive": len(cojo_in_win) > 0,
                    "cojo_standard_hits_in_window": int(len(cojo_in_win)),
                },
                "stature": {
                    "window_best_dom_d_neg_log10_p": float(neg_log10_p_from_z(st_dom[0])),
                    "window_best_dom_a_neg_log10_p": float(neg_log10_p_from_z(st_dom[1])),
                    "window_best_add_neg_log10_p": float(neg_log10_p_from_z(st_add[0])),
                },
                "bcs": {
                    "window_best_dom_d_neg_log10_p": float(neg_log10_p_from_z(bcs_dom[0])),
                    "window_best_dom_a_neg_log10_p": float(neg_log10_p_from_z(bcs_dom[1])),
                    "window_best_add_neg_log10_p": float(neg_log10_p_from_z(bcs_add[0])),
                },
            }
        )

    primary8 = [r for r in rows if r["is_primary_8"]]
    n_dom_detected = sum(r["bodyweight"]["detected_by_dominance_gws"] for r in primary8)
    n_cand_add_detected = sum(
        r["bodyweight"]["candidate_detected_by_additive_gws"] for r in primary8
    )
    n_cojo_add_detected = sum(
        r["bodyweight"]["detected_in_cojo_standard_additive"] for r in primary8
    )
    n_window_raw_add_detected = sum(
        r["bodyweight"]["window_best_add_neg_log10_p"] >= GWS_NEG_LOG10_P for r in primary8
    )

    fig1_reconciliation = {
        "bodyweight_chr14": {
            "paper_fig1_legend_p": "4.97e-442",
            "paper_fig1_legend_neg_log10_p": 441.3036,
            "other_1_1_dom_za_top_bp": 24974811,
            "other_1_1_dom_za_top_z": -45.802,
            "other_1_1_dom_za_top_neg_log10_p": float(neg_log10_p_from_z(-45.802)),
            "other_8_1_std_add_top_bp": 25015640,
            "other_8_1_std_add_top_z": 46.94,
            "other_8_1_std_add_top_neg_log10_p": float(neg_log10_p_from_z(46.94)),
            "other_8_1_matching_legend_bp": 25307116,
            "other_8_1_matching_legend_z": -44.991,
            "other_8_1_matching_legend_neg_log10_p": float(neg_log10_p_from_z(-44.991)),
            "other_8_1_matching_legend_p_from_3dp_z": "5.03e-442",
            "cause": (
                "Three factors explain the numbers: (1) 457.30 is the dominance-model genotypic "
                "additive term (Z_A=-45.802 at chr14:24,974,811 in other_1.1.txt), not the "
                "standard-additive model (other_8.1.txt); (2) in other_8.1.txt, 506 chr14 rows "
                "have P_A written as the literal string '<1e-308'; once parsed, the #1 variant is "
                "the reintroduced GGPv3 causal SNP chr14:25,015,640 (Z_A=46.940, -log10 P=480.22), "
                "while the #6 variant is BovineSNP50 SNP chr14:25,307,116 (Z_A=-44.991, two-sided "
                "-log10 P=441.299, P=5.03e-442), which matches the Fig. 1 legend's 4.97e-442 "
                "(different lead variant before/after custom GGPv3 site reintroduction); and "
                "(3) 3-d.p. Z rounding (-44.9913 rounds to -44.991, shifting P from 4.97e-442 "
                "to 5.03e-442). Tests are two-sided."
            ),
        },
        "bodyweight_chr5": {
            "paper_fig1_legend_p": "3.59e-127",
            "paper_fig1_legend_neg_log10_p": 126.4449,
            "other_8_1_std_add_top_bp": 106266836,
            "other_8_1_std_add_top_z": 24.45,
            "other_8_1_std_add_top_p_file": 5.05328e-132,
            "other_8_1_std_add_top_neg_log10_p": float(neg_log10_p_from_z(24.45)),
            "other_8_1_matching_legend_bp": 106250171,
            "other_8_1_matching_legend_z": -23.989,
            "other_8_1_matching_legend_p_file": 3.58619e-127,
            "other_8_1_matching_legend_neg_log10_p": float(neg_log10_p_from_z(-23.989)),
            "cause": (
                "Different lead variant in the same CCND2 peak: chr5:106,250,171 (A/G SNP, rank #7 "
                "on chr5 in other_8.1.txt) has file P_A=3.58619e-127, which rounds to the exact "
                "Fig. 1 legend value 3.59e-127 (-log10 P=126.441 from 3-d.p. Z_A=-23.989), whereas "
                "the #1 variant in other_8.1.txt (and in ExtFig2a.csv) is the 15-bp indel "
                "chr5:106,266,836 (Z_A=24.450, P_A=5.05e-132, -log10 P=131.30)."
            ),
        },
        "stature_chr5_and_chr14": {
            "chr5_paper_legend_p": "3.23e-97",
            "chr5_other_8_2_top_bp": 106264863,
            "chr5_other_8_2_top_z": 20.924,
            "chr5_other_8_2_top_p_file": 3.22659e-97,
            "chr5_other_8_2_top_neg_log10_p": float(neg_log10_p_from_z(20.924)),
            "chr14_paper_legend_p": "2.0e-306",
            "chr14_other_8_2_top_bp": 24973953,
            "chr14_other_8_2_top_z": 37.418,
            "chr14_other_8_2_top_p_file": 2.02645e-306,
            "chr14_other_8_2_top_neg_log10_p": float(neg_log10_p_from_z(37.418)),
            "cause": (
                "Exact match to the #1 variants in other_8.2.txt; third-digit difference between "
                "file P_A and log_ndtr(Z_A) is purely 3-d.p. rounding of Z_A."
            ),
        },
    }

    return {
        "gws_threshold_p": GWS_P,
        "gws_threshold_neg_log10_p": GWS_NEG_LOG10_P,
        "rows": rows,
        "fig1_legend_reconciliation": fig1_reconciliation,
        "summary": {
            "n_primary_non_additive_qtl": len(primary8),
            "detected_by_dominance_gws": int(n_dom_detected),
            "candidate_detected_by_standard_additive_gws": int(n_cand_add_detected),
            "detected_in_paper_cojo_standard_additive": int(n_cojo_add_detected),
            "invisible_to_standard_additive_at_candidate_and_cojo": int(
                len(primary8) - n_cojo_add_detected
            ),
            "window_1mb_raw_additive_gws_count": int(n_window_raw_add_detected),
            "note": (
                "7 of the 8 dominance QTL (all except PLAG1) are invisible to the "
                "standard-additive model at the candidate/lead variants (all candidate "
                "additive -log10(P) <= 2.34) and in the paper's GCTA-COJO standard-additive "
                "QTL set (ExtFig2a.csv). "
                "In the raw +-1 Mb window scan around DPF2/MUS81, an unlinked common SNP "
                "843 kb away (chr29:45,056,385, MAF=0.087) reaches -log10(P_A)=10.09, whereas "
                "DPF2 (chr29:44,213,160) has standard-additive -log10(P_A)=0.11 and MUS81 "
                "(chr29:44,645,469) has -log10(P_A)=0.13."
            ),
        },
    }


def fetch_ensembl_genes_ars(chrom: int, start_ars: int, end_ars: int) -> list[dict]:
    cache_path = WORK / "ensembl_genes_cache.json"
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    key = f"{chrom}:{start_ars}-{end_ars}"
    if key in cache:
        return cache[key]

    url = (
        f"https://rest.ensembl.org/overlap/region/bos_taurus/"
        f"{chrom}:{start_ars}-{end_ars}?feature=gene;content-type=application/json"
    )
    genes = []
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "lic-recessive-explorer/1.0"})
        with urllib.request.urlopen(req, timeout=20) as resp:
            raw = json.loads(resp.read().decode("utf-8"))
        for g in raw:
            if g.get("biotype") != "protein_coding":
                continue
            genes.append(
                {
                    "gene_id": g.get("id"),
                    "symbol": g.get("external_name") or g.get("id"),
                    "start_ars": int(g["start"]),
                    "end_ars": int(g["end"]),
                    "strand": int(g.get("strand", 1)),
                    "biotype": g.get("biotype"),
                }
            )
        cache[key] = genes
        cache_path.write_text(json.dumps(cache, indent=2))
    except Exception:
        genes = []
    return genes


def compute_regional_exports(con: duckdb.DuckDBPyConnection) -> list[dict]:
    o1_bw = parquet_path("other_1.1")
    o8_bw = parquet_path("other_8.1")
    o1_bcs = parquet_path("other_1.3")
    o8_bcs = parquet_path("other_8.3")
    o10 = parquet_path("other_10")

    ext3 = pd.read_csv(SOURCE / "source_data.ExtFig3.csv")

    manifest_entries = []
    index_rows = []
    for slug, label, gene, chrom, cand_bp, gen1_stem, cand_desc in REGIONAL_LOCI:
        bim_path = GENS / f"{gen1_stem}.bim"
        bed_path = GENS / f"{gen1_stem}.bed"
        bim = pd.read_csv(
            bim_path,
            sep=r"\s+",
            header=None,
            names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
        )
        lead_idx = int((bim["BP"] - cand_bp).abs().idxmin())
        lead_bp = int(bim.loc[lead_idx, "BP"])
        r2 = read_bed_chunked_r2(bed_path, lead_idx, chunk_size=2000)

        min_bp = int(bim["BP"].min())
        max_bp = int(bim["BP"].max())

        gwas_df = con.execute(
            f"""
            WITH m10 AS (
                SELECT "BP.UMD" AS BP,
                       FIRST("CHR.ARS") AS chr_ars,
                       FIRST("BP.ARS") AS bp_ars
                FROM '{o10}'
                WHERE "CHR.UMD" = {chrom} AND "BP.UMD" BETWEEN {min_bp} AND {max_bp}
                GROUP BY 1
            ),
            d1 AS (
                SELECT BP, FIRST(MAF) AS MAF, FIRST(B_A) AS B_A, FIRST(B_D) AS B_D,
                       FIRST(Z_A) AS Z_A, FIRST(Z_D) AS Z_D
                FROM '{o1_bw}'
                WHERE CHR = {chrom} AND BP BETWEEN {min_bp} AND {max_bp}
                GROUP BY 1
            ),
            a8 AS (
                SELECT BP, FIRST(B_A) AS ADD_B_A, FIRST(Z_A) AS ADD_Z_A
                FROM '{o8_bw}'
                WHERE CHR = {chrom} AND BP BETWEEN {min_bp} AND {max_bp}
                GROUP BY 1
            )
            SELECT d1.BP, d1.MAF, d1.B_A, d1.B_D, d1.Z_A, d1.Z_D,
                   a8.ADD_B_A, a8.ADD_Z_A,
                   m10.chr_ars, m10.bp_ars
            FROM d1
            LEFT JOIN a8 USING (BP)
            LEFT JOIN m10 USING (BP)
            """
        ).df()

        merged = pd.DataFrame(
            {
                "bp_umd": bim["BP"].astype(int),
                "snp": bim["SNP"].astype(str),
                "a1": bim["A1"].astype(str),
                "a2": bim["A2"].astype(str),
                "r2": np.round(r2, 4),
            }
        ).merge(gwas_df, left_on="bp_umd", right_on="BP", how="left")

        nlp_dom_d = neg_log10_p_from_z(merged["Z_D"].fillna(0.0).to_numpy())
        nlp_dom_a = neg_log10_p_from_z(merged["Z_A"].fillna(0.0).to_numpy())
        nlp_add = neg_log10_p_from_z(merged["ADD_Z_A"].fillna(0.0).to_numpy())

        # Lift gene annotations between ARS-UCD1.2 and UMD3.1 via other_10 anchors (with +-1.5 Mb
        # margin so genes near window boundaries are interpolated without edge clamping).
        wide_anchors = con.execute(
            f"""
            SELECT "BP.UMD" AS bp_umd, "BP.ARS" AS bp_ars
            FROM '{o10}'
            WHERE "CHR.UMD" = {chrom}
              AND "CHR.ARS" = '{chrom}'
              AND "BP.UMD" BETWEEN {max(1, min_bp - 1_500_000)} AND {max_bp + 1_500_000}
              AND "BP.ARS" IS NOT NULL
            ORDER BY "BP.ARS"
            """
        ).df()
        in_win_anchors = merged[
            (merged["chr_ars"].astype(str) == str(chrom)) & merged["bp_ars"].notna()
        ]
        genes_out = []
        if len(in_win_anchors) >= 10 and len(wide_anchors) >= 10:
            ars_min = int(in_win_anchors["bp_ars"].min())
            ars_max = int(in_win_anchors["bp_ars"].max())
            raw_genes = fetch_ensembl_genes_ars(chrom, ars_min, ars_max)
            xp = wide_anchors["bp_ars"].to_numpy(dtype=np.float64)
            fp = wide_anchors["bp_umd"].to_numpy(dtype=np.float64)
            order = np.argsort(xp)
            xp, fp = xp[order], fp[order]
            for g in raw_genes:
                s_umd = int(round(float(np.interp(g["start_ars"], xp, fp))))
                e_umd = int(round(float(np.interp(g["end_ars"], xp, fp))))
                if s_umd > e_umd:
                    s_umd, e_umd = e_umd, s_umd
                if e_umd < min_bp or s_umd > max_bp:
                    continue
                genes_out.append(
                    {
                        "gene_id": g["gene_id"],
                        "symbol": g["symbol"],
                        "assembly": "UMD3.1",
                        "start": s_umd,
                        "end": e_umd,
                        "start_umd": s_umd,
                        "end_umd": e_umd,
                        "start_ars": int(g["start_ars"]),
                        "end_ars": int(g["end_ars"]),
                        "strand": int(g["strand"]),
                        "biotype": g["biotype"],
                        "is_candidate_gene": (gene is not None and g["symbol"] == gene),
                    }
                )

        # In Ensembl ARS-UCD1.2, PLCD4 is merged into a CNOT9 readthrough transcript
        # (ENSBTAG00000003807), whereas in NCBI RefSeq (NM_001046509.1 / NP_001039954.1) and
        # Ensembl 86 (UMD3.1) PLCD4 is a distinct gene at chr2:107,298,023-107,321,320 UMD3.1.
        existing_symbols = {g["symbol"] for g in genes_out}
        if slug == "plcd4" and "PLCD4" not in existing_symbols:
            genes_out.append(
                {
                    "gene_id": "NCBI:509076",
                    "symbol": "PLCD4",
                    "assembly": "UMD3.1",
                    "start": 107298023,
                    "end": 107321320,
                    "start_umd": 107298023,
                    "end_umd": 107321320,
                    "start_ars": 106575326,
                    "end_ars": 106598623,
                    "strand": 1,
                    "biotype": "protein_coding",
                    "is_candidate_gene": True,
                    "note": (
                        "NCBI RefSeq NM_001046509.1 / Ensembl 86 UMD3.1 coordinates; "
                        "merged into CNOT9 readthrough in Ensembl ARS-UCD1.2."
                    ),
                }
            )

        # Keep all variants that have GWAS summary statistics in the regional window
        has_gwas = merged["Z_D"].notna() | (merged["bp_umd"] == lead_bp)
        sub = merged[has_gwas].copy().reset_index(drop=True)
        mask_arr = has_gwas.to_numpy()
        ld_available = slug != "plag1"

        columns = {
            "bp_umd": sub["bp_umd"].astype(int).tolist(),
            "bp_ars": sub["bp_ars"].fillna(0).astype(int).tolist(),
            "r2_x1000": (
                np.rint(sub["r2"].to_numpy() * 1000.0).astype(int).tolist()
                if ld_available
                else None
            ),
            "neg_log10_pd_x100": np.rint(nlp_dom_d[mask_arr] * 100.0).astype(int).tolist(),
            "neg_log10_pa_dom_x100": np.rint(nlp_dom_a[mask_arr] * 100.0).astype(int).tolist(),
            "neg_log10_pa_add_x100": np.rint(nlp_add[mask_arr] * 100.0).astype(int).tolist(),
        }

        extra_payload = {}
        if slug == "plag1":
            extra_payload["external_reference_gene"] = {
                "symbol": "PLAG1",
                "gene_id": "ENSBTAG00000004022",
                "assembly": "UMD3.1",
                "start_umd": 25007291,
                "end_umd": 25052403,
                "causal_variant_bp_umd": 25015640,
                "distance_outside_window_bp": int(min_bp - 25052403),
                "note": (
                    "PLAG1 and its published causal variant rs109815800 (chr14:25,015,640) lie "
                    "~263 kb upstream of the gen_1.6 window (25,315,124-26,814,744); "
                    "no LD colouring is shown for this window."
                ),
            }
        if slug == "myh1":
            # Include BCS dominance & additive stats for MYH1 + ExtFig3 comparison
            bcs_df = con.execute(
                f"""
                SELECT d.BP AS bp_umd, d.Z_D AS bcs_zd, a.Z_A AS bcs_add_za
                FROM '{o1_bcs}' d
                LEFT JOIN '{o8_bcs}' a USING (CHR, BP)
                WHERE d.CHR = 19 AND d.BP BETWEEN {min_bp} AND {max_bp}
                """
            ).df()
            sub_bcs = sub[["bp_umd", "r2"]].merge(bcs_df, on="bp_umd", how="left")
            bcs_pd = neg_log10_p_from_z(sub_bcs["bcs_zd"].fillna(0.0).to_numpy())
            bcs_add_pa = neg_log10_p_from_z(sub_bcs["bcs_add_za"].fillna(0.0).to_numpy())
            columns["bcs_neg_log10_pd_x100"] = np.rint(bcs_pd * 100.0).astype(int).tolist()
            columns["bcs_neg_log10_pa_add_x100"] = np.rint(bcs_add_pa * 100.0).astype(int).tolist()

            ext3_m = sub_bcs.merge(ext3, left_on="bp_umd", right_on="BP", how="inner")
            extra_payload["extfig3_comparison"] = {
                "n_extfig3_rows": int(len(ext3)),
                "n_matched": int(len(ext3_m)),
                "r2_correlation": float(ext3_m["r2"].corr(ext3_m["R2"])),
                "max_bcs_neg_log10_pd": float(np.max(bcs_pd)),
                "extfig3_min_p": float(ext3["P"].min()),
                "extfig3_max_neg_log10_p": float(-math.log10(ext3["P"].min())),
            }

        cand_pos_in_sub = sub.index[sub["bp_umd"] == lead_bp]
        cand_sub_idx = int(cand_pos_in_sub[0]) if len(cand_pos_in_sub) else None
        payload = {
            "locus": label,
            "slug": slug,
            "gene": gene,
            "chr": chrom,
            "x_axis_assembly": "UMD3.1",
            "gene_track_assembly": "UMD3.1",
            "ld_available": ld_available,
            "window_umd": [min_bp, max_bp],
            "candidate_variant": {
                "bp_umd": lead_bp,
                "bp_ars": (
                    int(sub.loc[cand_sub_idx, "bp_ars"])
                    if cand_sub_idx is not None and pd.notna(sub.loc[cand_sub_idx, "bp_ars"])
                    else None
                ),
                "snp": str(bim.loc[lead_idx, "SNP"]),
                "a1": str(bim.loc[lead_idx, "A1"]),
                "a2": str(bim.loc[lead_idx, "A2"]),
                "description": cand_desc,
                "index_in_columns": cand_sub_idx,
                "neg_log10_pd": (
                    float(columns["neg_log10_pd_x100"][cand_sub_idx] / 100.0)
                    if cand_sub_idx is not None
                    else None
                ),
                "neg_log10_pa_add": (
                    float(columns["neg_log10_pa_add_x100"][cand_sub_idx] / 100.0)
                    if cand_sub_idx is not None
                    else None
                ),
            },
            "n_variants_bim": int(len(bim)),
            "n_variants_exported": int(len(sub)),
            "n_high_ld_r2_gt_0_8": int((sub["r2"] >= 0.8).sum()) if ld_available else None,
            "annotation_source": (
                "Ensembl REST API bos_taurus (ARS-UCD1.2 protein_coding genes), "
                "lifted to UMD3.1 via piecewise-linear interpolation on other_10 anchors."
            ),
            "genes": genes_out,
            **extra_payload,
            "columns": columns,
        }
        entry = write_json(
            f"regional/{slug}.json",
            payload,
            [f"dryad/gens/{gen1_stem}.*", "dryad/other/other_1.1.txt", "dryad/other/other_8.1.txt"],
            compact=True,
        )
        manifest_entries.append(entry)
        index_rows.append(
            {
                "slug": slug,
                "locus": label,
                "gene": gene,
                "chr": chrom,
                "candidate_bp_umd": lead_bp,
                "file": f"regional/{slug}.json",
                "n_variants": int(len(sub)),
                "n_genes": len(genes_out),
            }
        )

    manifest_entries.append(
        write_json(
            "regional/index.json",
            {"loci": index_rows},
            ["dryad/gens/gen_1.*"],
        )
    )
    return manifest_entries


def compute_chip_blind_spot(con: duckdb.DuckDBPyConnection) -> dict:
    o15 = parquet_path("other_15")
    o1 = parquet_path("other_1.1")
    o10 = parquet_path("other_10")

    chip_df = con.execute(
        f"""
        WITH c15 AS (
            SELECT CHR, BP, MIN(P_D) AS chip_p_d
            FROM '{o15}'
            GROUP BY 1, 2
        ),
        s1 AS (
            SELECT a.CHR, a.BP,
                   FIRST(a.MAF) AS MAF, FIRST(a.Z_D) AS Z_D, FIRST(a.P_D) AS seq_p_d
            FROM '{o1}' a
            SEMI JOIN c15 USING (CHR, BP)
            GROUP BY 1, 2
        ),
        m10 AS (
            SELECT m."CHR.UMD" AS CHR, m."BP.UMD" AS BP,
                   FIRST(m."CHR.ARS") AS CHR_ARS, FIRST(m."BP.ARS") AS BP_ARS
            FROM '{o10}' m
            SEMI JOIN c15 ON m."CHR.UMD" = c15.CHR AND m."BP.UMD" = c15.BP
            GROUP BY 1, 2
        )
        SELECT c15.CHR, c15.BP, c15.chip_p_d, s1.MAF, s1.Z_D, s1.seq_p_d,
               m10.CHR_ARS, m10.BP_ARS
        FROM c15
        LEFT JOIN s1 USING (CHR, BP)
        LEFT JOIN m10 USING (CHR, BP)
        ORDER BY c15.CHR, c15.BP
        """
    ).df()

    nlp_z = neg_log10_p_from_z(chip_df["Z_D"].to_numpy(dtype=np.float64))
    nlp_file = -np.log10(np.clip(chip_df["chip_p_d"].to_numpy(dtype=np.float64), 1e-300, 1.0))
    chip_df["nlp"] = nlp_z

    locus_rows = []
    for q in NON_ADDITIVE_QTLS:
        chrom = q["chr"]
        lo, hi = q["window_start"], q["window_end"]
        seq_max_z = con.execute(
            f"SELECT MAX(ABS(Z_D)) FROM '{o1}' WHERE CHR={chrom} AND BP BETWEEN {lo} AND {hi}"
        ).fetchone()[0]
        seq_best_nlp = float(neg_log10_p_from_z(seq_max_z))

        sub_chip = chip_df[(chip_df["CHR"] == chrom) & (chip_df["BP"].between(lo, hi))]
        if len(sub_chip):
            best_idx = int(sub_chip["nlp"].idxmax())
            chip_best_nlp = float(sub_chip.loc[best_idx, "nlp"])
            chip_best_bp = int(sub_chip.loc[best_idx, "BP"])
            chip_n_gws = int((sub_chip["nlp"] >= GWS_NEG_LOG10_P).sum())
        else:
            chip_best_nlp = 0.0
            chip_best_bp = None
            chip_n_gws = 0

        locus_rows.append(
            {
                "qtl": q["qtl"],
                "gene": q["gene"],
                "chr": chrom,
                "candidate_bp": q["candidate_bp"],
                "window_start": lo,
                "window_end": hi,
                "is_primary_8": q["is_primary_8"],
                "is_major_recessive_syndrome": q["is_major_recessive_syndrome"],
                "sequence_best_neg_log10_pd": seq_best_nlp,
                "chip50k_best_neg_log10_pd": chip_best_nlp,
                "chip50k_best_bp": chip_best_bp,
                "chip50k_n_snps_in_window": int(len(sub_chip)),
                "chip50k_n_gws_in_window": chip_n_gws,
                "sequence_gws": seq_best_nlp >= GWS_NEG_LOG10_P,
                "chip50k_gws": chip_best_nlp >= GWS_NEG_LOG10_P,
            }
        )

    major_syndromes_kept = [
        r["qtl"]
        for r in locus_rows
        if r["is_major_recessive_syndrome"] and r["chip50k_gws"]
    ]
    all_kept = [r["qtl"] for r in locus_rows if r["is_primary_8"] and r["chip50k_gws"]]

    ars_chr = pd.to_numeric(chip_df["CHR_ARS"], errors="coerce").fillna(0).astype(int)
    ars_bp = chip_df["BP_ARS"].fillna(0).astype(int)

    return {
        "n_chip_distinct_positions": int(len(chip_df)),
        "fraction_of_chip_in_sequence": float(chip_df["Z_D"].notna().mean()),
        "fraction_of_sequence_on_chip": float(len(chip_df) / 16_128_339),
        "max_abs_diff_nlp_z_vs_file": float(np.max(np.abs(nlp_z - nlp_file))),
        "n_chip_gws_snps_genomewide": int((chip_df["nlp"] >= GWS_NEG_LOG10_P).sum()),
        "loci_comparison": locus_rows,
        "summary": {
            "primary_8_loci_significant_at_sequence": 8,
            "primary_8_loci_significant_on_50k_chip": all_kept,
            "major_recessive_syndromes_significant_on_50k_chip": major_syndromes_kept,
            "note": (
                "100% of the 45,865 distinct BovineSNP50 positions in other_15 are present in "
                "other_1.1 (0.284% of sequence variants). While the higher-MAF loci Chr2:22 "
                "(-log10 P=7.64) and PLAG1 (-log10 P=17.33) retain genome-wide significance on "
                "the 50k chip, only DPF2/MUS81 (-log10 P=14.76) remains genome-wide significant "
                "among the major-effect recessive syndromes; PLCD4 (3.10), FGD4 (3.95), "
                "ARSI (2.91), MTRF1 (5.42), and GALNT2 (3.61) are all missed at 50k-chip "
                "resolution."
            ),
        },
        "columns": {
            "chr": chip_df["CHR"].astype(int).tolist(),
            "bp_umd": chip_df["BP"].astype(int).tolist(),
            "chr_ars": ars_chr.tolist(),
            "bp_ars": ars_bp.tolist(),
            "neg_log10_pd_x100": np.rint(nlp_z * 100.0).astype(int).tolist(),
        },
    }


def compute_imputation_by_maf(con: duckdb.DuckDBPyConnection) -> dict:
    o9 = parquet_path("other_9")
    o1 = parquet_path("other_1.1")

    bin_meta = {
        "0.01-0.02": (0.01, 0.02),
        "0.02-0.03": (0.02, 0.03),
        "0.03-0.04": (0.03, 0.04),
        "0.04-0.05": (0.04, 0.05),
        "0.05-0.10": (0.05, 0.10),
        "0.10-0.20": (0.10, 0.20),
        "0.20-0.30": (0.20, 0.30),
        "0.30-0.40": (0.30, 0.40),
        "0.40-0.50": (0.40, 0.50),
    }
    agg_rows = con.execute(
        f"""
        SELECT
            CASE
                WHEN b.MAF < 0.02 THEN '0.01-0.02'
                WHEN b.MAF < 0.03 THEN '0.02-0.03'
                WHEN b.MAF < 0.04 THEN '0.03-0.04'
                WHEN b.MAF < 0.05 THEN '0.04-0.05'
                WHEN b.MAF < 0.10 THEN '0.05-0.10'
                WHEN b.MAF < 0.20 THEN '0.10-0.20'
                WHEN b.MAF < 0.30 THEN '0.20-0.30'
                WHEN b.MAF < 0.40 THEN '0.30-0.40'
                ELSE '0.40-0.50'
            END AS maf_bin,
            COUNT(*) AS n,
            AVG(a.AR2) AS mean_ar2,
            STDDEV_SAMP(a.AR2) AS sd_ar2,
            QUANTILE_CONT(a.AR2, [0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95]) AS qs,
            AVG(CASE WHEN a.AR2 >= 0.95 THEN 1.0 ELSE 0.0 END) AS frac_ge_095,
            AVG(CASE WHEN a.AR2 >= 0.80 THEN 1.0 ELSE 0.0 END) AS frac_ge_080
        FROM '{o1}' b
        JOIN '{o9}' a USING (CHR, BP)
        WHERE b.MAF >= 0.01 AND b.MAF <= 0.500001
        GROUP BY 1
        ORDER BY 1
        """
    ).fetchall()
    bin_rows = []
    for r in agg_rows:
        label = str(r[0])
        lo, hi = bin_meta[label]
        qs = r[4]
        bin_rows.append(
            {
                "maf_bin": label,
                "maf_min": lo,
                "maf_max": hi,
                "n_variants": int(r[1]),
                "mean_ar2": float(r[2]),
                "sd_ar2": float(r[3]),
                "quantiles": {
                    "p05": float(qs[0]),
                    "p10": float(qs[1]),
                    "p25": float(qs[2]),
                    "p50": float(qs[3]),
                    "p75": float(qs[4]),
                    "p90": float(qs[5]),
                    "p95": float(qs[6]),
                },
                "fraction_ar2_ge_0_95": float(r[5]),
                "fraction_ar2_ge_0_80": float(r[6]),
            }
        )

    candidates_spec = [
        ("Chr2:22", None, 2, 22836659, True),
        ("PLCD4", "PLCD4", 2, 107313998, True),
        ("FGD4", "FGD4", 5, 77632752, True),
        ("ARSI", "ARSI", 7, 63649071, True),
        ("MTRF1", "MTRF1", 12, 11463981, True),
        ("GALNT2", "GALNT2", 28, 1312334, True),
        ("DPF2", "DPF2", 29, 44213160, True),
        ("MUS81", "MUS81", 29, 44645469, True),
        ("MYH1", "MYH1", 19, 30114249, False),
        ("DGAT1", "DGAT1", 14, 1802265, False),
        ("PLAG1_causal", "PLAG1", 14, 25015640, False),
        ("PLAG1_dom_lead", "PLAG1", 14, 26064929, False),
    ]
    cand_rows = []
    for locus, gene, chrom, bp, is_bw_recessive in candidates_spec:
        r = con.execute(
            f"""
            SELECT a.AR2, b.MAF, b.B_A, b.B_D, b.Z_D
            FROM '{o9}' a
            LEFT JOIN '{o1}' b USING (CHR, BP)
            WHERE a.CHR = {chrom} AND a.BP = {bp}
            LIMIT 1
            """
        ).fetchone()
        ar2 = float(r[0]) if r else None
        cand_rows.append(
            {
                "locus": locus,
                "gene": gene,
                "chr": chrom,
                "bp": bp,
                "is_bw_recessive_table1": is_bw_recessive,
                "ar2": ar2,
                "maf": float(r[1]) if r and r[1] is not None else None,
                "abs_b_d_kg": abs(float(r[3])) if r and r[3] is not None else None,
                "hom_alt_derived_kg": -2.0 * float(r[2]) if r and r[2] is not None else None,
                "neg_log10_pd": (
                    float(neg_log10_p_from_z(r[4])) if r and r[4] is not None else None
                ),
                "exceeds_0_95": bool(ar2 is not None and ar2 > 0.95),
            }
        )

    bw_recessive_ar2 = [
        c["ar2"] for c in cand_rows if c["is_bw_recessive_table1"] and c["ar2"] is not None
    ]
    return {
        "maf_bins": bin_rows,
        "candidate_variants": cand_rows,
        "summary": {
            "min_ar2_bw_recessive_candidates": float(min(bw_recessive_ar2)),
            "all_bw_recessive_candidates_gt_0_95": all(v > 0.95 for v in bw_recessive_ar2),
            "note": (
                "All 8 body-weight recessive/dominance candidate variants in Table 1 have "
                "AR2 > 0.95 (range 0.979 to 1.000; PLCD4, FGD4, GALNT2 have AR2=1.000 as "
                "custom-chip physically typed content). The MYH1 left-flank proxy SNP "
                "(chr19:30,114,249, suggestive for BCS) has AR2 = 0.783."
            ),
        },
    }


def update_manifest(new_entries: list[dict]) -> None:
    manifest_path = OUT / "manifest.json"
    existing = {}
    if manifest_path.exists():
        data = json.loads(manifest_path.read_text())
        for item in data.get("files", []):
            existing[item["file"]] = item
    for item in new_entries:
        existing[item["file"]] = item
    # Ensure gzip_bytes is populated for all manifest entries
    files_list = []
    for fname in sorted(existing):
        entry = existing[fname]
        fpath = OUT / fname
        if fpath.exists():
            raw = fpath.read_bytes()
            entry["bytes"] = len(raw)
            entry["gzip_bytes"] = len(gzip.compress(raw, compresslevel=6))
        files_list.append(entry)
    manifest_payload = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": now_iso(),
        "code_version": CODE_VERSION,
        "files": files_list,
        "total_bytes": int(sum(m["bytes"] for m in files_list)),
        "total_gzip_bytes": int(sum(m.get("gzip_bytes", 0) for m in files_list)),
    }
    manifest_path.write_text(json.dumps(jsonable(manifest_payload), indent=2) + "\n")


def main() -> None:
    t0 = time.time()
    con = connect_duckdb()
    ensure_parquet_files(con)

    manifest_entries = []

    qtl_contrast = compute_qtl_contrast(con)
    manifest_entries.append(
        write_json(
            "qtl_contrast.json",
            qtl_contrast,
            [
                "dryad/other/other_1.1.txt",
                "dryad/other/other_8.1.txt",
                "dryad/other/other_1.2.txt",
                "dryad/other/other_8.2.txt",
                "dryad/source_data/source_data.ExtFig2a.csv",
            ],
        )
    )

    dom_payload, add_payload = compute_manhattan_exports(con)
    manifest_entries.append(
        write_json(
            "manhattan_dom.json",
            dom_payload,
            ["dryad/other/other_1.1.txt", "dryad/other/other_1.2.txt", "dryad/other/other_10.txt"],
            compact=True,
        )
    )
    manifest_entries.append(
        write_json(
            "manhattan_add.json",
            add_payload,
            ["dryad/other/other_8.1.txt", "dryad/other/other_8.2.txt", "dryad/other/other_10.txt"],
            compact=True,
        )
    )

    manifest_entries.extend(compute_regional_exports(con))

    chip_payload = compute_chip_blind_spot(con)
    manifest_entries.append(
        write_json(
            "chip_vs_sequence.json",
            chip_payload,
            ["dryad/other/other_15.txt", "dryad/other/other_1.1.txt", "dryad/other/other_10.txt"],
            compact=True,
        )
    )

    imputation_payload = compute_imputation_by_maf(con)
    manifest_entries.append(
        write_json(
            "imputation_by_maf.json",
            imputation_payload,
            ["dryad/other/other_9.txt", "dryad/other/other_1.1.txt"],
        )
    )

    update_manifest(manifest_entries)
    total_raw = sum(m["bytes"] for m in manifest_entries)
    total_gz = sum(m["gzip_bytes"] for m in manifest_entries)
    print(
        f"wrote {len(manifest_entries)} Phase 2 exports in {time.time() - t0:.1f}s "
        f"(raw={total_raw / 1e6:.2f} MB, gzip={total_gz / 1e6:.2f} MB)"
    )


if __name__ == "__main__":
    main()
