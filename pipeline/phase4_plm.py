"""Phase 4a: ESM-2 masked-marginal scoring of three missense candidates (ARSI dropped).

Method
- Model: facebook/esm2_t30_150M_UR50D (local cache only).
- For each position i, mask residue i and read the log-probabilities of the 20 standard amino
  acids. LLR(i, aa) = log p(aa | masked i) - log p(wt | masked i). Negative = less tolerated.
- Mutation percentile: 100 x fraction of all possible substitutions in the protein whose LLR is
  greater than or equal to the mutation's LLR (higher = more damaging than more substitutions).

Scope and limits
- n = 3 candidates, not validated against any phenotype beyond the anecdotal homozygote effects
  in effects.json (n = 4 homozygote-type comparisons; stated as anecdotal in the export).
- Sequences: NCBI RefSeq, stored in pipeline/data/refseq_proteins.fasta. Wild-type residues are
  asserted at each position before scoring.
- Not scorable this way: MUS81 p.Gly70* (stop-gain) and splice variants.
"""

from __future__ import annotations

import json
import os
import subprocess
from datetime import UTC, datetime
from pathlib import Path

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")

import numpy as np  # noqa: E402
import torch  # noqa: E402
from transformers import AutoTokenizer, EsmForMaskedLM  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web" / "public" / "data"
FASTA = ROOT / "pipeline" / "data" / "refseq_proteins.fasta"
EFFECTS = OUT / "effects.json"
MODEL_ID = "facebook/esm2_t30_150M_UR50D"
AA = "ACDEFGHIKLMNPQRSTVWY"

# gene, accession, 1-based position, wt, mut, paper p.-notation
CANDIDATES = [
    ("PLCD4", "NP_001039954.1", 326, "A", "T", "p.Ala326Thr"),
    ("MTRF1", "NP_001030185.1", 341, "R", "W", "p.Arg341Trp"),
    ("DPF2", "NP_001093826.1", 216, "K", "R", "p.Lys216Arg"),
]


def code_version() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        elif name is not None:
            seqs[name] += line.strip()
    return seqs


def masked_marginal_llr(model, tokenizer, seq: str, device: str) -> np.ndarray:
    """Return an (L, 20) array of LLR(i, aa) = log p(aa|masked i) - log p(wt_i|masked i)."""
    enc = tokenizer(seq, return_tensors="pt")
    ids = enc["input_ids"][0]  # [cls] seq [eos]
    L = len(seq)
    aa_ids = [tokenizer.convert_tokens_to_ids(a) for a in AA]
    mask_id = tokenizer.mask_token_id
    out = np.zeros((L, len(AA)), dtype=np.float64)
    batch = 32
    with torch.no_grad():
        for start in range(0, L, batch):
            idx = list(range(start, min(L, start + batch)))
            x = ids.unsqueeze(0).repeat(len(idx), 1)
            for row, i in enumerate(idx):
                x[row, i + 1] = mask_id
            logits = model(input_ids=x.to(device)).logits  # (b, T, V)
            logp = torch.log_softmax(logits.float(), dim=-1).cpu().numpy()
            for row, i in enumerate(idx):
                lp = logp[row, i + 1, aa_ids]
                wt_lp = lp[AA.index(seq[i])]
                out[i] = lp - wt_lp
    return out


def mutation_percentile(llr: np.ndarray, seq: str, pos0: int, mut: str) -> dict:
    mask = np.ones_like(llr, dtype=bool)
    for i, wt in enumerate(seq):
        mask[i, AA.index(wt)] = False  # exclude wt -> wt (LLR 0 by construction)
    all_subs = llr[mask]
    mut_llr = float(llr[pos0, AA.index(mut)])
    pct = 100.0 * float(np.mean(all_subs >= mut_llr))
    return {"llr": round(mut_llr, 3), "damage_percentile": round(pct, 2)}


def observed_effects() -> dict[str, dict]:
    rows = json.loads(EFFECTS.read_text())["rows"]
    out: dict[str, dict] = {}
    for r in rows:
        if r["trait"] != "bodyweight" or r["estimator"] not in ("A", "B"):
            continue
        if r["locus"] not in {"PLCD4", "MTRF1", "DPF2"}:
            continue
        out.setdefault(r["locus"], {})[r["estimator"]] = {
            "hom_alt_kg": round(float(r["hom_alt"]), 2),
            "se_kg": round(float(r["se_hom_alt"]), 2),
        }
    return out


def compute() -> dict:
    seqs = read_fasta(FASTA)
    tokenizer = AutoTokenizer.from_pretrained(MODEL_ID, local_files_only=True)
    model = EsmForMaskedLM.from_pretrained(MODEL_ID, local_files_only=True).eval()
    device = "mps" if torch.backends.mps.is_available() else "cpu"
    model = model.to(device)
    obs = observed_effects()

    genes = []
    for gene, acc, pos, wt, mut, notation in CANDIDATES:
        seq = seqs[acc]
        if seq[pos - 1] != wt:
            raise ValueError(f"{gene}: residue {pos} is {seq[pos - 1]}, expected {wt}")
        llr = masked_marginal_llr(model, tokenizer, seq, device)
        mm = mutation_percentile(llr, seq, pos - 1, mut)
        matrix = [
            [i + 1, seq[i], *[round(float(v), 2) for v in llr[i]]] for i in range(len(seq))
        ]
        genes.append(
            {
                "gene": gene,
                "accession": acc,
                "length": len(seq),
                "mutation": notation,
                "position": pos,
                "wt": wt,
                "mut": mut,
                **mm,
                "observed_bodyweight_hom_effect": obs.get(gene, {}),
                "matrix_columns": ["pos", "wt", *AA],
                "matrix": matrix,
            }
        )

    return {
        "schema_version": "1.0.0",
        "generated_at": datetime.now(UTC).isoformat(),
        "code_version": code_version(),
        "model": MODEL_ID,
        "device": device,
        "source_files": {
            "proteins": "pipeline/data/refseq_proteins.fasta",
            "effects": "web/public/data/effects.json",
        },
        "method": "ESM-2 masked-marginal LLR, one masked pass per position",
        "dropped": {
            "ARSI": "XP_002689338.1 residue 301 is N (G at 302); mismatch; "
            "dropped per user decision",
        },
        "genes": genes,
        "notes": [
            "damage_percentile = 100 x fraction of all possible substitutions with "
            "LLR >= the mutation LLR.",
            "n = 3 candidates; no validation beyond anecdotal homozygote effects (n = 4 "
            "homozygote comparisons). Not evidence of causality.",
            "MUS81 p.Gly70* (stop-gain) and splice variants are not scorable this way.",
        ],
    }


def main() -> None:
    payload = compute()
    (OUT / "plm_scores.json").write_text(json.dumps(payload) + "\n")


if __name__ == "__main__":
    main()
