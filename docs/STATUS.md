# Project Status — Recessive Variant Explorer

## Phase table

| Phase | Status | Commit | Note |
|---|---|---|---|
| 0 | done | 231380e | Gate A recorded; data inventory verified. |
| 1 | done | b5bfefb | Locus-level reproduction, population genetics, tests green. |
| 2 | done | 3b572c2 | Gate B fix-ups: Manhattan GGPv3 exclusion, national burden both definitions, known differences. |
| 3a | done | 8085809 | RNA splice/expression module, export, test green. |
| 3b | done | 8085809 | Farm profiles module, export, test green. |
| 3c | done | 8085809 | Inbreeding module, export, test green. |
| 4a | blocked | — | ESM-2 loaded via transformers OK; ARSI ref mismatch (G vs N) -> dropped; PLCD4/MTRF1/DPF2 scoring pending. |
| 4b | partial | 8085809, b900fea | AUC ~ 0.5, 101 pooled positives; homozygote means must match Phase 1 effects before modelling. |
| 4c | partial | 73c9962, b900fea | Detection summary module; MAF lookup for ~30 simulated QTL pending. |
| 4d | partial | 73c9962, b900fea, e52286c | Recessive power PASS (0.9477 vs 0.947); genotype-class 0.817 vs 0.911. |
| 5 | not started | — | Mating-risk simulator. |
| 6 | not started | — | Web app. |

## Open issues

(a) 4d genotype-class power is 0.817 vs ExtFig1 0.911 (suspected cause: the tail is computed with a 1-df Z-transform shortcut; fix: scipy.stats.ncx2.sf(crit, df=2, nc=lam)).

(b) 4b labels give AUC ~ 0.5 and only 101 positives (the homozygote means must match Phase 1 effects before modelling).

(c) 4a ARSI reference residue mismatch (dropped; score PLCD4, MTRF1 and DPF2 only).

(d) 4c MAF lookup for the ~30 simulated QTL.

(e) the non-additive fraction of inbreeding depression (paper: 23.5%).

## Resume commands

- lint: `uv run ruff check .`
- test: `uv run pytest`
- data: `make data`
- phase modules: `uv run python pipeline/phase<N>_<name>.py`
- phase4 TS: `node --test web/src/lib/power.test.ts`
- site (once built): `make web` / `cd web && npm run dev`

## Rules

- Nothing unvalidated goes into the site.
- Derived data only in git.
- No deploy before Gate D.
