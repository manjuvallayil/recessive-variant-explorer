# Project Status: Recessive Variant Explorer

Last verified: 2026-10-31, HEAD c01bf25. `ruff check .` clean, `uv run pytest` all passing,
`node --test web/src/lib/power.test.ts` 5/5.

## Phase table

| Phase | Status | Commit | Note |
|---|---|---|---|
| 0 | done | 231380e | Gate A recorded; data inventory verified. |
| 1 | done | b5bfefb | Table 1 reproduction. PLCD4 homozygote body weight -108.7 kg (paper -109) under the gen_1 key. |
| 2 | done | 3b572c2 | Gate B fix-ups (Manhattan, national burden, known differences). |
| 3a | done | 8085809 | RNA splice/expression module. Not re-verified this session. |
| 3b | done | 8085809 | Farm profiles module. Not re-verified this session. |
| 3c | done | c01bf25 | Unit bug fixed: -0.740 kg per 1% F (paper -0.77). Recessive attenuation 2.1%. The het+hom 23.6% is a population-structure artefact, not a reproduction of the paper's 23.5%. |
| 4a | blocked | - | ESM-2 loads. RefSeq residue check: PLCD4 A326, DPF2 K216, MTRF1 R341 match. ARSI N301 does not match (G is at 302). Decision needed (see issue c). No scoring code yet. |
| 4b | partial | 3b32ad9 | Key bug fixed. Body-weight flag AUC 0.678 [0.643, 0.713]; CV logistic 0.690; Mahalanobis 0.579; 170 pooled positives. Only DPF2/MUS81 (99 positives) is evaluated alone. Autoencoder, isolation forest, LightGBM and nested CV not done. |
| 4c | partial | 7d6271e | MAF and effects for the 30 truth QTL (MAF 0.96% to 5.06%). Dominance test detects 19/30 at GWS; additive detects 0/30. Recessive and genotype-class tests are not in ExtFig9b. Injected-effect scan not done. |
| 4d | done | b420f9b | Exact noncentral chi-square (ncx2 / Poisson mixture in TS). Cow recessive 0.9477 (paper 0.947), genotype-class 0.9122 (paper 0.911); max residual vs ExtFig1 0.0021. |
| 5 | not started | - | Mating-risk simulator. |
| 6 | not started | - | Web app. No `web/package.json` yet. |

## Verified findings (this session)

- Key semantics: genotype rows are shared between gen_1 and gen_3 (identical order), but each
  phenotype family has its own IID key. `phen_1` uses `gen_1.fam`; `phen_4` uses `gen_3.*.fam`.
  Evidence: PLAG1 (gen_3.7) vs body weight gives t = 216 under the gen_1 key and nothing under
  gen_3.7.fam; DGAT1-candidate (gen_3.6) vs milk/protein/fat gives t = 108, 71, -42 under
  gen_3.6.fam and nothing under gen_1. All gen_3 `.fam` files share one order.
- The earlier 4b AUC of about 0.5 came from joining phen_1 to gen_3 genotypes through the gen_3
  key. Fixed.
- ExtFig2c: its `ID` column does not match either genotype key (agreement at chance). Its
  within-row bw and genotype pairing reproduces the paper's homozygote effects. Not used as a key.
- phen_1 and phen_4 traits are uncorrelated when joined by ID. Expected if the families have
  separate keys. Not used for any inference.

## Open issues

(a) Closed (b420f9b).

(b) Closed for the core fix (3b32ad9). Still open: unsupervised autoencoder, isolation forest,
    LightGBM, nested CV. Sire-grouped CV is not possible (pedigree sires are outside the
    discovery ID space), so related animals can cross folds. Stated in ASSUMPTIONS item 36.

(c) **ARSI needs a decision (STOP, per brief).** XP_002689338.1 has N at residue 301 and G at 302
    (Met counted). PLCD4, DPF2 and MTRF1 numbering matches the brief with Met counted, so the
    ARSI mismatch is not a general numbering offset. Options: (1) you supply the ARSI accession or
    transcript the paper used (Supplementary Table 2); (2) drop ARSI from 4a and score
    PLCD4, DPF2 and MTRF1 only (the previous plan). I have not guessed a numbering.

(d) Closed (7d6271e).

(e) Inbreeding: the unit bug is fixed. The 23.5% attenuation is not reproduced. The primary
    recessive estimate is 2%. Open: a locus-consistent LOSO adjustment for a joint model (each
    locus has its own other_2 window), which may change the het-dummy results.

(f) New: Makefile `data` target runs only phases 0 to 2. Phases 3 and 4 are run by hand.

## Next steps (proposed, pending your decision on c)

1. 4a: score PLCD4, DPF2 and MTRF1 with ESM-2 (wild-type marginal, single pass). ARSI waits for (c).
2. 4b: add the autoencoder and isolation forest, and LightGBM with nested CV.
3. 4c: overlay the empirical detection points on the 4d power curve.
4. Phase 5: mating-risk simulator (TypeScript, unit-tested).
5. Phase 6: scaffold Next.js (`web/package.json`), then sections as each phase is complete.

## Resume commands

- lint: `uv run ruff check .`
- test: `uv run pytest`
- data: `make data` (phases 0 to 2 only)
- phase modules: `uv run python pipeline/phase<N>_<name>.py`
  - phase4 syndrome: `uv run python pipeline/phase4_syndrome.py` (about 95 s)
  - phase4 power: `uv run python pipeline/phase4_power.py`
  - phase4 detection: `uv run python pipeline/phase4_detection.py`
  - phase3 inbreeding: `uv run python pipeline/phase3_inbreeding.py`
- phase4 TS: `node --test web/src/lib/power.test.ts`
- site (once built): `make web` / `cd web && npm run dev`

## Rules

- Nothing unvalidated goes into the site.
- Derived data only in git.
- No deploy before Gate D.
