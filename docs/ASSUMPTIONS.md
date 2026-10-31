# Assumptions and Deviations (Phase 0)

Status: **Phase 0 complete; awaiting Stop Gate A go-ahead.** Every assumption below is either
verified in `work/phase0_summary.json` or explicitly marked as an open question.

## Data and assembly

1. **Assembly split.** DNA-based analyses (`gen_1`–`gen_5`, `gen_9`, `gen_10`, `other_1`,
   `other_8`, `other_9`, `other_15`) are UMD3.1. RNA-seq analyses (`gen_6.1`, `phen_7`) are
   ARS-UCD1.2, as stated in the paper. `other_10` provides the UMD3.1 → ARS-UCD1.2 remap and
   covers 99.553% of `other_1.1` keys (72,092 missing). Do not mix coordinates across assemblies.
2. **A1 is the mutant allele** at all nine lead variants. `bed-reader` counts A1 by default, and
   A1 frequencies reproduce the paper's alternative-allele frequencies (max delta 0.00043).
3. **Genotype coding** is 0 = reference homozygote, 1 = heterozygote, 2 = mutant homozygote for
   `gen_1`–`gen_9` and `gen_10`. This was checked against allele frequencies and class means.
4. **MYH1 is a left-flank proxy.** The paper's candidate is
   `g.30114250_30114327delins…`; the available PLINK variant is the left-flank SNP
   `30114249 C/CCCACCTGGG`, which the paper says was used for imputation. Report it as a proxy,
   not the exact frameshift.
5. **`gen_3` row order.** `gen_3.*` rows are identical to `gen_1.*` at all nine loci, but the
   `.fam` IID order is scrambled. Use row `i` of `gen_3` with row `i` of `gen_1.<locus>.fam`.
   Joining `gen_3` by its own `.fam` IID would silently misassign genotypes.
6. **Separate ID spaces.** Discovery (1–80,027), validation (1–31,580), depletion (1–49,115),
   carrier crosses (1–4,276), RNA (1–350), farm cohorts (1–12, 1–34), chip (1–98,002), and
   simulation (130,201–210,200) are separate. Never join across populations by numeric ID.

## Phenotypes and adjustment

7. **`other_2` trait mapping is swapped relative to the README.** Verified by exact ID-set and
   row-order equality: `other_2.1` = body weight (`phen_1.1`), `other_2.2` = BCS (`phen_1.3`),
   `other_2.3` = stature (`phen_1.2`). Use this mapping, not the README's order.
8. **LOSO adjustment.** Each `other_2.<trait>.<window>.tsv` has column 0 = phenotype ID and
   columns 1–500 = MCMC leave-one-segment-out population-effect samples. Adjusted phenotype =
   `yd_trait − mean(samples)`, then centered. The window nearest the locus under test is the
   relevant one (e.g. PLCD4 uses window 3 or 4; FGD4 uses window 5 or 6).
9. **`rinverse` weighting.** The brief asked to test unweighted, `rinverse`, and `1/rinverse`.
   Mean absolute error versus Table 1 is 1.357, 1.376, and 1.348 kg respectively — effectively
   tied. **Primary choice: unweighted**, because it makes no assumption about the unknown
   semantics of `rinverse` and matches the paper's mixed-model philosophy. `1/rinverse` is the
   sensitivity analysis. This is a documented deviation from “choose the closest” by 0.009 kg.
10. **`other_7` F is a fraction, not a percent.** F ranges 0–3.1 and 179 values exceed 1. Slope
    is −79.24 kg per unit F = −0.792 kg per 1% (paper −0.77). Values >1 are flagged, not clipped.
11. **`other_5` column order is offspring, sire, dam.** Chosen by Mendelian consistency: the two
    offspring-first permutations give 9/15,379 violations (0.0585%), versus 2.0–3.4% for the
    alternatives. `other_5` has 2,197 trios; 4,276 is the union of unique individuals.
12. **Pedigree zeros mean unknown**, not an animal with ID 0. `other_3` and `other_4` are
    headerless `ID SIRE DAM` files with 14–17% zero sires/dams.

## Variant and summary-statistic conventions

13. **Dominance parameterisation.** Class means `(mu0, mu1, mu2) = (mu, mu+a+d, mu+2a)`.
    Homozygote-alt effect = `2a`; heterozygote effect = `a+d`. A pure recessive has `d ≈ −a`.
14. **`other_1.1` sign convention.** `B_A` appears to be the additive effect of A2 (the
    reference allele), so hom-alt (A1/A1) ≈ `−2·B_A` and het ≈ `−B_A + B_D`. This reproduces
    Table 1 (PLCD4: `B_A = +55.675`, predicted hom-alt ≈ −111.4 kg vs observed −108.965 kg;
    predicted het ≈ −1.2 kg vs observed −0.354 kg). Treat this as a working hypothesis and
    confirm against the paper's Methods in Phase 1.
15. **P-values in summary statistics underflow** for strong hits. In Phase 2 compute `−log10 P`
    from Z with `scipy.special.log_ndtr`; do not clip. This is a Phase 2 note, not a Phase 0
    result.
16. **Duplicate keys.** `other_1.1` has 55 duplicate `(CHR,BP)` keys out of 16,128,394 rows.
    `other_9`, `other_10`, and `other_15` also have small duplicate counts. Join on `(CHR,BP)`
    and report coverage; do not assume uniqueness.
17. **`other_15` is a 50k-chip view, not sequence.** It covers 0.284% of `other_1.1` keys. Its
    role is the “chip blind spot” analysis in Phase 2, not a primary association dataset.

## RNA and candidate variants

18. **“4 RNA candidate variants” is unresolved.** The paper's RNA analyses focus on FGD4
    (splice donor), GALNT2 (splice acceptor), and MUS81 (premature stop); DPF2 is a missense
    candidate with expression data in `phen_7.2`. The brief's “4” may refer to FGD4, GALNT2,
    MUS81, and DPF2, or to a different set. Phase 3 should use FGD4/GALNT2/MUS81 for splicing
    and FGD4/GALNT2/MUS81 for expression, and state the ambiguity.
19. **Exact variant matching in `gen_6.1`.** Exact positions are present for PLCD4, FGD4, ARSI,
    DPF2, GALNT2, MYH1, Chr2:22, and the chr14 comparators. MTRF1 is absent; the nearest variant
    (`rs378062563`, chr12:11,464,059) is 78 bp away. MUS81 is absent; the nearest variant
    (`rs41652267`, chr29:44,645,366) is 103 bp away. Do not silently substitute these as exact
    candidates; report the distance and use them only as proxies if Phase 3 requires.
20. **`gen_9.1` README mismatch.** The README describes `gen_9.1` as DPF2/MUS81 chip genotypes;
    the file contains exactly chr29:44,213,160 (DPF2) and chr29:44,646,469 (MUS81), so the
    content matches the intent even though the README wording is loose. `gen_10.1` has only a
    `DPF2` column for the pair; treat it as the DPF2/MUS81 locus.

## Source data

21. **`ExtFig2c` cannot be row-joined.** Its IDs are renumbered 1..79,945 and are not the
    original animal IDs; the rows are not in `phen_1.1` order (body-weight correlation ≈ 0.002).
    Use it for genotype-class distributions and effect-size cross-checks only.
22. **Source-data images are not redistributed.** ExtFig6/7/8 images remain in gitignored
    `dryad/source_data/`; the web build links to the paper/Dryad record instead.
23. **`ExtFig9b` is the simulation truth** for Phase 4c: 624,323 rows, columns
    `CHR BP P_S P_D dominance_QTL`. Use it as primary; only fall back to injected effects if the
    file does not contain the needed MAF/effect information.

## Statistical and export assumptions

24. **Tests are two-sided** unless a one-sided depletion test is explicitly specified (Phase 1
    homozygous depletion). Class means use WLS; SEs are reported; no tuning to match the paper.
25. **Loci are independent** except DPF2/MUS81, which are in strong LD and treated as one locus
    in burden/simulator calculations.
26. **HWE** is assumed for expected homozygote counts and power calculations; observed vs
    expected is reported as a check.
27. **Privacy/export.** Only aggregates are exported. No animal or sire IDs, no individual-level
    rows. Suppress or merge cells with n < 5. Each JSON gets `schema_version` and a manifest
    entry. First-load payload < 1.5 MB gzipped; each file < ~500 KB gz where possible.

## Deviations from the brief

1. **V11 replaced** by the addendum's source-data inventory (10 CSVs + 16 images), as instructed.
2. **V4 weighting choice** is unweighted primary, `1/rinverse` sensitivity (see item 9).
3. **V10 exact-match caveat** for MTRF1 and MUS81 (see item 19); the “4 RNA candidates” phrase
   is documented as ambiguous rather than guessed.
4. **Table 1 values** in the brief are rounded; the paper's parsed values are used in
   `docs/DATA_INVENTORY.md`. No substantive mismatch.
5. **`other_5` row count** corrected to 2,197 trios (4,276 union individuals); the brief's
   “~2,197 rows” was right, but the V3 summary reports the union.
6. **`gen_3` mapping** is by row to `gen_1.fam`, not by its own IID (see item 5).
7. **`ExtFig2c`** is treated as unjoinable at row level (see item 21).

## Open questions for the user

1. Confirm the intended “4 RNA candidate variants” set for Phase 3 (FGD4/GALNT2/MUS81/DPF2?).
2. Confirm unweighted as the primary phenotype estimator, with `1/rinverse` as sensitivity.
3. Confirm that Table 1 parsed values are the authority where the brief's rounded values differ.
4. Confirm the Dryad licence terms for any reuse of source-data images (currently not
   redistributed).

## Stop Gate A checklist

- [x] V1–V11 pass; summary in `work/phase0_summary.json`.
- [x] `docs/DATA_INVENTORY.md` written.
- [x] `docs/ASSUMPTIONS.md` written.
- [x] `uv run ruff check .` passes.
- [x] `uv run pytest` passes.
- [x] `git ls-files` / `git status --ignored` show no PDF, no `dryad/` files, no source-data
      images.
- [ ] User go-ahead to start Phase 1.
