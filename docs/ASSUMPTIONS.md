# Assumptions and Deviations (Phase 0)

Status: **Phase 0 complete; Gate A approved; Phase 2 complete; Gate B not yet approved.** Every
assumption below is either verified in `work/phase0_summary.json`, updated through Phase 2, or
explicitly marked as an open question.

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
5. **Phenotype-specific `.fam` mapping for `gen_1`/`gen_3`.** `gen_3.*` genotype rows are
   identical to `gen_1.*` at all nine lead loci (identity = 1.0), but the `.fam` IID labels are
   permuted between the two sets. `phen_1.*` IDs align with `gen_1.*.fam`; `phen_4.*` IDs align
   with `gen_3.*.fam`. Join `phen_1` to `gen_1` by IID and `phen_4` to `gen_3` by IID. Do not use
   one `.fam` file for both phenotype families: doing so attenuates or reverses effects. Evidence:
   PLCD4 bodyweight raw effect is about −80.6 kg with `gen_1` mapping (B-adjusted −110.9 kg) but
   −7.9 kg with `gen_3` mapping; PLCD4 milk effect is about −2.42 with `gen_3` mapping (validation
   −2.85) but −0.25 with `gen_1` mapping. The earlier “always use the `gen_1` row mapping” note
   was wrong for `phen_4`.
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
   semantics of `rinverse` and matches the paper's mixed-model philosophy. `weights = rinverse`
   is kept as a labelled sensitivity run (Gate A decision 1).
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

18. **RNA targets are FGD4, GALNT2, and MUS81** (user decision, Gate A). Use `phen_7.1` splice
    junctions for FGD4 and GALNT2, and `phen_7.2` expression for FGD4, GALNT2, and MUS81. DPF2
    has no RNA phenotype. For MUS81, whose exact variant is absent from `gen_6.1`, use the exact
    DPF2 variant (`chr29:44213160`) as a clearly labelled proxy. The paper reports candidate
    DPF2↔MUS81 R² ≈ 0.98; in the 350 RNA animals the nearest MUS81 variant (`rs41652267`,
    103 bp away) has r² = 0.018 with DPF2, so the proxy is a caveat, not a validated tag.
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
28. **Stature ×5.0 cm factor was inferred** by scaling NZ TOP linear score units to cm
    (`×5.0 cm / score`) to match Table 1 across all nine loci; MAE = 0.058 cm. Both raw TOP-score
    and scaled cm values are preserved in `web/public/data/effects.json`.
29. **`other_8.1.txt` `P_A` underflow strings.** 506 `chr14` rows had `P_A` written as the literal
    string `"<1e-308"`, causing those rows to drop during DuckDB `read_csv` (DOUBLE type);
    `ensure_parquet_files` now reads it via
    `TRY_CAST(REPLACE(P_A, '<', '') AS DOUBLE) AS P_A` and unlinks/rebuilds `other_8.1.parquet` when
    the row count is below 16,128,364. After the rebuild, `chr14` standard-additive leads re-derive
    (e.g. GGPv3 causal SNP `chr14:25,015,640`, `Z_A = 46.94`).

## Regional plots, candidate-gene mapping, and Gate B (Phase 2 complete)

30. **PLAG1 chr14 regional window label.** The `plag1` regional plot is labelled
    `chr14 QTL (PLAG1 regulatory domain; gene ~263 kb outside the genotyped window)`, `gene` is
    `null`, LD colouring is disabled (`ld_available = false`, `r2_x1000 = null`), and the window
    must not be claimed to contain the PLAG1 gene.
31. **PLCD4 explicit gene entry.** Ensembl ARS-UCD1.2 merges PLCD4 into the CNOT9 readthrough
    locus; an explicit PLCD4 gene entry was added to the regional gene track (NCBI:509076).

## Stop Gate B checklist

- [ ] Full Gate B package re-posted to the user: full REPRODUCTION table (paper vs estimator A vs
      B vs `other_1.*` lead row, with n per genotype class), qtl_contrast table, export file sizes
      (raw and gzipped), chip_vs_sequence result, national_burden result vs paper (~1% and ~9,700),
      depletion and carrier-cross tests, and list of anything that disagrees with the paper.

## Known differences

- **Chr5 `-log10 P`:** the variant giving **-log10 P = 131.30** is the indel `chr5:106,266,836`
  (`Z_A = 24.45`, file P = 5.05e-132); the paper-matching lead is the SNP `chr5:106,250,171` A/G
  (file P = 3.58619e-127 -> rounds to `3.59e-127`, -log10 = 126.44). Different lead variants in the
  same CCND2 peak; file P is used where it does not underflow.
- **Chr14 `-log10 P`:** Manhattan y-axis follows the paper's definition (re-introduced GGPv3 sites
  excluded); exact value **-log10 P = 441.3** at BovineSNP50 SNP `chr14:25,307,116` (two-sided P =
  5.03e-442 ~ 4.97e-442); re-introduced GGPv3 causal SNP `chr14:25,015,640` gives **-log10 P = 480.2`;
  the 480.2 vs 441.3 difference is the GGPv3 re-introduction / lead change.

## Phase 4b syndrome-signal decisions

32. **Key check (verified, not assumed).** Genotype rows are shared across `gen_1.*` and `gen_3.*`,
    but each phenotype family has its own IID key. Checks: PLAG1 `gen_3.7` vs body weight
    (phen_1.1) gives slope 50 kg/allele, t = 216, under the `gen_1.fam` key and nothing under
    `gen_3.7.fam`. The DGAT1-candidate `gen_3.6` vs milk/protein/fat (phen_4.1-3) gives t = 108,
    71 and -42 under `gen_3.6.fam` and nothing under `gen_1.fam`. So `phen_1` -> `gen_1.fam` and
    `phen_4` -> `gen_3.fam` (item 5 stands). All `gen_3.*.fam` files share one order.
33. **Bug fixed in 4b.** Labels were built from `gen_3.*.fam` and joined to `phen_1` body weight,
    which mixes two keys. The AUC was about 0.5 for this reason. Labels now come from the
    `gen_1` windows. Check: PLCD4 homozygote body weight is -108.7 kg, against -109 in the paper.
34. **Labels.** Homozygous-mutant status (A1 = mutant, as in item 2). DPF2 and MUS81 form one
    locus (homozygous for either). The pooled "any" label covers PLCD4, FGD4, GALNT2 and
    DPF2/MUS81. Groups with fewer than 30 positives (PLCD4 20, FGD4 27, GALNT2 25 after
    dropping records without BCS) are not evaluated alone; they are only in "any".
35. **Scores.** B0 = -body weight (lower weight flagged; direction from Phase 1 effect signs,
    not tuned). B1 = Mahalanobis distance on bw/stature/BCS (unsupervised). S = logistic
    regression, 5-fold stratified out-of-fold predictions, scaler fit inside folds. AUC and AUPRC
    use bootstrap CIs with B = 1000.
36. **Known limitations.** Sire-grouped CV is not possible: pedigree sires are outside the
    discovery ID space (only dams partly match). Related animals may therefore cross folds. The
    unsupervised autoencoder, isolation forest, LightGBM and nested CV are not yet implemented.

## Phase 3c inbreeding decisions

37. **Unit bug fixed.** Per-1% slope was computed as slope x 100 (gave -7,402 kg). It is now
    slope / 100 = -0.740 kg per 1% F (95% bootstrap CI -0.764 to -0.714), against the paper's -0.77.
38. **F units assumed to be fractions.** F ranges 0 to 3.1; 179 animals exceed 1. Excluding them
    does not change the attenuation (23.8% vs 23.6%). Genomic cross-check of the scale is
    inconclusive (one 4,447-variant window; correlation 0.007 vs a -0.003 shuffled control).
39. **Attenuation: recessive estimate is homozygote-only (about 2%).** The het + hom model gives
    about 24% (CI 22-25%), but this comes from heterozygote dummies on raw yd absorbing population
    structure (PLCD4 het is +24 kg raw, -1 kg LOSO-adjusted). It is not evidence of a recessive
    effect and must not be presented as reproducing the paper's 23.5%. The paper's figure is a
    variance fraction from BayesC0; this is an OLS slope attenuation. Not directly comparable.
    LOSO adjustment is locus-specific and cannot be applied to a joint model without a decision;
    this is open.

## Phases 4a to 6 decisions

40. **4a protein model (one attempt).** ESM-2 150M masked-marginal LLR, one pass per position, for
    PLCD4 A326T, MTRF1 R341W and DPF2 K216R. Wild-type residues were checked against RefSeq at each
    position (all match). Percentiles: PLCD4 56, MTRF1 90, DPF2 66. Concordance with the observed
    homozygote effects failed: the model ordering (MTRF1 > DPF2 > PLCD4) is the reverse of the
    observed ordering (PLCD4 > DPF2 > MTRF1). With n = 3 this is not testable, so it is reported as
    failed and not used as evidence. ARSI dropped (residue 301 is N). Not scorable: MUS81 p.Gly70*, splicing.
41. **4b per-locus report.** Single-locus statistics use every genotyped animal with body weight
    (no stature or BCS needed). Homozygote counts: PLCD4 21, FGD4 28, GALNT2 29, DPF2/MUS81 110.
    The complete-case pooled "any" label is 170; with body weight only, 187.
    - Why the pooled count was about 101: the earlier rule pooled only loci with at least 30
      positives, which left only DPF2/MUS81 (99 complete-case). The rule now pools all four loci.
    - Why PLCD4 AUROC is 0.79, not 0.9+: under equal-variance normal classes the AUROC is set by the
      effect size. The shift from reference to homozygote is about 82 kg against a 63 kg SD, which
      caps AUROC near 0.82 (observed 0.79, CI 0.68 to 0.88). Reaching 0.9 would need a shift of about 114 kg.
    - Causes checked and ruled out as the main cause: genotype coding (A1 mutant, class means match
      Phase 1 and ExtFig2c homozygote means); missing phenotypes (dropped 1 PLCD4 and 4 GALNT2
      homozygotes at the complete-case step, now not used for per-locus statistics).
    - Depletion: observed homozygotes are 0.33 to 0.61 of the Hardy-Weinberg expectation (pooled
      q, sum 187 observed vs 429 expected). Population mixing would push observed above the pooled
      expectation, so this shortfall is not explained by structure alone.
42. **4c threshold and table.** Two-sided P < 5e-8 on ExtFig9b P_S (additive) or P_D (dominance).
    Dominance detects 19 of 30 truth QTL (MAF < 2%: 2 of 9; 2-3%: 4 of 6; 3-4%: 5 of 7; >= 4%: 8 of
    8). Additive detects 0 of 30; the closest additive P is 5.9e-8 (chr9). The MAF range is 0.96% to
    5.06%; the top bin is labelled >= 4%.
43. **Phase 2 determinism.** Manhattan thinning picked among tied |Z| values by DuckDB thread order.
    Winners are now chosen by a window with explicit tie-breaks (|Z|, then position, then sign), and
    the ARS mapping uses an explicit order. Imputation summaries differ only in the 15th significant
    digit (parallel float sums).
44. **QTL counts.** At candidate variants, 7 of 8 are below genome-wide significance in the additive
    model (PLAG1 is the exception, as the paper says). Within ±1 Mb, 6 of 8: the DPF2/MUS81 window
    holds an additive peak at 45.06 Mb, about 0.8 Mb from the candidate, which the candidate does not carry.
45. **Carrier union (hook).** Holstein-Friesian carriers of at least one of the seven variants:
    43.5% with DPF2 and MUS81 separate (as the paper counts them), 35.3% with them as one locus.
    Independence is assumed. The paper reports more than 40%.
46. **Phase 5 simulator.** Transmission 0, 0.5, 1 or the breed-weighted q for non-carrier, carrier,
    affected or unknown parents. Calf homozygote = product of transmissions; independence across loci;
    DPF2 proxy for DPF2/MUS81; purebred depletion ratios (primary breed) applied to homozygote
    probabilities. Inputs from validated exports only, plus one aggregate from other_11 (mean
    breed-weighted q, dam breed composition as proxy). ARSI excluded (no breed frequency).
    Herd: "avoid carrier x carrier" removes the (1,1) mating and renormalises; checked against the
    closed form in unit tests. The DPF2/MUS81 heterozygote effect is positive (+2.3 kg), so
    the net kilogram change can be positive; homozygote and net changes are shown separately.
47. **Calibration gap.** Random mating at national frequencies gives about 686 PLCD4 homozygous births
    per year (331 with depletion), against 1,034 in the national export. The export uses observed
    sire genotypes, and carrier sires are over-represented (sire amplification), so the export sits
    above a random-mating calculation. This is reported as a gap, not forced to match.
48. **Site scope (Phase 6).** Only validated exports are read; each is schema-checked at build time
    (zod) and checked for individual-identifier keys. The site shows the 3c slope (-0.74 kg per 1%)
    only, not the attenuation figure. Lazy loading per section; first-load JS and CSS about 256 KB gzipped.
    The Dryad dataset licence was not stated on the record text checked; the footer says to verify it.

## Deviations from the brief

1. **V11 replaced** by the addendum's source-data inventory (10 CSVs + 16 images), as instructed.
2. **V4 weighting choice** is unweighted primary, `weights = rinverse` labelled sensitivity (see item 9).
3. **V10 exact-match caveat** for MTRF1 and MUS81 (see item 19); MUS81 RNA analysis uses DPF2
   as a labelled proxy with a low-LD caveat (see item 18).
4. **Table 1 values** in the brief are rounded; the paper's parsed values are used in
   `docs/DATA_INVENTORY.md`. No substantive mismatch.
5. **`other_5` row count** corrected to 2,197 trios (4,276 union individuals); the brief's
   “~2,197 rows” was right, but the V3 summary reports the union.
6. **`other_15` coverage direction** clarified: 100% of `other_15` rows are in `other_1.1`;
   0.284% is the reverse (chip rows as a fraction of sequence rows).
7. **MUS81 RNA proxy** is DPF2 with an explicit low-LD caveat (r² = 0.018 in 350 animals).
8. **Phenotype-specific genotype mapping** follows item 5: `phen_1` → `gen_1.fam`; `phen_4` → `gen_3.fam`.
9. **`ExtFig2c`** is treated as unjoinable at row level (see item 21).

## Open questions for the user

1. ~~Confirm the intended “4 RNA candidate variants” set~~ — resolved at Gate A: FGD4, GALNT2,
   MUS81; DPF2 is the MUS81 proxy.
2. ~~Confirm unweighted as the primary phenotype estimator~~ — resolved at Gate A: unweighted
   primary, `rinverse` labelled sensitivity.
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
- [x] User go-ahead to start Phase 1.
