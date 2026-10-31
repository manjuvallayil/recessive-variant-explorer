# Data Inventory (Phase 0)

Status: **Phase 0 complete; V1–V11 pass.** Generated from `work/phase0_summary.json` by
`pipeline/phase0_verify.py` on 2026-10-31.

This is an inventory of the Dryad data used by the project, not a redistribution. Raw files
live under `dryad/` and are gitignored.

## Provenance and licence

- Paper: Reynolds et al., “Non-additive association analysis using proxy phenotypes identifies
  novel cattle syndromes”, *Nature Genetics* 53, 949–954 (2021).
  DOI: <https://doi.org/10.1038/s41588-021-00872-5>
- Data: Dryad, DOI: <https://doi.org/10.5061/dryad.kwh70rz1z>
- The paper PDF is in `paper/` (gitignored) and is not redistributed.
- Source-data figure images are in `dryad/source_data/` (gitignored) and are not copied into the
  web build. Link to the paper/Dryad record instead.
- Independent learning project; not affiliated with LIC. LIC holds patent applications on these
  mutations. No commercial framing.

## Directory layout

```text
dryad/
  gens/         PLINK genotypes and small genotype tables
  phens/        phenotype tables
  other/        summary statistics, pedigrees, breed/inbreeding tables
  source_data/  paper source data for Extended Data figures
pipeline/       Python verification and analysis code
work/           local derived data and summaries (gitignored)
docs/           this inventory, assumptions, and later reproduction notes
```

## Population and ID spaces

Numeric ID ranges overlap between populations. **Never join across populations by ID alone.**

| Population | Genotypes | Phenotypes | Animals | ID range | Notes |
|---|---|---|---:|---|---|
| Discovery | `gen_1.*` | `phen_1.*` | 80,027 genotyped; 79,945 body weight; 75,041 stature; 75,617 BCS | 1–80,027 | 82 genotyped animals lack body-weight phenotype. |
| Validation | `gen_2.*` | `phen_2.*` | 31,580 | 1–31,580 | First-lactation milk/protein/fat. |
| Pleiotropy / mixed model | `gen_3.*` | `phen_4.*` | 80,027 genotyped; 79,659 milk | 1–80,027 | 368 genotyped animals lack `phen_4.1`. |
| Homozygous depletion | `gen_4.*` + `other_6` | — | 49,115 purebreds | 1–49,115 | 29,771 HF, 19,344 Jersey (16/16ths). |
| Carrier crosses | `gen_5.*` + `other_5` | — | 4,276 individuals in 2,197 trios | 1–4,276 | Offspring/sire/dam; union of trio columns = 4,276. |
| RNA-seq | `gen_6.1` | `phen_7.*` | 350 | 1–350 | Mammary RNA-seq; DNA analyses UMD3.1, RNA ARS-UCD1.2. |
| GALNT2 farm cohort | `gen_7.1` | `phen_6.*` | 12 | 1–12 | 6 affected, 6 controls. |
| Research farm | `gen_8.1` | `phen_5.*` | 34 | 1–34 | PLCD4 8, FGD4 9, DPF2/MUS81 8, control 9. |
| Chip validation | `gen_9.1` | — | 98,002 | 1–98,002 | Two variants only: DPF2 and MUS81. |
| National burden | `gen_10.1` + `other_11` | — | 5,550 sires; 2,799,022 animals | sires 1–5,550 | Birth years 2010–2019. |
| Simulation | `gen_11.1`, `gen_12.1` | `phen_8.1` | 80,000 | 130,201–210,200 | `gen_11` 800,000 variants; `gen_12` 31,106. |
| Prospective calves | — | `phen_3.1` | 568 | — | Genotype group in `VARIANT`. |
| Inbreeding | — | `other_7` + `phen_1.1` | 68,578 | 1–68,578 | F and population-effect estimate. |

## Phenotypes

| File(s) | Rows | Columns | Notes |
|---|---:|---|---|
| `phen_1.1.txt` | 79,945 | `ID yd_trait rinverse` | Body weight. |
| `phen_1.2.txt` | 75,041 | `ID yd_trait rinverse` | Stature. |
| `phen_1.3.txt` | 75,617 | `ID yd_trait rinverse` | Body condition score (1–10). |
| `phen_2.1–2.3.txt` | 31,580 each | `ID yd_trait rinverse` | Milk, protein, fat. |
| `phen_3.1.txt` | 568 | `animal VARIANT HERD age breed_hf breed_j weight` | Prospectively genotyped calves. |
| `phen_4.1–4.16.txt` | 79,659 for 1–3; varies | `ID yd_trait rinverse` | 19 traits total (see trait map). |
| `phen_5.1.txt` | 34 | `ID Spine Tibia Height Chest Claw_length ToeAngle` | Research-farm anatomy; `NA` present. |
| `phen_5.2.txt`, `phen_5.3.txt` | 34 | biochemistry columns | February / June 2019; headers contain spaces/units. |
| `phen_5.4.txt` | 34 | weight time series | Date columns. |
| `phen_6.1.txt` | long format | `ID measure_date trait value` | GALNT2 cohort biochemistry. |
| `phen_6.2.txt` | 12 | weights at four dates + age | GALNT2 cohort. |
| `phen_7.1.txt` | 350 | `KEY` + 17 FGD4 junctions + 16 GALNT2 junctions | Logit splice-efficiency phenotypes. |
| `phen_7.2.txt` | 350 | `ID FGD4 GALNT2 MUS81` | Gene expression. |
| `phen_8.1.txt` | 80,000 | `anml_key yd_trait rinverse intercept` | Simulated phenotype, adjusted for 30 recessive QTL. |

### Trait map

- `phen_1`: 1 body weight, 2 stature, 3 BCS.
- `phen_2`: 1 milk, 2 protein, 3 fat.
- `phen_4`: 1 milk, 2 protein, 3 fat, 4 adaptability, 5 temperament, 6 farmer opinion,
  7 capacity, 8 rump angle, 9 rump width, 10 legs, 11 udder support, 12 front udder,
  13 front teat, 14 rear teat, 15 udder overall, 16 conformation.
- **`other_2` mapping is swapped relative to the Dryad README** (verified by exact ID-set and
  row-order equality):
  - `other_2.1.*` = body weight (`phen_1.1`, 79,945)
  - `other_2.2.*` = BCS (`phen_1.3`, 75,617)
  - `other_2.3.*` = stature (`phen_1.2`, 75,041)

### LOSO windows in `other_2`

| Window | Region | Window | Region |
|---:|---|---:|---|
| 1 | Chr2_20_30Mb | 8 | Chr7_55_65Mb |
| 2 | Chr2_15_25Mb | 9 | Chr12_5_15Mb |
| 3 | Chr2_100_110Mb | 10 | Chr14_20_30Mb |
| 4 | Chr2_105_115Mb | 11 | Chr19_25_35Mb |
| 5 | Chr5_70_80Mb | 12 | Chr28_0_10Mb |
| 6 | Chr5_75_85Mb | 13 | Chr29_40_50Mb |
| 7 | Chr7_60_70Mb | | |

Each `other_2.<trait>.<window>.tsv` is headerless, about 330 MB, 501 columns: column 0 is the
phenotype ID, columns 1–500 are MCMC leave-one-segment-out population-effect samples. Adjusted
phenotype = `yd_trait - mean(samples)`, then centered.

## Genotypes

| Family | Animals | Variants | Files | Notes |
|---|---:|---:|---|---|
| `gen_1.*` | 80,027 | 6,729–14,564 per window | 9 PLINK sets | Sequence-resolution windows around lead loci. |
| `gen_2.*` | 31,580 | 1 each | 12 PLINK sets | Prioritised mutations, validation population. |
| `gen_3.*` | 80,027 | 1 each | 12 PLINK sets | Same variants as `gen_2`, `phen_4` population. |
| `gen_4.*` | 49,115 | 1 each | 7 PLINK sets | Homozygous depletion; PLCD4, FGD4, MTRF1, MYH1, GALNT2, DPF2, MUS81. |
| `gen_5.*` | 4,276 | 1 each | 7 PLINK sets | Carrier crosses; same seven variants as `gen_4`. |
| `gen_6.1` | 350 | 17,180,153 | PLINK set | RNA-seq animals; stream/seek only. |
| `gen_7.1.txt` | 12 | 2 columns | `ID Status` | GALNT2 case/control. |
| `gen_8.1.txt` | 34 | 2 columns | headerless `ID group` | Research-farm cohort. |
| `gen_9.1` | 98,002 | 2 | PLINK set | DPF2 and MUS81 only. |
| `gen_10.1.txt` | 5,550 | 7 columns | `KEY PLCD4 FGD4 MTRF1 MYH1 GALNT2 DPF2` | Sire genotypes, 0/1/2. |
| `gen_11.1` | 80,000 | 800,000 | PLINK set | Simulated. |
| `gen_12.1` | 80,000 | 31,106 | PLINK set | Simulated marker subset. |

PLINK genotype coding: 0 = reference homozygote, 1 = heterozygote, 2 = mutant homozygote.
`bed-reader` counts A1 alleles by default, and A1 is the mutant allele at all nine lead variants.

### `gen_3` row-order caveat

`gen_3.*` row order is **identical to `gen_1.*` at all nine lead loci** (identity = 1.0), but the
`.fam` IID order differs from `gen_1`. To attach `gen_3` genotypes to discovery animals, map
row `i` to `gen_1.<locus>.fam` row `i`; do not join `gen_3` by its own `.fam` IID.

## Summary statistics and other tables

| File(s) | Rows | Columns | Notes |
|---|---:|---|---|
| `other_1.1.txt` | 16,128,394 | `CHR BP A1 A2 MAF B_A B_D Z_A Z_D P_A P_D` | Body-weight dominance GWAS, UMD3.1. |
| `other_1.2.txt`, `other_1.3.txt` | ~16.1M | same | Stature, BCS. |
| `other_8.1.txt` | 16,128,364 | `CHR BP A1 A2 MAF B_A Z_A P_A` | Body-weight standard-additive GWAS. |
| `other_8.2.txt`, `other_8.3.txt` | ~16.1M | same | Stature, BCS. |
| `other_9.txt` | 16,128,756 | `CHR BP AR2` | Imputation allelic R². |
| `other_10.txt` | 16,056,778 | `CHR.UMD BP.UMD CHR.ARS BP.ARS` | UMD3.1 → ARS-UCD1.2 remap. |
| `other_15.txt` | 45,915 | `CHR BP P_D` | BovineSNP50 positions and body-weight dominance P. |
| `other_2.*` | ~79,945 per trait/window | 501 | LOSO population-effect samples; 39 files. |
| `other_3.txt` | 172,572 | 3 headerless | Validation-population pedigree `ID SIRE DAM`; 0 = unknown. |
| `other_4.txt` | 314,677 | 3 headerless | `phen_4` population pedigree `ID SIRE DAM`; 0 = unknown. |
| `other_5.txt` | 2,197 | 3 headerless | Carrier-cross trios. Column order = offspring, sire, dam. |
| `other_6.txt` | 49,115 | `ID breed_hf breed_j` | Purebred fractions in sixteenths. |
| `other_7.txt` | 68,578 | 3 headerless | `ID F value`; F is a fraction, not percent. |
| `other_11.txt` | 2,799,022 | `IID SIRE_ID MGS_ID birth_year dam_HF_16 dam_J_16 mgd_HF_16 mgd_J_16` | National-burden population. |
| `other_12.*` | — | — | QMSim parameter/seed files. |
| `other_13.txt` | — | — | Simulated 31K marker effects; not used in Phase 0. |
| `other_14.txt` | — | — | ARS-UCD1.2 positions used for RNA-seq masking; not used in Phase 0. |

`other_5` column order was determined by Mendelian consistency across the seven variants:
offspring-first `(0,1,2)` and `(0,2,1)` give 9 violations / 15,379 trios (0.0585%), versus
2.0–3.4% for the other permutations. The 4,276 genotyped carrier-cross animals equal the union
of the 2,197 offspring, 456 sires, and 1,683 dams in `other_5`.

## Source data (`dryad/source_data/`)

| File | Rows | Columns | Contents | Figure |
|---|---:|---|---|---|
| `source_data.ExtFig1.csv` | 1,400 | `N model Power` | Power curves | Ext. Fig. 1 |
| `source_data.ExtFig2a.csv` | 100 | `Model CHR BP MAF BETA GENE` | Top additive effects | Ext. Fig. 2a |
| `source_data.ExtFig2c.csv` | 79,945 | 17 genotype/bodyweight columns | Genotype-class body weights | Ext. Fig. 2c |
| `source_data.ExtFig3.csv` | 7,862 | `CHR BP P R2` | Regional P / LD | Ext. Fig. 3 |
| `source_data.ExtFig4a.csv` | 25,330 | `chr BP Junction pvalue` | FGD4 splice junctions | Ext. Fig. 4a |
| `source_data.ExtFig4b.csv` | 56,040 | `chr BP Junction pvalue` | GALNT2 splice junctions | Ext. Fig. 4b |
| `source_data.ExtFig4c.csv` | 1,490 | `chr BP GE pvalue` | FGD4 expression | Ext. Fig. 4c |
| `source_data.ExtFig4d.csv` | 3,736 | `chr BP GE pvalue` | GALNT2 expression | Ext. Fig. 4d |
| `source_data.ExtFig4e.csv` | 856 | `chr BP GE pvalue` | MUS81 expression | Ext. Fig. 4e |
| `source_data.ExtFig9b.csv` | 624,323 | `CHR BP P_S P_D dominance_QTL` | Simulation truth | Ext. Fig. 9b |
| `source_data.ExtFig6a1–a8.jpg` | — | images | Phenotype images | Ext. Fig. 6 |
| `source_data.ExtFig7a–e.tif` | — | images | Histology | Ext. Fig. 7 |
| `source_data.ExtFig8a.jpg`, `8b1.jpg`, `8b2.jpg` | — | images | Phenotype images | Ext. Fig. 8 |

First three sample rows (whitespace-separated for readability):

```text
source_data.ExtFig1.csv: N model Power
500 Cow-Additive 0
1000 Cow-Additive 0
1500 Cow-Additive 0

source_data.ExtFig2a.csv: Model CHR BP MAF BETA GENE
Standard 10 16254330 0.285 1.608 .
Standard 10 81499004 0.276 1.605 .
Standard 11 3537255 0.299 1.601 .

source_data.ExtFig2c.csv: ID Bodyweight.CHR2_22MB Genotype.CHR2_22MB Bodyweight Genotype
Bodyweight.FGD4 Genotype.FGD4 Bodyweight.ARSI Genotype.ARSI Bodyweight.MTRF1 Genotype.MTRF1
Bodyweight.GALNT2 Genotype.GALNT2 Bodyweight.DPF2 Genotype.DPF2 Bodyweight.PLAG1 Genotype.PLAG1
1 -33.466 0/1 -30.999 0/0 -31.928 0/0 -33.095 0/0 -32.115 0/0 -30.424 0/0 -26.160 0/0 -28.762 0/1
2 1.175 0/0 -1.139 0/0 -0.610 0/0 0.295 0/0 -2.768 0/0 0.574 0/0 -0.991 0/0 -17.783 0/0
3 28.180 0/0 28.347 0/0 29.008 0/0 27.492 0/0 28.551 0/0 27.417 0/0 26.524 0/0 29.317 0/1

source_data.ExtFig3.csv: CHR BP P R2
19 29069945 0.444952 0.000458
19 29070574 0.851851 0.002172
19 29070843 0.330524 0.029013

source_data.ExtFig4a.csv: chr BP Junction pvalue
chr5 77508650 1 0.000569
chr5 77508676 1 0.934924
chr5 77508919 1 0.937357

source_data.ExtFig4b.csv: chr BP Junction pvalue
chr28 1027186 1 0.637812
chr28 1028814 1 0.802385
chr28 1028989 1 0.290769

source_data.ExtFig4c.csv: chr BP GE pvalue
chr5 77508650 FGD4 0.403875
chr5 77508676 FGD4 0.407580
chr5 77508919 FGD4 0.424268

source_data.ExtFig4d.csv: chr BP GE pvalue
chr28 1027186 GALNT2 0.504001
chr28 1028814 GALNT2 0.963174
chr28 1028989 GALNT2 0.085804

source_data.ExtFig4e.csv: chr BP GE pvalue
chr29 44645469 MUS81 2.383779e-20
chr29 44544915 MUS81 9.103718e-02
chr29 44545330 MUS81 2.244336e-02

source_data.ExtFig9b.csv: CHR BP P_S P_D dominance_QTL
1 1 0.515507 0.621710 False
1 2930 0.147199 0.847707 False
1 5840 0.348113 0.939974 False
```

**Important:** `ExtFig2c` IDs are renumbered `1..79,945`, not the original animal IDs. The rows
are not in `phen_1.1` order (body-weight correlation ≈ 0.002), so this file can be used for
distribution/class-mean comparisons only, not row-level joins.

## Locus and variant map

| Locus | QTL | Candidate mutation | Consequence | Discovery files | Depletion/cross files | `gen_1` window | A1/A2 | A1 freq | Paper AAF (all) |
|---|---|---|---|---|---|---:|---|---:|---:|
| Chr2:22 | 2:22 Mb | `rs211074296` (no candidate) | — | `gen_1.1`, `gen_2.1`, `gen_3.1` | — | 1.1 | A/G | 0.0736 | 0.074 |
| PLCD4 | 2:107 Mb | `g.107313998G>A` | p.Ala326Thr | `gen_1.2`, `gen_2.2`, `gen_3.2` | `gen_4.1`, `gen_5.1` | 1.2 | A/G | 0.0242 | 0.024 |
| FGD4 | 5:78 Mb | `g.77632752C>T` | splice donor | `gen_1.3`, `gen_2.3`, `gen_3.3` | `gen_4.2`, `gen_5.2` | 1.3 | T/C | 0.0241 | 0.024 |
| ARSI | 7:64 Mb | `g.63649071C>T` | p.Gly301Ser | `gen_1.4`, `gen_2.4`, `gen_3.4` | — | 1.4 | T/C | 0.0717 | 0.072 |
| MTRF1 | 12:11 Mb | `g.11463981G>A` | p.Arg341Trp | `gen_1.5`, `gen_2.5`, `gen_3.5` | `gen_4.3`, `gen_5.3` | 1.5 | A/G | 0.0449 | 0.045 |
| MYH1 | 19:30 Mb | `g.30114250_30114327delins…`; left-flank proxy `30114249 C/CCCACCTGGG` | frameshift | `gen_1.7`, `gen_2.9`, `gen_3.9` | `gen_4.4`, `gen_5.4` | 1.7 | C/CCCACCTGGG | 0.0422 | 0.042 |
| GALNT2 | 28:1 Mb | `g.1312334G>A` | splice acceptor | `gen_1.8`, `gen_2.10`, `gen_3.10` | `gen_4.5`, `gen_5.5` | 1.8 | A/G | 0.0332 | 0.033 |
| DPF2 | 29:44 Mb | `g.44213160A>G` | p.Lys216Arg | `gen_1.9`, `gen_2.11`, `gen_3.11` | `gen_4.6`, `gen_5.6` | 1.9 | G/A | 0.0392 | 0.039 |
| MUS81 | 29:44 Mb | `g.44645469G>T` | p.Gly70* | `gen_1.9`, `gen_2.12`, `gen_3.12` | `gen_4.7`, `gen_5.7` | 1.9 | T/G | 0.0396 | 0.040 |
| chr14:1802265 | comparator | `rs109234250` | — | `gen_1.6`, `gen_2.6`, `gen_3.6` | — | 1.6 | A/G | — | — |
| chr14:25015640 | PLAG1 comparator | `rs109815800` | — | `gen_1.6`, `gen_2.7`, `gen_3.7` | — | 1.6 | T/G | — | — |
| chr14:26064929 | comparator | `rs109686357` | — | `gen_1.6`, `gen_2.8`, `gen_3.8` | — | 1.6 | C/T | — | — |

`gen_9.1` contains chr29:44,213,160 (DPF2) and chr29:44,645,469 (MUS81) at exact positions.
`gen_10.1` has a single `DPF2` column that stands in for the DPF2/MUS81 pair.

## Paper cross-check: Table 1 vs brief

Table 1 was parsed from the paper PDF with `pdftotext -layout`. The brief's Section 3 values are
mostly rounded versions of the paper's values; no substantive mismatch was found. The paper is
the authority.

| Locus | Table 1 hom-alt BW (kg) | Brief | Delta | Table 1 het BW | Table 1 AAF (HF/J/all) | Flag |
|---|---:|---:|---:|---:|---|---|
| Chr2:22 | −11.8 | −11.8 | 0.0 | −0.007 | 0.124 / 0.002 / 0.074 | — |
| PLCD4 | −108.965 | −109.0 | +0.035 | −0.354 | 0.039 / 0.001 / 0.024 | rounding |
| FGD4 | −48.079 | −48.1 | +0.021 | +0.834 | 0.039 / 0.001 / 0.024 | rounding |
| ARSI | −11.981 | −12.0 | +0.019 | +0.564 | 0.124 / 0.002 / 0.072 | rounding |
| MTRF1 | −19.464 | −19.5 | +0.036 | +1.078 | 0.002 / 0.114 / 0.045 | rounding |
| GALNT2 | −112.008 | −112.0 | −0.008 | −0.837 | 0.055 / 0.001 / 0.033 | rounding |
| DPF2 | −53.092 | −53.1 | +0.008 | +4.059 | 0.066 / 0.001 / 0.039 | rounding |
| MUS81 | −51.911 | −51.9 | −0.011 | +4.126 | 0.066 / 0.001 / 0.040 | rounding |
| MYH1 | −4.643 | −4.6 | −0.043 | +0.798 | 0.006 / 0.109 / 0.042 | n.s.; rounding |

| Locus | Table 1 hom-alt stature (cm) | Brief | Delta |
|---|---:|---:|---:|
| PLCD4 | −5.421 | −5.4 | −0.021 |
| GALNT2 | −8.568 | −8.6 | +0.032 |
| FGD4 | −3.041 | −3.0 | −0.041 |
| DPF2 | −2.565 | −2.6 | +0.035 |

BCS Table 1 hom-alt values (paper): Chr2:22 −0.044; PLCD4 −0.284; FGD4 −0.238; ARSI −0.028;
MTRF1 −0.093; GALNT2 −0.079; DPF2 −0.050; MUS81 −0.054; MYH1 −0.133.

Additional corrections/notes:

- The brief's Table 1 values in Section 3 are rounded; the paper values above win.
- The MYH1 candidate is a delins (`g.30114250_30114327delins…`), not the SNP at 30114249.
  The available `gen_1.7`/`gen_2.9`/`gen_3.9`/`gen_4.4`/`gen_5.4` variant is the left-flank
  proxy `30114249 C/CCCACCTGGG`; the paper explicitly says imputed genotypes used the left
  flank.
- Table 1's "AAF" is the alternative (mutant) allele frequency and matches the A1 frequencies
  computed here.

## Verification results V1–V11

| Check | Result | Key numbers |
|---|---|---|
| V1 join key | PASS | IID-joined adjusted PLCD4 hom-alt = **−110.93 kg**; row-order raw = **+6.48 kg**. |
| V2 allele coding | PASS | A1 = mutant; max all-population frequency delta = 0.00043; max breed delta = 0.0099 (MYH1 J 0.1189 vs 0.109). |
| V3 ID spaces | PASS | `gen_1` 80,027 vs `phen_1.1` 79,945 (82 extra); `gen_3` 80,027 vs `phen_4.1` 79,659 (368 extra); `other_5` 2,197 trios / 4,276 union; `gen_3` row identity 1.0 but `.fam` order differs. |
| V4 weighting | PASS | Mean absolute error vs Table 1: unweighted 1.357, `rinverse` 1.376, `1/rinverse` 1.348 kg. Choose **unweighted** as primary; sensitivity = `1/rinverse`. |
| V5 `other_7` | PASS | 68,578 rows; F 0–3.1, 179 values >1; slope −79.24 kg per F = **−0.792 kg per 1%** (paper −0.77); col 3 r = 0.838 with `yd_trait`, 0.095 with adjusted. |
| V6 pedigree order | PASS | Best `(0,1,2)`/`(0,2,1)` = 9/15,379 violations (0.0585%); `other_3` 172,572 rows; `other_4` 314,677 rows. |
| V7 national burden | PASS | `gen_10` 5,550 sires; `other_11` 2,799,022 animals; birth years 2010–2019; 100% of SIRE_ID/MGS_ID in `gen_10.KEY`. |
| V8 LOSO files | PASS | Trait mapping 1→body weight, 2→BCS, 3→stature; 13 windows; 501 columns each. |
| V9 summary stats | PASS | `other_1.1` 16,128,394 rows / 16,128,339 distinct (CHR,BP) keys (55 dups); `other_9` 100% coverage; `other_10` 99.553% (72,092 missing); `other_15` 0.284%; chromosomes 1–29 only. |
| V10 RNA/chip IDs | PASS | `gen_6` 350 IDs = `phen_7.1.KEY` = `phen_7.2.ID`; 17,180,153 variants; `gen_7` = `phen_6` (12); `gen_8` = `phen_5` (34); exact candidate positions for 7/9; MTRF1 nearest 78 bp, MUS81 nearest 103 bp; DPF2 exact. |
| V11 source data | PASS | 10 CSVs + 16 images; see table above. |

## Data quirks and mismatches to remember

- `other_2` trait order is swapped relative to the README (verified).
- `gen_3` row order matches `gen_1`, but `.fam` IID order does not; use row mapping.
- `other_5` has 2,197 trios, not 4,276 rows; 4,276 is the union of offspring/sire/dam.
- `gen_9.1` README says DPF2/MUS81 chip; it contains exactly those two variants.
- MTRF1 and MUS81 exact mutations are **not** in `gen_6.1`; nearest variants are 78 bp and
  103 bp away. The phrase “4 RNA candidate variants” in the brief is ambiguous; the paper's RNA
  analyses focus on FGD4, GALNT2, and MUS81.
- MYH1 genotypes are the left-flank proxy SNP, not the full delins.
- `ExtFig2c` IDs are renumbered and cannot be row-joined to `phen_1.1`.
- `other_10` covers 99.55% of `other_1.1` keys; `other_15` covers 0.284%.
- `other_7` F is a fraction; 179 values exceed 1 and should be flagged, not silently clipped.
- Numeric ID ranges overlap across populations; never join across populations.
- `other_1.1` sign convention: `B_A` appears to be the A2/reference-allele additive effect;
  hom-alt (A1/A1) ≈ `−2·B_A`; het ≈ `−B_A + B_D`. Reproduces Table 1; confirm in Phase 1.
- `gen_10` genotype coding is 0/1/2 = reference hom / het / mutant hom; `DPF2` stands in for
  the DPF2/MUS81 pair.
- DNA analyses use UMD3.1; RNA-seq (`gen_6`, `phen_7`) uses ARS-UCD1.2.
