# Phase 1 Reproduction Report

Status: **Phase 1 complete.** Generated from `web/public/data/*.json` produced by `pipeline/phase1.py`.

## 1. Estimators and Weighting

Three estimators are reported for every Table 1 locus and trait:

- **Estimator A (unweighted `yd_trait`)**: plain OLS on `yd_trait` in `phen_1.{1,2,3}` without leave-one-segment-out (LOSO) population-structure adjustment.
- **Estimator A_rinverse (`weights = rinverse` sensitivity)**: WLS on `yd_trait` weighted by `rinverse` (Gate A Decision 1 sensitivity run).
- **Estimator B (LOSO-adjusted, paper-like)**: unweighted OLS on `yd_trait - mean(other_2.<trait>.<window>)`, centered, using the non-overlapping 10 Mb LOSO window containing each locus (`other_2.1` = body weight, `other_2.3` = stature, `other_2.2` = BCS).
- **Published `other_1.*` lead-variant row**: from the dominance GWAS summary statistics (`other_1.{1,2,3}.txt`), where `B_A` is the reference/A2 additive effect so `hom_alt = -2 * B_A` and `het = -B_A + B_D`.

**Unit note for Stature (`phen_1.2`, `other_2.3`, `other_1.2`):** Raw stature records and GWAS summary statistics are in New Zealand TOP (Traits Other than Production) linear score units (1–9 scale, where **1 TOP score unit = 5.0 cm**), whereas Table 1 in the paper reports stature in **cm**. Scaling stature by `5.0 cm / score` reproduces Table 1 across all nine loci (Section 2.2). Both cm and raw TOP-score values are preserved in `web/public/data/effects.json`.

### Estimator comparison summary

| Metric | Estimator A (unweighted `yd_trait`) | Estimator A_rinverse (`weights=rinverse`) | Estimator B (LOSO-adjusted) | Closest |
|---|---:|---:|---:|---|
| All 27 locus-trait pairs: Hom-Alt MAE | 12.963 | 12.944 | **0.474** | **B** |
| All 27 locus-trait pairs: Het MAE | 9.566 | 9.588 | **0.219** | **B** |
| Body weight (9 loci, kg): Hom-Alt MAE | 36.370 | 36.313 | **1.357** | **B** |
| Body weight (9 loci, kg): Het MAE | 26.875 | 26.942 | **0.629** | **B** |
| Large BW effects (\|Paper\| ≥ 20 kg): Hom-Alt MAE | 24.962 | 24.751 | **1.897** | **B** |
| Max \|Δ Hom-Alt\| across all 27 pairs | 61.912 | 62.183 | **3.355** | **B** |

**Why Estimator A diverges while Estimator B matches Table 1:** `yd_trait` is not pre-adjusted for breed/genomic structure (`M_α_` in Equation 1 of the paper). Because Holstein-Friesian (HF) cows are ~60 kg heavier and taller than Jersey (J) cows, HF-enriched alleles (`PLCD4`, `FGD4`, `GALNT2`, `DPF2`, `MUS81`, `Chr2:22`, `ARSI`) show a raw heterozygote shift of `+22` to `+30 kg` in Estimator A, while Jersey-enriched alleles (`MTRF1`, `MYH1`) show `-32` to `-34 kg`. Subtracting the LOSO polygenic background in Estimator B removes this breed confounding and matches both `other_1.*` and Table 1.

---

## 2. Table 1 Reproduction by Trait

### 2.1 Body Weight (`phen_1.1`, kg, N = 79,945)

| Locus | Class n (`0/1/2`) | Paper Het / Hom-Alt (kg) | Est. B Het / Hom-Alt (kg) | Δ Est. B vs Paper (Het / Hom) | `other_1.1` Derived Het / Hom-Alt | Est. A (unweighted) Het / Hom-Alt | Est. A (`rinverse`) Het / Hom-Alt |
|---|---|---:|---:|---:|---:|---:|---:|
| Chr2:22 | 68,624 / 10,869 / 452 | −0.007 / −11.800 | −0.476 / −12.601 | −0.469 / −0.801 | −0.412 / −12.534 | +24.955 / +27.573 | +25.139 / +27.647 |
| PLCD4 | 76,105 / 3,819 / 21 | −0.354 / −108.965 | −1.086 / −110.928 | −0.732 / −1.963 | −1.205 / −111.350 | +24.144 / −80.619 | +24.019 / −80.844 |
| FGD4 | 76,131 / 3,786 / 28 | +0.834 / −48.079 | +0.902 / −48.298 | +0.068 / −0.219 | +1.281 / −46.942 | +29.379 / −24.764 | +29.241 / −24.995 |
| ARSI | 68,930 / 10,559 / 456 | +0.564 / −11.981 | +0.196 / −12.886 | −0.368 / −0.905 | +0.078 / −12.678 | +30.051 / +31.695 | +29.847 / +31.728 |
| MTRF1 | 72,924 / 6,865 / 156 | +1.078 / −19.464 | +1.122 / −19.270 | +0.044 / +0.194 | +1.179 / −19.414 | −31.619 / −77.025 | −31.863 / −77.186 |
| GALNT2 | 74,665 / 5,251 / 29 | −0.837 / −112.008 | −0.836 / −112.644 | +0.001 / −0.636 | −0.808 / −112.054 | +22.111 / −89.057 | +22.249 / −89.013 |
| DPF2 | 73,787 / 6,052 / 106 | +4.059 / −53.092 | +2.325 / −56.447 | −1.734 / −3.355 | +2.286 / −56.552 | +26.151 / −28.649 | +26.320 / −28.916 |
| MUS81 | 73,725 / 6,111 / 109 | +4.126 / −51.911 | +2.385 / −55.221 | −1.741 / −3.310 | +2.360 / −55.474 | +26.029 / −26.157 | +26.193 / −26.534 |
| MYH1 | 73,374 / 6,393 / 178 | +0.798 / −4.643 | +0.289 / −5.475 | −0.509 / −0.832 | +0.380 / −5.424 | −33.944 / −66.555 | −34.117 / −66.826 |

**LOSO window sensitivity (body weight):** Four loci have a second overlapping 10 Mb window in `other_2.1.*`:
- `Chr2:22` window 2 (`Chr2_15_25Mb`): Het = `−0.050 kg`, Hom-Alt = `−11.796 kg` (Δ vs Table 1 = `−0.043 / +0.004 kg`).
- `PLCD4` window 4 (`Chr2_105_115Mb`): Het = `−0.950 kg`, Hom-Alt = `−110.591 kg` (Δ vs Table 1 = `−0.596 / −1.626 kg`).
- `FGD4` window 6 (`Chr5_75_85Mb`): Het = `+1.351 kg`, Hom-Alt = `−47.412 kg` (matches `ExtFig2c` and `other_1.1`).
- `ARSI` window 8 (`Chr7_55_65Mb`): Het = `+0.653 kg`, Hom-Alt = `−12.140 kg` (Δ vs Table 1 = `+0.089 / −0.159 kg`).

### 2.2 Stature (`phen_1.2`, cm = 5 × TOP score, N = 75,041)

| Locus | Class n (`0/1/2`) | Paper Het / Hom-Alt (cm) | Est. B Het / Hom-Alt (cm) | Δ Est. B vs Paper (Het / Hom) | `other_1.2` Derived Het / Hom-Alt (cm) | Est. A (unweighted, cm) | Est. A (`rinverse`, cm) |
|---|---|---:|---:|---:|---:|---:|---:|
| Chr2:22 | 64,596 / 10,032 / 413 | +0.025 / −0.428 | −0.011 / −0.518 | −0.036 / −0.090 | −0.015 / −0.530 | +1.598 / +1.944 | +1.598 / +1.944 |
| PLCD4 | 71,484 / 3,537 / 20 | −0.025 / −5.421 | −0.045 / −5.433 | −0.020 / −0.012 | −0.045 / −5.460 | +1.701 / −3.016 | +1.701 / −3.016 |
| FGD4 | 71,433 / 3,581 / 27 | +0.058 / −3.041 | +0.058 / −3.071 | 0.000 / −0.030 | +0.060 / −3.080 | +2.016 / −1.101 | +2.016 / −1.101 |
| ARSI | 64,688 / 9,925 / 428 | +0.055 / −1.152 | +0.045 / −1.177 | −0.010 / −0.025 | +0.045 / −1.170 | +2.096 / +1.839 | +2.096 / +1.839 |
| MTRF1 | 68,445 / 6,442 / 154 | +0.025 / −1.208 | +0.048 / −1.169 | +0.023 / +0.039 | +0.045 / −1.180 | −2.196 / −5.070 | −2.196 / −5.070 |
| GALNT2 | 70,156 / 4,860 / 25 | −0.134 / −8.568 | −0.134 / −8.606 | 0.000 / −0.038 | −0.135 / −8.630 | +1.307 / −7.382 | +1.307 / −7.382 |
| DPF2 | 69,283 / 5,662 / 96 | +0.132 / −2.565 | +0.072 / −2.695 | −0.060 / −0.130 | +0.075 / −2.700 | +1.568 / −0.888 | +1.568 / −0.888 |
| MUS81 | 69,222 / 5,721 / 98 | +0.136 / −2.487 | +0.080 / −2.622 | −0.056 / −0.135 | +0.080 / −2.620 | +1.563 / −0.748 | +1.563 / −0.748 |
| MYH1 | 68,840 / 6,033 / 168 | +0.086 / +0.285 | +0.088 / +0.307 | +0.002 / +0.022 | +0.085 / +0.300 | −2.237 / −3.826 | −2.237 / −3.826 |

### 2.3 Body Condition Score (`phen_1.3`, 1–10 scale, N = 75,617)

| Locus | Class n (`0/1/2`) | Paper Het / Hom-Alt | Est. B Het / Hom-Alt | Δ Est. B vs Paper (Het / Hom) | `other_1.3` Derived Het / Hom-Alt | Est. A (unweighted) | Est. A (`rinverse`) |
|---|---|---:|---:|---:|---:|---:|---:|
| Chr2:22 | 65,090 / 10,108 / 419 | +0.003 / −0.044 | +0.002 / −0.045 | −0.001 / −0.001 | +0.003 / −0.046 | +0.032 / +0.010 | +0.032 / +0.010 |
| PLCD4 | 72,041 / 3,556 / 20 | +0.005 / −0.284 | +0.009 / −0.275 | +0.004 / +0.009 | +0.009 / −0.276 | +0.031 / −0.252 | +0.031 / −0.252 |
| FGD4 | 71,983 / 3,606 / 28 | +0.012 / −0.238 | +0.012 / −0.234 | 0.000 / +0.004 | +0.014 / −0.232 | +0.034 / −0.237 | +0.034 / −0.237 |
| ARSI | 65,192 / 9,994 / 431 | −0.001 / −0.028 | 0.000 / −0.025 | +0.001 / +0.003 | 0.000 / −0.026 | +0.009 / −0.023 | +0.009 / −0.023 |
| MTRF1 | 68,961 / 6,502 / 154 | +0.009 / −0.093 | +0.010 / −0.093 | +0.001 / 0.000 | +0.010 / −0.094 | −0.011 / −0.126 | −0.011 / −0.126 |
| GALNT2 | 70,703 / 4,889 / 25 | +0.012 / −0.079 | +0.012 / −0.078 | 0.000 / +0.001 | +0.012 / −0.076 | +0.057 / +0.031 | +0.057 / +0.031 |
| DPF2 | 69,812 / 5,709 / 96 | +0.020 / −0.050 | +0.012 / −0.067 | −0.008 / −0.017 | +0.011 / −0.068 | +0.056 / −0.013 | +0.056 / −0.013 |
| MUS81 | 69,751 / 5,768 / 98 | +0.020 / −0.054 | +0.012 / −0.069 | −0.008 / −0.015 | +0.012 / −0.068 | +0.055 / −0.013 | +0.055 / −0.013 |
| MYH1 | 69,372 / 6,075 / 170 | −0.010 / −0.133 | −0.016 / −0.148 | −0.006 / −0.015 | −0.016 / −0.146 | −0.041 / −0.216 | −0.041 / −0.216 |

---

## 3. Aggregate Comparison with `source_data.ExtFig2c.csv`

Per Gate A Decision 4, `ExtFig2c.csv` is compared exclusively at the aggregate level (genotype-class counts, means, quantiles, and derived effects; never row-joined). Every single locus in `ExtFig2c.csv` matches our Estimator B on the corresponding LOSO window to `< 0.001 kg`:

| Locus | `ExtFig2c` Class n (`0/1/2`) | `ExtFig2c` Het / Hom-Alt (kg) | Est. B (matching window) Het / Hom-Alt (kg) | Matching Window | \|Δ\| (kg) |
|---|---|---:|---:|---:|---:|
| Chr2:22 | 68,624 / 10,869 / 452 | −0.476 / −12.601 | −0.476 / −12.601 | Window 1 (`Chr2_20_30Mb`) | < 0.001 |
| PLCD4 | 76,105 / 3,819 / 21 | −1.086 / −110.928 | −1.086 / −110.928 | Window 3 (`Chr2_100_110Mb`) | < 0.001 |
| FGD4 | 76,131 / 3,786 / 28 | +1.351 / −47.412 | +1.351 / −47.412 | Window 6 (`Chr5_75_85Mb`) | < 0.001 |
| ARSI | 68,930 / 10,559 / 456 | +0.196 / −12.886 | +0.196 / −12.886 | Window 7 (`Chr7_60_70Mb`) | < 0.001 |
| MTRF1 | 72,924 / 6,865 / 156 | +1.122 / −19.270 | +1.122 / −19.270 | Window 9 (`Chr12_5_15Mb`) | < 0.001 |
| GALNT2 | 74,665 / 5,251 / 29 | −0.836 / −112.644 | −0.836 / −112.644 | Window 12 (`Chr28_0_10Mb`) | < 0.001 |
| DPF2 | 73,787 / 6,052 / 106 | +2.325 / −56.447 | +2.325 / −56.447 | Window 13 (`Chr29_40_50Mb`) | < 0.001 |
| PLAG1 (`chr14:25015640`) | 25,524 / 33,586 / 20,835 | +21.285 / +37.880 | +21.285 / +37.880 | Window 10 (`Chr14_20_30Mb`) | < 0.001 |

---

## 4. Validation Cohort (`gen_2` + `phen_2`, N = 31,580) & Pleiotropy (`gen_3` × 19 Traits)

- **Validation (`validation_milk.json`)**: Across 9 loci × 3 first-lactation yield traits (milk, protein, fat = 27 locus-trait pairs), **26/27 (96.3%)** homozygote effects and **25/27 (92.6%)** heterozygote effects have the same sign in the independent validation population (`phen_2`) as in the discovery lactation cohort (`phen_4.1–4.3`). Minimum genotype-class size is `n = 6` (`PLCD4` hom-alt).
- **Pleiotropy (`pleiotropy.json`)**: 12 variants (9 Table 1 loci + 3 BTA14 comparators: `DGAT1` `chr14:1802265`, `PLAG1` causal `chr14:25015640`, and `PLAG1` GWAS lead `chr14:26064929`) × 19 traits = 228 cells (`phen_1` joined via `gen_1.fam`, `phen_4` joined via `gen_3.fam`). All classes have `n ≥ 20`.

---

## 5. Breed Frequencies (`gen_4` + `other_6`, N = 49,115 Purebreds)

Purebred counts: **29,771 HF** (`breed_hf = 16`) and **19,344 Jersey** (`breed_j = 16`).

| Locus | HF Alt Freq (mine) | Paper HF AAF | Δ HF | HF Carrier % | J Alt Freq (mine) | Paper J AAF | Δ J | J Carrier % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| PLCD4 | 0.0426 | 0.039 | +0.0036 | 8.42% | 0.0006 | 0.001 | −0.0004 | 0.11% |
| FGD4 | 0.0405 | 0.039 | +0.0015 | 8.02% | 0.0006 | 0.001 | −0.0004 | 0.12% |
| MTRF1 | 0.0018 | 0.002 | −0.0002 | 0.36% | 0.1173 | 0.114 | +0.0033 | 22.47% |
| MYH1 | 0.0040 | 0.006 | −0.0020 | 0.80% | 0.1189 | 0.109 | +0.0099 | 22.52% |
| GALNT2 | 0.0556 | 0.055 | +0.0006 | 10.92% | 0.0009 | 0.001 | −0.0001 | 0.18% |
| DPF2 | 0.0655 | 0.066 | −0.0005 | 12.72% | 0.0004 | 0.001 | −0.0006 | 0.09% |
| MUS81 | 0.0658 | 0.066 | −0.0002 | 12.78% | 0.0006 | 0.001 | −0.0004 | 0.11% |

**Cumulative carrier burden vs paper's ">40% of purebred HF and Jerseys":**
- Sum of carrier frequencies across the 6 distinct candidate loci (`PLCD4, FGD4, MTRF1, MYH1, GALNT2, DPF2`): **41.24% (HF)** and **45.50% (J)**.
- Empirical fraction carrying ≥1 mutant allele in `gen_4` (7 files, without `Chr2:22` or `ARSI`): **34.61% (HF)** and **39.64% (J)**.
- Including `Chr2:22` and `ARSI` at their Table 1 breed AAFs under independence: **60.04% (HF)** and **40.09% (J)**.

---

## 6. Homozygous Depletion (`gen_4`) and Carrier Crosses (`gen_5` + `other_5`)

### 6.1 Homozygous Depletion in Purebreds (`depletion.json`)

| Locus | Breed | N | Alt Freq | Observed Hom-Alt | Expected HWE Hom-Alt | Obs / Exp | One-sided LRT P | One-sided Binomial P |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| PLCD4 | HF | 29,771 | 0.0426 | 26 | 53.92 | **0.482** | 1.14 × 10⁻⁵ | 1.86 × 10⁻⁵ |
| FGD4 | HF | 29,771 | 0.0405 | 22 | 48.73 | **0.451** | 8.51 × 10⁻⁶ | 1.45 × 10⁻⁵ |
| MTRF1 | J | 19,344 | 0.1173 | 192 | 266.27 | **0.721** | 7.11 × 10⁻⁷ | 8.89 × 10⁻⁷ |
| MYH1 | J | 19,344 | 0.1189 | 245 | 273.59 | **0.896** | 0.0382 | 0.0418 |
| GALNT2 | HF | 29,771 | 0.0556 | 61 | 92.11 | **0.662** | 2.69 × 10⁻⁴ | 3.59 × 10⁻⁴ |
| DPF2 | HF | 29,771 | 0.0655 | 112 | 127.59 | **0.878** | 0.0788 | 0.0882 |
| MUS81 | HF | 29,771 | 0.0658 | 114 | 128.91 | **0.884** | 0.0897 | 0.1000 |

### 6.2 Carrier × Carrier Trios (`carrier_cross.json`, 2,197 trios in `other_5`)

Selected trios where `sire == 1` and `dam == 1` at the locus (`offspring` column 0, `sire` column 1, `dam` column 2):

| Locus | Het × Het Trios | Offspring `0 / 1 / 2` | Expected `1:2:1` | χ² (2 df) | P-value |
|---|---:|---|---:|---:|---:|
| PLCD4 | 140 | 34 / 80 / 26 | 35.0 / 70.0 / 35.0 | 3.771 | 0.1517 |
| FGD4 | 160 | 37 / 98 / 25 | 40.0 / 80.0 / 40.0 | 9.900 | **0.0071** |
| MTRF1 | 711 | 201 / 356 / 154 | 177.8 / 355.5 / 177.8 | 6.215 | **0.0447** |
| MYH1 | 662 | 173 / 345 / 144 | 165.5 / 331.0 / 165.5 | 3.725 | 0.1553 |
| GALNT2 | 224 | 68 / 120 / 36 | 56.0 / 112.0 / 56.0 | 10.286 | **0.0058** |
| DPF2 | 374 | 88 / 195 / 91 | 93.5 / 187.0 / 93.5 | 0.733 | 0.6933 |
| MUS81 | 376 | 89 / 197 / 90 | 94.0 / 188.0 / 94.0 | 0.867 | 0.6482 |

---

## 7. National Burden (`gen_10.1` + `other_11`, N = 2,799,022) & Sire Amplification

**Mendelian derivation:** For each animal with genotyped sire (`g_sire ∈ {0,1,2}`), genotyped maternal grandsire (`g_MGS ∈ {0,1,2}`), and ungenotyped maternal granddam with breed composition `(mgd_HF_16, mgd_J_16)`:
- `P(sire transmits mutant) = g_sire / 2`
- `P(dam transmits mutant) = 0.5 * (g_MGS / 2) + 0.5 * q_MGD`, where `q_MGD = (mgd_HF_16/16)*f_HF + (mgd_J_16/16)*f_J + remainder*f_all`
- `P(calf homozygous) = P(sire transmits) * P(dam transmits)`

| Mutation / Set | Expected Homozygotes (N = 2,799,022) | Rate per Animal | Scaled to 940,100 Females/yr | Paper Target |
|---|---:|---:|---:|---|
| PLCD4 | 3,079.3 | 0.110% | 1,034 / yr | — |
| FGD4 | 5,373.9 | 0.192% | 1,805 / yr | — |
| MTRF1 | 10,906.7 | 0.390% | 3,663 / yr | — |
| MYH1 | 8,568.7 | 0.306% | 2,878 / yr | — |
| GALNT2 | 2,508.7 | 0.090% | 843 / yr | — |
| DPF2/MUS81 | 6,653.8 | 0.238% | 2,235 / yr | — |
| **5 BW-significant loci (excl. MYH1)** | **28,415.7** | **1.015%** | **9,543 / yr** | **~1% (~9,700 / yr)** |
| **All 6 `gen_10` columns (incl. MYH1)** | **36,820.7** | **1.315%** | **12,367 / yr** | — |

**Explanation of the ~9,700 vs 9,543 / 12,367 gap:** In the paper, the `~1%` / `~9,700` headline appears in the section on the five genome-wide significant body-weight recessive syndromes (`PLCD4, FGD4, MTRF1, GALNT2, DPF2/MUS81`), *before* `MYH1` (a suggestive BCS locus, P = 7.6 × 10⁻⁸) is introduced in the following paragraph. Restricting to those 5 loci yields **1.015% (9,543 females/year)** — within 1.6% of the paper's ~9,700 figure. Including `MYH1` adds 0.306% to reach 1.315% (12,367/year).

---

## 8. Prospective Calves (`phen_3.1`, N = 568 across 51 Herds)

OLS model: `weight ~ VARIANT + age + HERD + breed_hf + breed_j` (`CONTROL` n = 526, mean unadjusted weight = 198.42 kg):

| Variant Group | n | Unadjusted Mean (kg) | Adjusted Effect vs Control (kg) | 95% CI (kg) | P-value |
|---|---:|---:|---:|---|---:|
| PLCD4 | 9 | 147.11 | **−55.01 ± 6.36** | [−67.51, −42.50] | 7.09 × 10⁻¹⁷ |
| FGD4 | 22 | 183.11 | **−10.39 ± 4.13** | [−18.50, −2.29] | 0.0121 |
| DPF2/MUS81 | 11 | 199.45 | **+3.33 ± 5.77** | [−8.01, +14.67] | 0.5641 |

---

## 9. Phase 2: Genome-Wide Dominance vs Additive Contrast, Chip Blind Spot, & Imputation

### 9.1 QTL Contrast Table (`qtl_contrast.json`, `chip_vs_sequence.json`, `imputation_by_maf.json`)

`-log10(P)` is computed directly from `Z` via `scipy.special.log_ndtr(-|Z|)` with zero clipping. Genome-wide significance threshold: `P < 5 × 10⁻⁸` (`−log10 P > 7.301`).

| QTL | Candidate (UMD3.1) | Window Best Dom `−log10(P_D)` | Candidate Dom `−log10(P_D)` | Candidate Std-Add `−log10(P_A)` | Window Best Std-Add `−log10(P_A)` | In COJO Std-Add (`ExtFig2a`) | 50k Chip Best `−log10(P_D)` | Candidate `AR2` |
|---|---|---:|---:|---:|---:|---|---:|---:|
| Chr2:22 | `2:22,836,659` | **8.91** | **8.91** | 1.64 | 6.65 | No | **7.64** | 0.984 |
| PLCD4 | `2:107,313,998` | **35.64** | **35.64** | 2.34 | 4.78 | No | 3.10 | 1.000 |
| FGD4 | `5:77,632,752` | **9.96** | **9.96** | 0.25 | 2.65 | No | 3.95 | 1.000 |
| ARSI | `7:63,649,071` | **13.22** | **8.53** | 0.77 | 3.26 | No | 2.91 | 0.988 |
| MTRF1 | `12:11,463,981` | **10.98** | **10.98** | 0.16 | 4.64 | No | 5.42 | 0.983 |
| PLAG1 | `14:26,064,929` | **28.43** | **28.43** | **311.06** | **323.10** | **Yes (7 hits)** | **17.33** | 1.000 (causal) / 0.455 (lead) |
| GALNT2 | `28:1,312,334` | **50.34** | **44.46** | 2.31 | 5.06 | No | 3.61 | 1.000 |
| DPF2/MUS81 | `29:44,213,160` / `44,645,469` | **60.15** | **50.86 / 52.20** | 0.11 / 0.13 | 10.09* | No | **14.76** | 0.979 / 0.986 |
| MYH1 (BCS) | `19:30,114,249` | 5.63 (BCS: 7.12) | 1.29 | 0.01 | 8.02 | No | 3.70 | 0.783 |

*\*Note on `DPF2/MUS81` window-best standard-additive `10.09`:* The candidate variants themselves are completely invisible to the standard-additive model (`−log10 P_A = 0.11` for `DPF2` and `0.13` for `MUS81`). An unlinked common SNP 843 kb away (`chr29:45,056,385`, `MAF = 0.087`) reaches `−log10 P_A = 10.09` in the raw `±1 Mb` window, but is eliminated by GCTA-COJO conditional analysis (`source_data.ExtFig2a.csv`, where the nearest standard-additive COJO QTL on BTA29 are at `37.76 Mb` and `47.10 Mb`). Thus **7 of 8 non-additive QTL are invisible to the standard-additive model (`PLAG1` is the exception)**.

### 9.2 Chip Blind Spot (`chip_vs_sequence.json`)

- **Coverage direction**: **100% (45,865 / 45,865)** of distinct BovineSNP50 positions in `other_15.txt` are present in `other_1.1.txt` (representing **0.284%** of the 16,128,339 distinct sequence positions).
- **Resolution loss**: At 50k-chip density, only 3 of the 8 dominance QTL retain a genome-wide significant SNP (`Chr2:22`, `PLAG1`, and `DPF2/MUS81`), and **only `DPF2/MUS81` (`−log10 P_D = 14.76`) remains genome-wide significant among the major-effect recessive syndromes** — reproducing Extended Data Fig. 9a.

### 9.3 Regional Plots (`regional/*.json`) & Imputation Accuracy (`imputation_by_maf.json`)

- **Regional LD & `ExtFig3` cross-check**: For `MYH1` (`gen_1.7`, BTA19), all 7,862 rows of `source_data.ExtFig3.csv` match our regional variants with `R²` correlation `r = 0.9949` and identical peak BCS dominance significance (`−log10 P = 7.12`, `P = 7.6 × 10⁻⁸`).
- **Imputation `AR2`**: Mean `AR2` increases monotonically from `0.931` (`MAF 0.01–0.02`, median `0.971`) to `0.974` (`MAF 0.40–0.50`, median `0.993`). All 8 body-weight recessive/dominance candidate variants in Table 1 exceed `AR2 > 0.95` (range `0.979–1.000`).
