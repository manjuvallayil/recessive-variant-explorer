// Zod schemas for the validated exports under web/public/data/.
// Used (1) at build time by scripts/validate-data.ts and (2) at runtime when a section lazy-loads
// its JSON. Objects are "loose" so extra export fields do not break validation; every field the
// site reads is typed here.

import { z } from "zod";

const num = z.number();
const str = z.string();
const bool = z.boolean();

const meta = z.looseObject({
  schema_version: str,
  code_version: str.optional(),
  generated_at: str.optional(),
});

export const FrequenciesSchema = z.looseObject({
  ...meta.shape,
  loci: z.record(
    str,
    z.looseObject({
      gene: str,
      breeds: z.record(
        str,
        z.looseObject({
          n: num,
          alt_freq: num,
          carrier_freq: num,
          hwe_p: num.optional(),
        }),
      ),
    }),
  ),
});

export const DepletionSchema = z.looseObject({
  ...meta.shape,
  rows: z.array(
    z.looseObject({
      locus: str,
      gene: str,
      breed: str,
      n: num,
      alt_freq: num,
      observed_hom_alt: num,
      expected_hom_alt: num,
      observed_over_expected: num,
      lrt_p_one_sided: num,
    }),
  ),
});

const burdenMutation = z.looseObject({
  freq_HF: num,
  freq_J: num,
  freq_all: num,
  expected_total: num,
  rate_per_animal: num,
  scaled_per_year: num,
});

export const NationalBurdenSchema = z.looseObject({
  ...meta.shape,
  n_animals: num,
  mutations: z.record(str, burdenMutation),
  any_mutation: z.looseObject({ rate_per_animal: num, scaled_per_year: num }),
  any_bw_significant_mutation: z.looseObject({
    rate_per_animal: num,
    scaled_per_year: num,
    mutations_included: z.array(str),
  }),
  paper_homozygous_rate: num,
  paper_scaled_homozygotes: num,
  females_born_per_year: num,
});

export const GenotypeClassSchema = z.looseObject({
  ...meta.shape,
  rows: z.array(
    z.looseObject({
      trait: str,
      locus: str,
      estimator: str,
      n_total: num,
      class_stats: z.record(
        str,
        z.looseObject({
          n: num,
          mean: num,
          sd: num,
          quantiles: z.looseObject({
            p05: num,
            p25: num,
            p50: num,
            p75: num,
            p95: num,
          }),
        }),
      ),
    }),
  ),
});

const manhattanTrait = z.looseObject({
  meta: z.looseObject({ n_thinned_points: num, max_points_cap: num }),
  columns: z.looseObject({
    chr: z.array(num),
    bp_umd: z.array(num),
    chr_ars: z.array(num),
    bp_ars: z.array(num),
    neg_log10_p_x100: z.array(num),
  }),
});

export const ManhattanSchema = z.looseObject({
  ...meta.shape,
  model: str,
  gws_threshold_neg_log10_p: num,
  traits: z.looseObject({ bodyweight: manhattanTrait, stature: manhattanTrait }),
});

export const ChipSchema = z.looseObject({
  ...meta.shape,
  summary: z.looseObject({
    primary_8_loci_significant_at_sequence: num,
    primary_8_loci_significant_on_50k_chip: z.array(str),
    major_recessive_syndromes_significant_on_50k_chip: z.array(str),
    note: str,
  }),
  columns: z.looseObject({
    chr: z.array(num),
    bp_umd: z.array(num),
    neg_log10_pd_x100: z.array(num),
  }),
});

const qtlWindow = z.looseObject({
  window_best_dom_d_neg_log10_p: num,
  window_best_dom_d_bp: num,
  window_best_add_neg_log10_p: num,
  window_best_add_bp: num,
  candidate_add_neg_log10_p: num,
  candidate_dom_d_neg_log10_p: num,
  candidate_detected_by_additive_gws: bool.optional(),
});

export const QtlContrastSchema = z.looseObject({
  ...meta.shape,
  gws_threshold_neg_log10_p: num,
  rows: z.array(
    z.looseObject({
      qtl: str,
      gene: z.string().nullable(),
      chr: num,
      candidate_bp: num,
      is_primary_8: bool,
      category: str,
      bodyweight: qtlWindow,
    }),
  ),
});

const regionalGene = z.looseObject({
  symbol: str.nullable(),
  start_umd: num,
  end_umd: num,
  strand: num,
  is_candidate_gene: bool,
});

export const RegionalSchema = z.looseObject({
  ...meta.shape,
  locus: str,
  slug: str,
  gene: z.string().nullable(),
  chr: num,
  ld_available: bool,
  window_umd: z.tuple([num, num]),
  candidate_variant: z.looseObject({
    bp_umd: num,
    snp: str,
    neg_log10_pd: num,
  }),
  genes: z.array(regionalGene),
  columns: z.looseObject({
    bp_umd: z.array(num),
    r2_x1000: z.array(num).nullable(),
    neg_log10_pd_x100: z.array(num),
  }),
});

export const RegionalIndexSchema = z.looseObject({
  loci: z.array(
    z.looseObject({
      slug: str,
      locus: str,
      candidate_bp_umd: num,
      file: str,
      n_variants: num,
    }),
  ),
});

export const PlmSchema = z.looseObject({
  ...meta.shape,
  model: str,
  device: str,
  genes: z.array(
    z.looseObject({
      gene: str,
      accession: str,
      length: num,
      mutation: str,
      position: num,
      wt: str,
      mut: str,
      llr: num,
      damage_percentile: num,
      observed_bodyweight_hom_effect: z.record(
        str,
        z.looseObject({ hom_alt_kg: num, se_kg: num }),
      ),
      matrix_columns: z.array(str),
      matrix: z.array(z.array(z.union([num, str]))),
    }),
  ),
  dropped: z.record(str, str),
});

export const RnaSchema = z.looseObject({
  rna: z.looseObject({
    splice: z.record(str, z.looseObject({ n_animals: num, mean_value: num })),
    expression: z.record(str, z.looseObject({ n_animals: num, mean_value: num })),
  }),
});

export const FarmSchema = z.looseObject({
  farm_profiles: z.record(
    str,
    z.record(str, z.looseObject({ n_animals: num })),
  ),
});

export const SireSchema = z.looseObject({
  per_mutation: z.record(
    str,
    z.looseObject({
      n_carrier_sires: num,
      total_progeny: num,
      carrier_share_of_all_progeny: num,
      carrier_share_by_birth_year: z.record(str, num),
    }),
  ),
});

export const ImputationSchema = z.looseObject({
  maf_bins: z.array(
    z.looseObject({
      maf_bin: str,
      maf_min: num,
      maf_max: num,
      n_variants: num,
      mean_ar2: num,
      quantiles: z.looseObject({
        p05: num,
        p10: num,
        p25: num,
        p50: num,
        p75: num,
        p90: num,
        p95: num,
      }),
      fraction_ar2_ge_0_95: num,
    }),
  ),
  candidate_variants: z.array(
    z.looseObject({
      locus: str,
      gene: z.string().nullable(),
      ar2: num,
      maf: num,
      exceeds_0_95: bool,
    }),
  ),
});

export const PowerSchema = z.looseObject({
  power: z.looseObject({
    variance_explained: num,
    sample_size: num,
    power: num,
    critical_chi2_df1: num,
  }),
  paper_claim: z.looseObject({
    maf: num,
    variance_explained: num,
    heritability: num,
    n_animals: num,
    power_recessive_class: num,
  }),
  validation: z.looseObject({
    passes: bool,
    max_abs_diff_vs_source: num,
    tolerance_abs: num,
    compared_models: z.array(str),
    source_ext_fig1: z.looseObject({
      cow_recessive_power: num,
      cow_genotype_class_power: num,
    }),
  }),
});

const syndromeScore = z.looseObject({
  auc: num,
  auc_ci: z.tuple([num, num]),
  auprc: num,
  auprc_ci: z.tuple([num, num]),
});

export const SyndromeSchema = z.looseObject({
  syndrome_signal: z.looseObject({
    n_records: num,
    n_positive_labels: num,
    base_rate: num,
    label_counts: z.record(str, num),
    scores: z.record(str, syndromeScore),
  }),
  per_group: z.record(
    str,
    z.looseObject({
      n_homozygous_with_bodyweight: num,
      n_homozygous_complete_case: num,
      n_expected_hwe_homozygous: num,
      observed_over_expected: num.nullable(),
      mean_bw_homozygous: num.nullable(),
      mean_bw_reference: num,
      auc_negated_bw_single_locus: num,
      auc_negated_bw_single_locus_ci: z.tuple([num, num]),
      auc_expected_from_effect_size: num,
    }),
  ),
  per_group_meta: z.looseObject({
    pooled_any_homozygous_with_bodyweight: num,
    pooled_any_homozygous_complete_case: num,
  }),
});

export const DetectionSchema = z.looseObject({
  detection: z.looseObject({
    n_total: num,
    n_detected: num,
    n_detected_additive: num,
    n_detected_dominance: num,
  }),
  qtl_truth: z.looseObject({
    threshold_p: num,
    threshold_definition: str,
    min_p_additive_truth: num,
    min_p_dominance_truth: num,
    n_detected_additive_truth: num,
    n_detected_dominance_truth: num,
    by_maf_bin: z.array(
      z.looseObject({
        bin: str,
        n_qtl: num,
        n_detected_additive: num,
        n_detected_dominance: num,
      }),
    ),
    by_effect_bin: z.array(
      z.looseObject({
        bin: str,
        n_qtl: num,
        n_detected_additive: num,
        n_detected_dominance: num,
      }),
    ),
  }),
});

export const InbreedingSchema = z.looseObject({
  inbreeding: z.looseObject({
    n_animals: num,
    slope_kg_per_1pct_F: num,
    slope_kg_per_1pct_F_ci: z.tuple([num, num]),
    F_n_above_1: num,
  }),
});

export const MatingInputsSchema = z.looseObject({
  loci: z.array(
    z.looseObject({
      id: str,
      label: str,
      in_body_weight_set: bool,
      q_HF: num,
      q_J: num,
      depletion_ratio_primary: num,
      depletion_primary_breed: str,
      national_mean_q: num,
      effect: z.looseObject({
        hom_kg: num,
        het_kg: num,
      }),
    }),
  ),
  national: z.looseObject({
    females_born_per_year: num,
    paper_homozygous_rate: num,
    paper_scaled_homozygotes: num,
    export_bw_significant_scaled: num,
    export_per_mutation_scaled: z.record(str, num),
  }),
});

export const ManifestSchema = z.looseObject({
  files: z.array(
    z.looseObject({
      file: str,
      bytes: num,
      gzip_bytes: num.optional(),
      source_files: z.array(str).optional(),
    }),
  ),
});

// Registry used by scripts/validate-data.ts (files relative to web/public/data).
export const DATA_FILES: { file: string; schema: z.ZodType }[] = [
  { file: "frequencies.json", schema: FrequenciesSchema },
  { file: "depletion.json", schema: DepletionSchema },
  { file: "national_burden.json", schema: NationalBurdenSchema },
  { file: "genotype_class_stats.json", schema: GenotypeClassSchema },
  { file: "manhattan_dom.json", schema: ManhattanSchema },
  { file: "manhattan_add.json", schema: ManhattanSchema },
  { file: "chip_vs_sequence.json", schema: ChipSchema },
  { file: "qtl_contrast.json", schema: QtlContrastSchema },
  { file: "regional/index.json", schema: RegionalIndexSchema },
  { file: "plm_scores.json", schema: PlmSchema },
  { file: "rna_phase3.json", schema: RnaSchema },
  { file: "farm_profiles.json", schema: FarmSchema },
  { file: "sire_amplification.json", schema: SireSchema },
  { file: "imputation_by_maf.json", schema: ImputationSchema },
  { file: "power_reference.json", schema: PowerSchema },
  { file: "syndrome_signal.json", schema: SyndromeSchema },
  { file: "detection_summary.json", schema: DetectionSchema },
  { file: "inbreeding.json", schema: InbreedingSchema },
  { file: "mating_inputs.json", schema: MatingInputsSchema },
  { file: "manifest.json", schema: ManifestSchema },
];

// Regional files are validated individually (one per locus slug from regional/index.json).
export const REGIONAL_SCHEMA = RegionalSchema;

export type Frequencies = z.infer<typeof FrequenciesSchema>;
export type Depletion = z.infer<typeof DepletionSchema>;
export type NationalBurden = z.infer<typeof NationalBurdenSchema>;
export type GenotypeClass = z.infer<typeof GenotypeClassSchema>;
export type Manhattan = z.infer<typeof ManhattanSchema>;
export type Chip = z.infer<typeof ChipSchema>;
export type QtlContrast = z.infer<typeof QtlContrastSchema>;
export type Regional = z.infer<typeof RegionalSchema>;
export type RegionalIndex = z.infer<typeof RegionalIndexSchema>;
export type Plm = z.infer<typeof PlmSchema>;
export type Rna = z.infer<typeof RnaSchema>;
export type Farm = z.infer<typeof FarmSchema>;
export type Sire = z.infer<typeof SireSchema>;
export type Imputation = z.infer<typeof ImputationSchema>;
export type Power = z.infer<typeof PowerSchema>;
export type Syndrome = z.infer<typeof SyndromeSchema>;
export type Detection = z.infer<typeof DetectionSchema>;
export type Inbreeding = z.infer<typeof InbreedingSchema>;
export type MatingInputs = z.infer<typeof MatingInputsSchema>;
