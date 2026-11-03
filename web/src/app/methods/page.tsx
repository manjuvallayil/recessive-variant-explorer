// Methods and provenance. A server component: it reads the manifest at build time (static export).
import { readFileSync } from "node:fs";
import { join } from "node:path";
import Link from "next/link";
import { ManifestSchema } from "@/lib/schemas";

export const metadata = {
  title: "Methods and data provenance: Recessive Variant Explorer",
};

type Row = { ui: string; value: string; source: string };

const UI_MAP: Row[] = [
  { ui: "Carrier frequencies (hook)", value: "PLCD4 HF 8.4%, DPF2 HF 12.7%, etc.", source: "frequencies.json loci.*.breeds.{HF,J}.carrier_freq" },
  { ui: "Carrier union, HF (hook)", value: "35.3% (chr29:44 once, DPF2 proxy); 43.5% (DPF2 and MUS81 separate)", source: "frequencies.json; independence assumed; computed in UI; paper treats chr29:44 as one locus" },
  { ui: "QTL and mutation sets (hook)", value: "6 novel BW QTL; 5 have candidates; 6 candidate mutations (DPF2 and MUS81 at chr29:44); MYH1 separate BCS-suggestive", source: "paper abstract/Table 1; qtl_contrast.json is_primary_8 and category" },
  { ui: "PLCD4 homozygote body weight (proxy)", value: "-110.9 kg vs reference 0.1", source: "genotype_class_stats.json bodyweight, estimator B, PLCD4 class 2 vs 0 mean" },
  { ui: "GALNT2 homozygote body weight (proxy)", value: "-112.6 kg", source: "genotype_class_stats.json bodyweight, estimator B, GALNT2 class 2 mean" },
  { ui: "Inbreeding slope (proxy)", value: "-0.74 kg per 1% F (CI -0.76 to -0.71)", source: "inbreeding.json inbreeding.slope_kg_per_1pct_F and _ci" },
  { ui: "GWS threshold (Manhattan)", value: "-log10 P 7.30 (P = 5e-8)", source: "qtl_contrast.json gws_threshold_neg_log10_p" },
  { ui: "Additive detection at candidate (Manhattan)", value: "7 of 8 below threshold", source: "qtl_contrast.json candidate_add_neg_log10_p (counted in UI)" },
  { ui: "Additive detection within ±1 Mb (Manhattan)", value: "6 of 8 below threshold", source: "qtl_contrast.json window_best_add_neg_log10_p (counted in UI)" },
  { ui: "PLAG1 additive -log10 P", value: "480 at reintroduced GGPv3 SNP; paper Fig. 1 441.3 at another lead", source: "qtl_contrast.json PLAG1 bodyweight.window_best_add_neg_log10_p and fig1_legend_reconciliation" },
  { ui: "Chip result", value: "Only DPF2/MUS81 GWS among major recessive syndromes", source: "chip_vs_sequence.json summary" },
  { ui: "Candidate -log10 P (zoom)", value: "e.g. PLCD4 35.6", source: "regional/*.json candidate_variant.neg_log10_pd" },
  { ui: "Protein model LLR and percentile", value: "PLCD4 -3.71 / 56; MTRF1 -7.48 / 90; DPF2 -2.11 / 66", source: "plm_scores.json genes.*.llr, damage_percentile" },
  { ui: "Protein model observed effects", value: "from effects.json (estimator B)", source: "plm_scores.json observed_bodyweight_hom_effect" },
  { ui: "RNA means", value: "splice and expression, n = 350 (export; paper 389)", source: "rna_phase3.json rna.*" },
  { ui: "Farm group sizes", value: "6 to 9 animals per group", source: "farm_profiles.json farm_profiles.*.*.n_animals" },
  { ui: "National homozygous births (bw set)", value: "1.02% of animals; about 9,548 per year", source: "national_burden.json any_bw_significant_mutation" },
  { ui: "National homozygous births (per mutation)", value: "PLCD4 1,034; FGD4 1,805; MTRF1 3,663; GALNT2 843; DPF2 2,235; MYH1 2,878", source: "national_burden.json mutations.*.scaled_per_year" },
  { ui: "Paper targets (national)", value: "1%; about 9,700 per year", source: "national_burden.json paper_homozygous_rate, paper_scaled_homozygotes (paper text)" },
  { ui: "Depletion ratio and p-value", value: "e.g. PLCD4 HF 0.48, p = 1.1e-5", source: "depletion.json rows.*.observed_over_expected, lrt_p_one_sided" },
  { ui: "Syndrome AUROC, pooled (signal)", value: "B0 0.678; logistic CV 0.690; Mahalanobis 0.579", source: "syndrome_signal.json syndrome_signal.scores.*" },
  { ui: "Per-locus AUROC and homozygote counts (signal)", value: "PLCD4 0.79 (n 21); GALNT2 0.86 (n 29); FGD4 0.63; DPF2/MUS81 0.63 (n 110)", source: "syndrome_signal.json per_group.*" },
  { ui: "Expected AUROC from effect size (signal)", value: "PLCD4 0.82; GALNT2 0.85", source: "syndrome_signal.json per_group.*.auc_expected_from_effect_size" },
  { ui: "Imputation AR² quantiles", value: "by MAF bin", source: "imputation_by_maf.json maf_bins.*.quantiles" },
  { ui: "Simulated-QTL detection (4c)", value: "dominance 19 of 30; additive 0 of 30 (P < 5e-8)", source: "detection_summary.json qtl_truth.*; min additive P 5.9e-8" },
  { ui: "Power (paper scenario)", value: "0.948 recessive (paper figure 0.947; text ~0.90)", source: "power_reference.json power.power; validation.max_abs_diff_vs_source 0.0021" },
  { ui: "Power curve (explorer)", value: "live noncentral chi-square", source: "web/src/lib/power.ts (tested against scipy in node --test)" },
  { ui: "Simulator inputs", value: "paper frequencies; site-re-estimated effects; depletion; national q", source: "mating_inputs.json (built by pipeline/phase5_inputs.py)" },
  { ui: "Simulator math", value: "transmission, HWE, carrier-avoidance", source: "web/src/lib/mating.ts (unit-tested in mating.test.ts)" },
];

export default function MethodsPage() {
  const manifest = ManifestSchema.parse(
    JSON.parse(readFileSync(join(process.cwd(), "public", "data", "manifest.json"), "utf8")),
  );
  return (
    <div className="mx-auto max-w-5xl space-y-10 px-4 py-10 sm:px-8">
      <p className="text-sm">
        <Link href="/">← Back to the story</Link>
      </p>
      <header>
        <h1 className="serif text-3xl font-semibold">Methods and data provenance</h1>
        <p className="mt-2 text-sm text-[var(--muted)]">
          Every number in the site traces to a derived export under <code>public/data/</code>, produced by the Python pipeline in
          this repository from the authors&apos; public data. Only aggregates are exported; no animal or sire identifiers appear in any
          file (checked at build time).
        </p>
      </header>

      <section aria-labelledby="map">
        <h2 id="map" className="serif text-xl font-semibold">UI number to source</h2>
        <div className="mt-3 overflow-x-auto">
          <table className="w-full text-left text-xs">
            <thead>
              <tr className="border-b border-[var(--rule)]">
                <th className="py-2 pr-3">Shown in the site</th>
                <th className="py-2 pr-3">Value</th>
                <th className="py-2">Source</th>
              </tr>
            </thead>
            <tbody>
              {UI_MAP.map((r) => (
                <tr key={r.ui} className="border-b border-[var(--rule)]/60 align-top">
                  <td className="py-2 pr-3">{r.ui}</td>
                  <td className="py-2 pr-3">{r.value}</td>
                  <td className="py-2 text-[var(--muted)]">{r.source}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </section>

      <section aria-labelledby="decisions">
        <h2 id="decisions" className="serif text-xl font-semibold">Decisions and corrections</h2>
        <ul className="mt-3 list-disc space-y-2 pl-5 text-sm">
          <li>
            <strong>Key semantics.</strong> phen_1 is keyed to gen_1.fam; phen_4 to gen_3.fam; genotype rows are shared. Checked with
            PLAG1 and DGAT1-candidate effects.
          </li>
          <li>
            <strong>ARSI dropped.</strong> RefSeq XP_002689338.1 has N at residue 301 (G at 302). Not guessed.
          </li>
          <li>
            <strong>Protein model.</strong> ESM-2 masked-marginal LLR for three candidates. Concordance with observed effects failed;
            reported as failed.
          </li>
          <li>
            <strong>Inbreeding.</strong> Per-1% slope corrected to -0.74 kg. The 24% het-plus-homozygote attenuation was a
            population-structure artefact; the homozygote-only estimate is 2%. Only the slope is shown on the site.
          </li>
          <li>
            <strong>Manhattan ties.</strong> Thinning is now deterministic (ties broken by position); the ARS mapping uses an explicit
            ordering.
          </li>
        </ul>
      </section>

      <section aria-labelledby="files">
        <h2 id="files" className="serif text-xl font-semibold">Export manifest</h2>
        <div className="mt-3 overflow-x-auto">
          <table className="w-full text-left text-xs">
            <thead>
              <tr className="border-b border-[var(--rule)]">
                <th className="py-2 pr-3">File</th>
                <th className="py-2 pr-3">Bytes</th>
                <th className="py-2 pr-3">Gzip bytes</th>
                <th className="py-2">Source files</th>
              </tr>
            </thead>
            <tbody>
              {manifest.files.map((f) => (
                <tr key={f.file} className="border-b border-[var(--rule)]/60 align-top">
                  <td className="py-1 pr-3">{f.file}</td>
                  <td className="py-1 pr-3">{f.bytes.toLocaleString("en-NZ")}</td>
                  <td className="py-1 pr-3">{f.gzip_bytes?.toLocaleString("en-NZ") ?? "—"}</td>
                  <td className="py-1 text-[var(--muted)]">{(f.source_files ?? []).join(", ")}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </section>
    </div>
  );
}
