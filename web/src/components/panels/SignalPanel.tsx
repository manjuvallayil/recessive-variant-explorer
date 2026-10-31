"use client";

import { DetectionSchema, ImputationSchema, SyndromeSchema, type Detection, type Imputation, type Syndrome } from "@/lib/schemas";
import { useData } from "@/lib/useData";
import { BarChart, IqrLines, OKABE, type BarItem } from "../charts";
import { Caption, PanelState, Source, num } from "./common";

export default function SignalPanel({ active }: { active: boolean }) {
  const syn = useData<Syndrome>("syndrome_signal.json", SyndromeSchema, active);
  const imp = useData<Imputation>("imputation_by_maf.json", ImputationSchema, active);
  const det = useData<Detection>("detection_summary.json", DetectionSchema, active);
  const error = syn.error ?? imp.error ?? det.error;
  const loading = active && !error && !(syn.data && imp.data && det.data);
  return (
    <PanelState loading={loading} error={error}>
      {syn.data && imp.data && det.data ? <SignalBody syn={syn.data} imp={imp.data} det={det.data} /> : null}
    </PanelState>
  );
}

function SignalBody({ syn, imp, det }: { syn: Syndrome; imp: Imputation; det: Detection }) {
  const s = syn.syndrome_signal;
  const overall: BarItem[] = [
    { label: "Body weight flag (lower weight)", value: s.scores.B0_body_weight_flag?.auc ?? 0.5, lo: s.scores.B0_body_weight_flag?.auc_ci[0], hi: s.scores.B0_body_weight_flag?.auc_ci[1], color: OKABE.blue },
    { label: "Logistic on 3 traits (5-fold CV)", value: s.scores.S_logistic_cv?.auc ?? 0.5, lo: s.scores.S_logistic_cv?.auc_ci[0], hi: s.scores.S_logistic_cv?.auc_ci[1], color: OKABE.green },
    { label: "Mahalanobis (unsupervised)", value: s.scores.B1_mahalanobis?.auc ?? 0.5, lo: s.scores.B1_mahalanobis?.auc_ci[0], hi: s.scores.B1_mahalanobis?.auc_ci[1], color: OKABE.grey },
  ];
  const perLocus: BarItem[] = Object.entries(syn.per_group).map(([g, v]) => ({
    label: `${g} (n=${v.n_homozygous_with_bodyweight}; exp. ${num(v.auc_expected_from_effect_size, 2)})`,
    value: v.auc_negated_bw_single_locus,
    lo: v.auc_negated_bw_single_locus_ci[0],
    hi: v.auc_negated_bw_single_locus_ci[1],
    color: OKABE.vermillion,
  }));
  const dt = det.qtl_truth;
  const imputeRows = imp.maf_bins.map((b) => ({
    x: (b.maf_min + b.maf_max) / 2,
    label: `${(b.maf_min * 100).toFixed(0)}–${(b.maf_max * 100).toFixed(0)}%`,
    p05: b.quantiles.p05,
    p25: b.quantiles.p25,
    p50: b.quantiles.p50,
    p75: b.quantiles.p75,
    p95: b.quantiles.p95,
  }));
  return (
    <div>
      <h3 className="text-sm font-semibold">Can body weight flag homozygous mutants? (AUROC, 95% bootstrap CI)</h3>
      <BarChart
        items={overall}
        domain={[0.5, 0.8]}
        format={(v) => v.toFixed(3)}
        ariaLabel="AUROC for the pooled any-homozygote label"
        xLabel="AUROC, pooled label (any of PLCD4, FGD4, GALNT2, DPF2/MUS81)"
        rowHeight={26}
      />
      <p className="mt-2 text-sm">
        Pooled label: <strong>{s.label_counts.any}</strong> homozygous animals with complete traits (base rate {num(s.base_rate * 100, 2)}%).
        Per-locus counts (all genotyped animals with body weight). Body weight here is the raw yield deviation, not the
        LOSO-adjusted phenotype used in the proxy panel:
      </p>
      <BarChart
        items={perLocus}
        domain={[0.5, 1]}
        format={(v) => v.toFixed(2)}
        ariaLabel="Single-locus AUROC of negated body weight, per mutant homozygote group"
        xLabel="AUROC of lower body weight for homozygous mutants"
        rowHeight={26}
      />
      <Caption>
        The expected AUROC from the effect size is the value for two equal-variance normal classes separated by the observed
        shift. With a 63 kg residual SD, a shift near 80 kg (PLCD4) caps AUROC near 0.82, so 0.9 is not reachable with this phenotype.
        DPF2/MUS81 has 110 homozygotes with a small shift (about 0.5 SD), which keeps its AUROC near 0.63.
      </Caption>
      <h3 className="mt-5 text-sm font-semibold">Imputation accuracy by minor allele frequency (AR², quantiles)</h3>
      <IqrLines
        rows={imputeRows}
        ariaLabel="Imputation AR-squared by minor allele frequency bin, median and interquartile range"
        yLabel="AR² (imputation R²)"
        xLabel="Minor allele frequency bin (equal spacing)"
        yDomain={[0.5, 1]}
      />
      <Caption>
        Candidate variants: AR² {imp.candidate_variants.filter((c) => c.exceeds_0_95).length} of {imp.candidate_variants.length} above 0.95.
        Low-frequency variants (below 2%) carry most of the imputation uncertainty.
      </Caption>
      <h3 className="mt-5 text-sm font-semibold">Simulated recessive QTL: detection at P &lt; 5 × 10⁻⁸</h3>
      <table className="mt-2 w-full text-left text-xs">
        <caption className="sr-only">Detections of 30 simulated recessive QTL by MAF bin</caption>
        <thead>
          <tr className="border-b border-[var(--rule)]">
            <th className="py-1">MAF bin</th>
            <th className="py-1">QTL</th>
            <th className="py-1">Additive detected</th>
            <th className="py-1">Dominance detected</th>
          </tr>
        </thead>
        <tbody>
          {dt.by_maf_bin.map((r) => (
            <tr key={r.bin} className="border-b border-[var(--rule)]/50">
              <td className="py-1">{r.bin}</td>
              <td className="py-1">{r.n_qtl}</td>
              <td className="py-1">{r.n_detected_additive}</td>
              <td className="py-1">{r.n_detected_dominance}</td>
            </tr>
          ))}
          <tr className="font-semibold">
            <td className="py-1">All</td>
            <td className="py-1">{dt.by_maf_bin.reduce((a, r) => a + r.n_qtl, 0)}</td>
            <td className="py-1">{dt.n_detected_additive_truth}</td>
            <td className="py-1">{dt.n_detected_dominance_truth}</td>
          </tr>
        </tbody>
      </table>
      <Caption>
        Threshold: two-sided P &lt; {dt.threshold_p.toExponential(0)} on the simulation&apos;s own P-values (ExtFig9b). The closest
        additive P-value among the 30 truth QTL is {dt.min_p_additive_truth.toExponential(1)}, just above the threshold, so the additive
        test detects none of them. The recessive and genotype-class tests are not in the source file.
      </Caption>
      <Source>syndrome_signal.json (4b), imputation_by_maf.json (other_9), detection_summary.json (4c, ExtFig9b and gen_11.1).</Source>
    </div>
  );
}

