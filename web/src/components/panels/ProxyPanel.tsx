"use client";

import { useState } from "react";
import { GenotypeClassSchema, InbreedingSchema, type GenotypeClass, type Inbreeding } from "@/lib/schemas";
import { useData } from "@/lib/useData";
import { BoxPlot, OKABE, type BoxStats } from "../charts";
import { LOCUS_LABELS, PanelState, Source, num } from "./common";

const LOCI = ["PLCD4", "GALNT2", "FGD4", "DPF2", "MUS81", "MTRF1", "ARSI", "Chr2:22", "MYH1"];

export default function ProxyPanel({ active }: { active: boolean }) {
  const { data, error } = useData<GenotypeClass>("genotype_class_stats.json", GenotypeClassSchema, active);
  const [locus, setLocus] = useState("PLCD4");
  const inb = useData<Inbreeding>("inbreeding.json", InbreedingSchema, active);
  return (
    <PanelState loading={active && !data && !error} error={error}>
      <label className="mb-2 flex items-center gap-2 text-sm">
        Locus
        <select
          value={locus}
          onChange={(e) => setLocus(e.target.value)}
          className="rounded border border-[var(--rule)] bg-[var(--bg)] px-2 py-1"
        >
          {LOCI.map((l) => (
            <option key={l} value={l}>
              {LOCUS_LABELS[l] ?? l}
            </option>
          ))}
        </select>
      </label>
      {data ? <ProxyBody data={data} locus={locus} /> : null}
      {inb.data ? (
        <p className="mt-4 text-sm">
          Inbreeding depression on body weight: <strong>{inb.data.inbreeding.slope_kg_per_1pct_F.toFixed(2)} kg per 1% inbreeding</strong>{" "}
          (95% bootstrap CI {inb.data.inbreeding.slope_kg_per_1pct_F_ci[0].toFixed(2)} to{" "}
          {inb.data.inbreeding.slope_kg_per_1pct_F_ci[1].toFixed(2)}; n = {inb.data.inbreeding.n_animals.toLocaleString("en-NZ")}). The paper reports -0.77.
        </p>
      ) : null}
      {inb.data ? <Source>inbreeding.json (3c; other_7 and phen_1.1). Slope on raw yield deviation; F assumed to be a fraction.</Source> : null}
    </PanelState>
  );
}

function ProxyBody({ data, locus }: { data: GenotypeClass; locus: string }) {
  const row = data.rows.find((r) => r.trait === "bodyweight" && r.locus === locus && r.estimator === "B");
  if (!row) return <p className="text-sm">No LOSO-adjusted class summary for this locus.</p>;
  const names = ["Reference (0)", "Carrier (1)", "Homozygous (2)"];
  const colors = [OKABE.grey, OKABE.sky, OKABE.vermillion];
  const stats: BoxStats[] = ["0", "1", "2"]
    .map((k, i) => {
      const c = row.class_stats[k];
      if (!c) return null;
      return {
        label: names[i] ?? k,
        n: c.n,
        mean: c.mean,
        p05: c.quantiles.p05,
        p25: c.quantiles.p25,
        p50: c.quantiles.p50,
        p75: c.quantiles.p75,
        p95: c.quantiles.p95,
        color: colors[i] ?? OKABE.grey,
      };
    })
    .filter((s): s is BoxStats => s !== null);
  const hom = row.class_stats["2"];
  const ref = row.class_stats["0"];
  return (
    <div>
      <BoxPlot
        stats={stats}
        ariaLabel={`Box plot of LOSO-adjusted body weight by ${LOCUS_LABELS[locus] ?? locus} genotype class`}
        yLabel="Body weight, LOSO-adjusted (kg)"
      />
      <p className="mt-2 text-sm">
        {hom ? (
          <>
            Homozygous mutants (n = {hom.n}) average <strong>{num(hom.mean, 1)} kg</strong> against{" "}
            {ref ? num(ref.mean, 1) : "—"} kg for reference animals. Box = interquartile range, line = median,
            dot = mean, whiskers = 5th to 95th percentile.
          </>
        ) : (
          "No homozygotes in this cohort."
        )}
      </p>
      <p className="mt-2 text-xs text-[var(--muted)]">
        Adjusted phenotype (estimator B) subtracts the LOSO population-structure term for the locus window, the paper-like
        analysis. Raw yield-deviation estimates are not shown here because raw het effects are confounded by structure.
      </p>
      <Source>genotype_class_stats.json (estimator B, bodyweight). Phase 1 class summaries.</Source>
    </div>
  );
}
