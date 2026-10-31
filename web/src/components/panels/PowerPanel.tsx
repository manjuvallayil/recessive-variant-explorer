"use client";

import { useMemo, useState } from "react";
import { PowerSchema, type Power } from "@/lib/schemas";
import { useData } from "@/lib/useData";
import { powerFor } from "@/lib/power";
import { PowerCurve } from "../charts";
import { Caption, PanelState, Source, num } from "./common";

const VARIANCES = [0.0005, 0.001, 0.002] as const;

export default function PowerPanel({ active }: { active: boolean }) {
  const { data, error } = useData<Power>("power_reference.json", PowerSchema, active);
  const [variance, setVariance] = useState<number>(0.001);
  const [n, setN] = useState<number>(50000);
  const [df, setDf] = useState<1 | 2>(1);
  const curve = useMemo(() => {
    const pts: { n: number; power: number }[] = [];
    for (let k = 1; k <= 60; k += 1) {
      const nn = Math.round((k * 100_000) / 60);
      pts.push({ n: nn, power: powerFor(variance, nn, df, 5e-8) });
    }
    return pts;
  }, [variance, df]);
  const current = powerFor(variance, n, df, 5e-8);
  return (
    <PanelState loading={active && !data && !error} error={error}>
      {data ? (
        <div>
          <div className="grid gap-3 text-sm sm:grid-cols-3">
            <label className="flex flex-col gap-1">
              Locus variance explained
              <select value={variance} onChange={(e) => setVariance(Number(e.target.value))} className="rounded border border-[var(--rule)] bg-[var(--bg)] px-2 py-1">
                {VARIANCES.map((v) => (
                  <option key={v} value={v}>
                    {(v * 100).toFixed(2)}%
                  </option>
                ))}
              </select>
            </label>
            <label className="flex flex-col gap-1">
              Test degrees of freedom
              <select value={df} onChange={(e) => setDf(Number(e.target.value) as 1 | 2)} className="rounded border border-[var(--rule)] bg-[var(--bg)] px-2 py-1">
                <option value={1}>1 (recessive / additive)</option>
                <option value={2}>2 (genotype class)</option>
              </select>
            </label>
            <label className="flex flex-col gap-1">
              Sample size: {n.toLocaleString("en-NZ")}
              <input type="range" min={1000} max={100000} step={500} value={n} onChange={(e) => setN(Number(e.target.value))} />
            </label>
          </div>
          <p className="mt-2 text-sm">
            Power at P &lt; 5 × 10⁻⁸: <strong>{num(current * 100, 1)}%</strong> with N = {n.toLocaleString("en-NZ")}, variance explained{" "}
            {(variance * 100).toFixed(2)}%, {df} df.
          </p>
          <PowerCurve
            points={curve}
            marker={{ n: 50000, power: powerFor(0.001, 50000, df, 5e-8), label: "Paper scenario (0.1%, N = 50,000)" }}
            ariaLabel="Power curve against sample size for the selected variance explained and degrees of freedom"
          />
          <Caption>
            Analytic power is the noncentral chi-square tail at the genome-wide threshold. The paper scenario is MAF 2.5%,
            heritability 0.25, 0.1% variance explained. Validation against the paper&apos;s Extended Data Fig. 1: maximum
            difference {data.validation.max_abs_diff_vs_source.toFixed(4)} (tolerance {data.validation.tolerance_abs}) for the cow
            recessive and genotype-class curves. The additive and sire-model curves are not reproduced by this closed form.
          </Caption>
          <Source>power_reference.json (4d; ExtFig1 reference). Live computation: web/src/lib/power.ts (tested against scipy).</Source>
        </div>
      ) : null}
    </PanelState>
  );
}
