"use client";

import { FrequenciesSchema, type Frequencies } from "@/lib/schemas";
import { useData } from "@/lib/useData";
import { BarChart, OKABE, type BarItem } from "../charts";
import { PanelState, Source, pct } from "./common";

const LOCI = ["PLCD4", "FGD4", "MTRF1", "MYH1", "GALNT2", "DPF2", "MUS81"];

export default function HookPanel({ active }: { active: boolean }) {
  const { data, error } = useData<Frequencies>("frequencies.json", FrequenciesSchema, active);
  return (
    <PanelState loading={active && !data && !error} error={error}>
      {data ? <HookBody data={data} /> : null}
    </PanelState>
  );
}

function HookBody({ data }: { data: Frequencies }) {
  const items: BarItem[] = [];
  for (const l of LOCI) {
    const b = data.loci[l]?.breeds;
    if (!b?.HF || !b.J) continue;
    items.push({ label: `${l} (Holstein-Friesian)`, value: b.HF.carrier_freq, color: OKABE.blue });
    items.push({ label: `${l} (Jersey)`, value: b.J.carrier_freq, color: OKABE.orange });
  }
  // Union of carriers under independence. DPF2 and MUS81 are shown separately as the paper counts them.
  const hf = (names: string[]) =>
    1 - names.reduce((p, n) => p * (1 - (data.loci[n]?.breeds.HF?.carrier_freq ?? 0)), 1);
  const unionSeparate = hf(["PLCD4", "FGD4", "MTRF1", "MYH1", "GALNT2", "DPF2", "MUS81"]);
  const unionOneLocus = hf(["PLCD4", "FGD4", "MTRF1", "MYH1", "GALNT2", "DPF2"]);
  return (
    <div>
      <p className="mb-2 text-sm">
        Share of animals that carry at least one copy of each candidate mutant allele (carriers, heterozygotes).
      </p>
      <BarChart
        items={items}
        domain={[0, Math.max(0.05, Math.ceil(Math.max(...items.map((i) => i.value)) * 10) / 10)]}
        format={(v) => pct(v, 1)}
        ariaLabel="Bar chart of carrier frequency by candidate variant and breed"
        xLabel="Carrier frequency (purebred animals)"
        rowHeight={24}
      />
      <p className="mt-3 text-sm">
        Holstein-Friesian animals carrying at least one of the seven variants (independence assumed):{" "}
        <strong>{pct(unionSeparate, 1)}</strong> with DPF2 and MUS81 counted separately, as in the paper, or{" "}
        <strong>{pct(unionOneLocus, 1)}</strong> with them counted as one locus. The paper reports more than 40%.
      </p>
      <Source>frequencies.json (gen_4 purebreds; other_6 breed fractions). Union calculation is in UI code and is not a pipeline export.</Source>
    </div>
  );
}
