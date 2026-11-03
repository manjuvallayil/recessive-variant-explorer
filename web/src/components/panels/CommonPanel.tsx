"use client";

import { DepletionSchema, FrequenciesSchema, NationalBurdenSchema, type Depletion, type Frequencies, type NationalBurden } from "@/lib/schemas";
import { useData } from "@/lib/useData";
import { BarChart, OKABE, type BarItem } from "../charts";
import { Caption, PanelState, Source, int, pct } from "./common";

export default function CommonPanel({ active }: { active: boolean }) {
  const nb = useData<NationalBurden>("national_burden.json", NationalBurdenSchema, active);
  const dep = useData<Depletion>("depletion.json", DepletionSchema, active);
  const freq = useData<Frequencies>("frequencies.json", FrequenciesSchema, active);
  const error = nb.error ?? dep.error ?? freq.error;
  const loading = active && !error && !(nb.data && dep.data && freq.data);
  return (
    <PanelState loading={loading} error={error}>
      {nb.data && dep.data && freq.data ? <CommonBody nb={nb.data} dep={dep.data} freq={freq.data} /> : null}
    </PanelState>
  );
}

function CommonBody({ nb, dep, freq }: { nb: NationalBurden; dep: Depletion; freq: Frequencies }) {
  const bwSet = nb.any_bw_significant_mutation;
  const six = nb.any_mutation;
  const perMut: BarItem[] = Object.entries(nb.mutations)
    .sort((a, b) => b[1].scaled_per_year - a[1].scaled_per_year)
    .map(([k, v]) => ({
      label: k === "MYH1" ? "MYH1 (BCS-suggestive)" : k,
      value: v.scaled_per_year,
      color: k === "MYH1" ? OKABE.grey : OKABE.blue,
    }));
  const depItems: BarItem[] = dep.rows.map((r) => ({
    label: `${r.gene} (${r.breed})`,
    value: r.observed_over_expected,
    color: r.lrt_p_one_sided < 0.05 ? OKABE.vermillion : OKABE.grey,
    note: `p=${r.lrt_p_one_sided.toExponential(1)}`,
  }));
  return (
    <div>
      <h3 className="text-sm font-semibold">Estimated mutant homozygous births per year (national export)</h3>
      <BarChart
        items={perMut}
        domain={[0, 4000]}
        format={(v) => int(v)}
        ariaLabel="Bar chart of estimated mutant homozygous births per year by mutation"
        xLabel="Homozygous births per year (2.8 million animals, scaled to 940,100 females)"
        rowHeight={24}
      />
      <p className="mt-3 text-sm">
        Body-weight-significant set (five variants; DPF2/MUS81 counted once):{" "}
        <strong>{pct(bwSet.rate_per_animal, 2)}</strong> of animals, about <strong>{int(bwSet.scaled_per_year)}</strong> females per
        year, against the paper&apos;s rounded {pct(nb.paper_homozygous_rate, 0)} and {int(nb.paper_scaled_homozygotes)}. Including
        MYH1 (BCS-suggestive only) gives {pct(six.rate_per_animal, 2)} and about {int(six.scaled_per_year)} per year.
      </p>
      <Caption>
        The national export uses sire and grandsire genotypes from the pedigree, so it includes elite-carrier sires. That is
        the main reason it sits above a random-mating calculation. The simulator&apos;s calibration panel shows the gap.
      </Caption>
      <h3 className="mt-5 text-sm font-semibold">Homozygote depletion (observed / expected, purebred animals)</h3>
      <BarChart
        items={depItems}
        domain={[0, 1.2]}
        format={(v) => v.toFixed(2)}
        ariaLabel="Bar chart of observed over expected mutant homozygotes by locus and breed"
        xLabel="Observed / expected homozygotes (1 = Hardy-Weinberg)"
        rowHeight={24}
      />
      <Caption>
        Red: one-sided likelihood-ratio test below 0.05. Values below 1 mean fewer homozygotes than Hardy-Weinberg predicts.
        Not every breed-locus combination has enough homozygotes to test.
      </Caption>
      <p className="mt-3 text-sm">
        Carrier frequency in Holstein-Friesian for PLCD4:{" "}
        <strong>{pct(freq.loci.PLCD4?.breeds.HF?.carrier_freq ?? 0, 1)}</strong>; in Jersey for MTRF1:{" "}
        <strong>{pct(freq.loci.MTRF1?.breeds.J?.carrier_freq ?? 0, 1)}</strong>.
      </p>
      <Source>national_burden.json (other_11, gen_10.1); depletion.json (gen_4 with other_6); frequencies.json.</Source>
    </div>
  );
}
