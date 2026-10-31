"use client";

import { useMemo, useState } from "react";
import { MatingInputsSchema, type MatingInputs } from "@/lib/schemas";
import { useData } from "@/lib/useData";
import {
  calfRisk,
  herdOutcome,
  nationalCalibration,
  type Breed,
  type LocusParams,
  type ParentState,
  type Status,
} from "@/lib/mating";
import { int, num, pct, PanelState, Caption, Source } from "./common";

const STATUSES: { v: Status; l: string }[] = [
  { v: "non", l: "Non-carrier" },
  { v: "carrier", l: "Carrier" },
  { v: "affected", l: "Affected (homozygous)" },
  { v: "unknown", l: "Unknown (population)" },
];

export default function SimulatorPanel({ active }: { active: boolean }) {
  const { data, error } = useData<MatingInputs>("mating_inputs.json", MatingInputsSchema, active);
  return (
    <PanelState loading={active && !data && !error} error={error}>
      {data ? <Simulator inputs={data} /> : null}
    </PanelState>
  );
}

function Simulator({ inputs }: { inputs: MatingInputs }) {
  const [hfShare, setHfShare] = useState(0.63);
  const [depletion, setDepletion] = useState(false);
  const [includeMyh1, setIncludeMyh1] = useState(false);
  const [avoid, setAvoid] = useState<"random" | "avoid">("random");
  const [calves, setCalves] = useState(1000);
  const [sireStatus, setSireStatus] = useState<Record<string, Status>>({});
  const [damStatus, setDamStatus] = useState<Record<string, Status>>({});

  const breed: Breed = { hf: hfShare, j: 1 - hfShare };
  const loci: LocusParams[] = useMemo(
    () =>
      inputs.loci
        .filter((l) => includeMyh1 || l.in_body_weight_set)
        .map((l) => ({
          id: l.id,
          label: l.label,
          qHF: l.q_HF,
          qJ: l.q_J,
          depletionRatio: l.depletion_ratio_primary,
          homKg: l.effect.hom_kg,
          hetKg: l.effect.het_kg,
          useHetEffect: l.id === "DPF2_MUS81",
        })),
    [inputs, includeMyh1],
  );
  const sire: ParentState = { status: sireStatus, breed };
  const dam: ParentState = { status: damStatus, breed };
  const calf = calfRisk(loci, sire, dam, depletion);
  const herd = herdOutcome(loci, breed, calves, avoid === "avoid", depletion);
  const herdRandom = herdOutcome(loci, breed, calves, false, depletion);
  const herdAvoid = herdOutcome(loci, breed, calves, true, depletion);

  const meanQ: Record<string, number> = {};
  for (const l of inputs.loci) meanQ[l.id] = l.national_mean_q;
  const calib = nationalCalibration(
    loci,
    { femalesPerYear: inputs.national.females_born_per_year, meanQ },
    false,
  );
  const calibDep = nationalCalibration(
    loci,
    { femalesPerYear: inputs.national.females_born_per_year, meanQ },
    true,
  );

  const setStatus = (set: (f: (s: Record<string, Status>) => Record<string, Status>) => void, id: string, v: Status) =>
    set((s) => ({ ...s, [id]: v }));

  return (
    <div className="space-y-4 text-sm">
      <div className="grid gap-3 sm:grid-cols-2">
        <label className="flex flex-col gap-1">
          Herd breed profile: Holstein-Friesian share {pct(hfShare, 0)} (Jersey {pct(1 - hfShare, 0)})
          <input type="range" min={0} max={1} step={0.01} value={hfShare} onChange={(e) => setHfShare(Number(e.target.value))} />
        </label>
        <label className="flex flex-col gap-1">
          Calves per year in herd: {int(calves)}
          <input type="range" min={100} max={5000} step={100} value={calves} onChange={(e) => setCalves(Number(e.target.value))} />
        </label>
        <label className="flex items-center gap-2">
          <input type="checkbox" checked={depletion} onChange={(e) => setDepletion(e.target.checked)} />
          Adjust for observed homozygote depletion
        </label>
        <label className="flex items-center gap-2">
          <input type="checkbox" checked={includeMyh1} onChange={(e) => setIncludeMyh1(e.target.checked)} />
          Include MYH1 (BCS-suggestive, body-weight effect not significant)
        </label>
      </div>

      <div className="grid gap-4 sm:grid-cols-2">
        <ParentTable title="Sire" loci={inputs.loci} status={sireStatus} onChange={(id, v) => setStatus(setSireStatus, id, v)} />
        <ParentTable title="Dam" loci={inputs.loci} status={damStatus} onChange={(id, v) => setStatus(setDamStatus, id, v)} />
      </div>

      <div>
        <h3 className="font-semibold">Calf outcome for this mating</h3>
        <table className="mt-1 w-full text-left text-xs">
          <thead>
            <tr className="border-b border-[var(--rule)]">
              <th className="py-1">Locus</th>
              <th className="py-1">P(affected calf)</th>
              <th className="py-1">P(carrier calf)</th>
              <th className="py-1">Expected kg</th>
            </tr>
          </thead>
          <tbody>
            {calf.perLocus.map((r) => (
              <tr key={r.id} className="border-b border-[var(--rule)]/50">
                <td className="py-1">{r.label}</td>
                <td className="py-1">{pct(r.pHom, 2)}</td>
                <td className="py-1">{pct(r.pHet, 2)}</td>
                <td className="py-1">{num(r.penaltyKg, 2)}</td>
              </tr>
            ))}
          </tbody>
        </table>
        <p className="mt-2">
          Probability the calf is affected at one or more loci (independence): <strong>{pct(calf.pAffectedAny, 2)}</strong>. Expected
          body-weight change: <strong>{num(calf.expectedPenaltyKg, 1)} kg</strong>.
        </p>
      </div>

      <div>
        <h3 className="font-semibold">Herd over one year: random mating vs avoiding carrier × carrier matings</h3>
        <fieldset className="mt-1 flex gap-4">
          <legend className="sr-only">Mating strategy</legend>
          <label className="flex items-center gap-1">
            <input type="radio" name="avoid" checked={avoid === "random"} onChange={() => setAvoid("random")} /> Random
          </label>
          <label className="flex items-center gap-1">
            <input type="radio" name="avoid" checked={avoid === "avoid"} onChange={() => setAvoid("avoid")} /> Avoid carrier × carrier
          </label>
        </fieldset>
        <table className="mt-2 w-full text-left text-xs">
          <thead>
            <tr className="border-b border-[var(--rule)]">
              <th className="py-1 pr-3 pr-3 align-top">Strategy</th>
              <th className="py-1 pr-3 pr-3 align-top">Affected calves per year</th>
              <th className="py-1 pr-3 pr-3 align-top">Homozygote kg change per year</th>
              <th className="py-1 pr-3 pr-3 align-top">Net kg change per year</th>
            </tr>
          </thead>
          <tbody>
            <tr className="border-b border-[var(--rule)]/50">
              <td className="py-1 pr-3 align-top">Random</td>
              <td className="py-1 pr-3 align-top">{num(herdRandom.affectedPerYear, 1)}</td>
              <td className="py-1 pr-3 align-top">{int(herdRandom.homKgChangePerYear)}</td>
              <td className="py-1 pr-3 align-top">{int(herdRandom.netKgChangePerYear)}</td>
            </tr>
            <tr className="border-b border-[var(--rule)]/50">
              <td className="py-1 pr-3 align-top">Avoid carrier × carrier</td>
              <td className="py-1 pr-3 align-top">{num(herdAvoid.affectedPerYear, 1)}</td>
              <td className="py-1 pr-3 align-top">{int(herdAvoid.homKgChangePerYear)}</td>
              <td className="py-1 pr-3 align-top">{int(herdAvoid.netKgChangePerYear)}</td>
            </tr>
          </tbody>
        </table>
        <p className="mt-1 text-xs">
          Selected strategy: {avoid === "avoid" ? "avoid" : "random"}; {num(herd.affectedPerYear, 1)} affected calves,{" "}
          {int(herd.homKgChangePerYear)} kg from homozygotes and {int(herd.netKgChangePerYear)} kg net per year (negative =
          loss). The DPF2/MUS81 heterozygote effect is positive (+2.3 kg), so net change can be smaller than the homozygote loss.
        </p>
      </div>

      <div>
        <h3 className="font-semibold">Calibration: random mating at national frequencies vs the national export</h3>
        <table className="mt-1 w-full text-left text-xs">
          <thead>
            <tr className="border-b border-[var(--rule)]">
              <th className="py-1">Locus</th>
              <th className="py-1">Simulated (no depletion)</th>
              <th className="py-1">Simulated (with depletion)</th>
              <th className="py-1">National export</th>
            </tr>
          </thead>
          <tbody>
            {calib.map((r, i) => (
              <tr key={r.id} className="border-b border-[var(--rule)]/50">
                <td className="py-1">{r.label}</td>
                <td className="py-1">{int(r.perYear)}</td>
                <td className="py-1">{int(calibDep[i]?.perYear ?? 0)}</td>
                <td className="py-1">{int(inputs.national.export_per_mutation_scaled[r.id === "DPF2_MUS81" ? "DPF2" : r.id] ?? 0)}</td>
              </tr>
            ))}
          </tbody>
        </table>
        <Caption>
          Units: females per year (940,100 born per year). The random-mating calibration ignores pedigree structure, so it sits below
          the export. The export draws sire genotypes from the pedigree, and carrier sires are over-represented. The paper reports
          about {pct(inputs.national.paper_homozygous_rate, 0)} and {int(inputs.national.paper_scaled_homozygotes)} per year.
        </Caption>
      </div>
      <Caption>
        Model: transmission 0, 0.5, 1 or the breed-weighted frequency for non-carrier, carrier, affected or unknown parents; calf
        homozygote = product of transmissions; loci combined assuming independence; DPF2 used for the DPF2/MUS81 locus (weak LD).
        Depletion ratios are purebred observed/expected and are applied to homozygote probabilities. Not a veterinary tool.
      </Caption>
      <Source>mating_inputs.json (frequencies, depletion, effects, national_burden, other_11 aggregate); math in web/src/lib/mating.ts (unit-tested).</Source>
    </div>
  );
}

function ParentTable({
  title,
  loci,
  status,
  onChange,
}: {
  title: string;
  loci: MatingInputs["loci"];
  status: Record<string, Status>;
  onChange: (id: string, v: Status) => void;
}) {
  return (
    <div>
      <h3 className="font-semibold">{title}</h3>
      <div className="mt-1 space-y-1">
        {loci.map((l) => (
          <label key={l.id} className="flex items-center justify-between gap-2">
            <span>{l.label}</span>
            <select
              value={status[l.id] ?? "unknown"}
              onChange={(e) => onChange(l.id, e.target.value as Status)}
              className="rounded border border-[var(--rule)] bg-[var(--bg)] px-2 py-1"
              aria-label={`${title} status for ${l.label}`}
            >
              {STATUSES.map((s) => (
                <option key={s.v} value={s.v}>
                  {s.l}
                </option>
              ))}
            </select>
          </label>
        ))}
      </div>
    </div>
  );
}
