"use client";

import { useState } from "react";
import { FarmSchema, PlmSchema, RnaSchema, type Farm, type Plm, type Rna } from "@/lib/schemas";
import { useData } from "@/lib/useData";
import { Heatmap } from "../charts";
import { Caption, PanelState, Source, num } from "./common";

type Tab = "plm" | "rna" | "farm";

export default function EvidencePanel({ active }: { active: boolean }) {
  const plm = useData<Plm>("plm_scores.json", PlmSchema, active);
  const rna = useData<Rna>("rna_phase3.json", RnaSchema, active);
  const farm = useData<Farm>("farm_profiles.json", FarmSchema, active);
  const [tab, setTab] = useState<Tab>("plm");
  const error = plm.error ?? rna.error ?? farm.error;
  const loading = active && !error && !(plm.data && rna.data && farm.data);
  return (
    <PanelState loading={loading} error={error}>
      <div role="tablist" aria-label="Evidence type" className="mb-3 flex flex-wrap gap-2 text-sm">
        {(
          [
            ["plm", "Protein model"],
            ["rna", "RNA splicing"],
            ["farm", "Research farm"],
          ] as [Tab, string][]
        ).map(([id, label]) => (
          <button
            key={id}
            role="tab"
            aria-selected={tab === id}
            onClick={() => setTab(id)}
            className={`rounded border px-3 py-1 ${tab === id ? "border-[var(--fg)] font-semibold" : "border-[var(--rule)] text-[var(--muted)]"}`}
          >
            {label}
          </button>
        ))}
      </div>
      {tab === "plm" && plm.data ? <PlmBody plm={plm.data} /> : null}
      {tab === "rna" && rna.data ? <RnaBody rna={rna.data} /> : null}
      {tab === "farm" && farm.data ? <FarmBody farm={farm.data} /> : null}
    </PanelState>
  );
}

function PlmBody({ plm }: { plm: Plm }) {
  const [gene, setGene] = useState("PLCD4");
  const g = plm.genes.find((x) => x.gene === gene) ?? plm.genes[0];
  if (!g) return null;
  const cols = g.matrix_columns.slice(2);
  const half = 12;
  const from = Math.max(0, g.position - 1 - half);
  const to = Math.min(g.matrix.length, g.position + half);
  const rows = g.matrix.slice(from, to);
  const values = rows.map((r) => r.slice(2).map((v) => (typeof v === "number" ? v : 0)));
  const labels = rows.map((r) => `${String(r[1])}${String(r[0])}`);
  const mutRow = g.position - 1 - from;
  const mutCol = cols.indexOf(g.mut);
  const hom = g.observed_bodyweight_hom_effect;
  return (
    <div>
      <div className="mb-2 flex flex-wrap items-center gap-2 text-sm">
        <label className="flex items-center gap-2">
          Candidate
          <select value={gene} onChange={(e) => setGene(e.target.value)} className="rounded border border-[var(--rule)] bg-[var(--bg)] px-2 py-1">
            {plm.genes.map((x) => (
              <option key={x.gene} value={x.gene}>
                {x.gene} {x.mutation}
              </option>
            ))}
          </select>
        </label>
      </div>
      <Heatmap
        rows={labels.map((l) => ({ label: l }))}
        columns={cols}
        values={values}
        highlight={mutRow >= 0 && mutCol >= 0 ? { row: mutRow, col: mutCol } : undefined}
        ariaLabel={`ESM-2 masked-marginal log-likelihood ratios around ${g.gene} ${g.mutation}`}
      />
      <Caption>
        Rows: residues around position {g.position}; columns: substitutions. Colour: log-likelihood ratio (darker = less tolerated).
        Outlined: {g.gene} {g.mutation} (wild type {g.wt}, mutant {g.mut}).
      </Caption>
      <p className="mt-2 text-sm">
        Mutation LLR <strong>{num(g.llr, 2)}</strong>, damage percentile <strong>{num(g.damage_percentile, 1)}</strong> (share of all
        possible substitutions in the protein that are less damaging). Observed homozygote effect (LOSO-adjusted):{" "}
        {hom?.B ? `${num(hom.B.hom_alt_kg, 1)} ± ${num(hom.B.se_kg, 1)} kg` : "not available"}.
      </p>
      <p className="mt-2 rounded border border-[var(--vermillion)]/50 p-2 text-xs">
        <strong>Concordance check failed.</strong> With three candidates, the protein-model ordering (MTRF1 &gt; DPF2 &gt; PLCD4) is the
        reverse of the observed homozygote effect ordering (PLCD4 &gt; DPF2 &gt; MTRF1). Three points cannot support a
        claim either way; this is recorded, not used as evidence. ARSI was dropped (residue mismatch). MUS81 p.Gly70* and splice
        variants are not scorable by this method.
      </p>
      <Source>plm_scores.json (facebook/esm2_t30_150M_UR50D, masked-marginal, RefSeq proteins; effects.json).</Source>
    </div>
  );
}

function RnaBody({ rna }: { rna: Rna }) {
  return (
    <div>
      <p className="text-sm">
        Splice-junction and expression values from 350 lactating cows. This export contains descriptive means only; no genotype-class
        comparison is included, so no effect size or p-value is shown.
      </p>
      <table className="mt-2 w-full text-left text-xs">
        <caption className="sr-only">Mean RNA values by gene</caption>
        <thead>
          <tr className="border-b border-[var(--rule)]">
            <th className="py-1">Measure</th>
            <th className="py-1">Gene</th>
            <th className="py-1">Mean</th>
            <th className="py-1">n</th>
          </tr>
        </thead>
        <tbody>
          {Object.entries(rna.rna.splice).map(([g, v]) => (
            <tr key={`s-${g}`} className="border-b border-[var(--rule)]/50">
              <td className="py-1">Splice junction (aggregate)</td>
              <td className="py-1">{g}</td>
              <td className="py-1">{num(v.mean_value, 3)}</td>
              <td className="py-1">{v.n_animals}</td>
            </tr>
          ))}
          {Object.entries(rna.rna.expression).map(([g, v]) => (
            <tr key={`e-${g}`} className="border-b border-[var(--rule)]/50">
              <td className="py-1">Expression</td>
              <td className="py-1">{g}</td>
              <td className="py-1">{num(v.mean_value, 2)}</td>
              <td className="py-1">{v.n_animals}</td>
            </tr>
          ))}
        </tbody>
      </table>
      <Caption>
        Carrier-versus-reference comparisons are the next step (Phase 3a output lacks class-level statistics). The paper reports
        splicing-based eQTL analysis for GALNT2 and FGD4; those comparisons are not reproduced here.
      </Caption>
      <Source>rna_phase3.json (phen_7.1, phen_7.2; 350 animals).</Source>
    </div>
  );
}

function FarmBody({ farm }: { farm: Farm }) {
  return (
    <div>
      <p className="text-sm">
        Research-farm cohorts (small groups, age, sire and breed matched in the paper). Only group sizes are exported so far.
      </p>
      <table className="mt-2 w-full text-left text-xs">
        <caption className="sr-only">Farm cohort sizes</caption>
        <thead>
          <tr className="border-b border-[var(--rule)]">
            <th className="py-1">Cohort</th>
            <th className="py-1">Group</th>
            <th className="py-1">n animals</th>
          </tr>
        </thead>
        <tbody>
          {Object.entries(farm.farm_profiles).flatMap(([cohort, groups]) =>
            Object.entries(groups).map(([g, v]) => (
              <tr key={`${cohort}-${g}`} className="border-b border-[var(--rule)]/50">
                <td className="py-1">{cohort}</td>
                <td className="py-1">{g}</td>
                <td className="py-1">{v.n_animals}</td>
              </tr>
            )),
          )}
        </tbody>
      </table>
      <Caption>
        With groups this small, no trait comparison is shown. Group differences with bootstrap intervals are the planned next step.
        Please read these as descriptive, not as evidence of a syndrome.
      </Caption>
      <Source>farm_profiles.json (gen_8.1 / phen_5.1, gen_7.1 / phen_6.1).</Source>
    </div>
  );
}
