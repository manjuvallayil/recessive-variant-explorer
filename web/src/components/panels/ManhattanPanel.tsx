"use client";

import { useMemo, useState } from "react";
import { ChipSchema, ManhattanSchema, QtlContrastSchema, type Chip, type Manhattan, type QtlContrast } from "@/lib/schemas";
import { useData } from "@/lib/useData";
import { OKABE, ScatterCanvas, type Pt } from "../charts";
import { Caption, PanelState, Source, num } from "./common";

const Y_CAP = 60;

export default function ManhattanPanel({ active }: { active: boolean }) {
  const dom = useData<Manhattan>("manhattan_dom.json", ManhattanSchema, active);
  const add = useData<Manhattan>("manhattan_add.json", ManhattanSchema, active);
  const chip = useData<Chip>("chip_vs_sequence.json", ChipSchema, active);
  const qtl = useData<QtlContrast>("qtl_contrast.json", QtlContrastSchema, active);
  const [model, setModel] = useState<"dom" | "add">("dom");
  const [trait, setTrait] = useState<"bodyweight" | "stature">("bodyweight");
  const [chipOn, setChipOn] = useState(false);
  const error = dom.error ?? add.error ?? chip.error ?? qtl.error;
  const loading = active && !error && !(dom.data && add.data && chip.data && qtl.data);
  return (
    <PanelState loading={loading} error={error}>
      {dom.data && add.data && chip.data && qtl.data ? (
        <ManhattanBody
          dom={dom.data}
          add={add.data}
          chip={chip.data}
          qtl={qtl.data}
          model={model}
          setModel={setModel}
          trait={trait}
          setTrait={setTrait}
          chipOn={chipOn}
          setChipOn={setChipOn}
        />
      ) : null}
    </PanelState>
  );
}

function ManhattanBody(props: {
  dom: Manhattan;
  add: Manhattan;
  chip: Chip;
  qtl: QtlContrast;
  model: "dom" | "add";
  setModel: (m: "dom" | "add") => void;
  trait: "bodyweight" | "stature";
  setTrait: (t: "bodyweight" | "stature") => void;
  chipOn: boolean;
  setChipOn: (v: boolean) => void;
}) {
  const { dom, add, chip, qtl, model, setModel, trait, setTrait, chipOn, setChipOn } = props;
  const src = model === "dom" ? dom : add;
  const cols = src.traits[trait].columns;
  const { points, offsets, total, mids } = useMemo(
    () => genomePoints(cols.chr, cols.bp_umd, cols.neg_log10_p_x100),
    [cols],
  );
  const chipPts = useMemo(() => {
    if (!chipOn || model !== "dom") return undefined;
    const c = chip.columns;
    return c.chr.map((ch, i) => ({
      x: (offsets[ch] ?? 0) + (c.bp_umd[i] ?? 0),
      y: (c.neg_log10_pd_x100[i] ?? 0) / 100,
      color: "#D55E00",
      size: 2,
    })) as Pt[];
  }, [chipOn, model, chip, offsets]);
  const threshold = src.gws_threshold_neg_log10_p;
  const maxY = Math.max(...points.map((p) => p.y));
  const summary = `Genome-wide ${model === "dom" ? "dominance" : "additive"} scan for ${trait}. ${src.traits[trait].meta.n_thinned_points} thinned points; dashed line at -log10 P = ${threshold.toFixed(2)} (P = 5e-8). Values above ${Y_CAP} are drawn at ${Y_CAP}. Maximum observed -log10 P = ${maxY.toFixed(1)}.`;
  return (
    <div>
      <div className="mb-2 flex flex-wrap items-center gap-3 text-sm">
        <fieldset className="flex gap-2">
          <legend className="sr-only">Model</legend>
          <label className="flex items-center gap-1">
            <input type="radio" name="model" checked={model === "dom"} onChange={() => setModel("dom")} /> Dominance
          </label>
          <label className="flex items-center gap-1">
            <input type="radio" name="model" checked={model === "add"} onChange={() => setModel("add")} /> Additive
          </label>
        </fieldset>
        <select value={trait} onChange={(e) => setTrait(e.target.value as "bodyweight" | "stature")} className="rounded border border-[var(--rule)] bg-[var(--bg)] px-2 py-1">
          <option value="bodyweight">Body weight</option>
          <option value="stature">Stature</option>
        </select>
        <label className="flex items-center gap-1" title="Dominance model only">
          <input type="checkbox" checked={chipOn} disabled={model !== "dom"} onChange={(e) => setChipOn(e.target.checked)} />
          Overlay 50k chip (BovineSNP50)
        </label>
      </div>
      <ScatterCanvas
        points={points.map((p) => ({ ...p, y: Math.min(p.y, Y_CAP) }))}
        overlay={chipPts}
        xDomain={[0, total]}
        yDomain={[0, Y_CAP]}
        threshold={threshold}
        height={300}
        ariaLabel={`Manhattan plot, ${trait}, ${model === "dom" ? "dominance" : "additive"} model`}
        xLabel="Chromosome (UMD3.1)"
        xTicks={mids}
        yLabel="-log10 P"
        summary={summary}
      />
      <Caption>{summary} Orange points: 50k chip SNPs (dominance model only). Chromosome offsets are cumulative UMD3.1 positions.</Caption>
      {model === "add" && trait === "bodyweight" ? (
        <Caption>
          PLAG1 (chr14) is the known additive locus. Its -log10 P is 480, above the axis cap. The additive model finds the PLAG1 locus and no other QTL in the list below.
        </Caption>
      ) : null}
      <QtlTable qtl={qtl} />
      <Source>manhattan_dom.json, manhattan_add.json, chip_vs_sequence.json, qtl_contrast.json (other_1.1, other_8.1, other_15). Thinning rule in each file&apos;s meta.</Source>
    </div>
  );
}

function genomePoints(chr: number[], bp: number[], nlp: number[]) {
  const maxBp = new Map<number, number>();
  chr.forEach((c, i) => maxBp.set(c, Math.max(maxBp.get(c) ?? 0, bp[i] ?? 0)));
  const offsets: Record<number, number> = {};
  let acc = 0;
  const mids: { x: number; label: string }[] = [];
  for (let c = 1; c <= 29; c += 1) {
    offsets[c] = acc;
    const len = maxBp.get(c) ?? 0;
    mids.push({ x: acc + len / 2, label: String(c) });
    acc += len + 5e6;
  }
  const points: Pt[] = chr.map((c, i) => ({
    x: (offsets[c] ?? 0) + (bp[i] ?? 0),
    y: (nlp[i] ?? 0) / 100,
    color: c % 2 === 0 ? OKABE.blue : OKABE.sky,
    size: 1.4,
  }));
  return { points, offsets, total: acc, mids };
}

function QtlTable({ qtl }: { qtl: QtlContrast }) {
  const rows = qtl.rows.filter((r) => r.is_primary_8);
  const t = qtl.gws_threshold_neg_log10_p;
  const candInvisible = rows.filter((r) => r.bodyweight.candidate_add_neg_log10_p < t).length;
  const winInvisible = rows.filter((r) => r.bodyweight.window_best_add_neg_log10_p < t).length;
  return (
    <div className="mt-4">
      <p className="text-sm">
        At the candidate variant, <strong>{candInvisible} of {rows.length}</strong> QTL are below genome-wide significance under
        the additive model. Within ±1 Mb of the candidate, <strong>{winInvisible} of {rows.length}</strong> are: the DPF2/MUS81
        window contains a nearby additive peak about 0.8 Mb from the candidate, which the candidate itself does not carry.
      </p>
      <div className="mt-2 overflow-x-auto">
        <table className="w-full text-left text-xs">
          <caption className="sr-only">Body weight: best -log10 P in ±1 Mb, additive vs dominance</caption>
          <thead>
            <tr className="border-b border-[var(--rule)]">
              <th className="py-1 pr-2">QTL</th>
              <th className="py-1 pr-2">Additive (candidate)</th>
              <th className="py-1 pr-2">Additive (±1 Mb best)</th>
              <th className="py-1 pr-2">Dominance (candidate)</th>
            </tr>
          </thead>
          <tbody>
            {rows.map((r) => (
              <tr key={r.qtl} className="border-b border-[var(--rule)]/50">
                <td className="py-1 pr-2">{r.qtl}</td>
                <td className="py-1 pr-2">{num(r.bodyweight.candidate_add_neg_log10_p, 2)}</td>
                <td className="py-1 pr-2">{num(r.bodyweight.window_best_add_neg_log10_p, 2)}</td>
                <td className="py-1 pr-2">{num(r.bodyweight.candidate_dom_d_neg_log10_p, 2)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      <p className="mt-1 text-[11px] text-[var(--muted)]">Threshold -log10 P = {num(t, 2)}. Values are -log10 P from Z (unclipped).</p>
    </div>
  );
}
