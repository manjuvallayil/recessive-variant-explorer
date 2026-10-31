"use client";

import { useState } from "react";
import { RegionalIndexSchema, RegionalSchema, type Regional, type RegionalIndex } from "@/lib/schemas";
import { useData } from "@/lib/useData";
import { OKABE, ScatterCanvas, type Pt } from "../charts";
import { Caption, PanelState, Source, num } from "./common";

export default function ZoomPanel({ active }: { active: boolean }) {
  const index = useData<RegionalIndex>("regional/index.json", RegionalIndexSchema, active);
  const [slug, setSlug] = useState<string>("plcd4");
  const entry = index.data?.loci.find((l) => l.slug === slug);
  const reg = useData<Regional>(entry ? `regional/${entry.slug}.json` : "regional/plcd4.json", RegionalSchema, active && !!index.data);
  return (
    <PanelState loading={active && !index.data && !index.error} error={index.error ?? reg.error}>
      {index.data ? (
        <label className="mb-2 flex items-center gap-2 text-sm">
          Locus
          <select value={slug} onChange={(e) => setSlug(e.target.value)} className="rounded border border-[var(--rule)] bg-[var(--bg)] px-2 py-1">
            {index.data.loci.map((l) => (
              <option key={l.slug} value={l.slug}>
                {l.locus}
              </option>
            ))}
          </select>
        </label>
      ) : null}
      {reg.data ? <RegionalBody reg={reg.data} /> : reg.error ? null : <p className="text-sm text-[var(--muted)]">Loading region…</p>}
    </PanelState>
  );
}

function RegionalBody({ reg }: { reg: Regional }) {
  const [lo, hi] = reg.window_umd;
  const cols = reg.columns;
  const pts: Pt[] = [];
  const colours: string[] = [];
  for (let i = 0; i < cols.bp_umd.length; i += 1) {
    const bp = cols.bp_umd[i] ?? 0;
    const y = (cols.neg_log10_pd_x100[i] ?? 0) / 100;
    const r2 = cols.r2_x1000 ? (cols.r2_x1000[i] ?? 0) / 1000 : null;
    pts.push({ x: bp, y, size: 2 });
    colours.push(r2 === null ? "#9a9a9a" : viridisHex(r2));
  }
  const indexed = pts.map((p, i) => ({ ...p, color: colours[i] }));
  const cand = reg.candidate_variant;
  const candPt: Pt = { x: cand.bp_umd, y: Math.min(cand.neg_log10_pd, 40), color: "#D55E00", size: 5 };
  const yMax = Math.min(40, Math.max(10, Math.max(...pts.map((p) => p.y))));
  const summary = `Regional plot for ${reg.locus} (${reg.x_axis_assembly}), ${cols.bp_umd.length} variants in window. Candidate ${cand.snp} at ${cand.bp_umd}, -log10 P (dominance) ${num(cand.neg_log10_pd, 1)}. ${reg.ld_available ? "Points coloured by r2 to the candidate." : "LD colouring not available for this window."}`;
  return (
    <div>
      <ScatterCanvas
        points={indexed}
        overlay={[candPt]}
        xFormat={(v) => `${(v / 1e6).toFixed(2)} Mb`}
        xDomain={[lo, hi]}
        yDomain={[0, yMax]}
        height={260}
        ariaLabel={`Regional association plot for ${reg.locus}`}
        xLabel={`chr${reg.chr} position (${reg.x_axis_assembly}); orange = candidate`}
        yLabel="-log10 P (dom.)"
        summary={summary}
      />
      <p className="mt-1 text-xs">
        Candidate <strong>{cand.snp}</strong> at {cand.bp_umd.toLocaleString("en-NZ")}, -log10 P = {num(cand.neg_log10_pd, 1)}.
        {reg.ld_available ? " Colour: r2 to the candidate (dark = low, yellow = high)." : " LD colouring not available for this window (grey points)."}
      </p>
      <GeneTrack reg={reg} />
      <Caption>
        Gene track: Ensembl ARS-UCD1.2 protein-coding genes lifted to UMD3.1 by piecewise-linear interpolation (annotation caveat; see Methods).
        {reg.gene ? null : " This window has no single candidate gene."}
      </Caption>
      <Source>regional/{reg.slug}.json (gen_1 window, other_1.1 dominance statistics, other_10 mapping)</Source>
    </div>
  );
}

function GeneTrack({ reg }: { reg: Regional }) {
  const [lo, hi] = reg.window_umd;
  const width = 600;
  const left = 10;
  const w = width - left * 2;
  const genes = reg.genes.filter((g) => g.end_umd >= lo && g.start_umd <= hi);
  const rows = Math.min(6, Math.max(1, genes.length));
  const height = 16 + rows * 16;
  return (
    <svg viewBox={`0 0 ${width} ${height}`} className="mt-2 w-full h-auto" role="img" aria-label={`Gene track: ${genes.length} genes in window`}>
      <desc>{genes.map((g) => g.symbol ?? "unnamed").join(", ")}</desc>
      {genes.slice(0, 40).map((g, i) => {
        const x0 = left + ((Math.max(g.start_umd, lo) - lo) / (hi - lo)) * w;
        const x1 = left + ((Math.min(g.end_umd, hi) - lo) / (hi - lo)) * w;
        const y = 8 + (i % rows) * 16;
        const fill = g.is_candidate_gene ? OKABE.vermillion : OKABE.grey;
        return (
          <g key={`${g.symbol ?? "g"}-${i}`}>
            <rect x={x0} y={y} width={Math.max(2, x1 - x0)} height={8} fill={fill} />
            <text x={x0} y={y - 1} className="fill-[var(--muted)] text-[9px]">
              {g.symbol ?? ""}
            </text>
          </g>
        );
      })}
    </svg>
  );
}

// Viridis-like sequential colour for r2 in [0,1] (colour-blind safe). Computed from d3's interpolateViridis.
function viridisHex(t: number): string {
  const clamp = Math.max(0, Math.min(1, t));
  const stops: number[][] = [
    [68, 1, 84],
    [59, 82, 139],
    [33, 145, 140],
    [94, 201, 98],
    [253, 231, 37],
  ];
  const seg = clamp * (stops.length - 1);
  const i = Math.min(stops.length - 2, Math.floor(seg));
  const f = seg - i;
  const a = stops[i] ?? [0, 0, 0];
  const b = stops[i + 1] ?? a;
  const c = a.map((v, k) => Math.round(v + ((b[k] ?? v) - v) * f));
  return `rgb(${c[0] ?? 0},${c[1] ?? 0},${c[2] ?? 0})`;
}
