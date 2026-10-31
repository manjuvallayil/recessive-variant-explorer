"use client";

// Shared chart primitives. D3 is used for scales, paths and colour; React renders SVG.
// Dense layers (Manhattan, regional) draw to <canvas>. All charts carry a text summary
// (aria-label and <desc>) so the information is available to screen readers.

import { useEffect, useRef } from "react";
import { interpolateBlues, interpolateViridis, line, max, scaleLinear, scaleSequential } from "d3";

export const OKABE = {
  blue: "#0072B2",
  orange: "#E69F00",
  green: "#009E73",
  vermillion: "#D55E00",
  sky: "#56B4E9",
  purple: "#CC79A7",
  grey: "#8a8a8a",
};

export interface BarItem {
  label: string;
  value: number;
  lo?: number;
  hi?: number;
  color?: string;
  note?: string;
}

export function BarChart({
  items,
  domain,
  format,
  ariaLabel,
  width = 600,
  rowHeight = 30,
  xLabel,
}: {
  items: BarItem[];
  domain: [number, number];
  format: (v: number) => string;
  ariaLabel: string;
  width?: number;
  rowHeight?: number;
  xLabel?: string;
}) {
  const left = 215;
  const right = 90;
  const top = 8;
  const bottom = xLabel ? 42 : 24;
  const height = top + bottom + items.length * rowHeight;
  const x = scaleLinear().domain(domain).range([left, width - right]);
  const ticks = x.ticks(5);
  const summary = items
    .map((d) => `${d.label}: ${format(d.value)}${d.lo !== undefined && d.hi !== undefined ? ` (${format(d.lo)} to ${format(d.hi)})` : ""}`)
    .join("; ");
  return (
    <svg viewBox={`0 0 ${width} ${height}`} role="img" aria-label={ariaLabel} className="w-full h-auto">
      <desc>{summary}</desc>
      {ticks.map((t) => (
        <g key={t}>
          <line x1={x(t)} x2={x(t)} y1={top} y2={height - bottom} className="stroke-[var(--grid)]" />
          <text x={x(t)} y={height - bottom + 16} textAnchor="middle" className="fill-[var(--muted)] text-[11px]">
            {t}
          </text>
        </g>
      ))}
      {items.map((d, i) => {
        const y = top + i * rowHeight;
        const baseline = domain[0] > 0 ? domain[0] : 0;
        const x0 = Math.min(x(baseline), x(d.value));
        const x1 = Math.max(x(baseline), x(d.value));
        const color = d.color ?? OKABE.blue;
        return (
          <g key={d.label}>
            <text x={left - 10} y={y + rowHeight / 2 + 4} textAnchor="end" className="fill-[var(--fg)] text-[12px]">
              {d.label}
            </text>
            <rect x={x0} width={Math.max(1, x1 - x0)} y={y + 6} height={rowHeight - 12} fill={color} rx={2} />
            {d.lo !== undefined && d.hi !== undefined ? (
              <line x1={x(d.lo)} x2={x(d.hi)} y1={y + rowHeight / 2} y2={y + rowHeight / 2} className="stroke-[var(--fg)]" strokeWidth={1.5} />
            ) : null}
            <text x={width - right + 6} y={y + rowHeight / 2 + 4} className="fill-[var(--fg)] text-[11px]">
              {format(d.value)}
              {d.note ? <tspan className="fill-[var(--muted)]"> {d.note}</tspan> : null}
            </text>
          </g>
        );
      })}
      {xLabel ? (
        <text x={(left + width - right) / 2} y={height - 6} textAnchor="middle" className="fill-[var(--muted)] text-[11px]">
          {xLabel}
        </text>
      ) : null}
    </svg>
  );
}

export interface BoxStats {
  label: string;
  n: number;
  mean: number;
  p05: number;
  p25: number;
  p50: number;
  p75: number;
  p95: number;
  color: string;
}

export function BoxPlot({
  stats,
  ariaLabel,
  yLabel,
  width = 560,
  height = 320,
}: {
  stats: BoxStats[];
  ariaLabel: string;
  yLabel: string;
  width?: number;
  height?: number;
}) {
  const left = 56;
  const right = 16;
  const top = 12;
  const bottom = 46;
  const lo = Math.min(...stats.map((s) => s.p05));
  const hi = Math.max(...stats.map((s) => s.p95));
  const pad = (hi - lo) * 0.08;
  const y = scaleLinear().domain([lo - pad, hi + pad]).range([height - bottom, top]);
  const slot = (width - left - right) / Math.max(1, stats.length);
  const summary = stats
    .map((s) => `${s.label} (n=${s.n}): mean ${s.mean.toFixed(1)}, median ${s.p50.toFixed(1)} kg`)
    .join("; ");
  const yTicks = y.ticks(6);
  return (
    <svg viewBox={`0 0 ${width} ${height}`} role="img" aria-label={ariaLabel} className="w-full h-auto">
      <desc>{summary}</desc>
      {yTicks.map((t) => (
        <g key={t}>
          <line x1={left} x2={width - right} y1={y(t)} y2={y(t)} className="stroke-[var(--grid)]" />
          <text x={left - 8} y={y(t) + 4} textAnchor="end" className="fill-[var(--muted)] text-[11px]">
            {t}
          </text>
        </g>
      ))}
      <text transform={`translate(14 ${(top + height - bottom) / 2}) rotate(-90)`} textAnchor="middle" className="fill-[var(--muted)] text-[11px]">
        {yLabel}
      </text>
      {stats.map((s, i) => {
        const cx = left + slot * (i + 0.5);
        const w = Math.min(70, slot * 0.5);
        return (
          <g key={s.label}>
            <line x1={cx} x2={cx} y1={y(s.p05)} y2={y(s.p95)} className="stroke-[var(--fg)]" />
            <rect x={cx - w / 2} width={w} y={y(s.p75)} height={Math.max(1, y(s.p25) - y(s.p75))} fill={s.color} fillOpacity={0.35} stroke={s.color} />
            <line x1={cx - w / 2} x2={cx + w / 2} y1={y(s.p50)} y2={y(s.p50)} className="stroke-[var(--fg)]" strokeWidth={2} />
            <circle cx={cx} cy={y(s.mean)} r={4} fill={s.color} stroke="var(--fg)" strokeWidth={1} />
            <text x={cx} y={height - bottom + 16} textAnchor="middle" className="fill-[var(--fg)] text-[11px]">
              {s.label}
            </text>
            <text x={cx} y={height - bottom + 30} textAnchor="middle" className="fill-[var(--muted)] text-[10px]">
              n={s.n}
            </text>
          </g>
        );
      })}
    </svg>
  );
}

export function IqrLines({
  rows,
  ariaLabel,
  yLabel,
  xLabel,
  yDomain,
  width = 600,
  height = 280,
}: {
  rows: { x: number; p05: number; p25: number; p50: number; p75: number; p95: number; label: string }[];
  ariaLabel: string;
  yLabel: string;
  xLabel: string;
  yDomain: [number, number];
  width?: number;
  height?: number;
}) {
  const left = 52;
  const right = 16;
  const top = 12;
  const bottom = 58;
  // Rows are equally spaced (bins are uneven in MAF); x is the row index.
  const x = scaleLinear().domain([0, Math.max(1, rows.length - 1)]).range([left + 10, width - right - 10]);
  const [y0, y1] = yDomain;
  const y = scaleLinear().domain([y0, y1]).range([height - bottom, top]);
  const med = line<{ x: number; p50: number }>().x((_, i) => x(i)).y((d) => y(d.p50));
  const band = rows
    .map((d, i) => `${x(i)},${y(d.p75)}`)
    .concat(rows.slice().reverse().map((d, j) => `${x(rows.length - 1 - j)},${y(d.p25)}`))
    .join(" ");
  const summary = rows.map((r) => `${r.label}: median ${r.p50.toFixed(3)} (IQR ${r.p25.toFixed(3)}-${r.p75.toFixed(3)})`).join("; ");
  return (
    <svg viewBox={`0 0 ${width} ${height}`} role="img" aria-label={ariaLabel} className="w-full h-auto">
      <desc>{summary}</desc>
      {y.ticks(5).map((t) => (
        <g key={t}>
          <line x1={left} x2={width - right} y1={y(t)} y2={y(t)} className="stroke-[var(--grid)]" />
          <text x={left - 8} y={y(t) + 4} textAnchor="end" className="fill-[var(--muted)] text-[11px]">
            {t.toFixed(2)}
          </text>
        </g>
      ))}
      <polygon points={band} fill={OKABE.sky} fillOpacity={0.35} />
      <path d={med(rows) ?? ""} fill="none" stroke={OKABE.blue} strokeWidth={2} />
      {rows.map((r, i) => (
        <g key={r.label}>
          <circle cx={x(i)} cy={y(r.p50)} r={3} fill={OKABE.blue} />
          <text x={x(i)} y={height - bottom + 14} textAnchor="end" transform={`rotate(-35 ${x(i)} ${height - bottom + 14})`} className="fill-[var(--muted)] text-[10px]">
            {r.label}
          </text>
        </g>
      ))}
      <text transform={`translate(14 ${(top + height - bottom) / 2}) rotate(-90)`} textAnchor="middle" className="fill-[var(--muted)] text-[11px]">
        {yLabel}
      </text>
      <text x={(left + width - right) / 2} y={height - 4} textAnchor="middle" className="fill-[var(--muted)] text-[11px]">
        {xLabel}
      </text>
    </svg>
  );
}

// Canvas helper: draw once per change, sized to the container's CSS width.
export function useCanvas(draw: (ctx: CanvasRenderingContext2D, w: number, h: number) => void, deps: unknown[], height: number) {
  const ref = useRef<HTMLCanvasElement | null>(null);
  useEffect(() => {
    const canvas = ref.current;
    if (!canvas) return;
    const dpr = window.devicePixelRatio || 1;
    const w = canvas.clientWidth || 600;
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(height * dpr);
    const ctx = canvas.getContext("2d");
    if (!ctx) return;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, height);
    draw(ctx, w, height);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, deps);
  return ref;
}

export interface Pt {
  x: number; // genomic position (or cumulative)
  y: number; // -log10 P
  color?: string;
  size?: number;
}

export function ScatterCanvas({
  points,
  xDomain,
  yDomain,
  height = 300,
  ariaLabel,
  threshold,
  overlay,
  colorFor,
  xLabel,
  yLabel,
  summary,
  xTicks,
  xFormat,
}: {
  points: Pt[];
  xDomain: [number, number];
  yDomain: [number, number];
  height?: number;
  ariaLabel: string;
  threshold?: number;
  overlay?: Pt[];
  colorFor?: (p: Pt) => string;
  xLabel: string;
  yLabel: string;
  summary: string;
  xTicks?: { x: number; label: string }[];
  xFormat?: (v: number) => string;
}) {
  const ref = useCanvas(
    (ctx, w, h) => {
      const left = 58;
      const right = 10;
      const top = 8;
      const bottom = 34;
      const x = scaleLinear().domain(xDomain).range([left, w - right]);
      const y = scaleLinear().domain(yDomain).range([h - bottom, top]);
      ctx.font = "11px system-ui, sans-serif";
      ctx.strokeStyle = "rgba(128,128,128,0.25)";
      ctx.fillStyle = "#777";
      ctx.textAlign = "right";
      for (const t of y.ticks(5)) {
        ctx.beginPath();
        ctx.moveTo(left, y(t));
        ctx.lineTo(w - right, y(t));
        ctx.stroke();
        ctx.fillText(String(t), left - 6, y(t) + 4);
      }
      ctx.textAlign = "center";
      if (xTicks) {
        for (const t of xTicks) {
          if (t.x < xDomain[0] || t.x > xDomain[1]) continue;
          ctx.fillText(t.label, x(t.x), h - 18);
        }
      } else {
        for (const t of x.ticks(6)) {
          ctx.fillText(xFormat ? xFormat(t) : String(Math.round(t)), x(t), h - 18);
        }
      }
      ctx.fillText(xLabel, (left + w) / 2, h - 3);
      ctx.save();
      ctx.translate(14, (top + h - bottom) / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.textAlign = "center";
      ctx.fillText(yLabel, 0, 0);
      ctx.restore();
      ctx.textAlign = "left";
      if (threshold !== undefined) {
        ctx.strokeStyle = "#D55E00";
        ctx.setLineDash([5, 4]);
        ctx.beginPath();
        ctx.moveTo(left, y(threshold));
        ctx.lineTo(w - right, y(threshold));
        ctx.stroke();
        ctx.setLineDash([]);
      }
      const draw = (arr: Pt[], fallback: string) => {
        for (const p of arr) {
          if (p.x < xDomain[0] || p.x > xDomain[1]) continue;
          ctx.fillStyle = p.color ?? colorFor?.(p) ?? fallback;
          const yy = Math.min(p.y, yDomain[1]);
          ctx.beginPath();
          ctx.arc(x(p.x), y(yy), p.size ?? 1.6, 0, Math.PI * 2);
          ctx.fill();
        }
      };
      draw(points, "#0072B2");
      if (overlay) draw(overlay, "#E69F00");
    },
    [points, overlay, xDomain[0], xDomain[1], yDomain[0], yDomain[1], threshold, xTicks, xFormat],
    height,
  );
  return (
    <div>
      <canvas ref={ref} role="img" aria-label={ariaLabel} style={{ width: "100%", height }} />
      <p className="sr-only">{summary}</p>
    </div>
  );
}

export function Heatmap({
  rows,
  columns,
  values,
  ariaLabel,
  highlight,
  cellW = 26,
  cellH = 16,
}: {
  rows: { label: string; mutantIndex?: number }[];
  columns: string[];
  values: number[][]; // [row][col] LLR
  ariaLabel: string;
  highlight?: { row: number; col: number };
  cellW?: number;
  cellH?: number;
}) {
  const left = 60;
  const top = 26;
  const width = left + columns.length * cellW + 10;
  const height = top + rows.length * cellH + 10;
  const color = scaleSequential(interpolateBlues).domain([-10, 0]);
  const summary = highlight
    ? `Mutation highlighted at row ${rows[highlight.row]?.label}, column ${columns[highlight.col]}. Values are log-likelihood ratios; more negative means less tolerated.`
    : "";
  return (
    <svg viewBox={`0 0 ${width} ${height}`} role="img" aria-label={ariaLabel} className="w-full h-auto">
      <desc>{summary}</desc>
      {columns.map((c, j) => (
        <text key={c} x={left + j * cellW + cellW / 2} y={top - 8} textAnchor="middle" className="fill-[var(--muted)] text-[11px]">
          {c}
        </text>
      ))}
      {rows.map((r, i) => (
        <g key={r.label + i}>
          <text x={left - 8} y={top + i * cellH + cellH / 2 + 4} textAnchor="end" className="fill-[var(--muted)] text-[10px]">
            {r.label}
          </text>
          {columns.map((c, j) => {
            const v = values[i]?.[j] ?? 0;
            return (
              <rect
                key={c}
                x={left + j * cellW}
                y={top + i * cellH}
                width={cellW - 1}
                height={cellH - 1}
                fill={color(Math.max(-10, Math.min(0, v)))}
              >
                <title>{`${r.label} ${c}: ${v.toFixed(2)}`}</title>
              </rect>
            );
          })}
          {highlight && highlight.row === i ? (
            <rect
              x={left + highlight.col * cellW - 1}
              y={top + i * cellH - 1}
              width={cellW + 1}
              height={cellH + 1}
              fill="none"
              stroke="#D55E00"
              strokeWidth={2}
            />
          ) : null}
        </g>
      ))}
    </svg>
  );
}

export function PowerCurve({
  points,
  marker,
  ariaLabel,
  width = 600,
  height = 260,
}: {
  points: { n: number; power: number }[];
  marker?: { n: number; power: number; label: string };
  ariaLabel: string;
  width?: number;
  height?: number;
}) {
  const left = 48;
  const right = 14;
  const top = 12;
  const bottom = 36;
  const x = scaleLinear().domain([0, max(points, (p) => p.n) ?? 1]).range([left, width - right]);
  const y = scaleLinear().domain([0, 1]).range([height - bottom, top]);
  const path = line<{ n: number; power: number }>().x((d) => x(d.n)).y((d) => y(d.power))(points) ?? "";
  const summary = marker ? `${marker.label}: power ${marker.power.toFixed(3)} at N=${marker.n}.` : "";
  return (
    <svg viewBox={`0 0 ${width} ${height}`} role="img" aria-label={ariaLabel} className="w-full h-auto">
      <desc>{summary}</desc>
      {[0, 0.25, 0.5, 0.75, 1].map((t) => (
        <g key={t}>
          <line x1={left} x2={width - right} y1={y(t)} y2={y(t)} className="stroke-[var(--grid)]" />
          <text x={left - 8} y={y(t) + 4} textAnchor="end" className="fill-[var(--muted)] text-[11px]">
            {t.toFixed(2)}
          </text>
        </g>
      ))}
      <path d={path} fill="none" stroke={OKABE.blue} strokeWidth={2} />
      {marker ? (
        <g>
          <circle cx={x(marker.n)} cy={y(marker.power)} r={5} fill={OKABE.vermillion} />
          <text x={x(marker.n) + 8} y={y(marker.power) - 8} className="fill-[var(--fg)] text-[11px]">
            {marker.label}
          </text>
        </g>
      ) : null}
      <text x={(left + width - right) / 2} y={height - 6} textAnchor="middle" className="fill-[var(--muted)] text-[11px]">
        Sample size (animals)
      </text>
    </svg>
  );
}

export const viridis = interpolateViridis;
