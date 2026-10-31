"use client";

import type { ReactNode } from "react";

export function PanelState({ loading, error, children }: { loading: boolean; error: string | null; children: ReactNode }) {
  if (error) return <p role="alert" className="text-sm text-[var(--vermillion)]">Could not load data: {error}</p>;
  if (loading) return <p className="text-sm text-[var(--muted)]" aria-live="polite">Loading data for this section…</p>;
  return <>{children}</>;
}

export function Caption({ children }: { children: ReactNode }) {
  return <p className="mt-2 text-xs leading-relaxed text-[var(--muted)]">{children}</p>;
}

export function Source({ children }: { children: ReactNode }) {
  return <p className="mt-1 text-[11px] text-[var(--muted)]">Source: {children}</p>;
}

export const pct = (v: number, d = 1) => `${(v * 100).toFixed(d)}%`;
export const num = (v: number, d = 1) => v.toFixed(d);
export const int = (v: number) => Math.round(v).toLocaleString("en-NZ");

export const LOCUS_LABELS: Record<string, string> = {
  "Chr2:22": "Chr2:22 Mb QTL",
  PLCD4: "PLCD4",
  FGD4: "FGD4",
  ARSI: "ARSI",
  MTRF1: "MTRF1",
  MYH1: "MYH1",
  GALNT2: "GALNT2",
  DPF2: "DPF2",
  MUS81: "MUS81",
};
