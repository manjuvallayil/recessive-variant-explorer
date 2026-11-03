const GENOTYPES = ["AA", "Aa", "Aa", "aa"];

function MiniPunnett() {
  return (
    <div aria-hidden className="grid grid-cols-2 gap-1 text-[11px]">
      {GENOTYPES.map((g, i) => (
        <span
          key={`${g}-${i}`}
          className={`flex h-7 w-8 items-center justify-center rounded-sm border ${
            g === "aa"
              ? "border-[var(--vermillion)] bg-[var(--vermillion)]/10 font-semibold"
              : "border-[var(--rule)] text-[var(--muted)]"
          }`}
        >
          {g}
        </span>
      ))}
    </div>
  );
}

export function SimulatorChip() {
  return (
    <a
      href="#simulator"
      style={{ textDecoration: "none" }}
      className="mt-4 inline-flex items-center gap-2 rounded border border-[var(--rule)] bg-[var(--surface)] px-3 py-2 text-xs font-semibold no-underline xl:hidden"
    >
      Try the mating-risk simulator <span aria-hidden>→</span>
    </a>
  );
}

export default function SimulatorTeaser() {
  return (
    <a
      href="#simulator"
      aria-label="Open the mating-risk simulator"
      style={{ textDecoration: "none" }}
      className="group block w-96 max-w-full rounded border border-[var(--rule)] bg-[var(--surface)] p-3 no-underline transition-colors hover:border-[var(--fg)]"
    >
      <p className="text-[10px] uppercase tracking-[0.14em] text-[var(--muted)]">Interactive · Instrument 09</p>
      <p className="serif mt-1 text-base font-semibold">Mating-risk simulator</p>
      <div className="mt-3 flex items-center gap-3">
        <MiniPunnett />
        <div>
          <p className="text-xs font-semibold">1 in 4 affected</p>
          <p className="text-[11px] text-[var(--muted)]">Aa × Aa example</p>
        </div>
      </div>
      <p className="mt-3 text-[11px] text-[var(--muted)]">Change sire, dam, locus and herd size.</p>
      <p className="mt-2 text-xs font-semibold">
        Open simulator{" "}
        <span aria-hidden className="inline-block transition-transform group-hover:translate-x-0.5">
          →
        </span>
      </p>
    </a>
  );
}
