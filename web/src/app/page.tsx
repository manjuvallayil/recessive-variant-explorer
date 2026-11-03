"use client";

import Link from "next/link";
import { useEffect, useRef, useState, useSyncExternalStore, type ReactNode } from "react";
import HookPanel from "@/components/panels/HookPanel";
import ProxyPanel from "@/components/panels/ProxyPanel";
import ManhattanPanel from "@/components/panels/ManhattanPanel";
import ZoomPanel from "@/components/panels/ZoomPanel";
import EvidencePanel from "@/components/panels/EvidencePanel";
import CommonPanel from "@/components/panels/CommonPanel";
import SignalPanel from "@/components/panels/SignalPanel";
import PowerPanel from "@/components/panels/PowerPanel";
import SimulatorPanel from "@/components/panels/SimulatorPanel";
import SimulatorTeaser, { SimulatorChip } from "@/components/SimulatorTeaser";

type Step = {
  id: string;
  title: string;
  body: ReactNode;
  Panel: (props: { active: boolean }) => ReactNode;
  panelTitle: string;
};

const STEPS: Step[] = [
  {
    id: "hook",
    title: "Rare alleles, hidden in plain sight",
    panelTitle: "Carrier frequencies by breed",
    Panel: HookPanel,
    body: (
      <>
        <p>
          The paper&apos;s abstract opens with a simple fact: mammalian species carry roughly 100 loss-of-function variants per
          individual. About one to five of them hit essential genes, and when homozygous they cause embryonic lethality or
          serious disease. The rest are harder to read, and are assumed to affect fitness in less obvious ways.
        </p>
        <p>
          The authors screened about 130,000 cattle and found eight body-weight QTL; seven were not detected by a standard
          additive scan. Six of the QTL were novel. Five of those six carry compelling candidate mutations, with two
          candidates at chr29:44 (DPF2 and MUS81), for six candidate mutations in total: missense (PLCD4, MTRF1, DPF2), a
          premature stop (MUS81), or a splice change (GALNT2, FGD4). The sixth novel locus (chr2:22) has no candidate. A
          separate BCS-suggestive frameshift in MYH1 is shown in the panel because the paper discusses it with these loci.
        </p>
        <p>
          The panel shows how common the carrier state is. A carrier is a healthy animal with one copy of the mutant allele; a
          mating of two carriers is how a recessive condition appears.
        </p>
      </>
    ),
  },
  {
    id: "proxy",
    title: "Routine data as a proxy phenotype",
    panelTitle: "Body weight by genotype class",
    Panel: ProxyPanel,
    body: (
      <>
        <p>
          Breeders already record body weight, stature and condition on large numbers of animals. The paper&apos;s idea was to treat
          weight as a proxy for a wide range of null alleles, which is how mouse studies detect them.
        </p>
        <p>
          Homozygous mutants weigh far less. For PLCD4, the homozygous class averages about 111 kg below reference animals after
          the population-structure adjustment (paper: 109 kg). GALNT2 is similar (about 113 kg; paper: 112 kg). The panel lets you
          switch loci.
        </p>
        <p>
          The lower body-weight classes are small, so the intervals are wide. Read the box plots as distributions, not as
          individual predictions.
        </p>
      </>
    ),
  },
  {
    id: "manhattan",
    title: "Additive models miss recessive QTL",
    panelTitle: "Genome-wide scan: dominance vs additive",
    Panel: ManhattanPanel,
    body: (
      <>
        <p>
          A standard additive GWAS looks for an allele that shifts the trait in proportion to the number of copies. A recessive
          allele only shows up when two copies are present, so an additive scan misses it.
        </p>
        <p>
          Fitting a dominance term recovers the recessive loci. The panel lists eight body-weight QTL, with PLAG1 included as the known
          additive comparator. At the candidate or lead variant, seven of the eight fall below genome-wide significance in the additive
          model; PLAG1 is the exception. Using the best additive SNP within ±1 Mb gives six of eight; the panel identifies the
          additional locus.
        </p>
        <p>
          The chip toggle overlays a 50k SNP array on the sequence-level scan. Of the major recessive syndromes, only DPF2/MUS81 stays
          genome-wide significant on the chip.
        </p>
      </>
    ),
  },
  {
    id: "zoom",
    title: "Zooming in: regional LD and candidate genes",
    panelTitle: "Regional association and LD",
    Panel: ZoomPanel,
    body: (
      <>
        <p>
          Each regional plot shows the dominance association in a window around the candidate, coloured by linkage disequilibrium
          (r²) with the candidate. Points that share the candidate&apos;s r² pattern are likely to be in the same haplotype.
        </p>
        <p>
          The gene track is a best-effort annotation. The coordinates are lifted from the ARS-UCD1.2 assembly to UMD3.1 by
          interpolation, so gene boundaries can be off by a few kilobases. The PLAG1 window has no LD data in this build, and is
          shown in grey.
        </p>
      </>
    ),
  },
  {
    id: "evidence",
    title: "Evidence, and what did not hold up",
    panelTitle: "Protein model, RNA and farm profiles",
    Panel: EvidencePanel,
    body: (
      <>
        <p>
          Three kinds of evidence are shown: a protein language model scoring each missense change, RNA splice and expression
          summaries, and small research-farm cohorts.
        </p>
        <p>
          The protein-model ranking runs opposite to the observed body-weight effects across three candidates. I report that as a
          failed check, not as support. The farm cohorts are a few animals per group, so they are descriptive only.
        </p>
      </>
    ),
  },
  {
    id: "common",
    title: "How common are these alleles?",
    panelTitle: "National estimates and depletion",
    Panel: CommonPanel,
    body: (
      <>
        <p>
          The national export uses the pedigree of about 2.8 million animals and the genotypes of their sires. It estimates that
          about 1.02% of animals are homozygous for at least one of five body-weight-significant variants (DPF2/MUS81 counted
          once); that is 9,548 per year in this export, close to the paper&apos;s rounded 1% and about 9,700 per year.
        </p>
        <p>
          Homozygotes are rarer than expected under random mating at most loci. That is consistent with selection against
          homozygotes, and the test is shown in the panel with its limits.
        </p>
      </>
    ),
  },
  {
    id: "signal",
    title: "Can data science find them?",
    panelTitle: "Syndrome signal, imputation and detection",
    Panel: SignalPanel,
    body: (
      <>
        <p>
          Could a routine record flag a homozygous mutant without genotyping? Body weight alone gives an AUROC of about 0.68 for the
          pooled label. A logistic model on three traits, with out-of-fold predictions, gives about 0.69. Neither is strong.
        </p>
        <p>
          Single-locus results depend on the effect size relative to the spread of body weight. PLCD4 and GALNT2 homozygotes are
          clearly lighter. In the site&apos;s combined DPF2/MUS81 homozygote group, the shift is smaller (about 0.5 SD); the paper
          reports per-mutation effects of about −53 kg. I explain the expected ceiling in the panel.
        </p>
        <p>
          Imputation accuracy falls with minor allele frequency, and the simulated-QTL scan shows which recessive effects a
          sequence-level additive test would miss.
        </p>
      </>
    ),
  },
  {
    id: "power",
    title: "Power to detect recessive effects",
    panelTitle: "Power explorer",
    Panel: PowerPanel,
    body: (
      <>
        <p>
          How many animals are needed to see a recessive locus that explains 0.1% of the variance at 2.5% minor allele
          frequency? The paper text rounds the answer to about 90% at about 50,000 animals; its Extended Data Fig. 1 reference
          point is about 94.7%. The analytic calculation here gives 94.8% for the recessive test, matching the figure.
        </p>
        <p>Move the sliders to see how power falls with a smaller locus, a lower sample size, or a two-degree-of-freedom test.</p>
      </>
    ),
  },
  {
    id: "simulator",
    title: "Try it: a mating-risk simulator",
    panelTitle: "Mating-risk simulator",
    Panel: SimulatorPanel,
    body: (
      <>
        <p>
          Set each parent&apos;s status for each locus. The simulator returns the chance of an affected calf, the expected body-weight
          change, and, for a herd, the difference between random mating and avoiding carrier-by-carrier matings.
        </p>
        <p>
          It is a teaching tool built from the paper&apos;s frequencies and the site&apos;s re-estimated effects, and it is not a
          veterinary or breeding decision tool. Calibration against the national export is shown in the panel.
        </p>
      </>
    ),
  },
  {
    id: "limits",
    title: "Limitations, and what I would do with LIC-scale data",
    panelTitle: "Limitations",
    Panel: LimitsPanel,
    body: (
      <>
        <p>
          <strong>Not done yet.</strong> The unsupervised autoencoder, isolation forest, LightGBM and nested cross-validation are not in
          the site. Sire-grouped cross-validation is not possible in this dataset.
        </p>
        <p>
          <strong>Dropped.</strong> ARSI, after a residue mismatch with the RefSeq record. MUS81 and splice variants are not scorable by the
          protein model.
        </p>
        <p>
          <strong>Single-locus effects are estimates.</strong> Homozygote classes are small (21 to 110 animals), and the population
          adjustment is locus-specific.
        </p>
        <p>
          <strong>Next with LIC-scale data:</strong> new phenotypes from routine and sensor data, pedigree-aware non-additive models,
          imputation against a current reference panel, and a prospective check of the syndrome signal in a held-out cohort.
        </p>
      </>
    ),
  },
];

function LimitsPanel() {
  return (
    <ul className="list-disc space-y-2 pl-5 text-sm">
      <li>Body-weight syndrome signal: AUROC 0.68 pooled; per-locus CIs are wide (FGD4 0.54–0.70; DPF2/MUS81 0.58–0.68) and exclude 0.5.</li>
      <li>Protein-model concordance with observed effects failed (n = 3).</li>
      <li>Farm cohorts: 6 to 9 animals per group, descriptive only.</li>
      <li>Simulator: independence across loci; DPF2 stands in for DPF2/MUS81; purebred depletion ratios applied to crossbreds.</li>
      <li>Not a breeding or veterinary tool.</li>
    </ul>
  );
}

export default function Page() {
  const [active, setActive] = useState<string>(STEPS[0]?.id ?? "hook");
  const sectionRefs = useRef<Record<string, HTMLElement | null>>({});

  const [visited, setVisited] = useState<Set<string>>(() => new Set([STEPS[0]?.id ?? "hook"]));

  useEffect(() => {
    const els = STEPS.map((s) => sectionRefs.current[s.id]).filter((e): e is HTMLElement => !!e);
    const io = new IntersectionObserver(
      (entries) => {
        for (const e of entries) {
          if (!e.isIntersecting) continue;
          const id = (e.target as HTMLElement).dataset.step;
          if (!id) continue;
          setActive(id);
          setVisited((v) => (v.has(id) ? v : new Set(v).add(id)));
        }
      },
      { rootMargin: "-45% 0px -45% 0px", threshold: 0 },
    );
    els.forEach((el) => io.observe(el));
    return () => io.disconnect();
  }, []);

  const current: Step = STEPS.find((s) => s.id === active) ?? (STEPS[0] as Step);

  return (
    <div className="min-h-screen">
      <header className="border-b border-[var(--rule)] px-4 py-6 sm:px-8">
        <div className="mx-auto grid max-w-7xl grid-cols-[minmax(0,1fr)] items-start gap-4 md:grid-cols-[minmax(0,1fr)_auto] xl:grid-cols-[minmax(0,1fr)_auto_auto]">
          <div className="min-w-0">
            <p className="text-xs uppercase tracking-wide text-[var(--muted)]">Independent project · not affiliated with LIC</p>
            <h1 className="serif mt-1 text-3xl font-semibold sm:text-4xl">Recessive Variant Explorer</h1>
            <p className="mt-2 max-w-2xl text-sm text-[var(--muted)]">
              A walkthrough of Reynolds et al. (2021), <em>Nature Genetics</em>: non-additive associations in cattle, reproduced
              from the authors&apos; public data.
            </p>
            <SimulatorChip />
          </div>
          <nav aria-label="Story steps" className="text-xs md:justify-self-end">
            <ol className="grid grid-cols-1 gap-x-4 gap-y-1 sm:grid-cols-2">
              {STEPS.map((s, i) => (
                <li key={s.id}>
                  <a
                    href={`#${s.id}`}
                    className={`inline-flex min-h-6 items-center py-1 ${active === s.id ? "font-semibold text-[var(--fg)]" : "text-[var(--muted)]"}`}
                    aria-current={active === s.id ? "step" : undefined}
                  >
                    {i + 1}. {s.title.split(":")[0]}
                  </a>
                </li>
              ))}
            </ol>
          </nav>
          <div className="hidden xl:block xl:justify-self-end">
            <SimulatorTeaser />
          </div>
        </div>
      </header>

      <div className="mx-auto grid max-w-7xl lg:grid-cols-[minmax(0,1fr)_minmax(0,1.05fr)]">
        <main className="px-4 sm:px-8">
          <div className="sticky top-0 z-20 -mx-4 h-[45vh] overflow-y-auto border-b border-[var(--rule)] bg-[var(--bg)] px-4 py-3 sm:-mx-8 sm:px-8 lg:hidden">
            <PanelFrame title={current.panelTitle} id={current.id}>
              {STEPS.filter((s) => visited.has(s.id)).map((s) => (
                <div key={s.id} hidden={s.id !== active}>
                  <s.Panel active={s.id === active} />
                </div>
              ))}
            </PanelFrame>
          </div>

          {STEPS.map((s) => (
            <section
              key={s.id}
              id={s.id}
              data-step={s.id}
              ref={(el) => {
                sectionRefs.current[s.id] = el;
              }}
              className="min-h-[70vh] border-b border-[var(--rule)] py-12 last:border-b-0 lg:min-h-[100vh]"
              aria-labelledby={`${s.id}-title`}
            >
              <h2 id={`${s.id}-title`} className="serif text-2xl font-semibold sm:text-3xl">
                {s.title}
              </h2>
              <div className="mt-4 space-y-4 text-[17px] leading-relaxed">{s.body}</div>
              <div className="mt-6 lg:hidden">
                <p className="text-xs text-[var(--muted)]">Panel above updates as you scroll.</p>
              </div>
            </section>
          ))}

          <Footer />
        </main>

        <aside className="hidden lg:block">
          <div className="sticky top-0 h-screen overflow-y-auto border-l border-[var(--rule)] p-6">
            <PanelFrame title={current.panelTitle} id={current.id}>
              {STEPS.filter((s) => visited.has(s.id)).map((s) => (
                <div key={s.id} hidden={s.id !== active}>
                  <s.Panel active={s.id === active} />
                </div>
              ))}
            </PanelFrame>
            <div className="mt-6">
              <ThemeToggle />
            </div>
          </div>
        </aside>
      </div>
    </div>
  );
}

function PanelFrame({ title, id, children }: { title: string; id: string; children: ReactNode }) {
  return (
    <div role="region" aria-label={`Interactive panel: ${title}`} data-panel={id}>
      <h2 className="serif mb-3 text-lg font-semibold">{title}</h2>
      {children}
    </div>
  );
}

const THEME_EVENT = "theme-change";

function subscribeTheme(cb: () => void) {
  window.addEventListener(THEME_EVENT, cb);
  return () => window.removeEventListener(THEME_EVENT, cb);
}

function readDark(): boolean {
  const stored = localStorage.getItem("theme");
  return stored === "dark" || (stored === null && window.matchMedia("(prefers-color-scheme: dark)").matches);
}

function ThemeToggle() {
  const dark = useSyncExternalStore(subscribeTheme, readDark, () => false);
  return (
    <button
      type="button"
      onClick={() => {
        const next = !dark;
        localStorage.setItem("theme", next ? "dark" : "light");
        document.documentElement.setAttribute("data-theme", next ? "dark" : "light");
        window.dispatchEvent(new Event(THEME_EVENT));
      }}
      className="rounded border border-[var(--rule)] px-3 py-1 text-xs"
      aria-pressed={dark}
    >
      {dark ? "Light theme" : "Dark theme"}
    </button>
  );
}

function Footer() {
  return (
    <footer className="mt-12 space-y-3 border-t border-[var(--rule)] py-8 text-sm text-[var(--muted)]">
      <p>
        <strong className="text-[var(--fg)]">Paper:</strong> Reynolds EGM, Neeley C, Lopdell TJ, et al. (2021). Non-additive
        association analysis using proxy phenotypes identifies novel cattle syndromes. <em>Nature Genetics</em> 53, 949–954.{" "}
        <a href="https://doi.org/10.1038/s41588-021-00872-5">doi:10.1038/s41588-021-00872-5</a>
      </p>
      <p>
        <strong className="text-[var(--fg)]">Data:</strong> Reynolds E, Littlejohn M (2021). Data from: Non-additive association
        analysis using proxy phenotypes identifies novel cattle syndromes. Dryad.{" "}
        <a href="https://doi.org/10.5061/dryad.kwh70rz1z">doi:10.5061/dryad.kwh70rz1z</a>. The dataset licence is not stated in the
        record text I checked; verify it on the Dryad record before reuse.
      </p>
      <p>
        This is an independent learning project, not affiliated with LIC. The site contains only derived aggregates. No animal
        or sire identifiers are published, and no paper figures or images are reproduced. LIC holds patent applications on these
        mutations, as noted in the project brief; nothing here is a commercial claim.
      </p>
      <p>
        <Link href="/methods/">Methods and data provenance</Link>
      </p>
    </footer>
  );
}
