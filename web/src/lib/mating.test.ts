import { test } from "node:test";
import assert from "node:assert/strict";
import {
  avoidCcHomClosedForm,
  breedQ,
  calfProbs,
  calfRisk,
  herdLocusRates,
  herdOutcome,
  nationalCalibration,
  transmit,
  type LocusParams,
  type ParentState,
} from "./mating.ts";

const first = <T,>(arr: T[]): T => {
  const v = arr[0];
  if (v === undefined) throw new Error("empty array");
  return v;
};
const close = (a: number, b: number, tol = 1e-12) =>
  assert.ok(Math.abs(a - b) <= tol, `${a} != ${b}`);

const PLCD4: LocusParams = {
  id: "PLCD4",
  label: "PLCD4",
  qHF: 0.04,
  qJ: 0.0006,
  depletionRatio: 0.5,
  homKg: -110,
  hetKg: 0,
  useHetEffect: false,
};
const DPF2: LocusParams = {
  id: "DPF2_MUS81",
  label: "DPF2/MUS81",
  qHF: 0.065,
  qJ: 0.0004,
  depletionRatio: 0.88,
  homKg: -56,
  hetKg: 2.3,
  useHetEffect: true,
};
const HF: ParentState = { status: {}, breed: { hf: 1, j: 0 } };
const withStatus = (s: Record<string, "non" | "carrier" | "affected" | "unknown">): ParentState => ({
  status: s,
  breed: { hf: 1, j: 0 },
});

test("transmission by status", () => {
  assert.equal(transmit("non", 0.3), 0);
  assert.equal(transmit("carrier", 0.3), 0.5);
  assert.equal(transmit("affected", 0.3), 1);
  assert.equal(transmit("unknown", 0.3), 0.3);
});

test("breed-weighted frequency is a linear mixture", () => {
  close(breedQ(PLCD4, { hf: 0.5, j: 0.5 }), 0.5 * 0.04 + 0.5 * 0.0006, 1e-15);
  close(breedQ(PLCD4, { hf: 1, j: 0 }), 0.04, 1e-15);
});

test("Mendelian outcomes for known parental status", () => {
  const cases: [string, string, number][] = [
    ["non", "non", 0],
    ["carrier", "non", 0],
    ["carrier", "carrier", 0.25],
    ["affected", "non", 0],
    ["affected", "carrier", 0.5],
    ["affected", "affected", 1],
  ];
  for (const [s, d, expected] of cases) {
    const risk = calfRisk(
      [PLCD4],
      withStatus({ PLCD4: s as "non" }),
      withStatus({ PLCD4: d as "non" }),
      false,
    );
    close(first(risk.perLocus).pHom, expected, 1e-15);
  }
});

test("unknown x unknown gives q squared (HWE random mating)", () => {
  const risk = calfRisk([PLCD4], HF, HF, false);
  close(first(risk.perLocus).pHom, 0.04 * 0.04, 1e-15);
});

test("depletion scales homozygote probability only", () => {
  const off = calfRisk([PLCD4], withStatus({ PLCD4: "carrier" }), withStatus({ PLCD4: "carrier" }), false);
  const on = calfRisk([PLCD4], withStatus({ PLCD4: "carrier" }), withStatus({ PLCD4: "carrier" }), true);
  close(first(on.perLocus).pHom, first(off.perLocus).pHom * PLCD4.depletionRatio, 1e-15);
  close(first(on.perLocus).pHet, first(off.perLocus).pHet, 1e-15);
});

test("calf probabilities are a valid split", () => {
  const { pHom, pHet } = calfProbs(0.3, 0.6, 1, false);
  assert.ok(pHom >= 0 && pHet >= 0 && pHom + pHet <= 1 + 1e-12);
  close(pHet, 0.3 * 0.4 + 0.7 * 0.6, 1e-15);
});

test("combined affected probability assumes independence across loci", () => {
  const loci = [PLCD4, { ...PLCD4, id: "X", label: "X", qHF: 0.02, depletionRatio: 1 }];
  const risk = calfRisk(loci, withStatus({ PLCD4: "affected", X: "carrier" }),
    withStatus({ PLCD4: "affected", X: "carrier" }), false);
  // PLCD4: both affected -> hom = 1, so P(any) = 1
  close(risk.pAffectedAny, 1, 1e-15);
});

test("DPF2/MUS81 contributes its het effect; PLCD4 does not", () => {
  const risk = calfRisk(
    [PLCD4, DPF2],
    withStatus({ PLCD4: "carrier", DPF2_MUS81: "carrier" }),
    withStatus({ PLCD4: "non", DPF2_MUS81: "non" }),
    false,
  );
  const dpf = risk.perLocus.find((l) => l.id === "DPF2_MUS81")!;
  const plcd = risk.perLocus.find((l) => l.id === "PLCD4")!;
  close(dpf.penaltyKg, dpf.pHom * DPF2.homKg + dpf.pHet * DPF2.hetKg, 1e-12);
  close(plcd.penaltyKg, plcd.pHom * PLCD4.homKg, 1e-12);
});

test("random-mating herd rates match HWE: P(hom)=q^2, P(het)=2q(1-q)", () => {
  const q = 0.07;
  const r = herdLocusRates(q, 1, false, false);
  close(r.pHom, q * q, 1e-14);
  close(r.pHet, 2 * q * (1 - q), 1e-14);
});

test("avoid carrier x carrier matches the closed form", () => {
  for (const q of [0.01, 0.04, 0.1, 0.3]) {
    const r = herdLocusRates(q, 1, true, false);
    close(r.pHom, avoidCcHomClosedForm(q), 1e-12);
  }
});

test("avoiding carrier x carrier reduces affected births", () => {
  for (const q of [0.01, 0.04, 0.1, 0.3]) {
    const random = herdLocusRates(q, 1, false, false).pHom;
    const avoid = herdLocusRates(q, 1, true, false).pHom;
    assert.ok(avoid < random, `q=${q}`);
  }
});

test("herd outcome: net kg = calves x expected penalty; homozygote kg excludes het effects", () => {
  const pop = { hf: 1, j: 0 };
  const out = herdOutcome([PLCD4, DPF2], pop, 1000, false, false);
  const sumPenalty = out.perLocus.reduce((s, l) => s + l.penaltyKg, 0);
  close(out.netKgChangePerYear, 1000 * sumPenalty, 1e-9);
  const homOnly = out.perLocus.reduce((s, l) => s + l.pHom * (l.id === "PLCD4" ? PLCD4.homKg : DPF2.homKg), 0);
  close(out.homKgChangePerYear, 1000 * homOnly, 1e-9);
  close(out.affectedPerYear, 1000 * out.pAffectedAny, 1e-9);
});

test("net change can be positive when DPF2 carrier gain outweighs rare homozygotes", () => {
  const pop = { hf: 1, j: 0 };
  const out = herdOutcome([DPF2], pop, 1000, true, false);
  assert.ok(out.homKgChangePerYear < 0);
});

test("depletion toggle lowers herd homozygote rate by the ratio", () => {
  const pop = { hf: 1, j: 0 };
  const off = herdOutcome([PLCD4], pop, 1000, false, false);
  const on = herdOutcome([PLCD4], pop, 1000, false, true);
  close(first(on.perLocus).pHom, first(off.perLocus).pHom * PLCD4.depletionRatio, 1e-15);
});

test("national calibration uses random mating at the national mean q", () => {
  const rows = nationalCalibration([PLCD4], { femalesPerYear: 940100, meanQ: { PLCD4: 0.02 } }, false);
  close(first(rows).perYear, 940100 * 0.02 * 0.02, 1e-9);
});
