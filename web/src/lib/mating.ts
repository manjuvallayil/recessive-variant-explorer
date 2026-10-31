// Mating-risk simulator math (Phase 5). Pure functions, no I/O.
//
// Model
// - Each locus is biallelic with a mutant allele. A parent transmits the mutant allele with
//   probability p: non-carrier 0, carrier 0.5, affected (homozygous) 1, unknown q (population
//   breed-weighted frequency). q = fHF * qHF + fJ * qJ (breed-weighted; crossbreds mixed linearly).
// - Calf homozygous probability P(hom) = pS * pD; P(het) = pS * (1 - pD) + (1 - pS) * pD.
// - Depletion toggle: P(hom) is multiplied by the observed/expected ratio from the purebred
//   depletion test (lost homozygotes are not renormalised; a simplification).
// - Loci are combined assuming independence, except DPF2/MUS81 which is one locus (DPF2 used as
//   proxy). DPF2-MUS81 LD is weak (r2 = 0.018 in 350 animals), so this is an approximation.
// - Penalty: P(hom) * homKg + P(het) * hetKg, with hetKg used only where useHetEffect is true
//   (DPF2/MUS81). Other het effects are not significant and are omitted.
// - Herd scenario (HWE mating within a population at a breed profile): genotype probabilities
//   P(0)=(1-q)^2, P(1)=2q(1-q), P(2)=q^2. Random mating: all sire-dam pairs allowed.
//   Avoid carrier x carrier: the (1,1) pairing is removed and the remaining pairs are renormalised.

export type Status = "non" | "carrier" | "affected" | "unknown";

export interface Breed {
  // Fractions of HF and J. Need not sum to 1 (remainder is treated as HF-equivalent via q_all).
  hf: number;
  j: number;
}

export interface LocusParams {
  id: string;
  label: string;
  qHF: number;
  qJ: number;
  depletionRatio: number; // observed/expected homozygotes (1 = no depletion)
  homKg: number; // homozygote effect (kg), negative = penalty
  hetKg: number; // heterozygote effect (kg)
  useHetEffect: boolean;
}

export interface ParentState {
  status: Record<string, Status>;
  breed: Breed;
}

export interface CalfLocusRisk {
  id: string;
  label: string;
  pHom: number;
  pHet: number;
  penaltyKg: number;
}

export interface CalfRisk {
  perLocus: CalfLocusRisk[];
  pAffectedAny: number; // P(homozygous for at least one locus), independence
  expectedPenaltyKg: number;
}

export function breedQ(locus: LocusParams, breed: Breed): number {
  const q = breed.hf * locus.qHF + breed.j * locus.qJ;
  return Math.min(1, Math.max(0, q));
}

export function transmit(status: Status, q: number): number {
  switch (status) {
    case "non":
      return 0;
    case "carrier":
      return 0.5;
    case "affected":
      return 1;
    case "unknown":
      return q;
  }
}

export function calfProbs(
  pSire: number,
  pDam: number,
  depletionRatio: number,
  applyDepletion: boolean,
): { pHom: number; pHet: number } {
  const hom0 = pSire * pDam;
  const het = pSire * (1 - pDam) + (1 - pSire) * pDam;
  const ratio = applyDepletion ? depletionRatio : 1;
  return { pHom: hom0 * ratio, pHet: het };
}

export function calfRisk(
  loci: LocusParams[],
  sire: ParentState,
  dam: ParentState,
  applyDepletion: boolean,
): CalfRisk {
  const perLocus: CalfLocusRisk[] = [];
  let noneAffected = 1;
  let expectedPenaltyKg = 0;
  for (const locus of loci) {
    const qS = breedQ(locus, sire.breed);
    const qD = breedQ(locus, dam.breed);
    const pS = transmit(sire.status[locus.id] ?? "unknown", qS);
    const pD = transmit(dam.status[locus.id] ?? "unknown", qD);
    const { pHom, pHet } = calfProbs(pS, pD, locus.depletionRatio, applyDepletion);
    const penaltyKg = pHom * locus.homKg + (locus.useHetEffect ? pHet * locus.hetKg : 0);
    perLocus.push({ id: locus.id, label: locus.label, pHom, pHet, penaltyKg });
    noneAffected *= 1 - pHom;
    expectedPenaltyKg += penaltyKg;
  }
  return { perLocus, pAffectedAny: 1 - noneAffected, expectedPenaltyKg };
}

// Genotype probabilities under HWE and the transmission probability for each genotype.
export function genotypeProbs(q: number): number[] {
  return [(1 - q) * (1 - q), 2 * q * (1 - q), q * q];
}

export const TRANSMIT_BY_GENOTYPE = [0, 0.5, 1];

// Herd rates for one locus at allele frequency q.
export function herdLocusRates(
  q: number,
  depletionRatio: number,
  avoidCarrierByCarrier: boolean,
  applyDepletion: boolean,
): { pHom: number; pHet: number } {
  const pg = genotypeProbs(q);
  let z = 0;
  let hom = 0;
  let het = 0;
  for (let gs = 0; gs < 3; gs += 1) {
    for (let gd = 0; gd < 3; gd += 1) {
      const allowed = !(avoidCarrierByCarrier && gs === 1 && gd === 1);
      if (!allowed) continue;
      const w = (pg[gs] ?? 0) * (pg[gd] ?? 0);
      const ps = TRANSMIT_BY_GENOTYPE[gs] ?? 0;
      const pd = TRANSMIT_BY_GENOTYPE[gd] ?? 0;
      z += w;
      hom += w * ps * pd;
      het += w * (ps * (1 - pd) + (1 - ps) * pd);
    }
  }
  const pHom0 = z > 0 ? hom / z : 0;
  const pHet = z > 0 ? het / z : 0;
  const ratio = applyDepletion ? depletionRatio : 1;
  return { pHom: pHom0 * ratio, pHet };
}

export interface HerdOutcome {
  affectedPerYear: number;
  // Change in kg per year (negative = loss). Homozygote-only, and net including het effects.
  homKgChangePerYear: number;
  netKgChangePerYear: number;
  perLocus: CalfLocusRisk[];
  pAffectedAny: number;
}

export function herdOutcome(
  loci: LocusParams[],
  population: Breed,
  calvesPerYear: number,
  avoidCarrierByCarrier: boolean,
  applyDepletion: boolean,
): HerdOutcome {
  const perLocus: CalfLocusRisk[] = [];
  let noneAffected = 1;
  let homPerCalf = 0;
  let netPerCalf = 0;
  for (const locus of loci) {
    const q = breedQ(locus, population);
    const { pHom, pHet } = herdLocusRates(
      q,
      locus.depletionRatio,
      avoidCarrierByCarrier,
      applyDepletion,
    );
    const homKg = pHom * locus.homKg;
    const hetKg = locus.useHetEffect ? pHet * locus.hetKg : 0;
    perLocus.push({ id: locus.id, label: locus.label, pHom, pHet, penaltyKg: homKg + hetKg });
    noneAffected *= 1 - pHom;
    homPerCalf += homKg;
    netPerCalf += homKg + hetKg;
  }
  const pAffectedAny = 1 - noneAffected;
  return {
    affectedPerYear: calvesPerYear * pAffectedAny,
    homKgChangePerYear: calvesPerYear * homPerCalf,
    netKgChangePerYear: calvesPerYear * netPerCalf,
    perLocus,
    pAffectedAny,
  };
}

// Closed form for the avoid carrier x carrier homozygous rate (random-mating HWE, no depletion).
// Used as a unit-test reference for herdLocusRates.
export function avoidCcHomClosedForm(q: number): number {
  const num = q * q - q * q * (1 - q) * (1 - q);
  const den = 1 - 4 * q * q * (1 - q) * (1 - q);
  return num / den;
}

export interface NationalInputs {
  femalesPerYear: number;
  meanQ: Record<string, number>; // locus id -> mean breed-weighted q over animals
}

// Calibration: random mating at the national mean q (unknown parents), no pedigree structure.
export function nationalCalibration(
  loci: LocusParams[],
  national: NationalInputs,
  applyDepletion: boolean,
): { id: string; label: string; perYear: number }[] {
  return loci.map((l) => {
    const q = national.meanQ[l.id] ?? 0;
    const { pHom } = calfProbs(q, q, l.depletionRatio, applyDepletion);
    return { id: l.id, label: l.label, perYear: national.femalesPerYear * pHom };
  });
}
