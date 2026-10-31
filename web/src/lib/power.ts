// Exact noncentral chi-square power for the recessive / genotype-class tests.
// Mirrors pipeline/phase4_power.py::power_for (scipy.stats.ncx2.sf vs chi2.ppf).
// Only erasable TypeScript syntax, so `node --test` can run this file directly.

const LANCZOS_G = 7;
const LANCZOS_COEF = [
  0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
  -176.61503916999385, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
  1.5056327351493116e-7,
];

export function logGamma(z: number): number {
  // Exact path for positive integers and half-integers (all shapes used here: df/2 + j).
  if (z > 0 && Number.isInteger(2 * z)) {
    if (Number.isInteger(z)) {
      let s = 0;
      for (let k = 2; k < z; k += 1) s += Math.log(k);
      return s; // log((z-1)!)
    }
    // Gamma(n + 1/2) = sqrt(pi) * prod_{k=1..n} (k - 1/2)
    let s = 0.5 * Math.log(Math.PI);
    for (let k = 1; k < z + 0.5; k += 1) s += Math.log(k - 0.5);
    return s;
  }
  if (z < 0.5) {
    return Math.log(Math.PI / Math.abs(Math.sin(Math.PI * z))) - logGamma(1 - z);
  }
  const x = z - 1;
  let a = LANCZOS_COEF[0] ?? 0;
  const t = x + LANCZOS_G + 0.5;
  for (let i = 1; i < LANCZOS_COEF.length; i += 1) {
    a += (LANCZOS_COEF[i] ?? 0) / (x + i);
  }
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
}

// Regularised upper incomplete gamma Q(a, x) = Gamma(a, x) / Gamma(a).
// Series for P when x < a + 1, Lentz continued fraction for Q otherwise
// (Numerical Recipes gammq). Accurate in the far upper tail, which matters at alpha = 5e-8.
export function gammaQ(a: number, x: number): number {
  if (x <= 0) return 1;
  const lnPrefactor = -x + a * Math.log(x) - logGamma(a);
  if (x < a + 1) {
    let sum = 1 / a;
    let term = sum;
    for (let n = 1; n < 10_000; n += 1) {
      term *= x / (a + n);
      sum += term;
      if (Math.abs(term) < Math.abs(sum) * 1e-16) break;
    }
    const p = sum * Math.exp(lnPrefactor);
    return 1 - p;
  }
  const tiny = 1e-300;
  let b = x + 1 - a;
  let c = 1 / tiny;
  let d = 1 / b;
  let h = d;
  for (let i = 1; i < 10_000; i += 1) {
    const an = -i * (i - a);
    b += 2;
    d = an * d + b;
    if (Math.abs(d) < tiny) d = tiny;
    c = b + an / c;
    if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d;
    const delta = d * c;
    h *= delta;
    if (Math.abs(delta - 1) < 1e-16) break;
  }
  return Math.exp(lnPrefactor) * h;
}

// Central chi-square survival function P(X > x), X ~ chi2(df).
export function chi2Sf(x: number, df: number): number {
  return gammaQ(df / 2, x / 2);
}

// Upper critical value: smallest x with P(X > x) <= alpha, by bisection on chi2Sf.
export function chi2UpperCritical(df: number, alpha: number): number {
  let lo = 0;
  let hi = 1000;
  for (let i = 0; i < 300; i += 1) {
    const mid = (lo + hi) / 2;
    if (chi2Sf(mid, df) > alpha) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return (lo + hi) / 2;
}

// Noncentral chi-square survival P(X > x), X ~ chi2(df, lambda), via the Poisson mixture
// sum_j Pois(j; lambda/2) * Q(df/2 + j, x/2). Weights are computed in log space.
export function ncChi2Sf(x: number, df: number, lambda: number): number {
  const mu = lambda / 2;
  const jMax = Math.ceil(mu + 12 * Math.sqrt(mu + 1) + 50);
  let total = 0;
  for (let j = 0; j <= jMax; j += 1) {
    const logW = -mu + (mu > 0 ? j * Math.log(mu) : j === 0 ? 0 : -Infinity) - logGamma(j + 1);
    if (logW === -Infinity) continue;
    total += Math.exp(logW) * gammaQ(df / 2 + j, x / 2);
  }
  return Math.min(1, Math.max(0, total));
}

export function powerFor(
  varianceExplained: number,
  sampleSize: number,
  df = 1,
  alpha = 5e-8,
): number {
  const lambda = (sampleSize * varianceExplained) / (1 - varianceExplained);
  const crit = chi2UpperCritical(df, alpha);
  return ncChi2Sf(crit, df, lambda);
}
