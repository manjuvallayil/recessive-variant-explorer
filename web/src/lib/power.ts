export function erf(z: number): number {
  const t = 1 / (1 + 0.3275911 * Math.abs(z));
  const t2 = t * t;
  const t3 = t * t2;
  const t4 = t * t3;
  const t5 = t * t4;
  const poly =
    0.254829592 * t +
    (-0.284496736) * t2 +
    1.421413741 * t3 +
    (-1.453152027) * t4 +
    1.061405429 * t5;
  const y = 1 - poly * Math.exp(-(z * z));
  return z >= 0 ? y : -y;
}

export function normalCdf(z: number): number {
  return 0.5 * (1 + erf(z / Math.SQRT2));
}

export function normalQuantile(p: number): number {
  let lo = -10;
  let hi = 10;
  for (let i = 0; i < 200; i += 1) {
    const mid = (lo + hi) / 2;
    if (normalCdf(mid) < p) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return (lo + hi) / 2;
}

export function chi2Cdf(x: number, df: number): number {
  if (df === 1) {
    const root = Math.sqrt(x);
    return 2 * normalCdf(root) - 1;
  }
  if (df === 2) {
    return 1 - Math.exp(-x / 2);
  }
  return NaN;
}

export function chi2Quantile(df: number, p: number): number {
  let lo = 0;
  let hi = 1000;
  for (let i = 0; i < 200; i += 1) {
    const mid = (lo + hi) / 2;
    if (chi2Cdf(mid, df) < p) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return (lo + hi) / 2;
}

export function powerFor(
  varianceExplained: number,
  sampleSize: number,
  df = 1,
  alpha = 5e-8,
): number {
  const lam = (sampleSize * varianceExplained) / (1 - varianceExplained);
  const crit = chi2Quantile(df, 1 - alpha);
  const rootLam = Math.sqrt(lam);
  const rootX = Math.sqrt(crit);
  const sd = Math.sqrt(df);
  const a = (rootX - rootLam) / sd;
  const b = -((rootX + rootLam) / sd);
  return 1 - normalCdf(a) + normalCdf(b);
}
