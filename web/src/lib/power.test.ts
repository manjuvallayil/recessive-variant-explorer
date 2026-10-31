import { test } from "node:test";
import assert from "node:assert/strict";
import { powerFor, gammaQ, chi2UpperCritical } from "./power.ts";

// Reference values: scipy.stats.ncx2.sf(chi2.ppf(1-5e-8, df), df, lambda)
// with lambda = 50000 * 0.001 / (1 - 0.001).
test("gammaQ matches closed forms", () => {
  // Q(1, x) = exp(-x)
  assert.ok(Math.abs(gammaQ(1, 1) - Math.exp(-1)) < 1e-12);
  assert.ok(Math.abs(gammaQ(1, 10) - Math.exp(-10)) < 1e-14);
  // Q(2, x) = (1 + x) exp(-x)
  assert.ok(Math.abs(gammaQ(2, 3) - 4 * Math.exp(-3)) < 1e-12);
});

test("central critical value at alpha 5e-8", () => {
  assert.ok(Math.abs(chi2UpperCritical(1, 5e-8) - 29.716785) < 1e-4);
  assert.ok(Math.abs(chi2UpperCritical(2, 5e-8) - 33.622485) < 1e-4);
});

test("powerFor recessive at paper variance matches scipy (0.94774)", () => {
  const p = powerFor(0.001, 50000, 1, 5e-8);
  assert.ok(Math.abs(p - 0.947737) < 1e-4, `got ${p}`);
});

test("powerFor genotype class at paper variance matches scipy (0.91219)", () => {
  const p = powerFor(0.001, 50000, 2, 5e-8);
  assert.ok(Math.abs(p - 0.912193) < 1e-4, `got ${p}`);
});

test("powerFor low power regime stays in range", () => {
  const p = powerFor(0.001, 500, 1, 5e-8);
  assert.ok(p >= 0 && p <= 1);
  assert.ok(p < 0.05);
});
