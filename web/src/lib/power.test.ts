import { test } from "node:test";
import assert from "node:assert/strict";
import { powerFor } from "./power.ts";

test("powerFor recessive at paper variance is ~0.948", () => {
  const p = powerFor(0.001, 50000, 1, 5e-8);
  assert.ok(Math.abs(p - 0.9477) < 0.005);
});

test("powerFor genotype class at paper variance is ~0.817", () => {
  const p = powerFor(0.001, 50000, 2, 5e-8);
  assert.ok(Math.abs(p - 0.8166) < 0.005);
});

test("powerFor low power regime stays in range", () => {
  const p = powerFor(0.001, 500, 1, 5e-8);
  assert.ok(p >= 0 && p <= 1);
});
