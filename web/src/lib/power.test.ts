import { test } from "node:test";
import assert from "node:assert/strict";
import { powerFor } from "./power.ts";

test("powerFor recessive at paper parameters is ~1", () => {
  const p = powerFor(2.2645540682891916, 50000, 1, 5e-8);
  assert.ok(p > 0.99);
});

test("powerFor genotype class at paper parameters is ~1", () => {
  const p = powerFor(2.2645540682891916, 50000, 2, 5e-8);
  assert.ok(p > 0.99);
});

test("powerFor low power regime", () => {
  const p = powerFor(0.5, 1000, 1, 5e-8);
  assert.ok(p >= 0 && p <= 1);
});
