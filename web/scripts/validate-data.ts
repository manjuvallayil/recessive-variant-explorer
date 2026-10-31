// Build-time validation of derived data (run by `npm run validate-data`, which `build` and `dev` call).
// 1. Every consumed export parses against its Zod schema.
// 2. No individual-identifier keys appear anywhere in any export (privacy rule).
// 3. Reports raw and gzipped sizes; warns above 500 KB gzipped (lazy-loaded sections).
import { readFileSync, readdirSync, statSync } from "node:fs";
import { join, dirname } from "node:path";
import { fileURLToPath } from "node:url";
import { gzipSync } from "node:zlib";
import { DATA_FILES, REGIONAL_SCHEMA, RegionalIndexSchema } from "../src/lib/schemas.ts";

const here = dirname(fileURLToPath(import.meta.url));
const DATA = join(here, "..", "public", "data");
const FORBIDDEN_KEYS = new Set(["ID", "IID", "iid", "anml_key", "SIRE_ID", "MGS_ID", "sire_id", "dam_id"]);

let failures = 0;
const warnings: string[] = [];

function walkKeys(node: unknown, path: string, hits: string[]): void {
  if (Array.isArray(node)) {
    node.forEach((v, i) => walkKeys(v, `${path}[${i}]`, hits));
  } else if (node && typeof node === "object") {
    for (const [k, v] of Object.entries(node)) {
      if (FORBIDDEN_KEYS.has(k)) hits.push(`${path}.${k}`);
      walkKeys(v, `${path}.${k}`, hits);
    }
  }
}

function checkFile(rel: string, schema: { safeParse: (x: unknown) => { success: boolean; error?: unknown } }): void {
  const abs = join(DATA, rel);
  const raw = readFileSync(abs);
  const json = JSON.parse(raw.toString("utf8")) as unknown;
  const parsed = schema.safeParse(json);
  if (!parsed.success) {
    failures += 1;
    console.error(`FAIL schema: ${rel}`);
    console.error(String(parsed.error).slice(0, 800));
    return;
  }
  const hits: string[] = [];
  walkKeys(json, "$", hits);
  if (hits.length > 0) {
    failures += 1;
    console.error(`FAIL privacy: ${rel} has individual-identifier keys: ${hits.slice(0, 5).join(", ")}`);
    return;
  }
  const gz = gzipSync(raw).length;
  const kb = (gz / 1024).toFixed(1);
  if (gz > 500 * 1024) warnings.push(`${rel}: ${kb} KB gzipped (> 500 KB target)`);
  console.log(`ok  ${rel.padEnd(28)} ${String(raw.length).padStart(9)} B  ${kb.padStart(8)} KB gz`);
}

for (const { file, schema } of DATA_FILES) {
  if (file === "regional/index.json") {
    checkFile(file, RegionalIndexSchema);
    continue;
  }
  checkFile(file, schema);
}

// Regional files: one per locus listed in regional/index.json.
const index = JSON.parse(readFileSync(join(DATA, "regional", "index.json"), "utf8")) as {
  loci: { file: string }[];
};
for (const entry of index.loci) {
  checkFile(entry.file, REGIONAL_SCHEMA);
}

// Any other JSON in data/ that is not in the registry is reported (not fatal).
const known = new Set([...DATA_FILES.map((d) => d.file), ...index.loci.map((l) => l.file)]);
const walk = (dir: string, prefix = ""): void => {
  for (const name of readdirSync(dir)) {
    const abs = join(dir, name);
    const rel = prefix ? `${prefix}/${name}` : name;
    if (statSync(abs).isDirectory()) walk(abs, rel);
    else if (name.endsWith(".json") && !known.has(rel)) warnings.push(`not validated: ${rel}`);
  }
};
walk(DATA);

for (const w of warnings) console.warn(`warn ${w}`);
if (failures > 0) {
  console.error(`${failures} validation failure(s)`);
  process.exit(1);
}
console.log("data validation passed");
