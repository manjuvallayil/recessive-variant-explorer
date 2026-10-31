"use client";

// Lazy per-section data loading. A panel calls useData(path, schema, enabled) and fetches only
// once it is enabled (its story step has been reached). Results are cached per path and parsed
// with the same Zod schema that the build-time validator uses.

import { useEffect, useState } from "react";
import type { z } from "zod";

const cache = new Map<string, Promise<unknown>>();

export function loadJson<T>(path: string, schema: z.ZodType<T>): Promise<T> {
  let p = cache.get(path);
  if (!p) {
    p = fetch(`/data/${path}`)
      .then((r) => {
        if (!r.ok) throw new Error(`HTTP ${r.status} for ${path}`);
        return r.json() as Promise<unknown>;
      })
      .then((json) => schema.parse(json));
    cache.set(path, p);
  }
  return p as Promise<T>;
}

export function useData<T>(path: string, schema: z.ZodType<T>, enabled: boolean) {
  const [state, setState] = useState<{ data: T | null; error: string | null }>({ data: null, error: null });
  useEffect(() => {
    if (!enabled) return;
    let alive = true;
    loadJson(path, schema)
      .then((data) => alive && setState({ data, error: null }))
      .catch((e: unknown) => alive && setState({ data: null, error: e instanceof Error ? e.message : String(e) }));
    return () => {
      alive = false;
    };
    // schema is a module-level constant; path controls the fetch.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [path, enabled]);
  return state;
}
