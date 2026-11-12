// Runs the browser's simulation engine (src/sim/) on one job from stdin and prints the result as JSON, so
// tools/check_js_engine.py can compare it with sim/engine.py.
//
//   echo '{"kind": "rules", "p": {...}}' | node tools/js_engine_cli.mjs

import { readFileSync } from 'node:fs';
import { binnedExact, pulseExact, pulseTotals } from '../src/sim/exact.js';
import { makeRng, runMolecules } from '../src/sim/mc.js';
import { buildRules, ruleArrays, siteDistances } from '../src/sim/rules.js';

const read = (name) => JSON.parse(readFileSync(new URL(`../samples/triplet/${name}`, import.meta.url)));
const flora = read('tta-flora-4.param.json');
const geometry = { flora, wire: read(flora.wire) };
const job = JSON.parse(readFileSync(0, 'utf8'));
const p = job.p;
const rules = buildRules(p, siteDistances(p.n, p.spacing, geometry));
const R = ruleArrays(rules);
const sites = p.n + 2;
const began = performance.now();
let out;
if (job.kind === 'rules') out = rules;
else if (job.kind === 'totals') out = Array.from(pulseTotals(R, p, rules), (v) => (Number.isNaN(v) ? null : v));
else if (job.kind === 'exact') {
  const result = pulseExact(R, p, job.window);
  out = { binned: binnedExact(result, job.edges).map((a) => Array.from(a)), photons: result.photons,
    excitations: result.excitations };
} else if (job.kind === 'mc') {
  const run = runMolecules(R, sites, p.duration, job.molecules, makeRng(job.seed), { photons: true, perMolecule: true });
  out = { counts: Array.from(run.counts), excitations: run.excitations, sums: run.sums,
    times: Array.from(run.photonTimes), sites: Array.from(run.photonSites) };
}
console.log(JSON.stringify({ out, seconds: (performance.now() - began) / 1000 }));
