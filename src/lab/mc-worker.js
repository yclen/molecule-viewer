// A Monte Carlo worker. Several run at once, one per spare CPU core.
//
// Live mode (the pulse page): { type: 'config', generation, R, sites, duration, target, seed } sets what to run and
// { type: 'on', on } starts and stops it. While on, it runs small batches of molecules (about 20 ms each, so new
// settings take effect quickly), paced to `target` molecules per second (0 = as fast as it can), and posts each
// batch's photons: { generation, molecules, excitations, times, sites }.
//
// Batch mode (the sweep page): { type: 'batch', id, R, sites, duration, molecules, seed } runs that many molecules
// and posts { id, molecules, counts, excitations, sums }.

import { makeRng, runMolecules } from '../sim/mc.js';

const BATCH_S = 0.02;
const MOST_MOLECULES = 2000;
const CATCH_UP_S = 0.25; // a worker that fell behind can make up this much of the schedule; the rest is lost

let cfg = null;
let on = false;
let running = false;
let rng = makeRng(1);
let speed = 0;
let seen = -1;
let nextStart = 0;

const now = () => performance.now() / 1000;
const channel = new MessageChannel();
let wake = null;
channel.port1.onmessage = () => wake?.();
const yieldNow = () => new Promise((resolve) => {
  wake = resolve;
  channel.port2.postMessage(0);
});
const sleep = (s) => (s > 0.001 ? new Promise((resolve) => setTimeout(resolve, s * 1000)) : yieldNow());

async function loop() {
  if (running) return;
  running = true;
  while (on && cfg) {
    const c = cfg;
    let size = seen === c.generation ? Math.min(Math.max(Math.round(BATCH_S * speed), 1), MOST_MOLECULES) : 1;
    if (c.target) size = Math.min(size, Math.max(1, Math.round(BATCH_S * c.target)));
    const start = nextStart;
    if (c.target) nextStart = Math.max(start, now() - CATCH_UP_S) + size / c.target;
    await sleep(c.target ? start - now() : 0); // also lets new messages in
    if (!on || cfg !== c) continue;
    const began = now();
    const run = runMolecules(c.R, c.sites, c.duration, size, rng, { photons: true });
    speed = size / Math.max(now() - began, 1e-6);
    seen = c.generation;
    postMessage({ generation: c.generation, molecules: size, excitations: run.excitations, times: run.photonTimes,
      sites: run.photonSites }, [run.photonTimes.buffer, run.photonSites.buffer]);
  }
  running = false;
}

onmessage = ({ data: m }) => {
  if (m.type === 'config') {
    if (!cfg || m.seed !== cfg.seed) rng = makeRng(m.seed);
    cfg = m;
    nextStart = now();
    loop();
  } else if (m.type === 'on') {
    on = m.on;
    nextStart = now();
    loop();
  } else if (m.type === 'batch') {
    const run = runMolecules(m.R, m.sites, m.duration, m.molecules, makeRng(m.seed), { perMolecule: true });
    postMessage({ id: m.id, molecules: m.molecules, counts: run.counts, excitations: run.excitations, sums: run.sums },
      [run.counts.buffer]);
  }
};
