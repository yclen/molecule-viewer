// The pulse Monte Carlo, ported from _rates, _pick and _pulse_chunk in sim/engine.py. Every molecule starts in all G,
// kExc is on until `duration`, and the molecule is followed until nothing else can happen (the Gillespie method: add
// up the rates of every rule that can fire, wait a random time, fire one).

import { EV, N_EVENTS } from './rules.js';

const EXCITATION = EV.excitation;
const PHOTON = EV.photon;

/** A seeded generator (xoshiro128**). random() is uniform on [0, 1) with 53 random bits. */
export function makeRng(seed) {
  let s0 = 0x9e3779b9 ^ seed, s1 = 0x243f6a88 ^ Math.imul(seed, 0x85ebca6b), s2 = 0xb7e15162, s3 = 0x13198a2e;
  const next = () => {
    const r = Math.imul(rotl(Math.imul(s1, 5), 7), 9);
    const t = s1 << 9;
    s2 ^= s0; s3 ^= s1; s1 ^= s2; s0 ^= s3; s2 ^= t; s3 = rotl(s3, 11);
    return r >>> 0;
  };
  for (let i = 0; i < 16; i++) next(); // let the seed spread through all four words
  return { random: () => ((next() >>> 5) * 67108864 + (next() >>> 6)) / 9007199254740992 };
}
const rotl = (x, k) => (x << k) | (x >>> (32 - k));

/**
 * Run `molecules` molecules one after the other.
 *
 * opts.photons: also return every photon's emission time and site (for the live histogram).
 * opts.perMolecule: also add up photons², excitations² and photons × excitations over the molecules (for the
 * sweep's error bars).
 * Returns { counts: Float64Array [event * sites + site], excitations, photonTimes?, photonSites?, sums? }.
 */
export function runMolecules(R, sites, duration, molecules, rng, opts = {}) {
  const { a, b, fa, fb, ta, tb, rate, event, where } = R;
  const k = R.count;
  const state = new Int32Array(sites);
  const weights = new Float64Array(k);
  const counts = new Float64Array(N_EVENTS * sites);
  let times = opts.photons ? new Float64Array(Math.max(256, 2 * molecules)) : null;
  let from = opts.photons ? new Int16Array(times.length) : null;
  let found = 0;
  let excitations = 0;
  const sums = opts.perMolecule ? [0, 0, 0] : null;
  const random = rng.random;

  for (let m = 0; m < molecules; m++) {
    state.fill(0);
    let t = 0;
    let molPhotons = 0;
    let molExc = 0;
    for (;;) {
      const lightOn = t < duration;
      let total = 0;
      for (let r = 0; r < k; r++) {
        const can = state[a[r]] === fa[r] && (b[r] < 0 || state[b[r]] === fb[r]) && (lightOn || event[r] !== EXCITATION);
        const w = can ? rate[r] : 0;
        weights[r] = w;
        total += w;
      }
      if (total === 0) break; // back to all G after the pulse: nothing else can happen
      const dt = -Math.log(1 - random()) / total;
      if (lightOn && t + dt >= duration) {
        t = duration; // the rates change here. Waits have no memory, so just draw a new one
        continue;
      }
      t += dt;
      // pick a rule, each with a chance in proportion to its rate
      const u = random() * total;
      let r = 0;
      let acc = weights[0];
      while ((acc < u || weights[r] === 0) && r < k - 1) acc += weights[++r];
      state[a[r]] = ta[r];
      if (b[r] >= 0) state[b[r]] = tb[r];
      const e = event[r];
      counts[e * sites + where[r]] += 1;
      if (e === EXCITATION) molExc++;
      else if (e === PHOTON) {
        molPhotons++;
        if (times) {
          if (found === times.length) {
            const bigger = new Float64Array(times.length * 2);
            bigger.set(times);
            times = bigger;
            const fromBigger = new Int16Array(bigger.length);
            fromBigger.set(from);
            from = fromBigger;
          }
          times[found] = t;
          from[found] = where[r];
          found++;
        }
      }
    }
    excitations += molExc;
    if (sums) {
      sums[0] += molPhotons * molPhotons;
      sums[1] += molExc * molExc;
      sums[2] += molPhotons * molExc;
    }
  }
  const out = { counts, excitations };
  if (times) {
    out.photonTimes = times.slice(0, found);
    out.photonSites = from.slice(0, found);
  }
  if (sums) out.sums = sums;
  return out;
}
