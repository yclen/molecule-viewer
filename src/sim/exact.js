// The exact (noise-free) solvers, ported from sim/engine.py: pulse_exact (the expected pulse histogram), pulse_totals
// (every event per molecule over one whole pulse) and fluorescence_exact (the two-state fluorophore, in closed form).
// They work on all 3^sites strings at once, so they are for n up to about 8.

import { EV, N_EVENTS } from './rules.js';
import { Coo, expmv, gmres, ilu0, matvec } from './sparse.js';

const EXCITATION = EV.excitation;
const PHOTON = EV.photon;
const HOP = EV.hop;

/** Every rule as { src, dst, rate, event, where } over all 3^sites strings. String k has letter (k / 3^s) % 3 at site s. */
export function transitions(R, sites) {
  const size = 3 ** sites;
  const pow = Array.from({ length: sites }, (_, s) => 3 ** s);
  const digit = new Uint8Array(size * sites);
  for (let k = 0; k < size; k++) for (let s = 0, v = k; s < sites; s++, v = (v - (v % 3)) / 3) digit[k * sites + s] = v % 3;
  const moves = [];
  for (let r = 0; r < R.count; r++) {
    const a = R.a[r], b = R.b[r], fa = R.fa[r], fb = R.fb[r];
    const shift = (R.ta[r] - fa) * pow[a] + (b >= 0 ? (R.tb[r] - fb) * pow[b] : 0);
    const src = [];
    for (let k = 0; k < size; k++) if (digit[k * sites + a] === fa && (b < 0 || digit[k * sites + b] === fb)) src.push(k);
    const s = Int32Array.from(src);
    moves.push({ src: s, dst: s.map((k) => k + shift), rate: R.rate[r], event: R.event[r], where: R.where[r] });
  }
  return { moves, digit, size };
}

// pulse_exact steps a stiff system with a 2-stage L-stable SDIRK; each stage is a GMRES solve with the same matrix.
const GAMMA = 1 - 2 ** -0.5;
const STEPS_PER_LIFETIME = 10;
const LEAST_STEPS = 10;
const MOST_STEPS = 2000;
const MOST_TAIL_STEPS = 2000;
export const QUIET = 1e-7; // the count after the window stops when the molecules that can still emit are this share

/** The time steps: [[light on, step length, number of steps]], up to the window. A step ends where the pulse does. */
export function exactPlan(p, window) {
  const step = Math.max(Math.min(p.tauS, p.tauSens) / STEPS_PER_LIFETIME, window / MOST_STEPS);
  const plan = [];
  for (const [light, start, end] of [[true, 0, Math.min(p.duration, window)], [false, p.duration, window]]) {
    if (end > start) {
      const count = Math.max(LEAST_STEPS, Math.ceil((end - start) / step - 1e-9));
      plan.push([light, (end - start) / count, count]);
    }
  }
  return plan;
}

/**
 * The expected pulse histogram for resolved parameters p, with kExc on until p.duration. `R` is ruleArrays(buildRules(p)).
 * progress(share done, what it's doing) is called as it goes. Returns a PulseExact.
 */
export function pulseExact(R, p, window, progress = () => {}) {
  const sites = p.n + 2;
  const { moves, size } = transitions(R, sites);
  const N = size + 3; // + perylene photons, DPA photons, excitations
  const counter = (event, where) => (event === PHOTON ? (where === sites - 1 ? size : size + 1) : event === EXCITATION ? size + 2 : -1);

  // A = I - γ h M, and M itself: the master equation plus rows that add up the counted events
  const matrices = (light, h) => {
    const M = new Coo(N, 4 * size);
    const A = new Coo(N, 4 * size);
    for (let i = 0; i < N; i++) A.push(i, i, 1);
    for (const { src, dst, rate, event, where } of moves) {
      if (event === EXCITATION && !light) continue;
      const row = counter(event, where);
      for (let i = 0; i < src.length; i++) {
        M.push(dst[i], src[i], rate);
        M.push(src[i], src[i], -rate);
        A.push(dst[i], src[i], -GAMMA * h * rate);
        A.push(src[i], src[i], GAMMA * h * rate);
        if (row >= 0) {
          M.push(row, src[i], rate);
          A.push(row, src[i], -GAMMA * h * rate);
        }
      }
    }
    return [M.toCsr(), A.toCsr()];
  };

  const plan = exactPlan(p, window);
  const planned = plan.reduce((s, [light, , count]) => s + count * (light ? 1 : 1.25), 0); // the tail is a guess
  let done = 0;
  const quiet = [0, ...Array.from({ length: sites }, (_, s) => 3 ** s)]; // all G, or one triplet and nothing else
  let x = new Float64Array(N);
  x[0] = 1;
  const times = [0];
  const counts = [[0, 0, 0]];
  const rates = [[0, 0]];
  let stillActive = NaN;
  const Mx = new Float64Array(N);
  for (const [light, h, count] of plan) {
    progress(done / planned, 'Preparing');
    const [M, A] = matrices(light, h);
    const precond = ilu0(A);
    let k1 = new Float64Array(N);
    let k2 = new Float64Array(N);
    const step = () => {
      matvec(M, x, Mx);
      k1 = gmres(A, Mx, k1, precond, { rtol: 1e-8 });
      const mid = new Float64Array(N);
      for (let i = 0; i < N; i++) mid[i] = x[i] + (1 - GAMMA) * h * k1[i];
      matvec(M, mid, Mx);
      k2 = gmres(A, Mx, k2, precond, { rtol: 1e-8 });
      for (let i = 0; i < N; i++) x[i] += h * ((1 - GAMMA) * k1[i] + GAMMA * k2[i]);
    };
    for (let s = 0; s < count; s++) {
      step();
      times.push(times[times.length - 1] + h);
      counts.push([x[size], x[size + 1], x[size + 2]]);
      matvec(M, x, Mx);
      rates.push([Mx[size], Mx[size + 1]]);
      progress(++done / planned, 'Solving');
    }
    if (!light) {
      for (let s = 0; s < MOST_TAIL_STEPS; s++) {
        let all = 0;
        for (let i = 0; i < size; i++) all += x[i];
        stillActive = all - quiet.reduce((t, q) => t + x[q], 0);
        if (stillActive < QUIET * x[size + 2]) break;
        step();
        progress(Math.min(++done / planned, 0.99), 'Counting the photons after the window');
      }
    }
  }
  const whole = p.duration <= window;
  return {
    times: Float64Array.from(times),
    perylene: Float64Array.from(counts, (c) => c[0]),
    dpa: Float64Array.from(counts, (c) => c[1]),
    peryleneRate: Float64Array.from(rates, (r) => r[0]),
    dpaRate: Float64Array.from(rates, (r) => r[1]),
    excitations: whole ? x[size + 2] : NaN,
    photons: whole ? x[size] + x[size + 1] : NaN,
    stillActive,
  };
}

/** A cubic Hermite curve through (times, totals) with the given slopes, evaluated at t (clipped to the times). */
function hermite(times, totals, slopes, t) {
  const last = times.length - 1;
  t = Math.min(Math.max(t, 0), times[last]);
  let lo = 0;
  let hi = last;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (times[mid] <= t) lo = mid;
    else hi = mid;
  }
  const h = times[hi] - times[lo];
  const s = (t - times[lo]) / h;
  const s2 = s * s;
  const s3 = s2 * s;
  return (2 * s3 - 3 * s2 + 1) * totals[lo] + (s3 - 2 * s2 + s) * h * slopes[lo] + (-2 * s3 + 3 * s2) * totals[hi]
    + (s3 - s2) * h * slopes[hi];
}

/** Expected [perylene, DPA] photons per molecule in each bin between `edges`, from a pulseExact result. */
export function binnedExact(result, edges) {
  if (result.fluorescence) return binnedFluorescence(result.p, edges);
  return [[result.perylene, result.peryleneRate], [result.dpa, result.dpaRate]].map(([totals, slopes]) => {
    const at = Array.from(edges, (t) => hermite(result.times, totals, slopes, t));
    return Float64Array.from({ length: edges.length - 1 }, (_, k) => Math.max(at[k + 1] - at[k], 0));
  });
}

/** pulse_exact's answer for the two-state fluorophore, in closed form. */
export function fluorescenceExact(p) {
  const excitations = fluorescenceDecays(p, Infinity);
  return { fluorescence: true, p, excitations, photons: p.fluorYield * excitations, stillActive: 0 };
}

function fluorescenceDecays(p, t) {
  const a = p.kExc;
  const k = 1 / p.tauFl;
  const s = a + k;
  const lit = Math.min(t, p.duration);
  const during = ((k * a) / s) * (lit + Math.expm1(-s * lit) / s);
  const atEnd = (-a / s) * Math.expm1(-s * p.duration);
  const dark = Math.max(t - p.duration, 0);
  return during + atEnd * -Math.expm1(-k * dark);
}

function binnedFluorescence(p, edges) {
  const at = Array.from(edges, (t) => p.fluorYield * fluorescenceDecays(p, t));
  return [Float64Array.from({ length: edges.length - 1 }, (_, k) => at[k + 1] - at[k]), new Float64Array(edges.length - 1)];
}

let after = { key: null }; // the last after-pulse matrix and preconditioner: they don't depend on kExc or the pulse

/**
 * Every event one molecule gives over one whole pulse, on average, with no random noise: a Float64Array
 * [event * sites + site], per molecule. Hops after the pulse aren't counted (NaN). See pulse_totals in engine.py.
 */
export function pulseTotals(R, p, rules) {
  const sites = p.n + 2;
  const { moves, digit, size } = transitions(R, sites);
  const excited = new Int32Array(size);
  const singlets = new Int32Array(size);
  for (let k = 0; k < size; k++)
    for (let s = 0; s < sites; s++) {
      const d = digit[k * sites + s];
      if (d > 0) excited[k]++;
      if (d === 2) singlets[k]++;
    }
  const lone = Uint8Array.from(excited, (e, k) => (e === 1 && singlets[k] === 0 ? 1 : 0));
  const active = Uint8Array.from(excited, (e, k) => (e > 0 && !lone[k] ? 1 : 0));

  let counts = new Float64Array(N_EVENTS * sites);
  let x = new Float64Array(size);
  x[0] = 1;
  if (p.duration > 0 && p.kExc > 0) {
    const N = size + N_EVENTS * sites;
    const A = new Coo(N, 6 * size);
    for (let i = 0; i < N; i++) A.push(i, i, 0); // every row has a diagonal entry
    for (const { src, dst, rate, event, where } of moves) {
      const row = size + event * sites + where;
      for (let i = 0; i < src.length; i++) {
        A.push(dst[i], src[i], rate);
        A.push(src[i], src[i], -rate);
        A.push(row, src[i], rate);
      }
    }
    const start = new Float64Array(N);
    start[0] = 1;
    const end = expmv(A.toCsr(), p.duration, start);
    x = end.slice(0, size);
    counts = end.slice(size);
  }

  const afterMoves = moves.filter((m) => m.event !== EXCITATION);
  const key = `${sites}|${JSON.stringify(rules.filter((r) => r[7] !== EXCITATION))}`;
  const keep = [];
  const at = new Int32Array(size).fill(-1);
  for (let k = 0; k < size; k++) if (active[k]) at[k] = keep.push(k) - 1;
  if (after.key !== key) {
    const Q = new Coo(keep.length, 4 * keep.length);
    for (let i = 0; i < keep.length; i++) Q.push(i, i, 0);
    for (const { src, dst, rate } of afterMoves)
      for (let i = 0; i < src.length; i++) {
        const c = at[src[i]];
        if (c < 0) continue;
        Q.push(c, c, rate); // this is -Q restricted to the active strings
        const r = at[dst[i]];
        if (r >= 0) Q.push(r, c, -rate);
      }
    const matrix = Q.toCsr();
    after = { key, matrix, precond: ilu0(matrix) };
  }
  const rhs = Float64Array.from(keep, (k) => x[k]);
  const tau = gmres(after.matrix, rhs, null, after.precond, { rtol: 1e-12, restart: 60, maxRestarts: 80 });
  const timeIn = new Float64Array(size); // ns spent in each string after the pulse, until all G or one triplet
  keep.forEach((k, i) => (timeIn[k] = tau[i]));
  let endsWithTriplet = 0;
  for (let k = 0; k < size; k++) if (lone[k]) endsWithTriplet += x[k];
  for (const { src, dst, rate, event, where } of afterMoves) {
    let spent = 0;
    for (let i = 0; i < src.length; i++) {
      spent += timeIn[src[i]];
      if (lone[dst[i]]) endsWithTriplet += rate * timeIn[src[i]];
    }
    counts[event * sites + where] += rate * spent;
  }
  counts[EV.triplet_lost * sites + sites - 1] += endsWithTriplet;
  for (let s = 0; s < sites; s++) counts[HOP * sites + s] = NaN;
  return counts;
}
