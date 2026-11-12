// The TTA Flora 4 model, ported line for line from sim/rules.py. Keep the two in step: sim/ is the reference, and
// tools/check_js_engine.py compares this engine with it.
//
// The system is a string of n + 2 letters, one per site: the sensitizer (site 0), DPA 1 to n (sites 1 to n) and the
// perylene (site n + 1). Each letter is G (ground), T (triplet) or S (singlet). Times are in ns and rates per ns.

export const EVENTS = ['excitation', 'isc', 'sens_singlet_lost', 'hop', 'triplet_lost', 'tta_singlet', 'tta_failed',
  'photon', 'back_transfer'];
export const N_EVENTS = EVENTS.length;
export const EV = Object.fromEntries(EVENTS.map((e, i) => [e, i]));
export const LETTERS = 'GTS';

export const DEFAULTS = {
  n: 8, // DPA mediators
  tauS: 5.0, // singlet lifetime on a DPA or the perylene, ns
  tauSens: 5.0, // sensitizer singlet decay lifetime, ns (rule 9), not counting ISC
  tripletYield: 0.5, // share of sensitizer singlets that become triplets. This sets kISC
  hopRatio: 100.0, // kHop = hopRatio / tauS
  excRatio: 0.1, // kExc = excRatio / tauSens
  kDecay: 1e-6, // triplet decay on a DPA or the perylene, /ns
  f: 0.4, // the share of TTA events that make a singlet
  fretRadius: 4.0, // nm where back transfer is 50% (0 = none)
  spacing: 0.0, // nm: 0 = distances from the real molecule; above 0 = site i at i * spacing
  block: 1, // 1 = a triplet can't hop onto a singlet
  kISC: null,
  kHop: null,
  kExc: null,
};
export const PULSE_DEFAULTS = { ...DEFAULTS, excRatio: 1.0, duration: 5.0 };
export const FLUORESCENCE_DEFAULTS = { tauFl: 5.0, fluorYield: 1.0 };

/** Fill in kISC, kHop and kExc from the yield and ratios, or the other way round if a rate was given directly. */
export function resolve(p) {
  p = { ...p };
  const decay = 1 / p.tauSens;
  if (p.kISC == null) {
    if (p.tripletYield >= 1) throw new Error('tripletYield must be below 1');
    p.kISC = (decay * p.tripletYield) / (1 - p.tripletYield);
  } else p.tripletYield = p.kISC / (p.kISC + decay);
  if (p.kHop == null) p.kHop = p.hopRatio / p.tauS;
  else p.hopRatio = p.kHop * p.tauS;
  if (p.kExc == null) p.kExc = p.excRatio * decay;
  else p.excRatio = p.kExc * p.tauSens;
  return p;
}

const centre = (rows, idx) => {
  const c = [0, 0, 0];
  for (const i of idx) for (let k = 0; k < 3; k++) c[k] += rows[i][k + 1] / idx.length;
  return c;
};
const mean = (rows) => centre(rows, rows.map((_, i) => i));
const dist = (u, v) => Math.hypot(u[0] - v[0], u[1] - v[1], u[2] - v[2]);

/**
 * Distance in nm from the sensitizer to sites 1 to n + 1, measured on the same geometry the viewer builds
 * (`geometry` = { flora, wire }: samples/triplet/tta-flora-4.param.json and the wire it names).
 */
export function siteDistances(n, spacing, geometry) {
  if (spacing > 0) return Array.from({ length: n + 1 }, (_, i) => spacing * (i + 1));
  const { flora, wire } = geometry;
  if (n === 0) return [dist(centre(flora.zero, flora.zeroRightCore), centre(flora.zero, flora.zeroLeftCore)) / 10];
  const t = wire.t;
  const sens = centre(flora.left, flora.leftCore);
  const dpa = mean(wire.dpa);
  const per = centre(flora.right, flora.rightCore).map((v, k) => v + (n - 1) * t[k]);
  const centres = [...Array.from({ length: n }, (_, i) => dpa.map((v, k) => v + i * t[k])), per];
  return centres.map((c) => dist(c, sens) / 10);
}

/** E = 1 / (1 + (r / R0)^6): the chance a singlet at distance r goes back to the sensitizer instead of emitting. */
export const fretEfficiency = (r, radius) => r.map((x) => (radius <= 0 ? 0 : 1 / (1 + (x / radius) ** 6)));

/** Every allowed hop as [from, to, rate]. Nothing hops onto the sensitizer or off the perylene. */
export function hops(n, kHop) {
  const moves = [[0, 1, kHop]];
  for (let site = 1; site <= n; site++) {
    const nb = [site - 1, site + 1].filter((s) => s >= 1);
    for (const to of nb) moves.push([site, to, kHop / nb.length]);
  }
  return moves;
}

const L = (c) => LETTERS.indexOf(c);

/**
 * The rules for parameters p, as [a, b (or -1), before a, before b, after a, after b, rate, event, where]. `distances`
 * is siteDistances(n, spacing, ...). Rules with rate 0 are left out.
 */
export function buildRules(p, distances) {
  p = resolve(p);
  const n = p.n;
  const E = fretEfficiency(distances, p.fretRadius);
  const rules = [];
  const one = (site, before, after, rate, event) => rules.push([site, -1, L(before), 0, L(after), 0, rate, EV[event], site]);
  const two = (a, b, before, after, rate, event, where) => {
    const [x, y] = before.split(' ');
    const [u, v] = after.split(' ');
    rules.push([a, b, L(x), L(y), L(u), L(v), rate, EV[event], where]);
  };
  one(0, 'G', 'S', p.kExc, 'excitation'); // rule 1
  one(0, 'S', 'T', p.kISC, 'isc'); // rule 2
  one(0, 'S', 'G', 1 / p.tauSens, 'sens_singlet_lost'); // rule 9
  const emit = 1 / p.tauS;
  for (let site = 1; site <= n + 1; site++) {
    one(site, 'T', 'G', p.kDecay, 'triplet_lost'); // rule 3
    const back = E[site - 1];
    two(0, site, 'G S', 'G G', emit * (1 - back), 'photon', site); // rules 4 and 10
    two(0, site, 'G S', 'S G', emit * back, 'back_transfer', site);
    two(0, site, 'S S', 'S G', emit, 'photon', site);
    two(0, site, 'T S', 'T G', emit, 'photon', site);
  }
  for (const [a, b, rate] of hops(n, p.kHop)) {
    two(a, b, 'T G', 'G T', rate, 'hop', b); // rule 5
    if (!p.block) two(a, b, 'T S', 'G T', rate, 'hop', b); // rule 6 is the absence of this
    two(a, b, 'T T', 'G S', p.f * rate, 'tta_singlet', b); // rule 7
    two(a, b, 'T T', 'G T', (1 - p.f) * rate, 'tta_failed', b); // rule 8
  }
  return rules.filter((r) => r[6] > 0);
}

/** A two-state fluorophore on site 0, excited at kExc and back to G at 1/tauFl, a photon with chance fluorYield. */
export function fluorescenceRules(p) {
  const decay = 1 / p.tauFl;
  return [
    [0, -1, 0, 0, 2, 0, p.kExc, EV.excitation, 0],
    [0, -1, 2, 0, 0, 0, decay * p.fluorYield, EV.photon, 0],
    [0, -1, 2, 0, 0, 0, decay * (1 - p.fluorYield), EV.sens_singlet_lost, 0],
  ].filter((r) => r[6] > 0);
}

/** The rules as typed arrays, the form the engines use (and can send to a worker). */
export function ruleArrays(rules) {
  const k = rules.length;
  const R = {
    count: k, a: new Int32Array(k), b: new Int32Array(k), fa: new Int32Array(k), fb: new Int32Array(k),
    ta: new Int32Array(k), tb: new Int32Array(k), rate: new Float64Array(k), event: new Int32Array(k),
    where: new Int32Array(k),
  };
  rules.forEach((r, i) => {
    [R.a[i], R.b[i], R.fa[i], R.fb[i], R.ta[i], R.tb[i], R.rate[i], R.event[i], R.where[i]] = r;
  });
  return R;
}
