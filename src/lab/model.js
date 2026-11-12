// What the pulse and sweep pages share: the molecule's geometry (for the FRET distances), rules for a parameter set,
// the laser-power calibration and the value syntax of the command line. Ported from sim/live.py and sim/tta_sim.py.

import flora from '../../samples/triplet/tta-flora-4.param.json';
import wire from '../../samples/triplet/dpa-wire.param.json';
import { buildRules, FLUORESCENCE_DEFAULTS, fluorescenceRules, PULSE_DEFAULTS, resolve, ruleArrays, siteDistances } from '../sim/rules.js';

const geometry = { flora, wire };
const distances = new Map();

/** Rules and their typed arrays for a parameter set p, the TTA molecule or (p.sample === 'fluorescence') the dye. */
export function rulesFor(p) {
  if (p.sample === 'fluorescence') {
    const rules = fluorescenceRules(p);
    return { rules, R: ruleArrays(rules), sites: 1 };
  }
  const key = `${p.n}|${p.spacing ?? 0}`;
  if (!distances.has(key)) distances.set(key, siteDistances(p.n, p.spacing ?? 0, geometry));
  const rules = buildRules(p, distances.get(key));
  return { rules, R: ruleArrays(rules), sites: p.n + 2 };
}

export const TTA_FIELDS = ['n', 'kHop', 'f', 'tripletYield', 'tauS', 'tauSens', 'kDecay', 'fretRadius'];
export const FLUOR_FIELDS = ['tauFl', 'fluorYield'];

// label, unit, lowest and highest allowed value (live.py FIELDS)
export const FIELDS = {
  duration: ['Pulse length', 'ns', 1e-3, 1e6],
  window: ['Time window', 'ns', 1e-2, 1e7],
  n: ['DPA mediators, n', '', 0, 40],
  kHop: ['Hop rate, kHop', '/ns', 1e-6, 1e7],
  f: ['Spin factor, f', '', 0, 1],
  tripletYield: ['Triplet yield', '', 0, 0.999],
  tauS: ['Singlet lifetime, τS', 'ns', 1e-3, 1e6],
  tauSens: ['Sensitizer lifetime, τSens', 'ns', 1e-3, 1e6],
  kDecay: ['Triplet decay, kDecay', '/ns', 0, 1e3],
  fretRadius: ['FRET radius, R0', 'nm', 0, 100],
  tauFl: ['Fluorescence lifetime, τF', 'ns', 1e-3, 1e6],
  fluorYield: ['Fluorescence yield, Φ_F', '', 0, 1],
};

/** The pulse defaults with the ratios turned into flat rates per ns: what both pages start with. */
export function defaultValues() {
  const p = resolve(PULSE_DEFAULTS);
  const values = Object.fromEntries([...TTA_FIELDS, 'duration'].map((k) => [k, p[k]]));
  return { ...values, ...FLUORESCENCE_DEFAULTS, window: 10 * p.tauS, kExc: p.kExc, sample: 'tta' };
}

/** The engine's parameters for the page's values: kHop and kExc given directly, kISC worked out from the yield. */
export function parameters(values) {
  const shared = { sample: values.sample, kExc: values.kExc, duration: values.duration, window: values.window };
  if (values.sample === 'fluorescence') return { ...shared, tauFl: values.tauFl, fluorYield: values.fluorYield };
  const p = resolve({ ...PULSE_DEFAULTS, ...shared, ...Object.fromEntries(TTA_FIELDS.map((k) => [k, values[k]])) });
  delete p.excRatio;
  delete p.hopRatio;
  return p;
}

// The assumed setup, used only to show the laser power that gives a kExc (live.py)
export const BEAM = { radiusUm: 1.0, extinction: 50_000, wavelengthNm: 635 };

/** kExc (/ns) per mW of laser power for a molecule at the centre of a Gaussian beam: 0.0389 /ns per mW here. */
export function kexcPerMw() {
  const crossSection = (Math.log(10) * 1000 * BEAM.extinction) / 6.02214076e23; // cm²
  const peak = (2 * 1e-3) / (Math.PI * (BEAM.radiusUm * 1e-4) ** 2); // W/cm² from 1 mW
  const photon = (6.62607015e-34 * 2.99792458e8) / (BEAM.wavelengthNm * 1e-9); // J
  return ((crossSection * peak) / photon) * 1e-9;
}

export function powerText(mw) {
  for (const [unit, size] of [['W', 1e3], ['mW', 1], ['µW', 1e-3]]) if (mw >= size) return `${fmt(mw / size, 3)} ${unit}`;
  return `${fmt(mw * 1e6, 3)} nW`;
}

/** '5' | '1,2,5' | '0:12' (whole numbers) | '1:300:6' (6 evenly spaced) | '1:300:log6' (6 spaced by ratio). */
export function parseValues(name, text) {
  text = String(text).trim();
  let values;
  if (text.includes(',')) values = text.split(',').map(Number);
  else if (text.includes(':')) {
    const parts = text.split(':');
    const lo = Number(parts[0]);
    const hi = Number(parts[1]);
    if (parts.length === 2) values = Array.from({ length: Math.floor(hi - lo) + 1 }, (_, i) => lo + i);
    else if (parts[2].startsWith('log')) {
      const count = parseInt(parts[2].slice(3), 10);
      values = Array.from({ length: count }, (_, i) => (count === 1 ? lo : lo * (hi / lo) ** (i / (count - 1))));
    } else {
      const count = parseInt(parts[2], 10);
      values = Array.from({ length: count }, (_, i) => (count === 1 ? lo : lo + ((hi - lo) * i) / (count - 1)));
    }
  } else values = [Number(text)];
  if (!values.length || values.some((v) => !Number.isFinite(v))) throw new Error(`can't read "${text}"`);
  return name === 'n' ? values.map((v) => Math.round(v)) : values;
}

/** A number with `digits` significant figures, like Python's :.3g. */
export function fmt(v, digits = 3) {
  if (!Number.isFinite(v)) return v == null || Number.isNaN(v) ? '–' : String(v);
  if (v === 0) return '0';
  const e = Math.floor(Math.log10(Math.abs(v)));
  if (e < -4 || e >= digits) return v.toExponential(digits - 1).replace(/\.?0+e/, 'e').replace('e+', 'e');
  return String(Number(v.toPrecision(digits)));
}

// How every plot's axes look: big, bold, near-white axis titles (each ending in its unit) and large tick numbers, so
// what an axis measures, and in what unit (ns, not µs), reads at a glance
export const PLOT_FONT = "'Segoe UI', system-ui, -apple-system, Roboto, sans-serif";
export const AXIS = {
  stroke: '#e8edf3',
  grid: { stroke: 'rgba(255,255,255,0.07)', width: 1 },
  ticks: { stroke: 'rgba(255,255,255,0.25)', width: 1 },
  font: `16px ${PLOT_FONT}`,
  labelFont: `700 20px ${PLOT_FONT}`,
  labelSize: 40,
  labelGap: 6,
};
export const Y_AXIS_SIZE = 86; // room for the tick numbers at that size

export const int = (v) => Math.round(v).toLocaleString('en-US');
export const percent = (v, d = 2) => (Number.isFinite(v) ? `${(100 * v).toFixed(d)}%` : '–');

export function howLong(seconds) {
  if (seconds < 1) return 'a second';
  if (seconds < 90) return `${Math.round(seconds)} s`;
  if (seconds < 90 * 60) return `${Math.round(seconds / 60)} min`;
  return `${Math.round(seconds / 3600)} h`;
}

/** Rough seconds for pulseExact, from timings of the JS engine (n = 5 in ~1 s, n = 8 in ~28 s at a 50 ns window). */
export function exactSeconds(p, window) {
  const size = 3 ** (p.n + 2);
  const step = Math.max(Math.min(p.tauS, p.tauSens) / 10, window / 2000);
  const steps = Math.max(10, Math.ceil(Math.min(p.duration, window) / step)) + 1.25 * Math.max(0, Math.ceil((window - p.duration) / step));
  return 0.2 + 4.5e-6 * size * steps * (1 + 0.3 * Math.log10(Math.max((p.kHop * step) / 5, 1)));
}
