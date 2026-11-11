// Parametric objects: molecules built in the app from a number the user sets (e.g. how many DPA units
// a molecular wire has). The geometry of each repeating piece is computed offline (tools/build_wire.py,
// saved as samples/<tab>/*.param.json); here we only stack copies, so any size builds instantly.

import { planAttachment, placeSensitizer } from './attach.js';

const rel = (path) => path.replace('../samples/', '');
const DATA = Object.fromEntries(
  Object.entries(import.meta.glob('../samples/*/*.param.json', { import: 'default', eager: true })).map(
    ([path, data]) => [rel(path), data],
  ),
);

const add = (p, v, k = 1) => [p[0] + k * v[0], p[1] + k * v[1], p[2] + k * v[2]];
const dist = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
const centroid = (pts) => [0, 1, 2].map((c) => pts.reduce((s, p) => s + p[c], 0) / pts.length);
const rows = (list, k = 0, t = [0, 0, 0]) => list.map(([s, ...p]) => [s, add(p, t, k)]);

const range = (start, n) => Array.from({ length: n }, (_, i) => start + i);

/** Sites (for the excited-state glow) of the DPA units in a wire core that starts at atom `offset`. */
function dpaSites(w, n, offset = 0) {
  const step = w.dpa.length + w.linker.length;
  return range(0, n).map((k) => ({
    label: `DPA ${k + 1}`, // numbered even when n = 1, so an excited unit stays excited as n changes
    kind: 'mediator',
    atoms: range(offset + k * step, w.dpa.length),
  }));
}

/** The DPA wire core: n DPA units with n - 1 linkers between them (no end caps). */
function wireCore(w, n) {
  const atoms = [];
  for (let k = 0; k < n; k++) {
    atoms.push(...rows(w.dpa, k, w.t));
    if (k < n - 1) atoms.push(...rows(w.linker, k, w.t));
  }
  return atoms;
}

const plans = new Map(); // "<sensitizer key>|wire" or "|zero" -> attachment plan (reused for every n)

/**
 * Replace the sensitizer atoms `sensIdx` of end-group rows with `swap` ({ key, atoms: [[symbol, [x, y, z]]] }).
 * Returns { rows, sens, core, map (old kept index -> new index), plan }; the new sensitizer comes first.
 */
function swapSensitizer(rows0, sensIdx, swap, context, variant) {
  const key = `${swap.key}|${variant}`;
  if (!plans.has(key)) plans.set(key, planAttachment(rows0, sensIdx, swap.atoms, context));
  const plan = plans.get(key);
  const placed = placeSensitizer(swap.atoms, plan);
  const old = new Set(sensIdx);
  const map = new Map();
  const kept = [];
  rows0.forEach((r, i) => {
    if (old.has(i)) return;
    map.set(i, placed.length + kept.length);
    kept.push(r);
  });
  return {
    rows: [...placed, ...kept],
    sens: range(0, placed.length),
    core: range(0, placed.length).filter((i) => placed[i][0] !== 'H'),
    map,
    plan,
  };
}

/** Rows [[symbol, x, y, z]] turned k steps about z (twist per step) and lifted k rises: one helical screw step. */
function screwRows(list, k, twist, rise) {
  const c = Math.cos(k * twist);
  const s = Math.sin(k * twist);
  return list.map(([sym, x, y, z]) => [sym, [x * c - y * s, x * s + y * c, z + k * rise]]);
}

const sub3 = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dot3 = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const unit3 = (a) => {
  const l = Math.hypot(...a);
  return [a[0] / l, a[1] / l, a[2] / l];
};

const BUILDERS = {
  /**
   * Cy3 | n base pairs of B-DNA | Cy5 (tools/build_fret.py). The base pair repeats by its helical screw; Cy3 is locked
   * to the first pair, Cy5 to the last, so the angle between their transition dipoles (along the N–N axis of each
   * cyanine) turns with n. Reports the distance, κ² and the FRET efficiency, against the freely tumbling κ² = 2/3.
   */
  helix(d, n) {
    const atoms = [];
    for (let k = 0; k < n; k++) atoms.push(...screwRows(d.unit, k, d.twist, d.rise));
    const cy3 = screwRows(d.cy3.atoms, 0, d.twist, d.rise);
    const cy5 = screwRows(d.cy5.atoms, n - 1, d.twist, d.rise);
    const start3 = atoms.length;
    const start5 = start3 + cy3.length;
    atoms.push(...cy3, ...cy5);
    const dipole = (rows, [i, j]) => ({ centre: centroid([rows[i][1], rows[j][1]]), dir: unit3(sub3(rows[j][1], rows[i][1])) });
    const a = dipole(cy3, d.cy3.dipole);
    const b = dipole(cy5, d.cy5.dipole);
    const r = sub3(b.centre, a.centre);
    const dist = Math.hypot(...r) / 10; // nm
    const rh = unit3(r);
    const kappa = dot3(a.dir, b.dir) - 3 * dot3(a.dir, rh) * dot3(b.dir, rh);
    const k2 = kappa * kappa;
    const angle = (Math.acos(Math.min(1, Math.abs(dot3(a.dir, b.dir)))) * 180) / Math.PI;
    const fret = (R0) => 1 / (1 + (dist / R0) ** 6);
    const e = fret(d.r0 * (k2 / (2 / 3)) ** (1 / 6));
    const eIso = fret(d.r0);
    const pct = (x) => `${(100 * x).toFixed(0)}%`;
    return {
      atoms,
      sites: [
        { label: 'Cy3', kind: 'dye', color: '#ffc23d', atoms: range(start3, cy3.length) },
        { label: 'Cy5', kind: 'dye', color: '#ff3d4d', atoms: range(start5, cy5.length) },
      ],
      comment: `Cy3 | ${n} bp B-DNA | Cy5 | ${dist.toFixed(2)} nm, κ² = ${k2.toFixed(2)}, E = ${pct(e)} | ${d.method}`,
      note:
        `Cy3 ↔ Cy5: ${dist.toFixed(2)} nm, dipoles ${angle.toFixed(0)}° apart, κ² = ${k2.toFixed(2)}. ` +
        `FRET ${pct(e)} stacked as drawn, ${pct(eIso)} if the dyes tumbled freely (κ² = 2/3; R₀ = ${d.r0} nm assumed).`,
    };
  },

  /**
   * One finished geometry per value of n, each built and optimised offline as a whole molecule
   * (e.g. dye –(C6H4)n– quencher from tools/build_selfheal.py). Each frame brings its own sites and notes.
   */
  frames(d, n) {
    const f = d.frames[String(n)];
    if (!f) throw new Error(`no geometry for n = ${n}`);
    return { atoms: rows(f.atoms), sites: f.sites, comment: `${f.comment} | ${d.method}`, note: f.note };
  },

  /**
   * n DPA units, anthracene C2 -> 1,4-phenylene -> C6 of the next one, H on both ends.
   * Each DPA + linker is the same set of atoms shifted by t (both pieces are centrosymmetric).
   */
  'dpa-wire'(d, n) {
    const atoms = wireCore(d, n);
    const capL = d.capLeft;
    const capR = add(d.capRight, d.t, n - 1);
    atoms.push(['H', capL], ['H', capR]);
    const length = dist(capL, capR);
    return {
      atoms,
      sites: dpaSites(d, n),
      comment: `${n} × DPA-2,6-diyl joined by 1,4-phenylene, H-capped | ${length.toFixed(1)} Å end to end | ${d.method}`,
      note: `${(length / 10).toFixed(2)} nm end to end (H to H)`,
    };
  },

  /**
   * Sensitizer -(C6H4)- [DPA -(C6H4)-] x n - perylene. For n >= 1: the DPA wire core plus two end groups
   * (left in the first DPA's frame, right in the last DPA's frame). n = 0 is a separately built molecule.
   */
  'tta-flora'(d, n, load, { sensitizer: swap } = {}) {
    let atoms, sens, pery, sites, plan;
    if (n === 0) {
      let zero = d.zero;
      let zSens = d.zeroSens;
      let zCore = d.zeroLeftCore;
      let remap = (i) => i;
      if (swap) {
        const s = swapSensitizer(d.zero, d.zeroSens, swap, [], 'zero');
        ({ rows: zero, sens: zSens, core: zCore, plan } = s);
        remap = (i) => s.map.get(i);
      }
      atoms = rows(zero);
      sens = zCore.map((i) => atoms[i][1]);
      pery = d.zeroRightCore.map((i) => atoms[remap(i)][1]);
      sites = [
        { label: 'Sensitizer', kind: 'sensitizer', atoms: zSens },
        { label: 'Perylene', kind: 'emitter', atoms: d.zeroEmit.map(remap) },
      ];
    } else {
      const w = load(d.wire);
      const core = wireCore(w, n);
      let leftRows = d.left;
      let lSens = d.leftSens;
      let lCore = d.leftCore;
      if (swap) ({ rows: leftRows, sens: lSens, core: lCore, plan } = swapSensitizer(d.left, d.leftSens, swap, core, 'wire'));
      const left = rows(leftRows);
      const right = rows(d.right, n - 1, w.t);
      sens = lCore.map((i) => left[i][1]);
      pery = d.rightCore.map((i) => right[i][1]);
      atoms = [...left, ...core, ...right];
      sites = [
        { label: 'Sensitizer', kind: 'sensitizer', atoms: lSens },
        ...dpaSites(w, n, left.length),
        { label: 'Perylene', kind: 'emitter', atoms: d.rightEmit.map((i) => left.length + core.length + i) },
      ];
    }
    const sep = dist(centroid(sens), centroid(pery));
    const what = n === 0 ? 'no mediator (one linker)' : `${n} DPA mediator${n > 1 ? 's' : ''}`;
    const who = swap ? swap.name : 'Sensitizer';
    const method = swap ? `${swap.name} placed rigidly on the phenylene (src/attach.js); rest: ${d.method}` : d.method;
    return {
      atoms,
      sites,
      comment: `${who}–perylene, ${what} | ${(sep / 10).toFixed(2)} nm core to core | ${method}`,
      note: `${who} ↔ perylene: ${(sep / 10).toFixed(2)} nm (core centre to centre)`,
      attach: plan && { atom: plan.attachC + 1, symbol: swap.atoms[plan.attachC][0] },
    };
  },
};

/**
 * Build a parametric example. Returns { text (XYZ), note, sites, attach } or throws.
 * options.sensitizer: { key, name, atoms: [[symbol, [x, y, z]]] } swaps the sensitizer (TTA Flora).
 */
export function buildParametric(tab, param, value, options = {}) {
  const load = (file) => {
    const data = DATA[`${tab}/${file}`];
    if (!data) throw new Error(`geometry data ${tab}/${file} is missing`);
    return data;
  };
  const data = load(param.data);
  const build = BUILDERS[data.builder];
  if (!build) throw new Error(`unknown builder "${data.builder}"`);
  const { atoms, comment, note, sites, attach } = build(data, value, load, options);
  const lines = atoms.map(([s, p]) => `${s.padEnd(2)} ${p.map((v) => v.toFixed(4).padStart(11)).join(' ')}`);
  return { text: `${atoms.length}\n${comment}\n${lines.join('\n')}\n`, note, sites, attach };
}

export const clampParam = (param, v) => Math.min(param.max, Math.max(param.min, Math.round(Number(v) || param.min)));
