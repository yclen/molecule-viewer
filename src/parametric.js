// Parametric objects: molecules built in the app from a number the user sets (e.g. how many DPA units
// a molecular wire has). The geometry of each repeating piece is computed offline (tools/build_wire.py,
// saved as samples/<tab>/*.param.json); here we only stack copies, so any size builds instantly.

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
    label: n === 1 ? 'DPA' : `DPA ${k + 1}`,
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

const BUILDERS = {
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
  'tta-flora'(d, n, load) {
    let atoms, sens, pery, sites;
    if (n === 0) {
      atoms = rows(d.zero);
      sens = d.zeroLeftCore.map((i) => atoms[i][1]);
      pery = d.zeroRightCore.map((i) => atoms[i][1]);
      sites = [
        { label: 'Sensitizer', kind: 'sensitizer', atoms: d.zeroSens },
        { label: 'Perylene', kind: 'emitter', atoms: d.zeroEmit },
      ];
    } else {
      const w = load(d.wire);
      const left = rows(d.left);
      const core = wireCore(w, n);
      const right = rows(d.right, n - 1, w.t);
      sens = d.leftCore.map((i) => left[i][1]);
      pery = d.rightCore.map((i) => right[i][1]);
      atoms = [...left, ...core, ...right];
      sites = [
        { label: 'Sensitizer', kind: 'sensitizer', atoms: d.leftSens },
        ...dpaSites(w, n, left.length),
        { label: 'Perylene', kind: 'emitter', atoms: d.rightEmit.map((i) => left.length + core.length + i) },
      ];
    }
    const sep = dist(centroid(sens), centroid(pery));
    const what = n === 0 ? 'no mediator (one linker)' : `${n} DPA mediator${n > 1 ? 's' : ''}`;
    return {
      atoms,
      sites,
      comment: `Sensitizer–perylene, ${what} | ${(sep / 10).toFixed(2)} nm core to core | ${d.method}`,
      note: `Sensitizer ↔ perylene: ${(sep / 10).toFixed(2)} nm (core centre to centre)`,
    };
  },
};

/** Build a parametric example. Returns { text (XYZ), note, sites } or throws. */
export function buildParametric(tab, param, value) {
  const load = (file) => {
    const data = DATA[`${tab}/${file}`];
    if (!data) throw new Error(`geometry data ${tab}/${file} is missing`);
    return data;
  };
  const data = load(param.data);
  const build = BUILDERS[data.builder];
  if (!build) throw new Error(`unknown builder "${data.builder}"`);
  const { atoms, comment, note, sites } = build(data, value, load);
  const lines = atoms.map(([s, p]) => `${s.padEnd(2)} ${p.map((v) => v.toFixed(4).padStart(11)).join(' ')}`);
  return { text: `${atoms.length}\n${comment}\n${lines.join('\n')}\n`, note, sites };
}

export const clampParam = (param, v) => Math.min(param.max, Math.max(param.min, Math.round(Number(v) || param.min)));
