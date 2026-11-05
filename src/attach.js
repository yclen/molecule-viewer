// Swap the sensitizer of an end group for another molecule (e.g. a different TADF sensitizer).
//
// The end group is [sensitizer + 1,4-phenylene] from tools/build_triplet.py. We drop the sensitizer atoms,
// which leaves the phenylene with one free ipso carbon, and bolt the new sensitizer on there:
//   1. candidate attachment points are the sensitizer's aromatic C–H groups (sp2 C with three neighbours);
//      the H is removed and the C bonds to the phenylene (aryl–aryl 1.48 Å, along the phenylene axis);
//   2. favoured are C–H groups on the outside of the molecule, pointing away from it, so the sensitizer
//      extends away from the wire like the original one does;
//   3. for each candidate the biaryl torsion is scanned, and the placement with the fewest close contacts
//      to the rest of the molecule (and a twist near 45°) wins.
// It's a rigid-body placement of the sensitizer's own optimised geometry (from its sample .xyz), not a new
// force-field optimisation, so it is meant for visualising, not for distances to better than ~0.5 Å.

import { detectBonds } from './bonds.js';

const BIARYL = 1.48; // Å, aryl–aryl single bond
const CANDIDATES = 8; // most outward C–H groups to try
const STEPS = 72; // torsion scan, 5° steps

const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const scale = (a, k) => [a[0] * k, a[1] * k, a[2] * k];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const unit = (a) => scale(a, 1 / Math.hypot(...a));
const dist2 = (a, b) => (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2;

/** Rotation matrix (rows) turning unit vector a onto unit vector b. */
function alignRotation(a, b) {
  const c = dot(a, b);
  if (c < -0.999999) {
    // opposite: half turn about any axis perpendicular to a
    const p = unit(Math.abs(a[0]) < 0.9 ? cross(a, [1, 0, 0]) : cross(a, [0, 1, 0]));
    return axisRotation(p, Math.PI);
  }
  const v = cross(a, b);
  const k = 1 / (1 + c);
  return [
    [v[0] * v[0] * k + c, v[0] * v[1] * k - v[2], v[0] * v[2] * k + v[1]],
    [v[1] * v[0] * k + v[2], v[1] * v[1] * k + c, v[1] * v[2] * k - v[0]],
    [v[2] * v[0] * k - v[1], v[2] * v[1] * k + v[0], v[2] * v[2] * k + c],
  ];
}

function axisRotation(u, t) {
  const c = Math.cos(t);
  const s = Math.sin(t);
  const C = 1 - c;
  const [x, y, z] = u;
  return [
    [c + x * x * C, x * y * C - z * s, x * z * C + y * s],
    [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
    [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
  ];
}

const apply = (R, p) => [dot(R[0], p), dot(R[1], p), dot(R[2], p)];
const mul = (A, B) => A.map((r) => [0, 1, 2].map((j) => r[0] * B[0][j] + r[1] * B[1][j] + r[2] * B[2][j]));

function neighbours(atoms) {
  const bonds = detectBonds(atoms.map(([symbol, [x, y, z]]) => ({ symbol, x, y, z })));
  const nb = atoms.map(() => []);
  for (let b = 0; b < bonds.length; b += 2) {
    nb[bonds[b]].push(bonds[b + 1]);
    nb[bonds[b + 1]].push(bonds[b]);
  }
  return nb;
}

/** Closest-contact penalty, Å²: heavy–heavy below 3.0 Å, heavy–H below 2.4 Å, H–H below 2.2 Å. */
function clash(placed, isH, skip, ctx) {
  let pen = 0;
  for (let i = 0; i < placed.length; i++) {
    const p = placed[i];
    for (const c of ctx) {
      if (skip[i] && c.near) continue; // atoms around the new bond: their distances are set by the bond itself
      const rmin = isH[i] ? (c.h ? 2.2 : 2.4) : c.h ? 2.4 : 3.0;
      const d2 = dist2(p, c.p);
      if (d2 < rmin * rmin) pen += (rmin - Math.sqrt(d2)) ** 2;
    }
  }
  return pen;
}

/**
 * Plan how to attach `sens` (rows [symbol, [x, y, z]]) in place of the atoms `oldSens` of `rows`
 * (rows [symbol, x, y, z], the end group), with `context` (rows [symbol, [x, y, z]], same frame) to avoid.
 * Returns a placement { R, from, to, drop, attachC } to reuse with placeSensitizer.
 */
export function planAttachment(rows, oldSens, sens, context) {
  const old = new Set(oldSens);
  const pos = (r) => [r[1], r[2], r[3]];
  const keep = rows.map((_, i) => i).filter((i) => !old.has(i));
  const keptC = keep.filter((i) => rows[i][0] === 'C');
  // the phenylene ipso carbon that held the old sensitizer: the kept C closest to any old heavy atom
  const oldHeavy = oldSens.filter((i) => rows[i][0] !== 'H').map((i) => pos(rows[i]));
  const near = (i) => Math.min(...oldHeavy.map((q) => dist2(pos(rows[i]), q)));
  const ipso = keptC.reduce((a, b) => (near(a) < near(b) ? a : b));
  const para = keptC
    .filter((i) => dist2(pos(rows[i]), pos(rows[ipso])) < 3.5 * 3.5) // same ring
    .reduce((a, b) => (dist2(pos(rows[a]), pos(rows[ipso])) > dist2(pos(rows[b]), pos(rows[ipso])) ? a : b));
  const Ci = pos(rows[ipso]);
  const u = unit(sub(Ci, pos(rows[para]))); // phenylene axis, pointing out of the ring at the free carbon
  const target = add(Ci, scale(u, BIARYL));
  const ortho = keptC.filter((i) => Math.abs(Math.sqrt(dist2(pos(rows[i]), Ci)) - 1.4) < 0.15);
  const ringNormal = unit(cross(sub(pos(rows[ortho[0]]), Ci), sub(pos(rows[ortho[1]]), Ci)));

  // atoms to stay clear of: the kept end group atoms plus the context
  const nearSet = new Set([ipso, ...ortho]);
  const ctx = [
    ...keep.map((i) => ({ p: pos(rows[i]), h: rows[i][0] === 'H', near: nearSet.has(i) })),
    ...context.map(([s, p]) => ({ p, h: s === 'H', near: false })),
  ];

  const nb = neighbours(sens);
  const isH = sens.map(([s]) => s === 'H');
  const heavy = sens.filter(([s]) => s !== 'H').map(([, p]) => p);
  const centroid = scale(heavy.reduce(add, [0, 0, 0]), 1 / heavy.length);
  const reach = Math.sqrt(Math.max(...sens.map(([, p]) => dist2(p, centroid)))) * 2 + 4;

  const cands = [];
  sens.forEach(([s, p], h) => {
    if (s !== 'H' || nb[h].length !== 1) return;
    const c = nb[h][0];
    if (sens[c][0] !== 'C' || nb[c].length !== 3) return;
    const v = unit(sub(p, sens[c][1]));
    cands.push({ h, c, v, out: dot(sub(p, centroid), v) });
  });
  if (!cands.length) throw new Error('no aromatic C–H to attach through');
  cands.sort((a, b) => b.out - a.out);
  const best = { score: Infinity };
  const outMax = cands[0].out;
  const localCtx = ctx.filter((c) => dist2(c.p, target) < reach * reach);

  for (const cand of cands.slice(0, CANDIDATES)) {
    const C0 = sens[cand.c][1];
    const R0 = alignRotation(cand.v, scale(u, -1)); // the removed H pointed at the phenylene
    const idx = sens.map((_, i) => i).filter((i) => i !== cand.h);
    const local = idx.map((i) => apply(R0, sub(sens[i][1], C0)));
    const lH = idx.map((i) => isH[i]);
    const skip = idx.map((i) => i === cand.c || nb[cand.c].includes(i));
    const nbr = nb[cand.c].filter((i) => i !== cand.h);
    for (let k = 0; k < STEPS; k++) {
      const Rt = axisRotation(u, (2 * Math.PI * k) / STEPS);
      const placed = local.map((p) => add(apply(Rt, p), target));
      // biaryl twist: angle between the phenylene plane and the plane at the sensitizer carbon
      const a = apply(Rt, apply(R0, sub(sens[nbr[0]][1], C0)));
      const b = apply(Rt, apply(R0, sub(sens[nbr[1]][1], C0)));
      const twist = Math.abs(dot(unit(cross(a, b)), ringNormal));
      // contacts dominate; between equally clear placements prefer outward C–H groups and a ~45° twist
      const score =
        4 * clash(placed, lH, skip, localCtx) + 0.3 * (outMax - cand.out) + 0.1 * (twist - Math.SQRT1_2) ** 2;
      if (score < best.score) Object.assign(best, { score, R: mul(Rt, R0), from: C0, to: target, drop: cand.h, attachC: cand.c });
    }
  }
  return best;
}

/** New sensitizer rows [symbol, x, y, z] for a placement from planAttachment. */
export function placeSensitizer(sens, plan) {
  return sens
    .filter((_, i) => i !== plan.drop)
    .map(([s, p]) => [s, ...add(apply(plan.R, sub(p, plan.from)), plan.to).map((x) => Math.round(x * 1e4) / 1e4)]);
}
