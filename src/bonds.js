import { getElement } from './elements.js';

/**
 * XYZ files carry no bond information, so infer it from distances:
 * two atoms are bonded when  d < r_cov(a) + r_cov(b) + tolerance.
 * (Same rule Open Babel / Avogadro use.)
 *
 * A uniform spatial grid keeps this O(n) instead of O(n²), so large
 * structures (tens of thousands of atoms) stay fast.
 *
 * Returns a flat Uint32Array of index pairs: [i0, j0, i1, j1, ...].
 */
export function detectBonds(atoms, tolerance = 0.4) {
  const n = atoms.length;
  if (n < 2) return new Uint32Array(0);

  const rad = new Float32Array(n);
  let maxR = 0;
  for (let i = 0; i < n; i++) {
    rad[i] = getElement(atoms[i].symbol).covalent;
    if (rad[i] > maxR) maxR = rad[i];
  }
  if (maxR <= 0) return new Uint32Array(0);

  const cell = 2 * maxR + tolerance; // longest possible bond
  const grid = new Map();
  const cx = new Int32Array(n);
  const cy = new Int32Array(n);
  const cz = new Int32Array(n);

  for (let i = 0; i < n; i++) {
    const a = atoms[i];
    cx[i] = Math.floor(a.x / cell);
    cy[i] = Math.floor(a.y / cell);
    cz[i] = Math.floor(a.z / cell);
    const key = `${cx[i]},${cy[i]},${cz[i]}`;
    let bucket = grid.get(key);
    if (!bucket) grid.set(key, (bucket = []));
    bucket.push(i);
  }

  const MIN_D2 = 0.1 * 0.1; // atoms this close are duplicates, not bonds
  const out = [];

  for (let i = 0; i < n; i++) {
    if (rad[i] === 0) continue;
    const a = atoms[i];
    for (let dx = -1; dx <= 1; dx++) {
      for (let dy = -1; dy <= 1; dy++) {
        for (let dz = -1; dz <= 1; dz++) {
          const bucket = grid.get(`${cx[i] + dx},${cy[i] + dy},${cz[i] + dz}`);
          if (!bucket) continue;
          for (const j of bucket) {
            if (j <= i || rad[j] === 0) continue;
            const b = atoms[j];
            const ex = a.x - b.x;
            const ey = a.y - b.y;
            const ez = a.z - b.z;
            const d2 = ex * ex + ey * ey + ez * ez;
            const cut = rad[i] + rad[j] + tolerance;
            if (d2 < cut * cut && d2 > MIN_D2) out.push(i, j);
          }
        }
      }
    }
  }

  return Uint32Array.from(out);
}
