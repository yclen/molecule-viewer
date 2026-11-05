// Where an excitation sits on a chromophore, for the excited-state glow.
//
// The glow weight of each atom is its share of the lowest π→π* excitation in simple Hückel theory:
//   w_i = (c_HOMO,i² + c_LUMO,i²) / 2      (average of where the hole and the excited electron are)
// Resonance integrals are scaled by the overlap of the two p orbitals, β_ij = k·β·|n_i · n_j|, using the local
// ring normals from the 3D structure. So twisted rings drop out by themselves: in DPA the phenyls (≈90°) carry
// almost nothing and the excitation sits on the anthracene, strongest at C9/C10 and C1/4/5/8, as in the real S1.
// It's a qualitative picture (no electron correlation, no charge transfer energetics), meant for visualisation.

// Hückel parameters (h: α_X = α + hβ, k: β_CX = kβ, e: π electrons), Streitwieser-style values.
// Atoms that can't be part of a π system (sp3 C, sulfone S, halogens, H) are left out.
function piParams(symbol, degree) {
  switch (symbol) {
    case 'C': return degree <= 3 ? { h: 0, k: 1, e: 1 } : null;
    case 'N': return degree === 3 ? { h: 1.5, k: 0.8, e: 2 } : degree <= 2 ? { h: 0.5, k: 1, e: 1 } : null;
    case 'O': return degree === 2 ? { h: 2, k: 0.8, e: 2 } : degree === 1 ? { h: 1, k: 1, e: 1 } : null;
    case 'S': return degree === 2 ? { h: 1, k: 0.6, e: 2 } : null;
    case 'B': return degree === 3 ? { h: -1, k: 0.7, e: 0 } : null;
    case 'Se': return degree === 2 ? { h: 0.8, k: 0.5, e: 2 } : null;
    default: return null;
  }
}

const MAX_PI = 160; // Jacobi is O(n³) per sweep; bigger sites fall back to an even glow

/**
 * Glow weights for one site.
 * atoms: [{symbol, x, y, z}], bonds: flat [i0, j0, i1, j1, ...], site: atom indices.
 * Returns { atoms: [index], weights: [0..1] } (max weight = 1) and the method used.
 */
export function excitationWeights(atoms, bonds, site) {
  const nbrs = atoms.map(() => []);
  for (let b = 0; b < bonds.length; b += 2) {
    nbrs[bonds[b]].push(bonds[b + 1]);
    nbrs[bonds[b + 1]].push(bonds[b]);
  }
  const inSite = new Set(site);
  const cand = new Map();
  for (const i of site) {
    const p = piParams(atoms[i].symbol, nbrs[i].length);
    if (p) cand.set(i, p);
  }
  // An atom only joins the π system through a π neighbour: drop isolated candidates (e.g. the O of a sulfone,
  // whose S is not conjugated), repeating until nothing changes.
  for (let changed = true; changed; ) {
    changed = false;
    for (const i of [...cand.keys()]) {
      if (!nbrs[i].some((j) => cand.has(j) && inSite.has(j))) {
        cand.delete(i);
        changed = true;
      }
    }
  }
  const pi = [...cand.keys()];
  const par = [...cand.values()];
  if (pi.length < 2 || pi.length > MAX_PI) return evenGlow(atoms, site);

  // p-orbital direction = normal of the plane through the atom and two of its neighbours
  const pos = (i) => [atoms[i].x, atoms[i].y, atoms[i].z];
  const normal = new Map();
  for (const i of pi) {
    if (nbrs[i].length < 2) continue;
    const [a, b] = nbrs[i];
    const o = pos(i);
    const u = sub(pos(a), o);
    const v = sub(pos(b), o);
    const n = cross(u, v);
    const len = Math.hypot(...n);
    if (len > 1e-6) normal.set(i, n.map((x) => x / len));
  }
  for (const i of pi) if (!normal.has(i) && nbrs[i].length === 1) normal.set(i, normal.get(nbrs[i][0]) ?? [0, 0, 1]);

  // Hamiltonian in units of β (β < 0, so H = -h on the diagonal, -k·cos on bonds): eigenvalues ascending = filling order
  const n = pi.length;
  const idx = new Map(pi.map((a, k) => [a, k]));
  const H = Array.from({ length: n }, () => new Float64Array(n));
  for (let k = 0; k < n; k++) H[k][k] = -par[k].h;
  for (const i of pi) {
    for (const j of nbrs[i]) {
      if (!inSite.has(j) || !idx.has(j) || j < i) continue;
      const a = idx.get(i);
      const b = idx.get(j);
      const overlap = Math.abs(dot(normal.get(i) ?? [0, 0, 1], normal.get(j) ?? [0, 0, 1]));
      const beta = -Math.min(par[a].k, par[b].k) * overlap;
      H[a][b] = H[b][a] = beta;
    }
  }
  const { values, vectors } = jacobi(H);
  const order = values.map((_, k) => k).sort((a, b) => values[a] - values[b]);
  const electrons = par.reduce((s, p) => s + p.e, 0);
  const homo = Math.ceil(electrons / 2) - 1;
  const lumo = homo + 1;
  if (homo < 0 || lumo >= n) return evenGlow(atoms, site);

  const w = new Float64Array(n);
  for (let k = 0; k < n; k++) w[k] = (vectors[k][order[homo]] ** 2 + vectors[k][order[lumo]] ** 2) / 2;
  const max = Math.max(...w);
  return {
    atoms: pi,
    weights: Array.from(w, (x) => x / max),
    method: 'Hückel HOMO→LUMO',
    gap: values[order[lumo]] - values[order[homo]], // in |β|
  };
}

function evenGlow(atoms, site) {
  const heavy = site.filter((i) => atoms[i].symbol !== 'H');
  return { atoms: heavy, weights: heavy.map(() => 1), method: 'even (no π system)' };
}

const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];

/** Symmetric eigen-decomposition (cyclic Jacobi). vectors[row][col]: column col is an eigenvector. */
function jacobi(A0) {
  const n = A0.length;
  const A = A0.map((r) => Float64Array.from(r));
  const V = Array.from({ length: n }, (_, i) => {
    const r = new Float64Array(n);
    r[i] = 1;
    return r;
  });
  for (let sweep = 0; sweep < 60; sweep++) {
    let off = 0;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += A[p][q] ** 2;
    if (off < 1e-20) break;
    for (let p = 0; p < n; p++) {
      for (let q = p + 1; q < n; q++) {
        if (Math.abs(A[p][q]) < 1e-14) continue;
        const theta = (A[q][q] - A[p][p]) / (2 * A[p][q]);
        const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const c = 1 / Math.sqrt(t * t + 1);
        const s = t * c;
        for (let k = 0; k < n; k++) {
          const akp = A[k][p];
          const akq = A[k][q];
          A[k][p] = c * akp - s * akq;
          A[k][q] = s * akp + c * akq;
        }
        for (let k = 0; k < n; k++) {
          const apk = A[p][k];
          const aqk = A[q][k];
          A[p][k] = c * apk - s * aqk;
          A[q][k] = s * apk + c * aqk;
        }
        for (let k = 0; k < n; k++) {
          const vkp = V[k][p];
          const vkq = V[k][q];
          V[k][p] = c * vkp - s * vkq;
          V[k][q] = s * vkp + c * vkq;
        }
      }
    }
  }
  return { values: A.map((r, i) => r[i]), vectors: V };
}
