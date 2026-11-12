// Sparse linear algebra for the exact solvers, standing in for what sim/engine.py takes from scipy: compressed-row
// matrices, an incomplete LU with no fill-in (ILU(0)) as the preconditioner, restarted GMRES, and the matrix
// exponential times a vector by uniformization.

/** Collects (row, col, value) entries. Duplicates are added together when it becomes a matrix. */
export class Coo {
  constructor(n, guess = 1024) {
    this.n = n;
    this.len = 0;
    this.rows = new Int32Array(guess);
    this.cols = new Int32Array(guess);
    this.vals = new Float64Array(guess);
  }

  push(r, c, v) {
    if (this.len === this.rows.length) {
      const grow = (old, Type) => {
        const bigger = new Type(old.length * 2);
        bigger.set(old);
        return bigger;
      };
      this.rows = grow(this.rows, Int32Array);
      this.cols = grow(this.cols, Int32Array);
      this.vals = grow(this.vals, Float64Array);
    }
    this.rows[this.len] = r;
    this.cols[this.len] = c;
    this.vals[this.len++] = v;
  }

  /** A compressed-row matrix with each row's columns sorted and duplicates summed. */
  toCsr() {
    const { n, len, rows, cols, vals } = this;
    const ptr = new Int32Array(n + 1);
    for (let i = 0; i < len; i++) ptr[rows[i] + 1]++;
    for (let i = 0; i < n; i++) ptr[i + 1] += ptr[i];
    const next = ptr.slice(0, n);
    const col = new Int32Array(len);
    const val = new Float64Array(len);
    for (let i = 0; i < len; i++) {
      const p = next[rows[i]]++;
      col[p] = cols[i];
      val[p] = vals[i];
    }
    // sort each row by column (rows are short) and merge duplicates, compacting as we go
    const outPtr = new Int32Array(n + 1);
    let w = 0;
    for (let i = 0; i < n; i++) {
      const start = ptr[i];
      const end = ptr[i + 1];
      if (end - start > 32) {
        // a long row (the rows that count events): sort an index, then rewrite the row in order
        const order = Array.from({ length: end - start }, (_, q) => start + q).sort((u, v) => col[u] - col[v]);
        const c = order.map((q) => col[q]);
        const v = order.map((q) => val[q]);
        for (let q = 0; q < c.length; q++) {
          col[start + q] = c[q];
          val[start + q] = v[q];
        }
      } else {
        for (let p = start + 1; p < end; p++) {
          const c = col[p];
          const v = val[p];
          let q = p - 1;
          while (q >= start && col[q] > c) {
            col[q + 1] = col[q];
            val[q + 1] = val[q];
            q--;
          }
          col[q + 1] = c;
          val[q + 1] = v;
        }
      }
      outPtr[i] = w;
      for (let p = start; p < end; p++) {
        if (w > outPtr[i] && col[w - 1] === col[p]) val[w - 1] += val[p];
        else {
          col[w] = col[p];
          val[w++] = val[p];
        }
      }
    }
    outPtr[n] = w;
    return { n, ptr: outPtr, col: col.slice(0, w), val: val.slice(0, w) };
  }
}

/** y = A x (y is overwritten). */
export function matvec(A, x, y = new Float64Array(A.n)) {
  const { n, ptr, col, val } = A;
  for (let i = 0; i < n; i++) {
    let s = 0;
    for (let p = ptr[i]; p < ptr[i + 1]; p++) s += val[p] * x[col[p]];
    y[i] = s;
  }
  return y;
}

/**
 * ILU(0): L and U with A's own sparsity, stored in one array. Every row needs its diagonal entry. Returns
 * solve(b, out), which applies (LU)^-1.
 */
export function ilu0(A) {
  const { n, ptr, col } = A;
  const lu = A.val.slice();
  const diag = new Int32Array(n);
  const iw = new Int32Array(n).fill(-1);
  for (let i = 0; i < n; i++) {
    for (let p = ptr[i]; p < ptr[i + 1]; p++) iw[col[p]] = p;
    let p = ptr[i];
    for (; p < ptr[i + 1] && col[p] < i; p++) {
      const k = col[p];
      const lik = (lu[p] /= lu[diag[k]]);
      for (let q = diag[k] + 1; q < ptr[k + 1]; q++) {
        const at = iw[col[q]];
        if (at !== -1) lu[at] -= lik * lu[q];
      }
    }
    if (p === ptr[i + 1] || col[p] !== i) throw new Error(`ILU(0): row ${i} has no diagonal entry`);
    diag[i] = p;
    if (lu[p] === 0) throw new Error(`ILU(0): zero pivot in row ${i}`);
    for (let q = ptr[i]; q < ptr[i + 1]; q++) iw[col[q]] = -1;
  }
  return (b, x = new Float64Array(n)) => {
    for (let i = 0; i < n; i++) {
      let s = b[i];
      for (let p = ptr[i]; p < diag[i]; p++) s -= lu[p] * x[col[p]];
      x[i] = s;
    }
    for (let i = n - 1; i >= 0; i--) {
      let s = x[i];
      for (let p = diag[i] + 1; p < ptr[i + 1]; p++) s -= lu[p] * x[col[p]];
      x[i] = s / lu[diag[i]];
    }
    return x;
  };
}

const dot = (u, v) => {
  let s = 0;
  for (let i = 0; i < u.length; i++) s += u[i] * v[i];
  return s;
};

/**
 * Solve A x = b by restarted GMRES, right-preconditioned with `precond` (so the residual it tests is A x - b itself).
 * Stops when |A x - b| <= rtol |b|. Throws if it hasn't converged after `maxRestarts` restarts.
 */
export function gmres(A, b, x0, precond, { rtol = 1e-8, restart = 40, maxRestarts = 50 } = {}) {
  const n = A.n;
  const x = x0 ? Float64Array.from(x0) : new Float64Array(n);
  const bnorm = Math.sqrt(dot(b, b));
  if (bnorm === 0) return new Float64Array(n);
  const tol = rtol * bnorm;
  const m = restart;
  const V = Array.from({ length: m + 1 }, () => new Float64Array(n));
  const H = Array.from({ length: m + 1 }, () => new Float64Array(m));
  const cs = new Float64Array(m);
  const sn = new Float64Array(m);
  const g = new Float64Array(m + 1);
  const z = new Float64Array(n);
  const w = new Float64Array(n);
  for (let outer = 0; outer <= maxRestarts; outer++) {
    matvec(A, x, w);
    const r = V[0];
    for (let i = 0; i < n; i++) r[i] = b[i] - w[i];
    const beta = Math.sqrt(dot(r, r));
    if (beta <= tol) return x;
    if (outer === maxRestarts) break;
    for (let i = 0; i < n; i++) r[i] /= beta;
    g.fill(0);
    g[0] = beta;
    let j = 0;
    for (; j < m; j++) {
      precond(V[j], z);
      matvec(A, z, w);
      for (let i = 0; i <= j; i++) {
        const h = dot(w, V[i]);
        H[i][j] = h;
        const Vi = V[i];
        for (let q = 0; q < n; q++) w[q] -= h * Vi[q];
      }
      const hn = Math.sqrt(dot(w, w));
      H[j + 1][j] = hn;
      if (hn > 0) {
        const Vn = V[j + 1];
        for (let q = 0; q < n; q++) Vn[q] = w[q] / hn;
      }
      for (let i = 0; i < j; i++) {
        const t = cs[i] * H[i][j] + sn[i] * H[i + 1][j];
        H[i + 1][j] = -sn[i] * H[i][j] + cs[i] * H[i + 1][j];
        H[i][j] = t;
      }
      const d = Math.hypot(H[j][j], H[j + 1][j]);
      cs[j] = H[j][j] / d;
      sn[j] = H[j + 1][j] / d;
      H[j][j] = d;
      H[j + 1][j] = 0;
      g[j + 1] = -sn[j] * g[j];
      g[j] *= cs[j];
      if (Math.abs(g[j + 1]) <= tol || hn === 0) {
        j++;
        break;
      }
    }
    // y = H^-1 g, then x += M^-1 (V y)
    const y = new Float64Array(j);
    for (let i = j - 1; i >= 0; i--) {
      let s = g[i];
      for (let k = i + 1; k < j; k++) s -= H[i][k] * y[k];
      y[i] = s / H[i][i];
    }
    w.fill(0);
    for (let i = 0; i < j; i++) {
      const Vi = V[i];
      for (let q = 0; q < n; q++) w[q] += y[i] * Vi[q];
    }
    precond(w, z);
    for (let q = 0; q < n; q++) x[q] += z[q];
  }
  throw new Error('the exact solve did not converge (GMRES)');
}

/**
 * exp(t A) v for a matrix whose off-diagonal entries are all >= 0 (a rate matrix, plus rows that count events), by
 * uniformization: with Λ at least every -A_ii, P = I + A / Λ has no negative entries and
 * exp(t A) v = Σ_k Poisson(k; Λt) P^k v. Nothing is subtracted, so even very small counts come out accurate.
 * Long times are split into pieces with Λ·dt <= 30, so the Poisson weights stay in range.
 */
export function expmv(A, t, v) {
  const { n, ptr, col, val } = A;
  let lambda = 0;
  for (let i = 0; i < n; i++)
    for (let p = ptr[i]; p < ptr[i + 1]; p++) if (col[p] === i && -val[p] > lambda) lambda = -val[p];
  let x = Float64Array.from(v);
  if (lambda === 0 || t === 0) return x;
  const pieces = Math.ceil((lambda * t) / 30);
  const L = (lambda * t) / pieces;
  const most = Math.ceil(L + 10 * Math.sqrt(L) + 20);
  const term = new Float64Array(n);
  const Aterm = new Float64Array(n);
  for (let piece = 0; piece < pieces; piece++) {
    term.set(x);
    let coef = Math.exp(-L);
    const out = new Float64Array(n);
    for (let i = 0; i < n; i++) out[i] = coef * term[i];
    for (let k = 1; k <= most; k++) {
      matvec(A, term, Aterm);
      for (let i = 0; i < n; i++) term[i] += Aterm[i] / lambda;
      coef *= L / k;
      for (let i = 0; i < n; i++) out[i] += coef * term[i];
      if (k > L && coef < 1e-18) break;
    }
    x = out;
  }
  return x;
}
