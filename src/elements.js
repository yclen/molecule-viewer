// Per-element data, indexed by atomic number (Z = row index + 1).
//   covalent: single-bond covalent radius in Å (Cordero et al. 2008) — used for bond detection
//   vdw:      van der Waals radius in Å (Bondi / Alvarez) — used for space-filling & ball sizes
//   color:    Jmol / CPK colour
const DATA = [
  // sym  covalent vdw   color
  ['H',  0.31, 1.20, 0xffffff],
  ['He', 0.28, 1.40, 0xd9ffff],
  ['Li', 1.28, 1.82, 0xcc80ff],
  ['Be', 0.96, 1.53, 0xc2ff00],
  ['B',  0.84, 1.92, 0xffb5b5],
  ['C',  0.76, 1.70, 0x909090],
  ['N',  0.71, 1.55, 0x3050f8],
  ['O',  0.66, 1.52, 0xff0d0d],
  ['F',  0.57, 1.47, 0x90e050],
  ['Ne', 0.58, 1.54, 0xb3e3f5],
  ['Na', 1.66, 2.27, 0xab5cf2],
  ['Mg', 1.41, 1.73, 0x8aff00],
  ['Al', 1.21, 1.84, 0xbfa6a6],
  ['Si', 1.11, 2.10, 0xf0c8a0],
  ['P',  1.07, 1.80, 0xff8000],
  ['S',  1.05, 1.80, 0xffff30],
  ['Cl', 1.02, 1.75, 0x1ff01f],
  ['Ar', 1.06, 1.88, 0x80d1e3],
  ['K',  2.03, 2.75, 0x8f40d4],
  ['Ca', 1.76, 2.31, 0x3dff00],
  ['Sc', 1.70, 2.15, 0xe6e6e6],
  ['Ti', 1.60, 2.11, 0xbfc2c7],
  ['V',  1.53, 2.07, 0xa6a6ab],
  ['Cr', 1.39, 2.06, 0x8a99c7],
  ['Mn', 1.39, 2.05, 0x9c7ac7],
  ['Fe', 1.32, 2.04, 0xe06633],
  ['Co', 1.26, 2.00, 0xf090a0],
  ['Ni', 1.24, 1.63, 0x50d050],
  ['Cu', 1.32, 1.40, 0xc88033],
  ['Zn', 1.22, 1.39, 0x7d80b0],
  ['Ga', 1.22, 1.87, 0xc28f8f],
  ['Ge', 1.20, 2.11, 0x668f8f],
  ['As', 1.19, 1.85, 0xbd80e3],
  ['Se', 1.20, 1.90, 0xffa100],
  ['Br', 1.20, 1.85, 0xa62929],
  ['Kr', 1.16, 2.02, 0x5cb8d1],
  ['Rb', 2.20, 3.03, 0x702eb0],
  ['Sr', 1.95, 2.49, 0x00ff00],
  ['Y',  1.90, 2.32, 0x94ffff],
  ['Zr', 1.75, 2.23, 0x94e0e0],
  ['Nb', 1.64, 2.18, 0x73c2c9],
  ['Mo', 1.54, 2.17, 0x54b5b5],
  ['Tc', 1.47, 2.16, 0x3b9e9e],
  ['Ru', 1.46, 2.13, 0x248f8f],
  ['Rh', 1.42, 2.10, 0x0a7d8c],
  ['Pd', 1.39, 1.63, 0x006985],
  ['Ag', 1.45, 1.72, 0xc0c0c0],
  ['Cd', 1.44, 1.58, 0xffd98f],
  ['In', 1.42, 1.93, 0xa67573],
  ['Sn', 1.39, 2.17, 0x668080],
  ['Sb', 1.39, 2.06, 0x9e63b5],
  ['Te', 1.38, 2.06, 0xd47a00],
  ['I',  1.39, 1.98, 0x940094],
  ['Xe', 1.40, 2.16, 0x429eb0],
  ['Cs', 2.44, 3.43, 0x57178f],
  ['Ba', 2.15, 2.68, 0x00c900],
  ['La', 2.07, 2.43, 0x70d4ff],
  ['Ce', 2.04, 2.42, 0xffffc7],
  ['Pr', 2.03, 2.40, 0xd9ffc7],
  ['Nd', 2.01, 2.39, 0xc7ffc7],
  ['Pm', 1.99, 2.38, 0xa3ffc7],
  ['Sm', 1.98, 2.36, 0x8fffc7],
  ['Eu', 1.98, 2.35, 0x61ffc7],
  ['Gd', 1.96, 2.34, 0x45ffc7],
  ['Tb', 1.94, 2.33, 0x30ffc7],
  ['Dy', 1.92, 2.31, 0x1fffc7],
  ['Ho', 1.92, 2.30, 0x00ff9c],
  ['Er', 1.89, 2.29, 0x00e675],
  ['Tm', 1.90, 2.27, 0x00d452],
  ['Yb', 1.87, 2.26, 0x00bf38],
  ['Lu', 1.87, 2.24, 0x00ab24],
  ['Hf', 1.75, 2.23, 0x4dc2ff],
  ['Ta', 1.70, 2.22, 0x4da6ff],
  ['W',  1.62, 2.18, 0x2194d6],
  ['Re', 1.51, 2.16, 0x267dab],
  ['Os', 1.44, 2.16, 0x266696],
  ['Ir', 1.41, 2.13, 0x175487],
  ['Pt', 1.36, 1.75, 0xd0d0e0],
  ['Au', 1.36, 1.66, 0xffd123],
  ['Hg', 1.32, 1.55, 0xb8b8d0],
  ['Tl', 1.45, 1.96, 0xa6544d],
  ['Pb', 1.46, 2.02, 0x575961],
  ['Bi', 1.48, 2.07, 0x9e4fb5],
  ['Po', 1.40, 1.97, 0xab5c00],
  ['At', 1.50, 2.02, 0x754f45],
  ['Rn', 1.50, 2.20, 0x428296],
];

export const ELEMENTS = {};
export const SYMBOLS = [null]; // SYMBOLS[Z] -> symbol

DATA.forEach(([symbol, covalent, vdw, color], i) => {
  ELEMENTS[symbol] = { Z: i + 1, symbol, covalent, vdw, color };
  SYMBOLS.push(symbol);
});

// Dummy / unknown atoms: drawn in hot pink, never bonded.
const UNKNOWN = { Z: 0, symbol: 'X', covalent: 0, vdw: 1.0, color: 0xff1493 };

export function getElement(symbol) {
  return ELEMENTS[symbol] ?? UNKNOWN;
}

/**
 * Turn whatever is in the first column of an XYZ line into an element symbol.
 * Handles "C", "c", "CL", "Cl", "C12" (labels), and atomic numbers like "6".
 */
export function normalizeSymbol(token) {
  if (/^\d+$/.test(token)) return SYMBOLS[Number(token)] ?? 'X';
  const letters = token.match(/^[A-Za-z]+/)?.[0];
  if (!letters) return 'X';
  if (letters.length >= 2) {
    const two = letters[0].toUpperCase() + letters[1].toLowerCase();
    if (ELEMENTS[two]) return two;
  }
  const one = letters[0].toUpperCase();
  return ELEMENTS[one] ? one : 'X';
}

/** Hill-order formula: C first, then H, then the rest alphabetically. Returns [[symbol, count], ...]. */
export function hillFormula(atoms) {
  const counts = new Map();
  for (const a of atoms) counts.set(a.symbol, (counts.get(a.symbol) ?? 0) + 1);
  const keys = [...counts.keys()].sort();
  const order = counts.has('C')
    ? ['C', ...(counts.has('H') ? ['H'] : []), ...keys.filter((k) => k !== 'C' && k !== 'H')]
    : keys;
  return order.map((k) => [k, counts.get(k)]);
}
