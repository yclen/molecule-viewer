import { normalizeSymbol } from './elements.js';

// Accept Fortran-style exponents too (1.0D-3).
const toNumber = (s) => Number.parseFloat(s.replace(/[dD]/, 'e'));

function parseAtomLine(line, lineNo) {
  const t = line.trim().split(/\s+/);
  const preview = line.trim().slice(0, 60);
  if (t.length < 4) throw new Error(`line ${lineNo}: expected "Element x y z", got "${preview}"`);
  const x = toNumber(t[1]);
  const y = toNumber(t[2]);
  const z = toNumber(t[3]);
  if (![x, y, z].every(Number.isFinite)) {
    throw new Error(`line ${lineNo}: couldn't read coordinates from "${preview}"`);
  }
  return { symbol: normalizeSymbol(t[0]), label: t[0], x, y, z };
}

/**
 * Parse an XYZ file. Supports multiple frames (trajectories) and extended XYZ
 * (extra columns are ignored). Coordinates are in Ångström.
 *
 * Returns [{ comment, atoms: [{ symbol, label, x, y, z }] }, ...]
 */
export function parseXYZ(text) {
  const lines = text.split(/\r?\n/);
  const frames = [];
  let i = 0;

  while (i < lines.length) {
    while (i < lines.length && lines[i].trim() === '') i++;
    if (i >= lines.length) break;

    const first = lines[i].trim().split(/\s+/)[0];
    if (!/^\d+$/.test(first)) {
      // Not a standard header. If it's the very start of the file, try a forgiving parse.
      if (frames.length === 0) return parseLoose(lines);
      throw new Error(`line ${i + 1}: expected an atom count to start frame ${frames.length + 1}`);
    }

    const n = Number(first);
    const comment = (lines[i + 1] ?? '').trim();
    const atoms = new Array(n);
    for (let k = 0; k < n; k++) {
      const ln = i + 2 + k;
      if (ln >= lines.length || lines[ln].trim() === '') {
        throw new Error(`frame ${frames.length + 1} says ${n} atoms but only ${k} were found`);
      }
      atoms[k] = parseAtomLine(lines[ln], ln + 1);
    }
    frames.push({ comment, atoms });
    i += 2 + n;
  }

  if (frames.length === 0 || frames[0].atoms.length === 0) throw new Error('no atoms found');
  return frames;
}

/** Fallback for files that are just "El x y z" lines with no header. */
function parseLoose(lines) {
  const atoms = [];
  lines.forEach((line, idx) => {
    if (!/^\s*[A-Za-z]{1,3}\w*\s+[-+.\d]/.test(line)) return;
    try {
      atoms.push(parseAtomLine(line, idx + 1));
    } catch {
      /* skip lines that don't parse */
    }
  });
  if (atoms.length === 0) throw new Error("this doesn't look like an XYZ file (no atom count or atom lines found)");
  return [{ comment: '', atoms }];
}
