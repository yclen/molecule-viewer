import './style.css';
import { MoleculeViewer, STYLES } from './viewer.js';
import { parseXYZ } from './parser.js';
import { detectBonds } from './bonds.js';
import { hillFormula } from './elements.js';
import { buildParametric, clampParam } from './parametric.js';
import { excitationWeights } from './excitation.js';
import { Color } from 'three';
import TABS from '../samples/tabs.json';
import faviconUrl from '../assets/favicon.png?inline';

// Inlined as a data URL so the single-file build keeps its icon.
document.head.append(Object.assign(document.createElement('link'), { rel: 'icon', href: faviconUrl }));

// Examples are bundled into the app (not fetched), so the built single-file
// version works when opened straight from disk with no server.
// samples/<tab>/index.json lists each tab's molecules; keys below are "<tab>/<file>".
const rel = (path) => path.replace('../samples/', '');
const SAMPLE_TEXT = Object.fromEntries(
  Object.entries(import.meta.glob('../samples/*/*.xyz', { query: '?raw', import: 'default', eager: true })).map(
    ([path, text]) => [rel(path), text],
  ),
);
const TAB_INDEX = Object.fromEntries(
  Object.entries(import.meta.glob('../samples/*/index.json', { import: 'default', eager: true })).map(
    ([path, list]) => [rel(path).split('/')[0], list],
  ),
);

const $ = (id) => document.getElementById(id);
const viewer = new MoleculeViewer($('viewport'));
if (import.meta.env.DEV) window.__viewer = viewer; // console access while developing

// Keep the molecule centred in the area not covered by the side panel.
const panel = $('panel');
const syncInset = () => viewer.setInsetLeft(panel.offsetLeft + panel.offsetWidth);
new ResizeObserver(syncInset).observe(panel);
window.addEventListener('resize', syncInset);
syncInset();

const DEFAULT_SAMPLE = 'xanthenes/fluorescein.xyz';
const PLAYBACK_FPS = 20;

const state = {
  name: '',
  frames: [], // [{ comment, atoms }]
  bonds: [], // bonds per frame, computed lazily
  frame: 0,
  timer: null,
};

const bondsFor = (i) => (state.bonds[i] ??= detectBonds(state.frames[i].atoms));

// ---------------------------------------------------------------------------
// Loading
// ---------------------------------------------------------------------------

function loadText(text, name, sites = null) {
  let frames;
  try {
    frames = parseXYZ(text);
  } catch (err) {
    toast(`Couldn't open ${name}: ${err.message}`);
    return false;
  }
  stopPlayback();
  Object.assign(state, { name, frames, bonds: new Array(frames.length), frame: 0 });
  setPicked([]);
  viewer.setMolecule(frames[0].atoms, bondsFor(0));
  setSites(sites ?? [{ label: 'Whole molecule', kind: 'molecule', atoms: frames[0].atoms.map((_, i) => i) }]);
  $('empty-state').hidden = true;
  $('info').hidden = false;
  setupFrames();
  updateInfo();
  return true;
}

async function loadFile(file) {
  try {
    if (loadText(await file.text(), file.name)) {
      setActiveExample(null);
      setParamObject(null);
    }
  } catch (err) {
    toast(`Couldn't read ${file.name}: ${err.message}`);
  }
}

// ---------------------------------------------------------------------------
// Example browser (tabs + list)
// ---------------------------------------------------------------------------

const examples = { tab: null, active: null }; // active = "<tab>/<file>"

function loadExample(tab, entry) {
  const key = `${tab}/${entry.file}`;
  if (entry.param) return loadParametric(tab, entry, paramValues.get(key) ?? entry.param.default);
  const text = SAMPLE_TEXT[key];
  if (text === undefined) return toast(`Example ${key} is missing.`);
  if (loadText(text, entry.name)) {
    setActiveExample(key);
    setParamObject(null);
  }
}

// ---------------------------------------------------------------------------
// Parametric objects (e.g. a molecular wire with an adjustable number of units)
// ---------------------------------------------------------------------------

const paramValues = new Map(); // remembers each object's last value, keyed "<tab>/<file>"
let paramObject = null; // { tab, entry } of the object on screen, if any

function loadParametric(tab, entry, value) {
  const key = `${tab}/${entry.file}`;
  const p = entry.param;
  value = clampParam(p, value);
  let built;
  try {
    built = buildParametric(tab, p, value);
  } catch (err) {
    return toast(`Couldn't build ${entry.name}: ${err.message}`);
  }
  if (!loadText(built.text, `${entry.name} (${p.symbol ?? 'n'} = ${value})`, built.sites)) return;
  paramValues.set(key, value);
  setActiveExample(key);
  setParamObject({ tab, entry }, value, built.note);
}

function setParamObject(obj, value, note) {
  paramObject = obj;
  $('params').hidden = !obj;
  if (!obj) return;
  const p = obj.entry.param;
  $('param-label').textContent = p.label;
  for (const el of [$('param-input'), $('param-slider')]) Object.assign(el, { min: p.min, max: p.max, value });
  $('param-dec').disabled = value <= p.min;
  $('param-inc').disabled = value >= p.max;
  $('param-note').textContent = note ?? '';
}

const curParam = () => (paramObject ? paramValues.get(`${paramObject.tab}/${paramObject.entry.file}`) : 0);

function setParamValue(v) {
  if (!paramObject) return;
  const next = clampParam(paramObject.entry.param, v);
  if (next !== curParam()) loadParametric(paramObject.tab, paramObject.entry, next);
  else $('param-input').value = next; // show the clamped value after out-of-range typing
}

$('param-dec').addEventListener('click', () => setParamValue(curParam() - 1));
$('param-inc').addEventListener('click', () => setParamValue(curParam() + 1));
$('param-slider').addEventListener('input', (e) => setParamValue(e.target.value));
$('param-input').addEventListener('change', (e) => setParamValue(e.target.value));

function setActiveExample(key) {
  examples.active = key;
  for (const li of $('example-list').children) li.classList.toggle('active', li.dataset.key === key);
}

function showTab(tab) {
  examples.tab = tab;
  for (const b of $('tab-bar').children) b.classList.toggle('active', b.dataset.tab === tab);
  const list = $('example-list');
  list.replaceChildren(
    ...(TAB_INDEX[tab] ?? []).map((entry) => {
      const li = document.createElement('li');
      li.dataset.key = `${tab}/${entry.file}`;
      li.tabIndex = 0;
      li.title = entry.description ?? '';
      const name = document.createElement('div');
      name.className = 'ex-name';
      name.textContent = entry.name;
      const desc = document.createElement('div');
      desc.className = 'ex-desc';
      desc.textContent = entry.description ?? '';
      li.append(name, desc);
      li.addEventListener('click', () => loadExample(tab, entry));
      li.addEventListener('keydown', (e) => {
        if (e.key === 'Enter') loadExample(tab, entry);
      });
      return li;
    }),
  );
  list.scrollTop = 0;
  setActiveExample(examples.active);
}

/** Step through the current tab's examples (↑ / ↓ keys). */
function stepExample(delta) {
  const tab = examples.tab;
  const entries = TAB_INDEX[tab] ?? [];
  if (!entries.length) return;
  const cur = entries.findIndex((e) => `${tab}/${e.file}` === examples.active);
  const next = entries[(cur + delta + entries.length) % entries.length];
  loadExample(tab, next);
  $('example-list').querySelector('.active')?.scrollIntoView({ block: 'nearest' });
}

function initExamples() {
  const bar = $('tab-bar');
  for (const { id, label } of TABS) {
    if (!TAB_INDEX[id]?.length) continue;
    const b = document.createElement('button');
    b.textContent = label;
    b.dataset.tab = id;
    b.addEventListener('click', () => showTab(id));
    bar.appendChild(b);
  }
  const [tab, file] = DEFAULT_SAMPLE.split('/');
  showTab(TAB_INDEX[tab] ? tab : TABS[0].id);
  const entry = TAB_INDEX[tab]?.find((e) => e.file === file);
  if (entry && !state.frames.length) loadExample(tab, entry);
}

// ---------------------------------------------------------------------------
// Info panel
// ---------------------------------------------------------------------------

function updateInfo() {
  const f = state.frames[state.frame];
  $('mol-name').textContent = state.name;
  // symbols come from our own element table, so this HTML is safe
  $('mol-formula').innerHTML = hillFormula(f.atoms)
    .map(([s, c]) => (c > 1 ? `${s}<sub>${c}</sub>` : s))
    .join('');
  $('mol-atoms').textContent = f.atoms.length.toLocaleString();
  $('mol-bonds').textContent = (bondsFor(state.frame).length / 2).toLocaleString();
  $('mol-comment').textContent = f.comment || '—';
  $('mol-comment').title = f.comment;
}

// ---------------------------------------------------------------------------
// Display controls
// ---------------------------------------------------------------------------

const styleBox = $('style-buttons');
Object.entries(STYLES).forEach(([key, s], idx) => {
  const b = document.createElement('button');
  b.textContent = s.label;
  b.dataset.style = key;
  b.title = `${s.label} (${idx + 1})`;
  b.addEventListener('click', () => setStyle(key));
  styleBox.appendChild(b);
});

function setStyle(key) {
  viewer.setStyle(key);
  for (const b of styleBox.children) b.classList.toggle('active', b.dataset.style === key);
}
setStyle('ballstick');

function setShowHydrogens(show) {
  $('show-h').checked = show;
  viewer.setShowHydrogens(show);
  const atoms = state.frames[state.frame]?.atoms ?? [];
  if (!show) setPicked(picked.filter((i) => atoms[i]?.symbol !== 'H')); // can't measure hidden atoms
}
$('show-h').addEventListener('change', (e) => setShowHydrogens(e.target.checked));
$('reset-view').addEventListener('click', () => viewer.fitView());

// ---------------------------------------------------------------------------
// Trajectories (multi-frame XYZ)
// ---------------------------------------------------------------------------

function setupFrames() {
  const n = state.frames.length;
  $('frames').hidden = n < 2;
  $('frame-slider').max = String(n - 1);
  updateFrameLabel();
}

function updateFrameLabel() {
  $('frame-slider').value = String(state.frame);
  $('frame-label').textContent = `· frame ${state.frame + 1} / ${state.frames.length}`;
}

function showFrame(i) {
  state.frame = i;
  const f = state.frames[i];
  viewer.setFrame(f.atoms, bondsFor(i));
  updateFrameLabel();
  updateInfo();
  updateMeasure(); // measurements follow the trajectory
  updateGlow(); // so does the glow
}

function startPlayback() {
  if (state.frames.length < 2 || state.timer) return;
  state.timer = setInterval(() => showFrame((state.frame + 1) % state.frames.length), 1000 / PLAYBACK_FPS);
  $('play-btn').textContent = '❚❚';
}

function stopPlayback() {
  clearInterval(state.timer);
  state.timer = null;
  $('play-btn').textContent = '▶';
}

$('play-btn').addEventListener('click', () => (state.timer ? stopPlayback() : startPlayback()));
$('frame-slider').addEventListener('input', (e) => {
  stopPlayback();
  showFrame(Number(e.target.value));
});

// ---------------------------------------------------------------------------
// Hover tooltip
// ---------------------------------------------------------------------------

const tip = $('tooltip');
viewer.onHover = (i, atom, x, y) => {
  if (i < 0) {
    tip.hidden = true;
    return;
  }
  const bonds = bondsFor(state.frame);
  let degree = 0;
  for (let b = 0; b < bonds.length; b++) if (bonds[b] === i) degree++;

  const fmt = (v) => v.toFixed(3).padStart(9);
  const title = document.createElement('b');
  title.textContent = `${atom.symbol} #${i + 1}${atom.label !== atom.symbol ? `  (${atom.label})` : ''}`;
  tip.replaceChildren(
    title,
    `\nx ${fmt(atom.x)} Å\ny ${fmt(atom.y)} Å\nz ${fmt(atom.z)} Å\n${degree} bond${degree === 1 ? '' : 's'}`,
  );
  tip.hidden = false;
  // keep it on-screen
  const w = tip.offsetWidth;
  const h = tip.offsetHeight;
  tip.style.left = `${Math.min(x + 16, window.innerWidth - w - 8)}px`;
  tip.style.top = `${Math.min(y + 16, window.innerHeight - h - 8)}px`;
};

// ---------------------------------------------------------------------------
// File open + drag & drop
// ---------------------------------------------------------------------------

const fileInput = $('file-input');
$('open-btn').addEventListener('click', () => fileInput.click());
fileInput.addEventListener('change', () => {
  if (fileInput.files[0]) loadFile(fileInput.files[0]);
  fileInput.value = '';
});

const overlay = $('drop-overlay');
const hasFiles = (e) => e.dataTransfer?.types?.includes('Files');
let dragDepth = 0;
window.addEventListener('dragenter', (e) => {
  if (!hasFiles(e)) return;
  e.preventDefault();
  dragDepth++;
  overlay.hidden = false;
});
window.addEventListener('dragover', (e) => {
  if (hasFiles(e)) e.preventDefault();
});
window.addEventListener('dragleave', () => {
  if (--dragDepth <= 0) {
    dragDepth = 0;
    overlay.hidden = true;
  }
});
window.addEventListener('drop', (e) => {
  e.preventDefault();
  dragDepth = 0;
  overlay.hidden = true;
  const file = e.dataTransfer?.files?.[0];
  if (file) loadFile(file);
});

// ---------------------------------------------------------------------------
// Measuring: click 2-4 atoms -> distance, angle, dihedral
// ---------------------------------------------------------------------------

const MAX_PICK = 4;
let picked = []; // atom indices, in click order

viewer.onAtomClick = (i) => {
  if (i < 0) return setPicked([]); // click on empty space clears
  if (picked.includes(i)) return setPicked(picked.filter((j) => j !== i)); // click again to unpick
  setPicked(picked.length >= MAX_PICK ? [i] : [...picked, i]);
};

function setPicked(list) {
  picked = list;
  viewer.setSelection(picked);
  updateMeasure();
}

const vec = (a, b) => [b.x - a.x, b.y - a.y, b.z - a.z];
const dot = (u, v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
const cross = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
const norm = (u) => Math.sqrt(dot(u, u));
const DEG = 180 / Math.PI;

function angleDeg(a, b, c) {
  const u = vec(b, a);
  const v = vec(b, c);
  return Math.acos(Math.min(1, Math.max(-1, dot(u, v) / (norm(u) * norm(v))))) * DEG;
}

function dihedralDeg(a, b, c, d) {
  const b1 = vec(a, b);
  const b2 = vec(b, c);
  const b3 = vec(c, d);
  const n1 = cross(b1, b2);
  const n2 = cross(b2, b3);
  const m1 = cross(n1, b2.map((x) => x / norm(b2)));
  return Math.atan2(dot(m1, n2), dot(n1, n2)) * DEG;
}

function updateMeasure() {
  const body = $('measure-body');
  const atoms = state.frames[state.frame]?.atoms ?? [];
  const sel = picked.filter((i) => i < atoms.length);
  if (sel.length < 2) {
    body.className = 'muted';
    body.textContent = sel.length
      ? `${tag(sel[0], atoms)} selected; click another atom.`
      : 'Click 2–4 atoms to measure distance, angle and dihedral.';
    return;
  }
  const A = sel.map((i) => atoms[i]);
  const rows = [];
  for (let k = 1; k < sel.length; k++) {
    rows.push(['Distance', `${tag(sel[k - 1], atoms)}–${tag(sel[k], atoms)}`,
      `${norm(vec(A[k - 1], A[k])).toFixed(3)} Å`]);
  }
  for (let k = 2; k < sel.length; k++) {
    rows.push(['Angle', sel.slice(k - 2, k + 1).map((i) => tag(i, atoms)).join('–'),
      `${angleDeg(A[k - 2], A[k - 1], A[k]).toFixed(1)}°`]);
  }
  if (sel.length === 4) {
    rows.push(['Dihedral', sel.map((i) => tag(i, atoms)).join('–'), `${dihedralDeg(...A).toFixed(1)}°`]);
  }
  body.className = 'measure-table';
  body.replaceChildren(
    ...rows.flatMap(([kind, what, value]) => {
      const k = Object.assign(document.createElement('span'), { className: 'm-kind', textContent: kind });
      const w = Object.assign(document.createElement('span'), { className: 'm-what', textContent: what });
      const v = Object.assign(document.createElement('span'), { className: 'm-value', textContent: value });
      return [k, w, v];
    }),
  );
}

const tag = (i, atoms) => `${atoms[i].symbol}${i + 1}`;

// ---------------------------------------------------------------------------
// Keyboard shortcuts
// ---------------------------------------------------------------------------

const styleKeys = Object.keys(STYLES);
// −/+ resize a parametric object; matched on the physical key so it works with or without Shift and on the numpad
const SIZE_KEYS = { Equal: 1, NumpadAdd: 1, Minus: -1, NumpadSubtract: -1 };
window.addEventListener('keydown', (e) => {
  if (e.ctrlKey || e.metaKey || e.altKey) return;
  if (e.target.closest?.('input[type=text], input[type=number], input[type=range], textarea, select')) return;
  const k = e.key.toLowerCase();
  if (k === ' ' && e.target.tagName === 'BUTTON') return; // the focused button handles Space itself
  if (k === 'escape') setPicked([]);
  else if (k === 'arrowdown' || k === 'arrowup') {
    e.preventDefault();
    stepExample(k === 'arrowdown' ? 1 : -1);
  } else if (k >= '1' && k <= String(styleKeys.length)) setStyle(styleKeys[Number(k) - 1]);
  else if (k === 'h') {
    setShowHydrogens(!$('show-h').checked);
  } else if (k === 'r') viewer.fitView();
  else if (k === 'o') fileInput.click();
  else if (k === 'e' && excite.sites.length) toggleDefaultSite();
  else if (paramObject && SIZE_KEYS[e.code] !== undefined) setParamValue(curParam() + SIZE_KEYS[e.code]);
  else if (k === ' ' && state.frames.length > 1) {
    e.preventDefault();
    state.timer ? stopPlayback() : startPlayback();
  }
});

// ---------------------------------------------------------------------------
// Excited state: a glowing Gaussian cloud on the excited chromophore(s)
// ---------------------------------------------------------------------------

const excite = {
  sites: [], // [{ label, kind, atoms }]
  on: new Set(), // labels of excited sites (kept when the molecule is rebuilt, e.g. a new n)
  weights: new Map(), // label -> { atoms, weights, method } for the current molecule
  color: new Color('#3d8bff'),
};

function setSites(sites) {
  excite.sites = sites;
  excite.weights.clear();
  const labels = new Set(sites.map((s) => s.label));
  for (const l of [...excite.on]) if (!labels.has(l)) excite.on.delete(l);
  $('site-chips').replaceChildren(
    ...sites.map((s) => {
      const b = document.createElement('button');
      b.textContent = s.label;
      b.dataset.label = s.label;
      b.addEventListener('click', () => toggleSite(s.label));
      return b;
    }),
  );
  $('excite').hidden = false;
  updateGlow();
}

function toggleSite(label) {
  if (excite.on.has(label)) excite.on.delete(label);
  else excite.on.add(label);
  updateGlow();
}

/** E key: switch everything off, or excite the first DPA (or the first site if there is no DPA). */
function toggleDefaultSite() {
  if (excite.on.size) excite.on.clear();
  else excite.on.add((excite.sites.find((s) => s.label.startsWith('DPA')) ?? excite.sites[0]).label);
  updateGlow();
}

function siteWeights(site) {
  if (!excite.weights.has(site.label)) {
    excite.weights.set(site.label, excitationWeights(state.frames[0].atoms, bondsFor(0), site.atoms));
  }
  return excite.weights.get(site.label);
}

function updateGlow() {
  const atoms = state.frames[state.frame]?.atoms ?? [];
  const sources = [];
  const methods = new Set();
  for (const site of excite.sites) {
    if (!excite.on.has(site.label)) continue;
    const w = siteWeights(site);
    w.atoms.forEach((i, k) => {
      if (w.weights[k] > 0.02) sources.push({ ...atoms[i], weight: w.weights[k], color: excite.color });
    });
    methods.add(w.method);
  }
  viewer.setGlow(sources);
  for (const b of $('site-chips').children) b.classList.toggle('excited', excite.on.has(b.dataset.label));
  $('excite-note').textContent = methods.size
    ? `Glow follows the π density of the lowest excitation (${[...methods].join(', ')}, twist-corrected).`
    : '';
}

function setGlowParams() {
  const gain = Number($('glow-gain').value);
  const sigma = Number($('glow-sigma').value);
  viewer.setGlowParams({ gain, sigma });
  $('glow-gain-val').textContent = `${gain.toFixed(2)}×`;
  $('glow-sigma-val').textContent = `${sigma.toFixed(1)} Å`;
}
$('glow-gain').addEventListener('input', setGlowParams);
$('glow-sigma').addEventListener('input', setGlowParams);
$('glow-color').addEventListener('input', (e) => {
  excite.color.set(e.target.value);
  updateGlow();
});
setGlowParams();

// ---------------------------------------------------------------------------
// Toast
// ---------------------------------------------------------------------------

let toastTimer;
function toast(msg) {
  const el = $('toast');
  el.textContent = msg;
  el.hidden = false;
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => (el.hidden = true), 6000);
}

initExamples();
