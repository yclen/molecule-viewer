import './style.css';
import { MoleculeViewer, STYLES } from './viewer.js';
import { parseXYZ } from './parser.js';
import { detectBonds } from './bonds.js';
import { hillFormula } from './elements.js';
import { buildParametric, clampParam } from './parametric.js';
import { excitationWeights } from './excitation.js';
import { ExcitedStates, STATES, STATE_NAME, STATE_SYMBOL, copyState } from './excited-states.js';
import { Simulation, parseSliders } from './simulation.js';
import { Color } from 'three';
import TABS from '../samples/tabs.json';
import faviconUrl from '../assets/favicon.png?inline';
import { initTabs } from './lab/tabs.js';

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

/**
 * objectKey: what excited states are remembered under (an example's key, or the file name).
 * view: { surroundings, keepAxes } for the viewer (see MoleculeViewer.setMolecule).
 */
function loadText(text, name, sites = null, objectKey = name, view = {}) {
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
  viewer.setMolecule(frames[0].atoms, bondsFor(0), view);
  setSites(objectKey, sites ?? [{ label: 'Whole molecule', kind: 'molecule', atoms: frames[0].atoms.map((_, i) => i) }]);
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
      showAbout(null);
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
  const view = { surroundings: entry.surroundings ?? 0, keepAxes: Boolean(entry.keepAxes) };
  if (loadText(text, entry.name, namedSites(entry, text), key, view)) {
    setActiveExample(key);
    setParamObject(null);
    showAbout(entry);
  }
}

/** Sites an example names in its index.json: each is `last` atoms from the file's end, less the final `skip`. Null if none. */
function namedSites(entry, text) {
  const n = Number.parseInt(text, 10);
  return (
    entry.sites?.map((s) => {
      const first = n - s.last - (s.skip ?? 0);
      return { ...s, atoms: Array.from({ length: s.last }, (_, k) => first + k) };
    }) ?? null
  );
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
  const sens = entry.sensitizers ? sensitizerFor(key, entry) : null;
  let built;
  try {
    built = buildParametric(tab, p, value, { sensitizer: sens });
  } catch (err) {
    return toast(`Couldn't build ${entry.name}${sens ? ` with ${sens.name}` : ''}: ${err.message}`);
  }
  const name = `${entry.name} (${p.symbol ?? 'n'} = ${value}${sens ? `, ${sens.name}` : ''})`;
  if (!loadText(built.text, name, built.sites, key)) return;
  paramValues.set(key, value);
  setActiveExample(key);
  setParamObject({ tab, entry }, value, built.note);
  showSensitizers(key, entry, built.attach);
  showAbout(entry);
}

// ---------------------------------------------------------------------------
// Sensitizer swap: objects with a "sensitizers" entry (TTA Flora) can take any sensitizer from the sample
// library in place of their own. sensitizers.from lists sample tabs ("uc-sens") or files ("tta/ptoep.xyz").
// ---------------------------------------------------------------------------

const sensChoice = new Map(); // "<tab>/<file>" of the object -> chosen sample key ('' = the original)
const sensAtoms = new Map(); // sample key -> [[symbol, [x, y, z]]]

function sensitizerOptions(entry) {
  const groups = [];
  for (const src of entry.sensitizers.from) {
    const [tab, file] = src.split('/');
    const label = TABS.find((t) => t.id === tab)?.label ?? tab;
    const items = (TAB_INDEX[tab] ?? [])
      .filter((e) => !e.param && (!file || e.file === file) && SAMPLE_TEXT[`${tab}/${e.file}`] !== undefined)
      .map((e) => ({ key: `${tab}/${e.file}`, name: e.name.replace(/\s*\(sensitizer\)$/i, ''), desc: e.description }));
    const g = groups.find((x) => x.label === label);
    if (g) g.items.push(...items);
    else if (items.length) groups.push({ label, items });
  }
  return groups;
}

function sensitizerFor(objKey, entry) {
  const key = sensChoice.get(objKey);
  if (!key) return null;
  if (!sensAtoms.has(key)) {
    sensAtoms.set(key, parseXYZ(SAMPLE_TEXT[key])[0].atoms.map((a) => [a.symbol, [a.x, a.y, a.z]]));
  }
  const name = sensitizerOptions(entry).flatMap((g) => g.items).find((i) => i.key === key)?.name ?? key;
  return { key, name, atoms: sensAtoms.get(key) };
}

let sensFor = null; // object key the dropdown is showing
function showSensitizers(objKey, entry, attach) {
  const box = $('sensitizer');
  box.hidden = !entry?.sensitizers;
  if (box.hidden) return;
  const sel = $('sens-select');
  if (sensFor !== objKey) {
    sensFor = objKey;
    const orig = Object.assign(document.createElement('option'), { value: '', textContent: entry.sensitizers.original });
    sel.replaceChildren(
      orig,
      ...sensitizerOptions(entry).map((g) => {
        const og = Object.assign(document.createElement('optgroup'), { label: g.label });
        og.append(...g.items.map((i) => Object.assign(document.createElement('option'), { value: i.key, textContent: i.name, title: i.desc ?? '' })));
        return og;
      }),
    );
  }
  sel.value = sensChoice.get(objKey) ?? '';
  $('sens-note').textContent = attach
    ? `Linked through ${attach.symbol}${attach.atom} of its sample structure (the outermost free aromatic C–H) ` +
      'to the phenylene, torsion set to avoid clashes. Rigid placement, not re-optimised.'
    : 'As built in tools/build_triplet.py (conformer search + MMFF94s).';
}

$('sens-select').addEventListener('change', (e) => {
  if (!paramObject) return;
  sensChoice.set(sensFor, e.target.value);
  loadParametric(paramObject.tab, paramObject.entry, curParam());
  e.target.blur(); // so arrow keys step examples again rather than the dropdown
});

function setParamObject(obj, value, note) {
  paramObject = obj;
  $('params').hidden = !obj;
  if (!obj) $('sensitizer').hidden = true;
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

// ---------------------------------------------------------------------------
// About this molecule: an example's `info` text (from tools/catalog.py), a small subset of Markdown:
// blank-line-separated blocks; "### Heading"; lines starting "- " make a list; **bold** and *italic* inline;
// a block starting "Sources:" is set small. Everything else is escaped text.
// ---------------------------------------------------------------------------

const escapeHtml = (s) => s.replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' })[c]);
const inline = (s) =>
  escapeHtml(s)
    .replace(/\*\*(.+?)\*\*/g, '<strong>$1</strong>')
    .replace(/(^|[^*])\*(?!\s)(.+?)\*/g, '$1<em>$2</em>');

/**
 * Within a block: "### " starts a heading line; "- " at the start of a line starts a list item, and the indented
 * lines after it continue that item; any other unindented line is paragraph text (so a sentence can lead into a
 * list, and text can follow one, without a blank line).
 */
function renderInfo(text) {
  const out = [];
  for (const block of text.replace(/\r/g, '').trim().split(/\n\s*\n/)) {
    let para = [];
    let items = null;
    const flushPara = () => {
      if (!para.length) return;
      const cls = para[0].startsWith('Sources:') ? ' class="src"' : '';
      out.push(`<p${cls}>${inline(para.join(' '))}</p>`);
      para = [];
    };
    const flushList = () => {
      if (!items) return;
      out.push(`<ul>${items.map((it) => `<li>${inline(it.join(' '))}</li>`).join('')}</ul>`);
      items = null;
    };
    for (const raw of block.split('\n')) {
      const line = raw.trim();
      if (!line) continue;
      if (line.startsWith('### ')) {
        flushPara();
        flushList();
        out.push(`<h4>${inline(line.slice(4))}</h4>`);
      } else if (line.startsWith('- ')) {
        flushPara();
        (items ??= []).push([line.slice(2)]);
      } else if (items && /^\s/.test(raw)) {
        items[items.length - 1].push(line);
      } else {
        flushList();
        para.push(line);
      }
    }
    flushPara();
    flushList();
  }
  return out.join('');
}

/** Show (or hide, for null / no info) the About panel for an example entry. */
function showAbout(entry) {
  const info = entry?.info;
  $('about').hidden = !info;
  $('about-body').innerHTML = info ? renderInfo(info) : '';
  $('about-body').scrollTop = 0;
}

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
const onViewerTab = () => (document.body.dataset.page ?? 'viewer') === 'viewer';
window.addEventListener('dragenter', (e) => {
  if (!hasFiles(e) || !onViewerTab()) return;
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
  if (!onViewerTab()) return; // the other tabs handle their own drops
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
  const site = excite.sites.find((s) => s.atoms.includes(i));
  if (site) selectSite(site.label); // clicking a molecule also opens it in the Excited states editor
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
  if (e.ctrlKey || e.metaKey || e.altKey || !onViewerTab()) return;
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
  else if (k === 'e' && excite.sites.length) toggleSelectedSite('S');
  else if (k === 't' && excite.sites.length) toggleSelectedSite('T');
  else if (paramObject && SIZE_KEYS[e.code] !== undefined) setParamValue(curParam() + SIZE_KEYS[e.code]);
  else if (k === ' ' && simMode) {
    e.preventDefault();
    simPlayPause();
  } else if (k === ' ' && state.frames.length > 1) {
    e.preventDefault();
    state.timer ? stopPlayback() : startPlayback();
  }
});

// ---------------------------------------------------------------------------
// Excited states: a glowing Gaussian cloud on the excited chromophore(s)
// ---------------------------------------------------------------------------

// Each site (molecule / chromophore) has two excited states, the singlet S₁ and the triplet T₁, each with its own
// { on, color, gain, sigma, fade } (see excited-states.js). One site is selected at a time; the editor below the list
// changes only that one (both of its states), and every other site keeps what was set for it.
//
// Every ground ⇄ excited change goes through setOn(), whether from a click or a simulation script, so the
// simulation hears about clicks too. The glow of each (site, state) fades toward its state (a level 0..1 per
// glow group, see GlowLayer.setLevels). A state with an afterglow (fade > 0) fades, and even a very short
// excitation of it shows as a flash; with fade = 0 the glow follows the state in the same frame. Group of
// (site i, state) = 2i for S, 2i + 1 for T.
const excite = new ExcitedStates();
excite.key = null; // the current object's key
excite.selected = null; // label of the site being edited
excite.weights = new Map(); // label -> { atoms, weights, method } for the current molecule
excite.levels = new Map(); // 'label|S' -> displayed glow level, 0..1
excite.litAt = new Map(); // 'label|S' -> real time (s) it last went ground → excited
excite.built = new Set(); // 'label|S' keys whose sprites are in the glow layer right now
if (import.meta.env.DEV) window.__excite = excite; // e.g. __setOn('DPA 1', true, 'T')

const FADE_IN = 0.04; // s (real time), time constant of the glow coming on (going off: each state's own fade)
const MIN_FLASH = 0.12; // s: with a fade, an excitation shorter than this still shows for this long

let realTime = 0;
let glowDirty = true; // colours, sizes or atom positions changed: rebuild the glow sprites
let listDirty = true; // on/off changed: refresh the site list

const glowKey = (label, s) => `${label}|${s}`;

/** The one place a state is switched on or off. state: 'S' (singlet, default) or 'T' (triplet). */
function setOn(label, on, state = 'S') {
  const st = excite.get(label, state);
  if (!st || st.on === on) return;
  st.on = on;
  if (on) excite.litAt.set(glowKey(label, state), realTime);
  if (label === excite.selected) $(`edit-on-${state}`).checked = on;
  listDirty = true;
  sim.notify(label, on, state);
}
if (import.meta.env.DEV) window.__setOn = setOn;

function setSites(objectKey, sites) {
  excite.setSites(objectKey, sites);
  excite.key = objectKey;
  excite.weights.clear();
  for (const s of sites) {
    for (const k of STATES) excite.levels.set(glowKey(s.label, k), excite.get(s.label, k).on ? 1 : 0); // no fade on load
  }
  $('site-list').replaceChildren(
    ...sites.map((s) => {
      const li = document.createElement('li');
      li.dataset.label = s.label;
      li.tabIndex = 0;
      const switches = STATES.map((k) => {
        const sw = document.createElement('button');
        sw.className = `e-switch e-switch-${k}`;
        sw.textContent = k;
        sw.title = `${STATE_NAME[k][0].toUpperCase()}${STATE_NAME[k].slice(1)} ${STATE_SYMBOL[k]}: excite / de-excite`;
        sw.addEventListener('click', (e) => {
          e.stopPropagation();
          selectSite(s.label);
          toggleSite(s.label, k);
        });
        return sw;
      });
      const name = Object.assign(document.createElement('span'), { className: 'e-name', textContent: s.label });
      const st = Object.assign(document.createElement('span'), { className: 'e-state' });
      li.append(...switches, name, st);
      li.addEventListener('click', () => selectSite(s.label));
      li.addEventListener('keydown', (e) => {
        if (e.key === 'Enter') selectSite(s.label);
      });
      return li;
    }),
  );
  $('excite').hidden = false;
  // keep editing the same site across rebuilds (e.g. a new n); otherwise start on the first excited one
  const keep = sites.some((s) => s.label === excite.selected);
  selectSite(keep ? excite.selected : (excite.excited()[0] ?? sites[0]).label);
  // a running simulation's script was written against the old sites
  if (sim.started) {
    simReset();
    simStatus('Reset: the molecule changed.');
  }
}

function selectSite(label) {
  excite.selected = label;
  $('edit-name').textContent = label;
  for (const k of STATES) {
    const st = excite.get(label, k);
    $(`edit-on-${k}`).checked = st.on;
    $(`glow-color-${k}`).value = st.color;
    $(`glow-gain-${k}`).value = st.gain;
    $(`glow-sigma-${k}`).value = st.sigma;
    $(`glow-fade-${k}`).value = st.fade;
  }
  updateGlow();
}

/** Change colour / brightness / spread / afterglow of one state of a site (not on/off: that's setOn). */
function setProps(label, patch, state = 'S') {
  const st = excite.get(label, state);
  if (!st) return;
  for (const k of ['color', 'gain', 'sigma', 'fade']) if (patch[k] !== undefined) st[k] = patch[k];
  if (label === excite.selected) selectSite(label);
  else updateGlow();
}

function editSelected(state, patch) {
  if (excite.selected) setProps(excite.selected, patch, state);
}

function toggleSite(label, state = 'S') {
  setOn(label, !excite.get(label, state).on, state);
}

/** E / T keys: switch the selected site's singlet / triplet. */
function toggleSelectedSite(state) {
  if (excite.selected) toggleSite(excite.selected, state);
}

function siteWeights(site) {
  if (!excite.weights.has(site.label)) {
    excite.weights.set(site.label, excitationWeights(state.frames[0].atoms, bondsFor(0), site.atoms));
  }
  return excite.weights.get(site.label);
}

/** Something about the glow changed (colour, size, frame): rebuild it on the next frame. */
function updateGlow() {
  glowDirty = true;
  listDirty = true;
  if (!excite.selected) return;
  for (const k of STATES) {
    const st = excite.get(excite.selected, k);
    $(`glow-gain-${k}-val`).textContent = `${st.gain.toFixed(2)}×`;
    $(`glow-sigma-${k}-val`).textContent = `${st.sigma.toFixed(1)} Å`;
    $(`glow-fade-${k}-val`).textContent = st.fade > 0 ? `${st.fade.toFixed(2)} s` : 'instant';
  }
}
if (import.meta.env.DEV) window.__updateGlow = updateGlow;

/** Sprites for the (site, state) keys in `visible`, each its own level group. */
function rebuildGlow(visible) {
  const atoms = state.frames[state.frame]?.atoms ?? [];
  const sources = [];
  const methods = new Set();
  excite.sites.forEach((site, i) => {
    STATES.forEach((k, j) => {
      if (!visible.has(glowKey(site.label, k))) return;
      const st = excite.get(site.label, k);
      const color = new Color(st.color);
      const group = 2 * i + j;
      // both states glow where the lowest π→π* excitation sits: S₁ and T₁ share that HOMO/LUMO pair
      const w = siteWeights(site);
      w.atoms.forEach((a, n) => {
        if (w.weights[n] > 0.02) sources.push({ ...atoms[a], weight: w.weights[n] * st.gain, sigma: st.sigma, color, group });
      });
      methods.add(w.method);
    });
  });
  viewer.setGlow(sources);
  excite.built = visible;
  glowDirty = false;
  $('excite-note').textContent = methods.size
    ? `Glow follows the π density of the lowest excitation (${[...methods].join(', ')}, twist-corrected).`
    : '';
}

/** Every frame: fade each glow toward its state, rebuilding the sprites only when needed. */
function stepGlow(dt) {
  const levels = [];
  const visible = new Set();
  excite.sites.forEach((site, i) => {
    STATES.forEach((k, j) => {
      const key = glowKey(site.label, k);
      const st = excite.get(site.label, k);
      const held = st.fade > 0 && realTime - (excite.litAt.get(key) ?? -Infinity) < MIN_FLASH;
      const target = st.on || held ? 1 : 0;
      const prev = excite.levels.get(key) ?? 0;
      let lvl = target; // no afterglow: the glow is exactly the state
      if (st.fade > 0) {
        lvl = prev + (target - prev) * (1 - Math.exp(-dt / (target > prev ? FADE_IN : st.fade)));
        if (Math.abs(target - lvl) < 0.004) lvl = target;
      }
      excite.levels.set(key, lvl);
      levels[2 * i + j] = lvl;
      if (target > 0 || lvl > 0) visible.add(key);
    });
  });
  const same = visible.size === excite.built.size && [...visible].every((l) => excite.built.has(l));
  if (glowDirty || !same) rebuildGlow(visible);
  viewer.setGlowLevels(levels);
}

function updateSiteList() {
  for (const li of $('site-list').children) {
    const { S, T } = excite.get(li.dataset.label);
    li.classList.toggle('selected', li.dataset.label === excite.selected);
    li.classList.toggle('on-S', S.on);
    li.classList.toggle('on-T', T.on);
    li.style.setProperty('--s-color', S.color);
    li.style.setProperty('--t-color', T.color);
    li.style.setProperty('--e-color', T.on && !S.on ? T.color : S.color);
    const text = [S.on && STATE_SYMBOL.S, T.on && STATE_SYMBOL.T].filter(Boolean).join(' + ');
    li.querySelector('.e-state').textContent = text || 'ground';
  }
  listDirty = false;
}

for (const k of STATES) {
  $(`edit-on-${k}`).addEventListener('change', (e) => excite.selected && setOn(excite.selected, e.target.checked, k));
  $(`glow-color-${k}`).addEventListener('input', (e) => editSelected(k, { color: e.target.value }));
  $(`glow-gain-${k}`).addEventListener('input', (e) => editSelected(k, { gain: Number(e.target.value) }));
  $(`glow-sigma-${k}`).addEventListener('input', (e) => editSelected(k, { sigma: Number(e.target.value) }));
  $(`glow-fade-${k}`).addEventListener('input', (e) => editSelected(k, { fade: Number(e.target.value) }));
  $(`excite-all-${k}`).addEventListener('click', () => {
    for (const s of excite.sites) setOn(s.label, true, k);
  });
}
$('excite-all-off').addEventListener('click', () => {
  for (const s of excite.sites) for (const k of STATES) setOn(s.label, false, k);
});

// ---------------------------------------------------------------------------
// Simulation mode: scripts that switch excited states on and off over time (see simulation.js)
// ---------------------------------------------------------------------------

const sim = new Simulation({
  sites: () => excite.sites,
  isOn: (label, s) => excite.get(label, s).on,
  setOn,
  setProps,
});
let simMode = false;
let simUiDirty = true;
let simSnapshot = null; // { key, states } taken when a run starts, put back by Reset
let simNote = ''; // a message from the app (not the script), e.g. why the run was reset
sim.onChange = () => (simUiDirty = true);
if (import.meta.env.DEV) window.__sim = sim;

// Preset scripts: samples/sims/*.js, bundled as text. The first line (a // comment) is the name.
const CUSTOM_ID = 'custom';
const CUSTOM_TEMPLATE = `// My script
// Runs once at t = 0; every / after / onExcite keep running in the background.
// The full API is listed at the top of src/simulation.js.

every(1, () => {
  for (const s of sites()) toggle(s);
});
`;
const SIM_PRESETS = Object.entries(
  import.meta.glob('../samples/sims/*.js', { query: '?raw', import: 'default', eager: true }),
)
  .sort(([a], [b]) => a.localeCompare(b))
  .map(([path, code]) => {
    const id = path.split('/').pop();
    return { id, name: code.match(/^\s*\/\/\s*(.+)/)?.[1].trim() ?? id, code };
  });
SIM_PRESETS.push({ id: CUSTOM_ID, name: 'My script', code: CUSTOM_TEMPLATE });

const store = {
  get(k) {
    try {
      return localStorage.getItem(k);
    } catch {
      return null;
    }
  },
  set(k, v) {
    try {
      if (v === null) localStorage.removeItem(k);
      else localStorage.setItem(k, v);
    } catch {
      // storage unavailable: edits just aren't remembered
    }
  },
};
const codeKey = (id) => `molview:sim-code:${id}`;
const sliderKey = (id) => `molview:sim-sliders:${id}`;
const preset = (id) => SIM_PRESETS.find((p) => p.id === id) ?? SIM_PRESETS[0];
let simScript = preset(store.get('molview:sim-script')).id;

function initSim() {
  $('sim-script').replaceChildren(
    ...SIM_PRESETS.map((p) => Object.assign(document.createElement('option'), { value: p.id, textContent: p.name })),
  );
  selectScript(simScript);
  setSpeed(Number(store.get('molview:sim-speed') ?? 0));
  setMode(store.get('molview:mode') === 'sim' ? 'sim' : 'view');
}

function selectScript(id) {
  simScript = preset(id).id;
  store.set('molview:sim-script', simScript);
  $('sim-script').value = simScript;
  $('sim-code').value = store.get(codeKey(simScript)) ?? preset(simScript).code;
  loadSliders();
  updateRevert();
  if (sim.started) simStatus('Script changed: Reset, then Run to use it.');
}

function updateRevert() {
  $('sim-revert').disabled = $('sim-code').value === preset(simScript).code && !sim.sliderValues.size;
}

/** The current script's moved sliders, as saved in the browser. */
function loadSliders() {
  let saved = {};
  try {
    saved = JSON.parse(store.get(sliderKey(simScript))) ?? {};
  } catch {
    // unreadable: the sliders start from the script's values
  }
  sim.sliderValues = new Map(Object.entries(saved));
  renderSliders();
}

function setSlider(name, value) {
  sim.sliderValues.set(name, value);
  store.set(sliderKey(simScript), JSON.stringify(Object.fromEntries(sim.sliderValues)));
  updateRevert();
}

/** One row per slider() in the script's text. A running script reads the values live (see simulation.js). */
function renderSliders() {
  const sliders = parseSliders($('sim-code').value);
  $('sim-sliders').hidden = !sliders.length;
  $('sim-sliders').replaceChildren(
    ...sliders.flatMap((s) => {
      const toPos = (v) => (s.log ? Math.log10(v) : v);
      const lo = toPos(s.min);
      const hi = toPos(s.max);
      const input = Object.assign(document.createElement('input'), { type: 'range', min: lo, max: hi, step: (hi - lo) / 200 });
      const val = Object.assign(document.createElement('span'), { className: 'val' });
      const show = (v) => (val.textContent = `${v} ${s.unit}`.trim());
      const value = Math.min(s.max, Math.max(s.min, sim.sliderValues.get(s.name) ?? s.value));
      input.value = toPos(value);
      show(value);
      input.addEventListener('input', () => {
        const p = Number(input.value);
        const v = Number((s.log ? 10 ** p : p).toPrecision(3));
        setSlider(s.name, v);
        show(v);
      });
      return [Object.assign(document.createElement('span'), { textContent: s.name }), input, val];
    }),
  );
}

function setMode(mode) {
  simMode = mode === 'sim';
  store.set('molview:mode', mode);
  for (const b of $('mode-buttons').children) b.classList.toggle('active', b.dataset.mode === mode);
  $('sim').hidden = !simMode;
  $('panel').classList.toggle('wide', simMode);
  if (!simMode) sim.pause(); // leaving the mode freezes the run; Resume carries on from there
  simUiDirty = true;
}

/** v: log10 of the speed, so the slider is even from 0.1× to 20×. */
function setSpeed(v) {
  sim.speed = 10 ** v;
  store.set('molview:sim-speed', String(v));
  $('sim-speed').value = v;
  const x = sim.speed;
  $('sim-speed-val').textContent = `${x < 1 ? x.toFixed(2) : x < 10 ? x.toFixed(1) : x.toFixed(0)}×`;
}

function simStatus(msg) {
  simNote = msg;
  simUiDirty = true;
}

function simStart() {
  simSnapshot = { key: excite.key, states: new Map([...excite.states].map(([l, st]) => [l, copyState(st)])) };
  for (const s of excite.sites) siteWeights(s); // so the first excitation of a big site doesn't stall the clock
  simNote = '';
  sim.start($('sim-code').value);
}

/** Stop the run and put the excited states back as they were when it started. */
function simReset() {
  sim.reset();
  if (simSnapshot) {
    const here = simSnapshot.key === excite.key; // or the run's molecule was swapped out
    const states = excite.objects.get(simSnapshot.key);
    for (const [label, saved] of simSnapshot.states) {
      const st = states?.get(label);
      if (!st) continue;
      for (const k of STATES) {
        const { on, ...props } = saved[k];
        Object.assign(st[k], props);
        if (here) setOn(label, on, k);
        else st[k].on = on;
      }
    }
    if (here && excite.selected) selectSite(excite.selected);
  }
  simSnapshot = null;
  simStatus('');
}

function simPlayPause() {
  if (sim.status === 'running') sim.pause();
  else if (sim.status === 'paused') sim.resume();
  else {
    if (sim.status === 'error') simReset();
    simStart();
  }
}

const SIM_STEP = 0.1; // sim s per Step

function simStep() {
  if (sim.status === 'error') simReset();
  if (!sim.started) simStart();
  sim.pause();
  sim.step(SIM_STEP);
}

function updateSimUi() {
  simUiDirty = false;
  const s = sim.status;
  $('sim-play').textContent = s === 'running' ? '❚❚ Pause' : s === 'paused' ? '▶ Resume' : '▶ Run';
  const err = sim.error ? `Error: ${sim.error.message}` : '';
  $('sim-status').textContent = err || simNote;
  $('sim-status').classList.toggle('error', Boolean(err));
  $('sim-counters').replaceChildren(
    ...[...sim.counters].flatMap(([name, n]) => [
      Object.assign(document.createElement('dt'), { textContent: name }),
      Object.assign(document.createElement('dd'), { textContent: typeof n === 'number' && !Number.isInteger(n) ? n.toFixed(3) : String(n) }),
    ]),
  );
  $('sim-counters').hidden = !sim.counters.size;
  $('sim-log').textContent = sim.logLines.slice(-6).join('\n');
  $('sim-log').hidden = !sim.logLines.length;
}

$('mode-buttons').addEventListener('click', (e) => {
  const b = e.target.closest('button[data-mode]');
  if (b) setMode(b.dataset.mode);
});
$('sim-script').addEventListener('change', (e) => selectScript(e.target.value));
$('sim-code').addEventListener('input', (e) => {
  const code = e.target.value;
  store.set(codeKey(simScript), code === preset(simScript).code ? null : code);
  renderSliders();
  updateRevert();
  if (sim.started) simStatus('Edited: Reset, then Run to use the new code.');
});
$('sim-code').addEventListener('keydown', (e) => {
  if (e.key !== 'Tab') return; // Tab indents instead of leaving the editor
  e.preventDefault();
  e.target.setRangeText('  ', e.target.selectionStart, e.target.selectionEnd, 'end');
  e.target.dispatchEvent(new Event('input'));
});
$('sim-revert').addEventListener('click', () => {
  store.set(codeKey(simScript), null);
  store.set(sliderKey(simScript), null);
  selectScript(simScript);
});
$('sim-play').addEventListener('click', simPlayPause);
$('sim-step').addEventListener('click', simStep);
$('sim-reset').addEventListener('click', simReset);
$('sim-speed').addEventListener('input', (e) => setSpeed(Number(e.target.value)));

viewer.onFrame = (dt) => {
  realTime += dt;
  if (simMode) sim.advance(dt * sim.speed);
  stepGlow(dt);
  if (listDirty) updateSiteList();
  if (simMode) {
    $('sim-clock').textContent = `t = ${sim.t.toFixed(2)} s`;
    if (simUiDirty) updateSimUi();
  }
};

initSim();

// ---------------------------------------------------------------------------
// Collapsible panel sections: click a heading. Which ones are closed is remembered.
// ---------------------------------------------------------------------------

function initFolds() {
  let closed = [];
  try {
    closed = JSON.parse(store.get('molview:collapsed')) ?? [];
  } catch {
    // unreadable: everything starts open
  }
  const folds = new Set(closed);
  for (const section of document.querySelectorAll('#panel > section[id]')) {
    const head = section.querySelector(':scope > .label');
    if (!head) continue;
    const apply = () => {
      const shut = folds.has(section.id);
      section.classList.toggle('collapsed', shut);
      head.setAttribute('aria-expanded', String(!shut));
    };
    const toggle = () => {
      if (!folds.delete(section.id)) folds.add(section.id);
      store.set('molview:collapsed', JSON.stringify([...folds]));
      apply();
    };
    head.classList.add('fold');
    head.tabIndex = 0;
    head.setAttribute('role', 'button');
    head.title = 'Show / hide this section';
    head.addEventListener('click', (e) => {
      toggle();
      if (e.detail) head.blur(); // a mouse click must not keep the focus: Space still runs the simulation
    });
    head.addEventListener('keydown', (e) => {
      if (e.key !== 'Enter' && e.key !== ' ') return;
      e.preventDefault();
      e.stopPropagation();
      toggle();
    });
    apply();
  }
}
initFolds();

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

// The app tabs along the top: the viewer, the pulse experiment and the intensity sweep. The viewer stops drawing
// while it's hidden.
initTabs((page) => viewer.setPaused(page !== 'viewer'));
