// The intensity sweep, ported from sim/sweep_gui.py, the intensity command in sim/tta_sim.py and sim/plot_intensity.py:
// photons and excitations per molecule over one pulse, against kExc. Exact points are solved in a worker
// (pulseTotals); Monte Carlo points run in a pool of workers until the precision or the time limit. Runs are kept in
// this browser and can be saved as files.

import uPlot from 'uplot';
import 'uplot/dist/uPlot.min.css';
import { EVENTS, N_EVENTS, PULSE_DEFAULTS, resolve } from '../sim/rules.js';
import ExactWorker from './exact-worker.js?worker&inline';
import McWorker from './mc-worker.js?worker&inline';
import { defaultValues, FIELDS, fmt, int, kexcPerMw, parseValues, rulesFor, TTA_FIELDS } from './model.js';
import { deleteRun, escapeHtml, fromCsv, listRuns, putRun, saveLocation, saveRun, stamp, toCsv } from './store.js';

const EXACT_MAX_N = 8;
const THREADS = Math.min(8, Math.max(1, Math.floor((navigator.hardwareConcurrency || 2) / 2)));
const BATCH = 100_000; // Monte Carlo molecules per round, shared between the workers
const METHODS = { exact: 'Exact solve', mc: 'Monte Carlo', both: 'Both at every point' };
const KINDS = {
  counts: 'Photons per molecule, linear',
  counts_log: 'Photons per molecule, log, with slope',
  yield: 'Quantum yield, linear',
  yield_log: 'Quantum yield, log kExc',
};
const LOG_KEXC = new Set(['counts_log', 'yield_log']);
const CONTEXT = ['n', 'kHop', 'f', 'tripletYield', 'tauS', 'tauSens', 'kDecay', 'fretRadius', 'duration'];
const ALWAYS_TOGETHER = new Set(['kExc', 'excRatio']);
const SHADES = ['#cfe3fb', '#8fbcef', '#5598e7', '#2f74c8', '#1f5596']; // one hue, light to dark (on a dark page)
const CATEGORICAL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4'];
const PERYLENE = '#11d8f2';
const DPA = '#2a78d6';
const MOST_CURVES = 5;
const PER_MOLECULE = {
  excitations_per_molecule: 'excitation', isc_per_molecule: 'isc', sens_singlets_lost_per_molecule: 'sens_singlet_lost',
  back_transfers_per_molecule: 'back_transfer', tta_singlet_per_molecule: 'tta_singlet',
  tta_failed_per_molecule: 'tta_failed', triplets_decayed_per_molecule: 'triplet_lost',
};
const MC_ONLY = ['molecules', 'mean_photons_sq', 'mean_excitations_sq', 'mean_photons_x_excitations', 'photons_relative_error'];

/** The raw counts of a point, all per molecule, from an [event * sites + site] table (intensity_counts in tta_sim.py). */
function intensityCounts(counts, sites, molecules) {
  const at = (event, s) => counts[EVENTS.indexOf(event) * sites + s] / molecules;
  const sum = (event) => Array.from({ length: sites }, (_, s) => at(event, s)).reduce((a, b) => a + b, 0);
  const photons = Array.from({ length: sites }, (_, s) => at('photon', s));
  return {
    photons_per_molecule: photons.reduce((a, b) => a + b, 0),
    perylene_photons_per_molecule: photons[sites - 1],
    dpa_photons_per_molecule: photons.slice(0, -1).reduce((a, b) => a + b, 0),
    ...Object.fromEntries(Object.entries(PER_MOLECULE).map(([name, event]) => [name, sum(event)])),
  };
}

const box = (name, label, unit, value, hint = '') =>
  `<label for="sweep-${name}" title="${hint}">${label}</label><input type="text" id="sweep-${name}" data-box="${name}"
    value="${value}" autocomplete="off" title="${hint}"><span class="unit">${unit}</span>`;

const TEMPLATE = (v) => `
<header class="lab-head">
  <div><div class="lab-title">Intensity sweep</div>
  <div class="lab-sub">Photons and excitations per molecule, per pulse, against kExc</div></div>
  <div class="spacer"></div>
  <div class="lamp" data-id="lamp"></div>
</header>
<div class="tabs2" role="tablist">
  <button role="tab" data-tab="run" aria-selected="true">Run</button>
  <button role="tab" data-tab="plot" aria-selected="false">Plot</button>
</div>
<div class="lab-body" data-pane="run" style="padding-top:12px">
  <aside class="lab-side">
    <section class="panel">
      <div class="section">Sweep of kExc</div>
      <div class="grid">
        ${box('from', 'From', '/ns', '0.001')}
        ${box('to', 'To', '/ns', '100')}
        ${box('points', 'Points', '', '31')}
        <label for="sweep-spacing">Spacing</label>
        <select id="sweep-spacing" data-box="spacing" style="grid-column: span 2"><option value="log">spaced by ratio (log)</option><option value="linear">evenly spaced (linear)</option></select>
        ${box('duration', 'Pulse length', 'ns', fmt(v.duration, 6), 'One value, or several (5,50) for one curve each')}
      </div>
      <div class="caption mono" data-id="rangeText" style="color:var(--text)"></div>
      <div class="note">Log spacing suits the log plot and the slope. For the linear plots, an even spacing fills the axis better.</div>
    </section>
    <section class="panel">
      <div class="section">Method</div>
      <div class="grid wide">
        <label for="sweep-method">Method</label>
        <select id="sweep-method" data-box="method" style="grid-column: span 2">${Object.entries(METHODS).map(([k, t]) => `<option value="${k}">${t}</option>`).join('')}</select>
        ${box('checks', 'Monte Carlo checks', 'pts', '3', 'With the exact solve: also run Monte Carlo at this many points, as a check')}
        ${box('precision', 'Monte Carlo precision', '%', '1', 'Each Monte Carlo point runs until photons per molecule is known to this relative standard error')}
        ${box('max_seconds', 'Most time per point', 's', '600')}
        ${box('name', 'Run name', '', 'kexc-sweep')}
      </div>
      <div class="caption mono" data-id="methodText" style="color:var(--text)"></div>
      <div class="row" style="margin-top:12px"><button class="fire" data-id="run" style="flex:3">RUN SWEEP</button><button data-id="stop" disabled>STOP</button></div>
      <div class="note">A stopped sweep keeps every point it finished.</div>
    </section>
    <section class="panel">
      <div class="section">Sample <span class="hint">several values give one curve each</span></div>
      <div class="grid wide">${TTA_FIELDS.map((k) => box(k, FIELDS[k][0], FIELDS[k][1], fmt(v[k], 6), 'One value, or several (a,b or a:b:count) for one curve each')).join('')}</div>
    </section>
  </aside>
  <div class="lab-main">
    <section class="panel">
      <div class="section">Progress</div>
      <div class="mono" data-id="progress" style="color:var(--text)">Not running.</div>
      <div class="bar"><div data-id="bar"></div></div>
      <div class="note mono" data-id="command" style="margin-top:8px"></div>
    </section>
    <section class="panel plot-panel"><div class="plot-box" data-id="preview"><div class="empty">Photons per molecule appear here as the sweep runs.</div></div></section>
    <section class="panel" style="padding:8px 14px"><div class="log" data-id="log" style="margin:0"></div></section>
  </div>
</div>
<div class="lab-body" data-pane="plot" hidden style="padding-top:12px; grid-template-columns: minmax(0, 1fr) 330px">
  <div class="lab-main">
    <section class="panel" style="display:flex; flex-direction:column; flex: 0 1 34%; min-height: 160px">
      <div class="section">Intensity runs <span class="hint">click to pick · Ctrl or Shift for several</span></div>
      <div class="runs-box"><table class="runs"><thead><tr><th>Run</th><th>Date</th><th>Rows</th><th>Swept</th><th>Method</th><th>n</th><th>Pulse, ns</th></tr></thead><tbody data-id="runs"></tbody></table></div>
    </section>
    <section class="panel plot-panel"><div class="plot-title" data-id="plotTitle"></div><div class="plots" data-id="plots"><div class="empty">Pick a run, then a plot.</div></div></section>
  </div>
  <aside class="lab-side">
    <section class="panel stack">
      <div class="section" style="margin:0">Plots</div>
      ${Object.entries(KINDS).map(([k, t]) => `<button data-kind="${k}">${t}</button>`).join('')}
      <label class="check"><input type="checkbox" data-id="logx"> Log x axis on the linear plots</label>
      <label class="check" style="margin-top:0"><input type="checkbox" data-id="perMw"> Laser power on the x axis, not kExc</label>
      <div class="grid" style="grid-template-columns: 1fr 100px 30px"><label>kExc per mW</label><input type="text" data-id="kpm" value="${fmt(kexcPerMw(), 4)}"><span class="unit">/ns</span></div>
      <div class="note">The default is the pulse page's, from its beam settings (1 µm beam, ε 50,000 /M/cm at 635 nm).</div>
    </section>
    <section class="panel stack">
      <div class="section" style="margin:0">Files</div>
      <button data-id="saveRun">Save the run</button>
      <button data-id="import">Import a results.csv…</button>
      <input type="file" data-id="file" accept=".csv,.json" multiple hidden>
      <button data-id="delete">Delete from this browser</button>
      <div data-id="where"></div>
      <div class="note">Runs live in this browser until you delete them. Saving writes results.csv and params.json, like
        sim/data. Import takes a results.csv from here or from the Python app (with its params.json, if you pick both).</div>
    </section>
  </aside>
</div>`;

export function mount(root) {
  const values = defaultValues();
  root.innerHTML = TEMPLATE(values);
  const $ = (id) => root.querySelector(`[data-id="${id}"]`);
  const boxes = Object.fromEntries([...root.querySelectorAll('[data-box]')].map((el) => [el.dataset.box, el]));
  $('where').append(saveLocation());

  // ---- tabs ----
  const showTab = (tab) => {
    for (const b of root.querySelectorAll('[data-tab]')) b.setAttribute('aria-selected', b.dataset.tab === tab);
    for (const pane of root.querySelectorAll('[data-pane]')) pane.hidden = pane.dataset.pane !== tab;
    if (tab === 'plot') refreshRuns().then(drawPicked);
    else previewChart?.resize();
  };
  for (const b of root.querySelectorAll('[data-tab]')) b.addEventListener('click', () => showTab(b.dataset.tab));

  // ---- reading the boxes ----
  function read() {
    const b = Object.fromEntries(Object.entries(boxes).map(([k, el]) => [k, el.value.trim()]));
    const low = Number(b.from);
    const high = Number(b.to);
    const count = Number(b.points);
    if (!Number.isFinite(low)) throw new Error('from');
    if (!Number.isFinite(high) || !Number.isInteger(count)) throw new Error('to');
    const log = b.spacing === 'log';
    if (!(0 < low && low < high && count >= 2) && !(low === high && count === 1)) throw new Error('to');
    const kExc = parseValues('kExc', `${low}:${high}:${log ? 'log' : ''}${count}`);
    const lists = { kExc };
    for (const name of ['duration', ...TTA_FIELDS]) {
      let found;
      try {
        found = parseValues(name, b[name]);
      } catch {
        throw new Error(name);
      }
      const [, , lo, hi] = FIELDS[name];
      if (found.some((x) => x < lo || x > hi)) throw new Error(name);
      lists[name] = found;
    }
    const method = b.method;
    const checks = method === 'exact' ? Number(b.checks) : 0;
    if (method === 'exact' && !(Number.isInteger(checks) && checks >= 0)) throw new Error('checks');
    const precision = Number(b.precision) / 100;
    const maxSeconds = Number(b.max_seconds);
    if ((method !== 'exact' || checks) && !(precision > 0 && precision < 1 && maxSeconds > 0)) throw new Error('precision');
    const name = b.name.replace(/[^\w.-]+/g, '-') || 'intensity';
    const curves = Object.entries(lists).filter(([k, l]) => k !== 'kExc' && l.length > 1).map(([k]) => k);
    const points = Object.values(lists).reduce((s, l) => s * l.length, 1);
    return { lists, method, checks, precision, maxSeconds, name, curves, points, mostN: Math.max(...lists.n),
      command: `python sim/tta_sim.py intensity --kExc ${low}:${high}:${log ? 'log' : ''}${count} `
        + ['duration', ...TTA_FIELDS].map((k) => `--${k} ${b[k]}`).join(' ') + ` --method ${method}`
        + (method === 'exact' ? ` --checks ${checks}` : '')
        + (method !== 'exact' || checks ? ` --precision ${precision} --max-seconds ${maxSeconds}` : '') + ` --name ${name}` };
  }

  let job = null; // the sweep that's running
  function check() {
    for (const el of Object.values(boxes)) el.classList.remove('bad');
    let r;
    try {
      r = read();
    } catch (wrong) {
      boxes[wrong.message]?.classList.add('bad');
      $('rangeText').textContent = 'Check the box marked in red.';
      $('run').disabled = true;
      return null;
    }
    const tooBig = r.method !== 'mc' && r.mostN > EXACT_MAX_N;
    boxes.checks.disabled = r.method !== 'exact';
    const usesMc = r.method !== 'exact' || r.checks > 0;
    boxes.precision.disabled = !usesMc;
    boxes.max_seconds.disabled = !usesMc;
    $('rangeText').textContent = `${r.points} point${r.points > 1 ? 's' : ''}, ${r.curves.length ? `one curve per value of ${r.curves.join(', ')}` : 'one curve'}.`;
    const text = $('methodText');
    text.classList.toggle('warn', tooBig);
    if (tooBig) text.textContent = `The exact solve only goes up to n = ${EXACT_MAX_N}. Pick Monte Carlo.`;
    else if (r.method === 'mc') text.textContent = 'Monte Carlo: each point runs until the precision or the time limit. Low kExc needs many molecules, since a photon needs two absorptions in one pulse.';
    else text.textContent = 'Exact: no noise at any kExc. A few seconds a point at n = 8 (longer for the first, and for long, bright pulses); well under a second at n of 6 or less.';
    $('run').disabled = !!job || tooBig;
    if (!job) $('command').textContent = r.command;
    return r;
  }
  for (const el of Object.values(boxes)) el.addEventListener('input', check);
  for (const el of Object.values(boxes)) el.addEventListener('change', check);
  check();

  // ---- running a sweep ----
  let exactWorker = null;
  let pool = [];
  const mcPool = () => {
    if (!pool.length) pool = Array.from({ length: THREADS }, () => new McWorker());
    return pool;
  };
  const stopWorkers = () => {
    exactWorker?.terminate();
    exactWorker = null;
    pool.forEach((w) => w.terminate());
    pool = [];
  };

  const solveExact = (p, R, rules) => new Promise((resolve, reject) => {
    exactWorker ??= new ExactWorker();
    const key = Math.random();
    exactWorker.onmessage = ({ data }) => {
      if (data.key !== key) return;
      if (data.type === 'totals') resolve(data.counts);
      else if (data.type === 'failed') reject(new Error(data.error));
    };
    exactWorker.postMessage({ type: 'totals', key, p, R, rules });
  });

  /** Monte Carlo for one point: rounds of BATCH molecules until the precision or the time limit (intensity_mc). */
  async function monteCarlo(p, R, sites, cfg, seed, report) {
    const started = performance.now();
    let molecules = 0;
    const counts = new Float64Array(N_EVENTS * sites);
    const sums = [0, 0, 0];
    let error = Infinity;
    const workers = mcPool();
    let round = 0;
    for (;;) {
      const each = Math.ceil(BATCH / workers.length);
      const results = await Promise.all(workers.map((w, i) => new Promise((resolve) => {
        const id = `${seed}-${round}-${i}`;
        w.onmessage = ({ data }) => data.id === id && resolve(data);
        w.postMessage({ type: 'batch', id, R, sites, duration: p.duration, molecules: each,
          seed: (seed + 7919 * round + 104729 * i) >>> 0 });
      })));
      if (!job || job.stopped) throw new Error('stopped');
      round++;
      for (const r of results) {
        molecules += r.molecules;
        for (let k = 0; k < counts.length; k++) counts[k] += r.counts[k];
        for (let k = 0; k < 3; k++) sums[k] += r.sums[k];
      }
      let photons = 0;
      for (let s = 0; s < sites; s++) photons += counts[EVENTS.indexOf('photon') * sites + s];
      const mean = photons / molecules;
      if (mean > 0) error = Math.sqrt(Math.max(sums[0] / molecules - mean * mean, 0) / molecules) / mean;
      report(`${int(molecules)} molecules, error ${Number.isFinite(error) ? `${(100 * error).toFixed(2)}%` : '–'}`);
      if (error <= cfg.precision || (performance.now() - started) / 1000 > cfg.maxSeconds) break;
    }
    return { ...intensityCounts(counts, sites, molecules), molecules, mean_photons_sq: sums[0] / molecules,
      mean_excitations_sq: sums[1] / molecules, mean_photons_x_excitations: sums[2] / molecules,
      photons_relative_error: error };
  }

  async function run() {
    const cfg = check();
    if (!cfg || job) return;
    const { lists } = cfg;
    // every combination, kExc changing fastest so each curve fills in turn (and the after-pulse solve is reused)
    const others = ['duration', ...TTA_FIELDS];
    let combos = [{}];
    for (const name of others) combos = combos.flatMap((c) => lists[name].map((v) => ({ ...c, [name]: v })));
    const points = combos.flatMap((c) => lists.kExc.map((kExc) => {
      const p = resolve({ ...PULSE_DEFAULTS, ...c, kExc });
      return Object.fromEntries(Object.keys(PULSE_DEFAULTS).map((k) => [k, p[k]]));
    }));
    let checks = new Set();
    if (cfg.method === 'exact' && cfg.checks) {
      checks = new Set(Array.from({ length: cfg.checks }, (_, i) => Math.round(((i + 1) * (points.length - 1)) / (cfg.checks + 1))));
    }
    const swept = Object.fromEntries(Object.entries(lists).filter(([, l]) => l.length > 1).map(([k, l]) => [k, [...new Set(l)].sort((a, b) => a - b)]));
    const { day, time } = stamp();
    const runs = await listRuns();
    let id = `${day}_${cfg.name}`;
    for (let k = 2; runs.some((r) => r.id === id); k++) id = `${day}_${cfg.name}_${k}`;
    const seed = crypto.getRandomValues(new Uint32Array(1))[0] % 1_000_000_000;
    const record = {
      id, date: `${day} ${time}`, rows: [],
      params: { mode: 'intensity', command: cfg.command, date: `${day} ${time}`, units: 'times in ns, rates per ns',
        defaults: PULSE_DEFAULTS, swept, points: points.length, method: cfg.method, checks: [...checks].sort((a, b) => a - b),
        precision: cfg.precision, max_seconds: cfg.maxSeconds, batch: BATCH, seed, threads: THREADS,
        source: 'Molecule Viewer website (src/sim/, a port of sim/)' },
    };
    job = { stopped: false, record };
    $('run').disabled = true;
    $('stop').disabled = false;
    $('log').textContent = '';
    setLamp('#f2a93b', 'RUNNING');
    const started = performance.now();
    await putRun(record);
    try {
      for (let i = 0; i < points.length; i++) {
        const p = points[i];
        const { rules, R, sites } = rulesFor(p);
        const hows = { exact: ['exact'], mc: ['mc'], both: ['exact', 'mc'] }[cfg.method].concat(checks.has(i) ? ['mc'] : []);
        for (const how of hows) {
          const began = performance.now();
          const head = `Point ${i + 1} of ${points.length}, ${how === 'exact' ? 'exact' : 'Monte Carlo'}, kExc ${fmt(p.kExc, 4)} /ns`;
          $('progress').textContent = `${head}…`;
          let found;
          if (how === 'exact') {
            const counts = await solveExact(p, R, rules);
            if (job.stopped) throw new Error('stopped');
            found = { ...intensityCounts(counts, sites, 1), ...Object.fromEntries(MC_ONLY.map((k) => [k, NaN])) };
          } else {
            found = await monteCarlo(p, R, sites, cfg, seed + 1_000_003 * i, (s) => ($('progress').textContent = `${head}: ${s}`));
          }
          const row = { point: i, method: how, ...p, ...found, wall_s: (performance.now() - began) / 1000 };
          record.rows.push(row);
          await putRun(record); // a stopped sweep keeps what it finished
          $('log').textContent += `point ${i + 1}/${points.length} ${how} kExc ${fmt(p.kExc, 4)} photons ${fmt(row.photons_per_molecule, 4)} in ${row.wall_s.toFixed(1)} s\n`;
          $('log').scrollTop = 1e9;
          $('bar').style.width = `${(100 * (i + 1)) / points.length}%`;
          drawPreview(record);
        }
      }
      $('progress').textContent = `Done: ${points.length} point${points.length > 1 ? 's' : ''} in ${((performance.now() - started) / 1000).toFixed(1)} s. Saved in this browser as ${id}; see the Plot tab.`;
      setLamp('#2fd27a', 'DONE');
    } catch (error) {
      if (error.message === 'stopped') {
        $('progress').textContent = `Stopped. ${record.rows.length} row${record.rows.length === 1 ? '' : 's'} kept as ${id}.`;
        setLamp(null, 'STOPPED');
      } else {
        $('progress').textContent = `Failed: ${error.message}`;
        setLamp('#ff2100', 'FAILED');
      }
    }
    if (!record.rows.length) await deleteRun(id);
    job = null;
    $('stop').disabled = true;
    check();
  }
  $('run').addEventListener('click', run);
  $('stop').addEventListener('click', () => {
    if (!job) return;
    job.stopped = true;
    stopWorkers(); // a solve or batch can't be interrupted, so end the workers
    $('progress').textContent = 'Stopping…';
  });

  const setLamp = (colour, text) => {
    const el = $('lamp');
    el.innerHTML = `<i style="color:${colour || 'var(--dim)'}">●</i>${text}`;
    el.classList.toggle('lit', !!colour);
  };
  setLamp(null, 'IDLE');

  // ---- the preview while running ----
  let previewChart = null;
  function drawPreview(record) {
    previewChart?.destroy();
    const holder = $('preview');
    holder.innerHTML = '';
    try {
      previewChart = chart(holder, [record], 'preview', {});
    } catch {
      holder.innerHTML = '<div class="empty">Photons per molecule appear here as the sweep runs.</div>';
    }
  }

  // ---- the plot tab ----
  let runs = [];
  const picked = new Set();
  let kind = null;
  let charts = [];
  let lastClicked = null;

  async function refreshRuns() {
    runs = (await listRuns()).filter((r) => r.rows.length).sort((a, b) => (b.date || '').localeCompare(a.date || ''));
    for (const id of [...picked]) if (!runs.some((r) => r.id === id)) picked.delete(id);
    const body = $('runs');
    body.innerHTML = runs.length ? '' : '<tr><td colspan="7" style="color:var(--dim)">No runs yet. Run a sweep, or import a results.csv.</td></tr>';
    runs.forEach((r, i) => {
      const rows = r.rows;
      const tr = document.createElement('tr');
      const shared = (k) => {
        const s = new Set(rows.map((x) => x[k]));
        return s.size === 1 ? fmt([...s][0], 4) : [...s].map((v) => fmt(v, 3)).join(', ');
      };
      const methods = [...new Set(rows.map((x) => x.method))].join(' + ');
      tr.innerHTML = `<td>${escapeHtml(r.id)}</td><td>${escapeHtml((r.date || '').slice(0, 16))}</td><td>${rows.length}</td>
        <td>${escapeHtml(Object.keys(r.params?.swept ?? {}).join(' '))}</td><td>${escapeHtml(methods)}</td><td>${shared('n')}</td><td>${shared('duration')}</td>`;
      tr.classList.toggle('sel', picked.has(r.id));
      tr.addEventListener('click', (e) => {
        if (e.shiftKey && lastClicked != null) {
          const [a, b] = [lastClicked, i].sort((x, y) => x - y);
          for (let k = a; k <= b; k++) picked.add(runs[k].id);
        } else if (e.ctrlKey || e.metaKey) {
          if (picked.has(r.id)) picked.delete(r.id);
          else picked.add(r.id);
        } else {
          picked.clear();
          picked.add(r.id);
        }
        lastClicked = i;
        for (const [k, row] of [...body.children].entries()) row.classList.toggle('sel', picked.has(runs[k]?.id));
        kind ??= 'counts_log';
        drawPicked();
      });
      tr.addEventListener('dblclick', () => {
        kind = 'counts_log';
        drawPicked();
      });
      body.append(tr);
    });
  }

  function drawPicked() {
    charts.forEach((c) => c.destroy());
    charts = [];
    const holder = $('plots');
    holder.innerHTML = '';
    for (const b of root.querySelectorAll('[data-kind]')) b.setAttribute('aria-pressed', b.dataset.kind === kind);
    const chosen = runs.filter((r) => picked.has(r.id));
    if (!chosen.length || !kind) {
      $('plotTitle').textContent = '';
      holder.innerHTML = `<div class="empty">${runs.length ? 'Pick a run, then a plot.' : 'No runs yet.'}</div>`;
      return;
    }
    const perMw = $('perMw').checked ? Number($('kpm').value) : null;
    try {
      charts = [chart(holder, chosen, kind, { logx: $('logx').checked, perMw: perMw > 0 ? perMw : null, title: $('plotTitle') })];
    } catch (error) {
      $('plotTitle').textContent = '';
      holder.innerHTML = `<div class="empty">${escapeHtml(error.message)}</div>`;
    }
  }
  for (const b of root.querySelectorAll('[data-kind]'))
    b.addEventListener('click', () => {
      kind = b.dataset.kind;
      drawPicked();
    });
  for (const id of ['logx', 'perMw', 'kpm']) $(id).addEventListener('change', drawPicked);

  $('saveRun').addEventListener('click', async () => {
    const chosen = runs.filter((r) => picked.has(r.id));
    if (!chosen.length) return;
    for (const r of chosen) {
      try {
        await saveRun(r.id, [{ name: 'results.csv', data: toCsv(r.rows) }, { name: 'params.json', data: JSON.stringify(r.params, null, 2) }]);
      } catch {
        return;
      }
    }
  });
  $('delete').addEventListener('click', async () => {
    const chosen = runs.filter((r) => picked.has(r.id));
    if (!chosen.length || !confirm(`Delete ${chosen.length === 1 ? chosen[0].id : `${chosen.length} runs`} from this browser? Saved files aren't touched.`)) return;
    for (const r of chosen) await deleteRun(r.id);
    picked.clear();
    await refreshRuns();
    drawPicked();
  });
  $('import').addEventListener('click', () => $('file').click());
  $('file').addEventListener('change', async () => {
    const files = [...$('file').files];
    $('file').value = '';
    const csv = files.find((f) => f.name.endsWith('.csv'));
    if (!csv) return;
    const json = files.find((f) => f.name.endsWith('.json'));
    const rows = fromCsv(await csv.text()).map((r) => ({ ...r, method: String(r.method ?? 'exact') }));
    let params = json ? JSON.parse(await json.text()) : null;
    if (!params) {
      // no params.json: call swept whatever varies
      const swept = {};
      for (const k of ['kExc', ...CONTEXT]) {
        const s = [...new Set(rows.map((r) => r[k]))].filter((v) => Number.isFinite(v));
        if (s.length > 1) swept[k] = s.sort((a, b) => a - b);
      }
      params = { mode: 'intensity', swept, imported: csv.name };
    }
    const { day, time } = stamp();
    const base = csv.name.replace(/\.csv$/i, '') === 'results' ? `imported_${day}` : csv.name.replace(/\.csv$/i, '');
    let id = base;
    for (let k = 2; runs.some((r) => r.id === id); k++) id = `${base}_${k}`;
    await putRun({ id, date: params.date || `${day} ${time}`, rows, params });
    picked.clear();
    picked.add(id);
    kind ??= 'counts_log';
    await refreshRuns();
    drawPicked();
  });

  new ResizeObserver(() => charts.forEach((c) => c.resize())).observe($('plots'));
  new ResizeObserver(() => previewChart?.resize()).observe($('preview'));

  return {
    shown(on) {
      if (on) {
        previewChart?.resize();
        charts.forEach((c) => c.resize());
      }
    },
  };
}

// ---------------------------------------------------------------------------
// The figures: a port of plot_intensity.draw, drawn with uPlot on the dark page
// ---------------------------------------------------------------------------

const column = (rows, name) => rows.map((r) => r[name]);

/** [{ label, exact, mc, rows, colour }] for the picked runs, one per value of whatever each swept besides kExc. */
function curves(runs) {
  const found = [];
  let ordered = false;
  for (const run of runs) {
    const swept = Object.keys(run.params?.swept ?? {}).filter((k) => !ALWAYS_TOGETHER.has(k));
    const keys = [...new Set(run.rows.map((r) => JSON.stringify(swept.map((k) => r[k]))))].sort();
    for (const key of keys) {
      const vals = JSON.parse(key);
      const mine = run.rows.filter((r) => JSON.stringify(swept.map((k) => r[k])) === key);
      const parts = swept.map((k, i) => `${k} = ${fmt(vals[i], 4)}`);
      if (runs.length > 1) parts.unshift(run.id);
      const byK = (a, b) => a.kExc - b.kExc;
      found.push({ label: parts.join(', '), exact: mine.filter((r) => r.method === 'exact').sort(byK),
        mc: mine.filter((r) => r.method === 'mc').sort(byK), rows: mine });
    }
    ordered = runs.length === 1 && swept.length === 1;
  }
  if (!found.length) throw new Error('No rows to plot yet.');
  if (found.length > MOST_CURVES) throw new Error(`${found.length} curves: at most ${MOST_CURVES} can be told apart. Pick fewer runs.`);
  const colours = found.length === 1 ? [CATEGORICAL[0]]
    : ordered ? found.map((_, i) => SHADES[Math.round((i * (SHADES.length - 1)) / (found.length - 1))]) : CATEGORICAL;
  found.forEach((c, i) => (c.colour = colours[i]));
  return found;
}

function yieldSe(rows) {
  return rows.map((r) => {
    const N = r.molecules;
    const P = r.photons_per_molecule;
    const E = r.excitations_per_molecule;
    const Q = P / E;
    const varP = r.mean_photons_sq - P * P;
    const varE = r.mean_excitations_sq - E * E;
    const cov = r.mean_photons_x_excitations - P * E;
    return Math.sqrt(Math.max(varP - 2 * Q * cov + Q * Q * varE, 0) / N) / E;
  });
}
const photonsSe = (rows) => rows.map((r) => Math.sqrt(Math.max(r.mean_photons_sq - r.photons_per_molecule ** 2, 0) / r.molecules));

/** d log(photons) / d log(kExc), by numpy.gradient's second-order differences. */
function localSlope(rows) {
  const pts = rows.filter((r) => r.photons_per_molecule > 0).map((r) => [Math.log(r.kExc), Math.log(r.photons_per_molecule), r.kExc]);
  if (pts.length < 3) return null;
  const x = pts.map((p) => p[0]);
  const y = pts.map((p) => p[1]);
  const n = x.length;
  const g = x.map((_, i) => {
    if (i === 0) return (y[1] - y[0]) / (x[1] - x[0]);
    if (i === n - 1) return (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
    const h1 = x[i] - x[i - 1];
    const h2 = x[i + 1] - x[i];
    return (h1 * h1 * y[i + 1] - h2 * h2 * y[i - 1] + (h2 * h2 - h1 * h1) * y[i]) / (h1 * h2 * (h1 + h2));
  });
  return { x: pts.map((p) => p[2]), slope: g };
}

const shared = (found, name) => {
  const s = new Set(found.flatMap((c) => c.rows.map((r) => r[name])));
  return s.size === 1 ? [...s][0] : null;
};

const SUPER = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
const isDecade = (v) => v > 0 && Math.abs(Math.log10(v) - Math.round(Math.log10(v))) < 1e-6;
const decadeLabel = (v) => {
  const e = Math.round(Math.log10(v));
  return Math.abs(e) <= 3 ? String(Number(v.toPrecision(1))) : `10${String(e).split('').map((c) => SUPER[c]).join('')}`;
};
/** Props for a log axis: ticks, grid and labels at whole decades only (all splits if it spans under two decades). */
const logAxis = (scaleBy = 1) => ({
  filter: (u, splits) => (splits.filter(isDecade).length >= 2 ? splits.map((v) => (isDecade(v) ? v : null)) : splits),
  values: (u, splits) => splits.map((v) => (v == null ? '' : isDecade(v) ? decadeLabel(v * scaleBy) : fmt(v * scaleBy, 2))),
});

const AXIS = { stroke: '#7c8796', grid: { stroke: 'rgba(255,255,255,0.05)', width: 1 }, ticks: { stroke: 'rgba(255,255,255,0.12)', width: 1 },
  font: '12px Segoe UI, system-ui, sans-serif', labelFont: '12px Segoe UI, system-ui, sans-serif' };

/**
 * Draw `kind` for `runs` into `holder`. Returns { destroy, resize }. kind 'preview' is the run tab's log-log of photons,
 * one colour per curve.
 */
function chart(holder, runs, kind, { logx = false, perMw = null, title = null }) {
  const found = curves(runs);
  if (kind === 'preview') found.forEach((c, i) => (c.colour = [PERYLENE, ...CATEGORICAL.slice(1)][i % 5]));
  const isLogK = kind === 'preview' || LOG_KEXC.has(kind);
  logx = logx || isLogK;
  const split = found.length === 1 && kind !== 'preview';
  const scale = perMw ? 1 / perMw : 1;
  const xlabel = perMw ? 'Laser power, mW' : 'kExc, /ns';
  const isYield = kind === 'yield' || kind === 'yield_log';
  const logy = kind === 'counts_log' || kind === 'preview';
  const yOf = isYield ? (rows) => rows.map((r) => (100 * r.photons_per_molecule) / r.excitations_per_molecule) : (rows) => column(rows, 'photons_per_molecule');
  const seOf = isYield ? (rows) => yieldSe(rows).map((v) => 100 * v) : photonsSe;
  const ok = (v) => (Number.isFinite(v) && (!logy || v > 0) ? v : null);

  // one table per line: exact, Monte Carlo, and (one curve) the perylene and DPA parts
  const tables = [];
  const series = [{ label: perMw ? 'Power, mW' : 'kExc, /ns', value: (u, v) => (v == null ? '–' : fmt(v, 4)) }];
  const errorBars = [];
  for (const c of found) {
    for (const how of ['exact', 'mc']) {
      const rows = c[how];
      if (!rows.length) continue;
      const method = how === 'exact' ? 'Exact solve' : 'Monte Carlo';
      const label = !c.label ? method : c.label + (how === 'exact' || !c.exact.length ? '' : ', Monte Carlo');
      tables.push([column(rows, 'kExc').map((k) => k * scale), yOf(rows).map(ok)]);
      if (how === 'exact') series.push({ label, stroke: c.colour, width: 2, spanGaps: true, points: { show: true, size: 5, fill: c.colour }, value: (u, v) => fmt(v, 4) });
      else {
        const lone = !c.exact.length;
        series.push({ label, stroke: c.colour, width: lone ? 2 : 0, spanGaps: true, paths: lone ? undefined : () => null,
          points: { show: true, size: 9, width: 1.6, fill: '#0d1117', stroke: c.colour }, value: (u, v) => fmt(v, 4) });
        if (kind !== 'preview') errorBars.push({ series: series.length - 1, se: seOf(rows), colour: c.colour, x: column(rows, 'kExc').map((k) => k * scale) });
      }
      if (split && how === 'exact' && !isYield) {
        tables.push([column(rows, 'kExc').map((k) => k * scale), column(rows, 'perylene_photons_per_molecule').map(ok)]);
        series.push({ label: 'from the perylene', stroke: PERYLENE, width: 1.5, dash: [8, 4], spanGaps: true, points: { show: false }, value: (u, v) => fmt(v, 4) });
        tables.push([column(rows, 'kExc').map((k) => k * scale), column(rows, 'dpa_photons_per_molecule').map(ok)]);
        series.push({ label: 'from the DPAs', stroke: DPA, width: 1.5, dash: [2, 4], spanGaps: true, points: { show: false }, value: (u, v) => fmt(v, 4) });
      }
    }
  }
  const data = uPlot.join(tables);

  const duration = shared(found, 'duration');
  let topAxis = null;
  const topTicks = (by) => (logx ? { filter: logAxis().filter, values: (u, sp) => sp.map((v) => (v == null ? '' : fmt(v * by, 3))) }
    : { values: (u, sp) => sp.map((v) => fmt(v * by, 3)) });
  if (perMw) topAxis = { ...AXIS, scale: 'x', side: 0, grid: { show: false }, label: `kExc, /ns, at ${fmt(perMw, 4)} /ns per mW`, ...topTicks(perMw) };
  else if (duration && kind !== 'preview') topAxis = { ...AXIS, scale: 'x', side: 0, grid: { show: false }, label: 'kExc × pulse: absorptions per pulse if the sensitizer were always free', ...topTicks(duration) };

  const triplet = shared(found, 'tripletYield');
  const f = shared(found, 'f');
  const limit = isYield && triplet != null && f != null ? (100 * triplet * f) / (1 + f) : null;

  const hooks = {
    draw: [(u) => {
      const { ctx } = u;
      const dpr = devicePixelRatio;
      ctx.save();
      const { left, top, width, height } = u.bbox;
      ctx.beginPath();
      ctx.rect(left, top, width, height);
      ctx.clip();
      for (const bar of errorBars) {
        if (!u.series[bar.series].show) continue;
        ctx.strokeStyle = bar.colour;
        ctx.lineWidth = 1.2 * dpr;
        const ys = data[bar.series];
        bar.x.forEach((x, i) => {
          const k = data[0].indexOf(x);
          const y = ys[k];
          if (y == null || !Number.isFinite(bar.se[i])) return;
          const px = u.valToPos(x, 'x', true);
          const lo = logy ? Math.max(y - bar.se[i], y * 1e-3) : y - bar.se[i];
          const y0 = u.valToPos(lo, 'y', true);
          const y1 = u.valToPos(y + bar.se[i], 'y', true);
          ctx.beginPath();
          ctx.moveTo(px, y0);
          ctx.lineTo(px, y1);
          ctx.moveTo(px - 3 * dpr, y0);
          ctx.lineTo(px + 3 * dpr, y0);
          ctx.moveTo(px - 3 * dpr, y1);
          ctx.lineTo(px + 3 * dpr, y1);
          ctx.stroke();
        });
      }
      if (limit != null) {
        const py = u.valToPos(limit, 'y', true);
        ctx.setLineDash([4 * dpr, 3 * dpr]);
        ctx.strokeStyle = '#7c8796';
        ctx.lineWidth = dpr;
        ctx.beginPath();
        ctx.moveTo(left, py);
        ctx.lineTo(left + width, py);
        ctx.stroke();
        ctx.setLineDash([]);
        ctx.fillStyle = '#7c8796';
        ctx.font = `${11 * dpr}px Segoe UI, system-ui, sans-serif`;
        ctx.textAlign = 'right';
        ctx.fillText(`Limit ${fmt(limit, 3)}%: triplet yield × f/(1+f)`, left + width - 6 * dpr, py - 5 * dpr);
      }
      ctx.restore();
    }],
  };

  const yLabel = kind === 'counts_log' ? 'Photons per molecule (log scale)' : kind === 'counts' ? 'Photons per molecule, per pulse'
    : kind === 'preview' ? 'Photons per molecule' : 'Quantum yield (%), photons per excitation';
  let yMax = null;
  if (limit != null) yMax = limit * 1.12;
  const xs = data[0].filter((v) => v > 0);
  const xRange = logx ? [Math.min(...xs) / 1.3, Math.max(...xs) * 1.3] : [0, Math.max(...data[0]) * 1.03];

  const parts = kind === 'counts_log' ? [0.72, 0.28] : [1];
  holder.innerHTML = '';
  Object.assign(holder.style, { display: 'flex', flexDirection: 'column' });
  const boxes = parts.map(() => {
    const div = document.createElement('div');
    div.className = 'plot-box';
    holder.append(div);
    return div;
  });
  boxes.forEach((b, i) => (b.style.flex = String(parts[i])));
  const sizeOf = (b) => ({ width: Math.max(b.clientWidth, 200), height: Math.max(b.clientHeight - (b === boxes[0] ? 30 : 0), 100) });
  const sync = uPlot.sync(`sweep-${Math.random()}`);
  const plots = [];
  plots.push(new uPlot({
    ...sizeOf(boxes[0]),
    legend: { show: kind !== 'preview' || found.length > 1, live: true },
    cursor: { sync: { key: sync.key }, points: { size: 8 } },
    scales: {
      x: { time: false, distr: logx ? 3 : 1, auto: false, range: () => xRange },
      y: { distr: logy ? 3 : 1, range: (u, min, max) => (logy ? [min / 1.6, max * 1.6] : [0, yMax ?? max * 1.08]) },
    },
    axes: [
      { ...AXIS, ...(logx ? logAxis() : {}), label: parts.length > 1 ? undefined : `${xlabel}${logx ? ' (log scale)' : ''}`,
        labelSize: parts.length > 1 ? 0 : 22, ...(parts.length > 1 ? { values: (u, sp) => sp.map(() => ''), size: 10 } : {}) },
      { ...AXIS, label: yLabel, labelSize: 22, size: 70, ...(logy ? logAxis() : {}) },
      ...(topAxis ? [topAxis] : []),
    ],
    series,
    hooks,
  }, data, boxes[0]));

  if (kind === 'counts_log') {
    const slopes = found.map((c) => ({ c, s: localSlope(c.exact.length ? c.exact : c.mc) })).filter((x) => x.s);
    const t2 = slopes.map(({ s }) => [s.x.map((k) => k * scale), s.slope]);
    const d2 = t2.length ? uPlot.join(t2) : [[xRange[0], xRange[1]]];
    plots.push(new uPlot({
      ...sizeOf(boxes[1]),
      legend: { show: false },
      cursor: { sync: { key: sync.key }, points: { size: 6 } },
      scales: { x: { time: false, distr: 3, auto: false, range: () => xRange }, y: { auto: false, range: () => [-0.15, 2.3] } },
      axes: [{ ...AXIS, ...logAxis(), label: `${xlabel} (log scale)`, labelSize: 22 }, { ...AXIS, label: 'Local slope', labelSize: 22, size: 70, splits: () => [0, 1, 2] }],
      series: [{}, ...slopes.map(({ c }) => ({ label: `${c.label || 'Local'} slope`, stroke: c.colour, width: 2, spanGaps: true, points: { show: true, size: 5, fill: c.colour } }))],
      hooks: { drawClear: [(u) => {
        const { ctx } = u;
        const dpr = devicePixelRatio;
        ctx.save();
        ctx.setLineDash([4 * dpr, 3 * dpr]);
        ctx.strokeStyle = 'rgba(124,135,150,0.7)';
        ctx.fillStyle = '#7c8796';
        ctx.font = `${10 * dpr}px Segoe UI, system-ui, sans-serif`;
        ctx.textAlign = 'right';
        for (const [level, text] of [[2, '2: needs two absorptions'], [1, '1'], [0, '0: saturated']]) {
          const py = u.valToPos(level, 'y', true);
          ctx.beginPath();
          ctx.moveTo(u.bbox.left, py);
          ctx.lineTo(u.bbox.left + u.bbox.width, py);
          ctx.stroke();
          ctx.fillText(text, u.bbox.left + u.bbox.width - 4 * dpr, py - 3 * dpr);
        }
        ctx.restore();
      }] },
    }, d2, boxes[1]));
  }

  if (title) {
    const context = CONTEXT.map((k) => [k, shared(found, k)]).filter(([, v]) => v != null).map(([k, v]) => `${k} ${fmt(v, 4)}`).join(', ');
    title.innerHTML = `<b>${KINDS[kind]}</b> &nbsp; ${escapeHtml(context)}`;
  }
  return {
    destroy: () => plots.forEach((u) => u.destroy()),
    resize: () => plots.forEach((u, i) => boxes[i].clientWidth && u.setSize(sizeOf(boxes[i]))),
  };
}
