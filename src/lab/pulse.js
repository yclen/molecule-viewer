// The live pulse experiment, ported from sim/live.py: fresh molecules are flashed all the time, and each photon lands in
// a histogram by when it came out. The Monte Carlo runs in workers, one per spare CPU core; the exact overlay is solved
// in another worker. Everything is in ns and rates per ns, except the wall-clock rates (per real second).

import uPlot from 'uplot';
import 'uplot/dist/uPlot.min.css';
import { binnedExact, fluorescenceExact, QUIET } from '../sim/exact.js';
import ExactWorker from './exact-worker.js?worker&inline';
import McWorker from './mc-worker.js?worker&inline';
import {
  AXIS, defaultValues, exactSeconds, FIELDS, FLUOR_FIELDS, fmt, howLong, int, kexcPerMw, parameters, percent, powerText,
  rulesFor, TTA_FIELDS, Y_AXIS_SIZE,
} from './model.js';
import { saveLocation, saveRun, stamp, toCsv } from './store.js';

const KEXC_LOW = 1e-4;
const KEXC_HIGH = 1e2;
const SLIDER_STEPS = 600; // 100 per decade
const TARGET = 100_000; // molecules per second to aim for. 0 means as fast as the machine can go
const CELLS_PER_LIFETIME = 100; // counts are kept this fine, so bins can change without losing anything
const MOST_CELLS = 200_000;
const BIN_CELLS = [1, 2, 5, 10, 20, 50, 100];
const DEFAULT_BIN_CELLS = 10; // a tenth of the singlet lifetime
const THREADS = Math.min(8, Math.max(1, Math.floor((navigator.hardwareConcurrency || 2) / 2)));
const FRAME_MS = 50;
const FRAMES_PER_TEXT = 5;
const RATE_OVER_S = 2;
const EMPTY_BIN = 0.5; // where an empty bin sits on a log axis, which can't show 0
const PAST_WARNING = 0.01;
const EXACT_MOST_N = 9;
const EXACT_SLOW_S = 60;
const EXACT_MOST_S = 600;
const EXACT_SPEED = { alone: 0.7, counting: 1.8 }; // measured at n = 8 against exactSeconds: the laser on slows it
const COLOURS = { perylene: '#11d8f2', dpa: '#2a78d6', fluorescence: '#c792ff' };
const EXACT_COLOURS = { perylene: '#b5f4fc', dpa: '#9cc3f5', fluorescence: '#e6d0ff' };
const PILES = { tta: ['perylene', 'dpa'], fluorescence: ['fluorescence'] };
const NAMES = { perylene: 'Perylene photons', dpa: 'DPA photons', fluorescence: 'Fluorescence photons' };
const HINT = 'Hover over the histogram to read a bin.';

const rateToStep = (k) => Math.round((SLIDER_STEPS * Math.log10(k / KEXC_LOW)) / Math.log10(KEXC_HIGH / KEXC_LOW));
const stepToRate = (s) => Number((KEXC_LOW * (KEXC_HIGH / KEXC_LOW) ** (s / SLIDER_STEPS)).toPrecision(3));
const pad = (v) => String(Math.floor(v)).padStart(2, '0');
const clockText = (s) => `${pad(s / 3600)}:${pad((s / 60) % 60)}:${pad(s % 60)}`;
const fluorescent = (p) => p.sample === 'fluorescence';

/** Everything counted for one set of parameters, in fine cells. `main` is the emitter's pile, `other` the DPAs'. */
class Measurement {
  constructor(p) {
    this.p = p;
    this.mainSite = fluorescent(p) ? 0 : p.n + 1;
    const lifetime = fluorescent(p) ? p.tauFl : p.tauS;
    this.cell = Math.max(lifetime / CELLS_PER_LIFETIME, p.window / MOST_CELLS);
    this.cells = Math.ceil(p.window / this.cell);
    this.main = new Float64Array(this.cells);
    this.other = new Float64Array(this.cells);
    this.totals = { molecules: 0, excitations: 0, photons: 0, past_window: 0 };
  }

  add({ molecules, excitations, times, sites }) {
    let past = 0;
    for (let i = 0; i < times.length; i++) {
      const c = Math.floor(times[i] / this.cell);
      if (c < this.cells) (sites[i] === this.mainSite ? this.main : this.other)[c]++;
      else past++;
    }
    this.totals.molecules += molecules;
    this.totals.excitations += excitations;
    this.totals.photons += times.length;
    this.totals.past_window += past;
  }

  /** The bin edges, then the main and other photons in each bin. */
  binned(cellsPerBin) {
    const bins = Math.ceil(this.cells / cellsPerBin);
    const edges = Float64Array.from({ length: bins + 1 }, (_, k) => this.cell * cellsPerBin * k);
    const main = new Float64Array(bins);
    const other = new Float64Array(bins);
    for (let c = 0; c < this.cells; c++) {
      main[(c / cellsPerBin) | 0] += this.main[c];
      other[(c / cellsPerBin) | 0] += this.other[c];
    }
    return { edges, main, other };
  }
}

const field = (name) => {
  const [label, unit] = FIELDS[name];
  return `<label for="pulse-${name}">${label}</label><input type="text" id="pulse-${name}" data-field="${name}"
    inputmode="decimal" autocomplete="off"><span class="unit">${unit}</span>`;
};

const TEMPLATE = `
<header class="lab-head">
  <div class="lab-title">Live pulse experiment</div>
  <div class="lab-sub">time-resolved emission, Monte Carlo · times in ns, rates per ns</div>
  <div class="spacer"></div>
  <div class="clock"><small>ACQ</small><span data-id="clock">00:00:00</span></div>
  <div class="lamp" data-id="laserLamp"></div>
  <div class="lamp" data-id="exactLamp"></div>
</header>
<div class="lab-body">
  <aside class="lab-side">
    <section class="panel">
      <div class="section">Laser · 635 nm</div>
      <button class="fire" data-id="fire"></button>
      <div class="section" style="margin-top:14px">Excitation rate, kExc</div>
      <div class="power"><input type="text" data-id="kexc" inputmode="decimal" autocomplete="off" aria-label="Excitation rate, kExc, per ns"><span>/ns</span></div>
      <input type="range" class="laser" data-id="slider" min="0" max="${SLIDER_STEPS}" step="1" aria-label="Laser power">
      <div class="ends"><span>0.0001 /ns</span><span>100 /ns</span></div>
      <div class="row steps" style="margin-top:8px">
        <button data-factor="0.1">÷10</button><button data-factor="0.5">÷2</button>
        <button data-factor="2">×2</button><button data-factor="10">×10</button>
      </div>
      <div class="caption mono" data-id="power" style="color:var(--text)"></div>
      <div class="note">w0 1 µm · ε 50,000 /M/cm · 635 nm · ${fmt(kexcPerMw(), 3)} /ns per mW</div>
      <div class="grid" style="margin-top:12px">${field('duration')}</div>
    </section>
    <section class="panel">
      <div class="section">Acquisition</div>
      <div class="grid">
        ${field('window')}
        <label for="pulse-bin">Bin width</label><select id="pulse-bin" data-id="bin"></select><span class="unit"></span>
        <label for="pulse-target">Throughput</label><input type="text" id="pulse-target" data-id="target" inputmode="numeric"><span class="unit">mol/s</span>
      </div>
      <label class="check"><input type="checkbox" data-id="log"> Log counts</label>
      <button class="exact-btn" data-id="exact"></button>
      <div class="bar" data-id="exactBar" hidden><div></div></div>
      <div class="caption" data-id="exactNote"></div>
      <div class="row" style="margin-top:12px"><button data-id="clear">CLEAR</button><button data-id="save">SAVE</button></div>
      <div data-id="where"></div>
    </section>
    <section class="panel">
      <div class="section">Sample</div>
      <div class="segmented"><button data-sample="tta">TTA</button><button data-sample="fluorescence">FLUORESCENCE</button></div>
      <div class="grid wide" data-id="tta" style="margin-top:12px">${TTA_FIELDS.map(field).join('')}</div>
      <div class="grid wide" data-id="fluorescence" style="margin-top:12px" hidden>${FLUOR_FIELDS.map(field).join('')}</div>
      <div class="note" style="margin-top:10px">Changing a value clears the histogram. Bin width, the log axis and
        throughput don't. TTA: sensitizer, n DPA mediators, perylene. Fluorescence: a plain dye, for comparison.</div>
    </section>
  </aside>
  <div class="lab-main">
    <div class="tiles">
      <div class="tile"><div class="name">Quantum yield</div><div class="value hero" data-tile="yield">–</div><div class="sub" data-sub="yield">&nbsp;</div></div>
      <div class="tile"><div class="name">Count rate</div><div class="value" data-tile="rate">0</div><div class="sub" data-sub="rate">photons / s</div></div>
      <div class="tile"><div class="name">Photons</div><div class="value" data-tile="photons">0</div><div class="sub" data-sub="photons">in the histogram</div></div>
      <div class="tile"><div class="name">Molecules</div><div class="value" data-tile="molecules">0</div><div class="sub" data-sub="molecules">0 /s</div></div>
      <div class="tile"><div class="name">Absorbed / molecule</div><div class="value" data-tile="absorbed">–</div><div class="sub" data-sub="absorbed">photons</div></div>
    </div>
    <section class="panel plot-panel">
      <div class="plot-head">
        <div class="section">Histogram <span class="hint">red band: pulse on · drag to zoom · double-click to reset · click a legend entry to hide it</span></div>
        <div class="spacer"></div>
        <button data-id="autoY" aria-pressed="true">↕ AUTO Y</button><button data-id="home">⌂ HOME</button>
      </div>
      <div class="plot-title" data-id="title"></div>
      <div class="plot-box" data-id="plotBox"><div class="legend" data-id="legend"></div></div>
    </section>
    <div class="status" data-id="status">Press START to flash molecules. ${THREADS} worker${THREADS > 1 ? 's' : ''} on this computer.</div>
  </div>
</div>`;

export function mount(root) {
  root.innerHTML = TEMPLATE;
  const $ = (id) => root.querySelector(`[data-id="${id}"]`);
  const fields = Object.fromEntries([...root.querySelectorAll('[data-field]')].map((el) => [el.dataset.field, el]));
  $('where').append(saveLocation());

  const values = defaultValues();
  let p = parameters(values);
  let measurement = new Measurement(p);
  let generation = 0;
  let running = false;
  let visible = true;
  let cellsPerBin = DEFAULT_BIN_CELLS;
  let target = TARGET;
  let autoY = true;
  let acquired = 0;
  let lastTick = performance.now();
  let recent = []; // [time, molecules, photons], over the last RATE_OVER_S
  let frame = 0;
  let lastDrawn = -1;
  let drawn = null; // what's on screen: { edges, main, other, expected }
  const hidden = new Set(); // piles hidden from the legend
  const exact = { on: false, key: null, result: null, doing: null, share: 0, began: 0, worker: null };
  const seedBase = crypto.getRandomValues(new Uint32Array(1))[0] % 1_000_000;

  // ---- the Monte Carlo workers ----
  const workers = Array.from({ length: THREADS }, (_, i) => {
    const w = new McWorker();
    w.onmessage = ({ data }) => {
      if (data.generation === generation) measurement.add(data);
    };
    w.index = i;
    return w;
  });
  // While an exact solve runs, only two workers count: on a laptop with fast and slow cores, a full set of workers
  // makes the solve several times slower (n = 8: 25 s alone, over 2 min beside 8 workers). Two still keep up 100,000/s.
  const active = () => (exact.doing ? Math.min(2, THREADS) : THREADS);
  const configure = () => {
    const { R, sites } = rulesFor(p);
    const n = active();
    for (const w of workers) {
      w.postMessage({ type: 'config', generation, R, sites, duration: p.duration, target: target / n,
        seed: seedBase + 7919 * w.index + 104729 * generation });
      w.postMessage({ type: 'on', on: running && w.index < n });
    }
  };

  // ---- the controls ----
  const showValues = () => {
    for (const [name, el] of Object.entries(fields)) el.value = fmt(values[name], 6);
    $('kexc').value = fmt(values.kExc, 3);
    $('slider').value = rateToStep(values.kExc);
    $('power').textContent = `laser power ≈ ${powerText(values.kExc / kexcPerMw())}`;
    $('target').value = String(target);
    for (const b of root.querySelectorAll('[data-sample]')) b.setAttribute('aria-pressed', b.dataset.sample === values.sample);
    $('tta').hidden = values.sample !== 'tta';
    $('fluorescence').hidden = values.sample !== 'fluorescence';
  };

  const showBinChoices = () => {
    const select = $('bin');
    select.innerHTML = BIN_CELLS.map((c) => `<option value="${c}">${fmt(c * measurement.cell, 3)} ns</option>`).join('');
    select.value = cellsPerBin;
  };

  /** New parameters: a fresh histogram. */
  function restart() {
    p = parameters(values);
    generation++;
    measurement = new Measurement(p);
    recent = [];
    acquired = 0;
    lastDrawn = -1;
    configure();
    showBinChoices();
    showValues();
    updateExact();
    makePlot();
    draw();
    showNumbers();
  }

  for (const [name, el] of Object.entries(fields)) {
    el.addEventListener('change', () => {
      const [, , lo, hi] = FIELDS[name];
      let v = Number(el.value.trim());
      if (name === 'n') v = Math.round(v);
      if (!Number.isFinite(v) || v < lo || v > hi) {
        el.classList.add('bad');
        el.title = `Between ${lo} and ${hi}`;
        setTimeout(() => {
          el.classList.remove('bad');
          el.value = fmt(values[name], 6);
        }, 900);
        return;
      }
      el.title = '';
      if (v !== values[name]) {
        values[name] = v;
        restart();
      }
    });
    el.addEventListener('keydown', (e) => e.key === 'Enter' && el.blur());
  }

  const setRate = (k) => {
    k = Math.min(Math.max(k, KEXC_LOW), KEXC_HIGH);
    if (k === values.kExc) return showValues();
    values.kExc = k;
    restart();
  };
  $('slider').addEventListener('input', () => {
    const k = stepToRate(Number($('slider').value));
    $('kexc').value = fmt(k, 3);
    $('power').textContent = `laser power ≈ ${powerText(k / kexcPerMw())}`;
  });
  $('slider').addEventListener('change', () => setRate(stepToRate(Number($('slider').value))));
  $('kexc').addEventListener('change', () => {
    const k = Number($('kexc').value.trim());
    if (Number.isFinite(k) && k > 0) setRate(k);
    else showValues();
  });
  $('kexc').addEventListener('keydown', (e) => e.key === 'Enter' && $('kexc').blur());
  for (const b of root.querySelectorAll('[data-factor]'))
    b.addEventListener('click', () => setRate(Number((values.kExc * Number(b.dataset.factor)).toPrecision(3))));

  for (const b of root.querySelectorAll('[data-sample]'))
    b.addEventListener('click', () => {
      if (values.sample === b.dataset.sample) return;
      values.sample = b.dataset.sample;
      restart();
    });

  $('target').addEventListener('change', () => {
    const v = Math.round(Number($('target').value));
    if (Number.isFinite(v) && v >= 0) {
      target = v;
      configure(); // the same histogram: only the pace changes
    }
    showValues();
  });
  $('bin').addEventListener('change', () => {
    cellsPerBin = Number($('bin').value);
    lastDrawn = -1;
    draw();
  });
  $('log').addEventListener('change', () => {
    makePlot();
    draw();
  });
  $('autoY').addEventListener('click', () => {
    autoY = !autoY;
    $('autoY').setAttribute('aria-pressed', autoY);
    if (autoY) fitY();
  });
  $('home').addEventListener('click', home);
  $('clear').addEventListener('click', restart);
  $('save').addEventListener('click', save);

  const fire = $('fire');
  const showFire = () => {
    fire.innerHTML = running ? '■&nbsp;&nbsp; STOP' : '▶&nbsp;&nbsp; START';
    fire.classList.toggle('on', running);
  };
  fire.addEventListener('click', () => {
    running = !running;
    for (const w of workers) w.postMessage({ type: 'on', on: running && w.index < active() });
    if (running) recent = [];
    showFire();
    updateExact(); // its time estimate depends on the laser
    updateLamps();
  });

  // ---- the exact overlay ----
  function exactVerdict() {
    if (fluorescent(p)) return [true, 'Instant: a two-state molecule has a closed-form answer.'];
    if (p.n > EXACT_MOST_N) return [false, `n = ${p.n} is too large to solve exactly: it has ${int(3 ** (p.n + 2))} strings.`];
    const seconds = exactSeconds(p, p.window) * EXACT_SPEED[running ? 'counting' : 'alone'];
    if (seconds > EXACT_MOST_S) return [false, `Too slow to solve exactly at these settings: about ${howLong(seconds)}.`];
    return [true, `Takes about ${howLong(seconds)}${seconds > EXACT_SLOW_S ? ', which is slow' : ''}.`];
  }

  const stopSolve = () => {
    exact.worker?.terminate();
    exact.worker = null;
    if (exact.doing) {
      exact.doing = null;
      configure(); // all the workers again
    }
  };

  function updateExact() {
    const [allowed, verdict] = exactVerdict();
    const btn = $('exact');
    btn.disabled = !allowed;
    const key = JSON.stringify(p);
    if (!(allowed && exact.on)) {
      if (exact.doing) stopSolve();
      if (exact.key !== key) {
        exact.key = null;
        exact.result = null;
      }
      $('exactNote').textContent = verdict;
      $('exactBar').hidden = true;
    } else if (key !== exact.key) {
      stopSolve();
      exact.key = key;
      exact.result = null;
      if (fluorescent(p)) {
        exact.result = fluorescenceExact(p);
        $('exactNote').textContent = verdict;
      } else {
        exact.doing = 'Preparing';
        exact.share = 0;
        exact.began = performance.now();
        $('exactBar').hidden = false;
        $('exactNote').textContent = 'Preparing.';
        const worker = new ExactWorker();
        exact.worker = worker;
        worker.onmessage = ({ data }) => readExact(data, key);
        worker.postMessage({ type: 'pulseExact', key, p, R: rulesFor(p).R });
        configure(); // fewer Monte Carlo workers while it solves
      }
    }
    showExactButton();
    updateLamps();
  }

  function readExact(data, key) {
    if (key !== exact.key) return; // from a solve that was dropped
    const spent = (performance.now() - exact.began) / 1000;
    if (data.type === 'progress') {
      exact.share = data.share;
      exact.doing = data.doing;
      const left = data.share > 0.05 ? (spent / data.share) * (1 - data.share) : null;
      $('exactBar').firstElementChild.style.width = `${(100 * data.share).toFixed(1)}%`;
      $('exactNote').textContent = `${data.doing}${left ? `, about ${howLong(left)} left.` : '.'}`;
    } else {
      exact.doing = null;
      exact.worker?.terminate();
      exact.worker = null;
      configure(); // all the workers again
      $('exactBar').hidden = true;
      if (data.type === 'done') {
        exact.result = data.result;
        $('exactNote').textContent = `Exact solve done in ${howLong(spent)}.`;
        lastDrawn = -1;
        draw();
        showNumbers();
      } else $('exactNote').textContent = `The exact solve failed: ${data.error}`;
    }
    showExactButton();
    updateLamps();
  }

  const exactShown = () => exact.on && exact.result && exact.key === JSON.stringify(p);

  function showExactButton() {
    const btn = $('exact');
    let state = 'off';
    let text = '◇&nbsp;&nbsp; EXACT OVERLAY · OFF';
    if (btn.disabled) text = '◇&nbsp;&nbsp; EXACT OVERLAY · N/A';
    else if (exact.on && exact.doing) {
      state = 'solving';
      text = `◆&nbsp;&nbsp; EXACT OVERLAY · SOLVING ${Math.round(100 * exact.share)}%`;
    } else if (exact.on) {
      state = 'on';
      text = '◆&nbsp;&nbsp; EXACT OVERLAY · ON';
    }
    btn.innerHTML = text;
    btn.className = `exact-btn ${state === 'off' ? '' : state}`;
  }
  $('exact').addEventListener('click', () => {
    exact.on = !exact.on;
    updateExact();
    lastDrawn = -1;
    draw();
  });

  function updateLamps() {
    const lamp = (el, colour, text) => {
      el.innerHTML = `<i style="color:${colour || 'var(--dim)'}">●</i>${text}`;
      el.classList.toggle('lit', !!colour);
    };
    lamp($('laserLamp'), running && '#ff2100', running ? `LASER ON · kExc ${fmt(values.kExc, 3)} /ns` : 'LASER OFF');
    if (exact.on && exact.doing) lamp($('exactLamp'), '#f2a93b', `EXACT · SOLVING ${Math.round(100 * exact.share)}%`);
    else if (exactShown()) lamp($('exactLamp'), '#11d8f2', 'EXACT · READY');
    else lamp($('exactLamp'), null, 'EXACT');
  }

  // ---- the plot ----
  let plot = null;
  const box = $('plotBox');
  const legend = $('legend');

  function makePlot() {
    plot?.destroy();
    const log = $('log').checked;
    const piles = PILES[p.sample];
    plot = new uPlot({
      width: Math.max(box.clientWidth, 200),
      height: Math.max(box.clientHeight, 150),
      legend: { show: false },
      cursor: { drag: { x: true, y: false, setScale: false }, points: { show: false }, y: false,
        bind: { dblclick: () => null } },
      scales: { x: { time: false, auto: false, min: 0, max: p.window }, y: { distr: log ? 3 : 1, auto: false, min: log ? EMPTY_BIN : 0, max: 10 } },
      axes: [
        { ...AXIS, label: 'Time after the pulse starts (ns)' },
        { ...AXIS, label: 'Photons per bin', size: Y_AXIS_SIZE },
      ],
      series: [{}, ...piles.map((pile, i) => ({
        stroke: COLOURS[pile], width: 1.5, fill: `${COLOURS[pile]}${i === 0 ? '30' : '40'}`,
        paths: uPlot.paths.stepped({ align: 1 }), points: { show: false }, fillTo: () => (log ? EMPTY_BIN : 0),
        show: !hidden.has(pile),
      }))],
      hooks: {
        drawClear: [(u) => {
          const { ctx } = u;
          const x0 = u.valToPos(0, 'x', true);
          const x1 = u.valToPos(Math.min(p.duration, u.scales.x.max), 'x', true);
          const { left, top, width, height } = u.bbox;
          const a = Math.max(x0, left);
          const b = Math.min(x1, left + width);
          if (b > a) {
            ctx.fillStyle = 'rgba(255, 33, 0, 0.17)';
            ctx.fillRect(a, top, b - a, height);
          }
        }],
        draw: [drawExactCurves],
        setCursor: [(u) => {
          const left = u.cursor.left;
          if (left == null || left < 0) {
            $('status').textContent = HINT;
            return;
          }
          showBin(u.posToVal(left, 'x'));
        }],
        setSelect: [(u) => {
          if (u.select.width > 4) {
            u.setScale('x', { min: u.posToVal(u.select.left, 'x'), max: u.posToVal(u.select.left + u.select.width, 'x') });
            if (autoY) fitY();
          }
          u.setSelect({ width: 0, height: 0 }, false);
        }],
      },
    }, [[0, 1], ...piles.map(() => [0, 0])], box);
    plot.over.addEventListener('dblclick', home);
    plot.root.querySelector('.u-over').title = '';
    legend.innerHTML = '';
    piles.forEach((pile, i) => {
      const item = document.createElement('div');
      item.innerHTML = `<i style="background:${COLOURS[pile]}"></i>${NAMES[pile]}`;
      item.classList.toggle('off', hidden.has(pile));
      item.onclick = () => {
        if (hidden.has(pile)) hidden.delete(pile);
        else hidden.add(pile);
        item.classList.toggle('off', hidden.has(pile));
        plot.setSeries(i + 1, { show: !hidden.has(pile) });
        if (autoY) fitY();
      };
      legend.append(item);
    });
    const band = document.createElement('div');
    band.style.cursor = 'default';
    band.innerHTML = '<i style="background:rgba(255,33,0,0.55); height:10px"></i>Pulse on (laser)';
    legend.append(band);
    const ex = document.createElement('div');
    ex.className = 'exact-legend';
    ex.style.cursor = 'default';
    ex.innerHTML = `<i class="dash" style="border-color:${EXACT_COLOURS[piles[0]]}"></i>Exact`;
    legend.append(ex);
    legend.lastChild.hidden = true;
  }

  function drawExactCurves(u) {
    if (!drawn?.expected) return;
    const { ctx } = u;
    const { edges, expected } = drawn;
    const piles = PILES[p.sample];
    const log = $('log').checked;
    ctx.save();
    const { left, top, width, height } = u.bbox;
    ctx.beginPath();
    ctx.rect(left, top, width, height);
    ctx.clip();
    ctx.setLineDash([6 * devicePixelRatio, 4 * devicePixelRatio]);
    ctx.lineWidth = 2 * devicePixelRatio;
    piles.forEach((pile, i) => {
      if (hidden.has(pile)) return;
      ctx.strokeStyle = EXACT_COLOURS[pile];
      ctx.beginPath();
      let pen = false;
      for (let k = 0; k < expected[i].length; k++) {
        const y = expected[i][k];
        if (log && !(y > 0)) {
          pen = false;
          continue;
        }
        const px = u.valToPos((edges[k] + edges[k + 1]) / 2, 'x', true);
        const py = u.valToPos(y, 'y', true);
        if (pen) ctx.lineTo(px, py);
        else ctx.moveTo(px, py);
        pen = true;
      }
      ctx.stroke();
    });
    ctx.restore();
  }

  function draw() {
    if (!plot) return;
    const { edges, main, other } = measurement.binned(cellsPerBin);
    const log = $('log').checked;
    const floor = log ? EMPTY_BIN : 0;
    const steps = (counts) => {
      const y = new Float64Array(counts.length + 1);
      for (let k = 0; k < counts.length; k++) y[k] = Math.max(counts[k], floor);
      y[counts.length] = y[counts.length - 1] ?? floor;
      return y;
    };
    let expected = null;
    if (exactShown()) {
      expected = binnedExact(exact.result, edges).map((share) => share.map((v) => v * measurement.totals.molecules));
    }
    drawn = { edges, main, other, expected };
    legend.lastChild.hidden = !expected;
    const piles = PILES[p.sample];
    plot.setData([Array.from(edges), ...[main, other].slice(0, piles.length).map((c) => Array.from(steps(c)))], false);
    if (autoY) fitY();
    else plot.redraw(false);
  }

  function fitY() {
    if (!drawn) return;
    const { edges, main, other, expected } = drawn;
    const { min, max } = plot.scales.x;
    const piles = PILES[p.sample];
    let top = 0;
    const lists = [main, other].slice(0, piles.length);
    lists.forEach((counts, i) => {
      if (hidden.has(piles[i])) return;
      for (let k = 0; k < counts.length; k++) {
        if (edges[k + 1] < min || edges[k] > max) continue;
        top = Math.max(top, counts[k], expected ? expected[i][k] : 0);
      }
    });
    const log = $('log').checked;
    if (log) plot.setScale('y', { min: EMPTY_BIN, max: Math.max(top * 2, 10) });
    else plot.setScale('y', { min: 0, max: Math.max(top * 1.08, 5) });
  }

  function home() {
    plot.setScale('x', { min: 0, max: p.window });
    autoY = true;
    $('autoY').setAttribute('aria-pressed', true);
    fitY();
  }

  function showBin(t) {
    if (!drawn) return;
    const { edges, main, other, expected } = drawn;
    let k = 0;
    while (k < edges.length - 1 && edges[k + 1] <= t) k++;
    if (t < 0 || k >= main.length) return;
    const piles = PILES[p.sample];
    const parts = [`${fmt(edges[k], 4)} to ${fmt(edges[k + 1], 4)} ns`];
    piles.forEach((pile, i) => {
      const counts = i === 0 ? main : other;
      parts.push(`${pile} ${int(counts[k])}${expected ? ` (exact ${expected[i][k].toFixed(1)})` : ''}`);
    });
    $('status').textContent = parts.join('   ·   ');
  }

  new ResizeObserver(() => {
    if (plot && box.clientWidth > 0) plot.setSize({ width: box.clientWidth, height: Math.max(box.clientHeight, 150) });
  }).observe(box);

  // ---- every frame ----
  function showNumbers() {
    const { molecules, photons, excitations, past_window: pastWindow } = measurement.totals;
    const spent = recent.length > 1 ? recent[recent.length - 1][0] - recent[0][0] : 0;
    const perS = [1, 2].map((k) => (spent ? (recent[recent.length - 1][k] - recent[0][k]) / spent : 0));
    const past = photons ? pastWindow / photons : 0;
    const held = running || !perS[0] ? '' : ' · held';
    const tiles = { yield: '–', rate: int(perS[1]), photons: int(photons), molecules: int(molecules),
      absorbed: molecules ? fmt(excitations / molecules, 3) : '–' };
    const subs = { yield: ' ', rate: `photons / s${held}`, molecules: `${int(perS[0])} /s${held}`, absorbed: 'photons',
      photons: past >= PAST_WARNING ? `${Math.round(100 * past)}% past the window` : 'in the histogram' };
    if (excitations) {
      tiles.yield = percent(photons / excitations);
      subs.yield = `± ${percent(Math.sqrt(photons) / excitations)}`;
    }
    if (exactShown() && Number.isFinite(exact.result.excitations)) {
      const unfinished = exact.result.stillActive >= QUIET * exact.result.excitations;
      subs.yield += `  ·  exact ${unfinished ? '≥ ' : ''}${percent(exact.result.photons / exact.result.excitations)}`;
    }
    for (const [name, text] of Object.entries(tiles)) root.querySelector(`[data-tile="${name}"]`).textContent = text;
    for (const [name, text] of Object.entries(subs)) root.querySelector(`[data-sub="${name}"]`).textContent = text;
    root.querySelector('[data-sub="photons"]').classList.toggle('warn', past >= PAST_WARNING);
    $('clock').textContent = clockText(acquired);
    updateLamps();
    showExactButton();
    const width = measurement.cell * cellsPerBin;
    const sample = fluorescent(p) ? `fluorophore, τF = ${fmt(p.tauFl, 4)} ns, Φ_F = ${fmt(p.fluorYield, 4)}`
      : `TTA, kHop = ${fmt(p.kHop, 3)} /ns, n = ${p.n}`;
    $('title').innerHTML = `<b>Quantum yield ${tiles.yield}</b> &nbsp; ${sample}, kExc = ${fmt(p.kExc, 3)} /ns, pulse = `
      + `${fmt(p.duration, 4)} ns, ${int(molecules)} molecules, ${fmt(width, 3)} ns bins`;
  }

  setInterval(() => {
    const now = performance.now() / 1000;
    const { molecules, photons } = measurement.totals;
    if (running) {
      acquired += now - lastTick;
      recent.push([now, molecules, photons]);
      while (now - recent[0][0] > RATE_OVER_S) recent.shift();
      if (!molecules) $('status').textContent = 'Starting…';
      else if ($('status').textContent === 'Starting…') $('status').textContent = HINT;
    }
    lastTick = now;
    if (!visible) return;
    frame++;
    if (lastDrawn !== molecules) {
      lastDrawn = molecules;
      draw();
    }
    if (frame % FRAMES_PER_TEXT === 0) showNumbers();
  }, FRAME_MS);

  // ---- saving ----
  async function save() {
    const m = measurement;
    const { edges, main, other } = m.binned(cellsPerBin);
    const piles = PILES[m.p.sample];
    const counted = [main, other].slice(0, piles.length);
    const totals = { ...m.totals };
    if (!totals.molecules) {
      $('status').textContent = 'Nothing to save yet.';
      return;
    }
    const exactBins = exactShown() ? binnedExact(exact.result, edges).slice(0, piles.length) : null;
    const rows = [];
    for (let k = 0; k < main.length; k++) {
      const total = main[k] + other[k];
      const scale = 1 / (totals.molecules * (edges[k + 1] - edges[k]));
      const row = { point: 0, bin_start_ns: edges[k], bin_end_ns: edges[k + 1], pulse_on: edges[k] < m.p.duration ? 1 : 0,
        photons: total };
      piles.forEach((pile, i) => (row[`${pile}_photons`] = counted[i][k]));
      row.photons_per_molecule_per_ns = total * scale;
      row.error = Math.sqrt(total) * scale;
      if (exactBins) piles.forEach((pile, i) => (row[`exact_${pile}_photons`] = exactBins[i][k] * totals.molecules));
      rows.push(row);
    }
    const summary = { ...m.p, ...totals, photons_per_excitation: totals.photons / Math.max(totals.excitations, 1) };
    if (exactBins) summary.exact_photons_per_excitation = exact.result.photons / exact.result.excitations;
    const { day, time } = stamp();
    const files = [
      { name: 'params.json', data: JSON.stringify({ mode: 'live', date: `${day} ${time}`, units: 'times in ns, rates per ns',
        parameters: m.p, bin_ns: edges[1], source: 'Molecule Viewer website (src/sim/, a port of sim/)' }, null, 2) },
      { name: 'histogram.csv', data: toCsv(rows) },
      { name: 'summary.csv', data: toCsv([summary]) },
    ];
    const png = await plotPng();
    if (png) files.push({ name: 'histogram.png', data: png });
    try {
      const where = await saveRun(`${day}_${fluorescent(m.p) ? 'live-fluorescence' : 'live'}`, files);
      $('status').textContent = `Saved to ${where}.`;
    } catch {
      $('status').textContent = 'Not saved.';
    }
  }

  async function plotPng() {
    const src = plot.ctx.canvas;
    const out = document.createElement('canvas');
    out.width = src.width;
    out.height = src.height;
    const ctx = out.getContext('2d');
    ctx.fillStyle = '#151b23';
    ctx.fillRect(0, 0, out.width, out.height);
    ctx.drawImage(src, 0, 0);
    const blob = await new Promise((resolve) => out.toBlob(resolve, 'image/png'));
    return blob ? new Uint8Array(await blob.arrayBuffer()) : null;
  }

  showFire();
  restart();
  $('status').textContent = `Press START to flash molecules. ${THREADS} worker${THREADS > 1 ? 's' : ''} on this computer.`;

  return {
    shown(on) {
      visible = on;
      if (on && plot) {
        plot.setSize({ width: box.clientWidth, height: Math.max(box.clientHeight, 150) });
        lastDrawn = -1;
      }
    },
  };
}
