// Simulation mode: a clock, a timer queue and a small script API that switches excited states on and off.
//
// A script is plain JavaScript run once when the simulation starts. It sets up background loops (every, after,
// onExcite, ...) that the clock then runs. Time is in "sim seconds": the clock advances by real time × speed, and
// is frozen while paused. Random numbers come from a seeded generator, so seed(n) makes a run repeatable.
//
// Rates are per sim second. chance(rate, dt) uses p = 1 - exp(-rate·dt), so results don't depend on the step size.
//
// Script API (all names are plain globals inside a script):
//   now()                      sim time, s
//   sites(kind?)               labels of the sites, in molecule order ('sensitizer', 'mediator', 'emitter', 'molecule')
//   kindOf(label)              a site's kind
//   neighbors(label)           the sites before and after it in molecule order
//   isOn(label)                is it excited
//   on(label) / off(label)     excite / de-excite; returns true if that changed anything
//   toggle(label)
//   transfer(from, to)         off(from), then on(to)
//   set(label, { color, gain, sigma })   change how its glow looks
//   every(interval, fn(dt))    run fn every interval s; returns a handle with .stop()
//   after(delay, fn)           run fn once, delay s from now; returns a handle with .stop()
//   at(time, fn)               run fn once at sim time `time`
//   onExcite(who, fn(label))   run fn whenever a site goes ground → excited. who: a label, a kind, an array, or '*'
//   onDecay(who, fn(label))    ... excited → ground
//   decay(label, tau)          switch it off after a random lifetime (exponential, mean tau), unless it was already
//                              switched off (or off and on again) before then
//   chance(rate, dt)           true with probability 1 - exp(-rate·dt)
//   random()                   uniform in [0, 1)
//   exponential(tau)           random lifetime with mean tau
//   pick(array)                random element
//   seed(n)                    restart the random generator from n
//   count(name, n = 1)         add to a counter shown under the controls
//   countOf(name)              a counter's value (0 if never counted)
//   show(name, value)          show a value under the controls (like a counter, but set rather than added to)
//   log(...values)             show a line under the controls
//   pause()                    pause the simulation

export const MAX_EVENTS_PER_FRAME = 200000;
const MAX_DEPTH = 64; // on/off handlers calling on/off calling handlers...

/** Small, fast, seedable generator (mulberry32). */
function rng(seed) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

export class Simulation {
  /**
   * host: {
   *   sites(): [{ label, kind }]          current sites, in molecule order
   *   isOn(label), setOn(label, on)        read / write the excited state (setOn must call back sim.notify)
   *   setProps(label, patch)               color / gain / sigma
   * }
   */
  constructor(host) {
    this.host = host;
    this.status = 'reset'; // reset | running | paused | error
    this.speed = 1;
    this.onChange = null; // called after status, counters or log change
    this._clear();
  }

  _clear() {
    this.t = 0;
    this.queue = []; // sorted by time, then seq
    this.seq = 0;
    this.handlers = []; // { type: 'excite' | 'decay', match(label), fn }
    this.gen = new Map(); // label -> excitation generation (for decay())
    this.counters = new Map();
    this.logLines = [];
    this.error = null;
    this.depth = 0;
    this.random = rng((Math.random() * 2 ** 32) >>> 0);
  }

  get started() {
    return this.status !== 'reset';
  }

  /** Compile and run a script from t = 0. Throws on a syntax error or an error in the script's top level. */
  start(code) {
    this._clear();
    const api = this._api();
    let fn;
    try {
      fn = new Function(...Object.keys(api), `"use strict";\n${code}\n`);
    } catch (e) {
      this._fail(e);
      return false;
    }
    this.status = 'running';
    if (!this._guard(() => fn(...Object.values(api)))) return false;
    this._changed();
    return true;
  }

  pause() {
    if (this.status === 'running') this.status = 'paused';
    this._changed();
  }

  resume() {
    if (this.status === 'paused') this.status = 'running';
    this._changed();
  }

  /** Stop and forget everything the script set up (the host restores the excited states itself). */
  reset() {
    this._clear();
    this.status = 'reset';
    this._changed();
  }

  /** Advance the clock by dt sim seconds, running every timer that falls due on the way. */
  advance(dt) {
    if (this.status !== 'running') return;
    const end = this.t + dt;
    let n = 0;
    while (this.queue.length && this.queue[0].time <= end) {
      const ev = this.queue.shift();
      if (ev.stopped) continue;
      if (++n > MAX_EVENTS_PER_FRAME) {
        this._fail(new Error(`more than ${MAX_EVENTS_PER_FRAME} timer calls in one frame: is an interval 0?`));
        return;
      }
      this.t = ev.time;
      if (ev.interval) {
        ev.time += ev.interval;
        this._insert(ev);
      }
      if (!this._guard(() => ev.fn(ev.interval || 0))) return;
      if (this.status !== 'running') return; // the script paused itself: stop right here
    }
    this.t = end;
  }

  /** While paused: advance by dt, then stay paused. */
  step(dt) {
    if (this.status !== 'paused') return;
    this.status = 'running';
    this.advance(dt);
    if (this.status === 'running') this.status = 'paused';
    this._changed();
  }

  /** The host calls this on every ground ⇄ excited change, from the script or from a click. */
  notify(label, on) {
    if (on) this.gen.set(label, (this.gen.get(label) ?? 0) + 1);
    if (!this.started || this.status === 'error') return;
    const type = on ? 'excite' : 'decay';
    for (const h of this.handlers) {
      if (h.type === type && h.match(label)) if (!this._guard(() => h.fn(label))) return;
    }
  }

  // ---------------------------------------------------------------------------

  _changed() {
    this.onChange?.();
  }

  _fail(e) {
    this.status = 'error';
    this.error = e;
    this._changed();
  }

  /** Run script code; on an error, stop the simulation and report it. */
  _guard(fn) {
    if (++this.depth > MAX_DEPTH) {
      this.depth = 0;
      this._fail(new Error('on/off handlers kept triggering each other (more than 64 deep)'));
      return false;
    }
    try {
      fn();
      return true;
    } catch (e) {
      this._fail(e);
      return false;
    } finally {
      this.depth = Math.max(0, this.depth - 1);
    }
  }

  _insert(ev) {
    ev.seq = this.seq++;
    const q = this.queue;
    let lo = 0;
    let hi = q.length;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (q[mid].time < ev.time || (q[mid].time === ev.time && q[mid].seq < ev.seq)) lo = mid + 1;
      else hi = mid;
    }
    q.splice(lo, 0, ev);
  }

  _schedule(time, fn, interval = 0) {
    if (typeof fn !== 'function') throw new TypeError('expected a function');
    const ev = { time, fn, interval, stopped: false };
    this._insert(ev);
    return { stop: () => (ev.stopped = true) };
  }

  _label(label) {
    if (!this.host.sites().some((s) => s.label === label)) {
      const known = this.host.sites().map((s) => `'${s.label}'`).join(', ');
      throw new Error(`no site '${label}' on this molecule (sites: ${known})`);
    }
    return label;
  }

  _matcher(who) {
    if (who === '*') return () => true;
    const list = Array.isArray(who) ? who : [who];
    return (label) => {
      const kind = this.host.sites().find((s) => s.label === label)?.kind;
      return list.includes(label) || list.includes(kind);
    };
  }

  _api() {
    const sim = this;
    const host = this.host;
    const isOn = (label) => host.isOn(sim._label(label));
    const on = (label) => !isOn(label) && (host.setOn(label, true), true);
    const off = (label) => isOn(label) && (host.setOn(label, false), true);
    const exponential = (tau) => -tau * Math.log(1 - sim.random());
    const positive = (x, what) => {
      if (!(x > 0)) throw new RangeError(`${what} must be > 0 (got ${x})`);
      return x;
    };
    return {
      now: () => sim.t,
      sites: (kind) => host.sites().filter((s) => kind === undefined || s.kind === kind).map((s) => s.label),
      kindOf: (label) => host.sites().find((s) => s.label === sim._label(label)).kind,
      neighbors: (label) => {
        const all = host.sites();
        const i = all.findIndex((s) => s.label === sim._label(label));
        return [all[i - 1], all[i + 1]].filter(Boolean).map((s) => s.label);
      },
      isOn,
      on,
      off,
      toggle: (label) => (isOn(label) ? off(label) : on(label), isOn(label)),
      transfer: (from, to) => {
        off(from);
        return on(to);
      },
      set: (label, patch) => host.setProps(sim._label(label), patch),
      every: (interval, fn) => sim._schedule(sim.t + positive(interval, 'every() interval'), fn, interval),
      after: (delay, fn) => sim._schedule(sim.t + Math.max(0, delay), fn),
      at: (time, fn) => sim._schedule(Math.max(sim.t, time), fn),
      onExcite: (who, fn) => sim.handlers.push({ type: 'excite', match: sim._matcher(who), fn }),
      onDecay: (who, fn) => sim.handlers.push({ type: 'decay', match: sim._matcher(who), fn }),
      decay: (label, tau) => {
        if (!isOn(label)) return null;
        const g = sim.gen.get(label);
        return sim._schedule(sim.t + exponential(positive(tau, 'decay() lifetime')), () => {
          if (sim.gen.get(label) === g) off(label);
        });
      },
      chance: (rate, dt) => sim.random() < 1 - Math.exp(-rate * dt),
      random: () => sim.random(),
      exponential: (tau) => exponential(positive(tau, 'exponential() mean')),
      pick: (arr) => arr[Math.floor(sim.random() * arr.length)],
      seed: (n) => (sim.random = rng(n)),
      count: (name, n = 1) => {
        sim.counters.set(name, (sim.counters.get(name) ?? 0) + n);
        sim._changed();
      },
      countOf: (name) => sim.counters.get(name) ?? 0,
      show: (name, value) => {
        sim.counters.set(name, value);
        sim._changed();
      },
      log: (...values) => {
        const text = values.map((v) => (typeof v === 'string' ? v : JSON.stringify(v))).join(' ');
        sim.logLines.push(`${sim.t.toFixed(2)} s  ${text}`);
        if (sim.logLines.length > 50) sim.logLines.shift();
        sim._changed();
      },
      pause: () => sim.pause(),
    };
  }
}
