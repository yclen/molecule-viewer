// Simulation mode: a clock, a timer queue and a small script API that switches excited states on and off.
//
// A script is plain JavaScript run once when the simulation starts. It sets up background loops (every, after,
// onExcite, ...) that the clock then runs. Time is in "sim seconds": the clock advances by real time × speed, and
// is frozen while paused. Random numbers come from a seeded generator, so seed(n) makes a run repeatable.
//
// Rates are per sim second. chance(rate, dt) uses p = 1 - exp(-rate·dt), so results don't depend on the step size.
//
// Every site has two excited states: the singlet S₁ ('S') and the triplet T₁ ('T'). They are independent (a site
// can have either, both or neither on). Every call that takes a `state` defaults to 'S', so a script that never
// mentions states drives the singlets. 'singlet' / 'triplet' work too, and '*' means "either / both" where noted.
//
// Script API (all names are plain globals inside a script):
//   now()                      sim time, s
//   sites(kind?)               labels of the sites, in molecule order ('sensitizer', 'mediator', 'emitter', 'molecule')
//   kindOf(label)              a site's kind
//   neighbors(label)           the sites before and after it in molecule order
//   isOn(label, state?)        is that state on ('*': is any state on)
//   on(label, state?)          switch a state on; returns true if that changed anything
//   off(label, state?)         switch a state off ('*': both); returns true if that changed anything
//   toggle(label, state?)
//   transfer(from, to, state?, toState?)   off(from, state), then on(to, toState); toState defaults to state
//   set(label, { color, gain, sigma }, state?)   change how a state's glow looks ('*': both)
//   every(interval, fn(dt))    run fn every interval s; returns a handle with .stop()
//   after(delay, fn)           run fn once, delay s from now; returns a handle with .stop()
//   at(time, fn)               run fn once at sim time `time`
//   onExcite(who, fn(label, state), state?)   run fn whenever that state of a site switches on. who: a label, a
//                              kind, an array, or '*'; state '*' = either state
//   onDecay(who, fn(label, state), state?)    ... switches off
//   decay(label, tau, state?)  switch that state off after a random lifetime (exponential, mean tau), unless it was
//                              already switched off (or off and on again) before then
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

import { STATES, stateKey } from './excited-states.js';

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
   *   isOn(label, state), setOn(label, on, state)   read / write a state, 'S' or 'T' (setOn must call back
   *                                        sim.notify)
   *   setProps(label, patch, state)        color / gain / sigma
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
    this.handlers = []; // { type: 'excite' | 'decay', state: 'S' | 'T' | '*', match(label), fn }
    this.gen = new Map(); // 'label|S' -> excitation generation (for decay())
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

  /** The host calls this on every on ⇄ off change of a state ('S' / 'T'), from the script or from a click. */
  notify(label, on, state = 'S') {
    const key = `${label}|${state}`;
    if (on) this.gen.set(key, (this.gen.get(key) ?? 0) + 1);
    if (!this.started || this.status === 'error') return;
    const type = on ? 'excite' : 'decay';
    for (const h of this.handlers) {
      if (h.type !== type || (h.state !== '*' && h.state !== state) || !h.match(label)) continue;
      if (!this._guard(() => h.fn(label, state))) return;
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
    const key = (state) => stateKey(state ?? 'S'); // 'S' | 'T', or throws
    const keys = (state) => (state === '*' ? STATES : [key(state)]); // '*' = both
    const isOn = (label, state) => keys(state).some((k) => host.isOn(sim._label(label), k));
    const on = (label, state) => {
      const k = key(state);
      return !host.isOn(sim._label(label), k) && (host.setOn(label, true, k), true);
    };
    const off = (label, state) => {
      let changed = false;
      for (const k of keys(state)) {
        if (host.isOn(sim._label(label), k)) {
          host.setOn(label, false, k);
          changed = true;
        }
      }
      return changed;
    };
    const handler = (type, who, fn, state) => {
      if (typeof fn !== 'function') throw new TypeError(`on${type === 'excite' ? 'Excite' : 'Decay'}() needs a function`);
      sim.handlers.push({ type, state: state === '*' ? '*' : key(state), match: sim._matcher(who), fn });
    };
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
      toggle: (label, state) => {
        const k = key(state);
        if (isOn(label, k)) off(label, k);
        else on(label, k);
        return isOn(label, k);
      },
      transfer: (from, to, state, toState) => {
        sim._label(to); // fail before switching anything off
        off(from, key(state));
        return on(to, toState ?? state);
      },
      set: (label, patch, state) => {
        for (const k of keys(state)) host.setProps(sim._label(label), patch, k);
      },
      every: (interval, fn) => sim._schedule(sim.t + positive(interval, 'every() interval'), fn, interval),
      after: (delay, fn) => sim._schedule(sim.t + Math.max(0, delay), fn),
      at: (time, fn) => sim._schedule(Math.max(sim.t, time), fn),
      onExcite: (who, fn, state) => handler('excite', who, fn, state),
      onDecay: (who, fn, state) => handler('decay', who, fn, state),
      decay: (label, tau, state) => {
        const k = key(state);
        if (!isOn(label, k)) return null;
        const g = sim.gen.get(`${label}|${k}`);
        return sim._schedule(sim.t + exponential(positive(tau, 'decay() lifetime')), () => {
          if (sim.gen.get(`${label}|${k}`) === g) off(label, k);
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
