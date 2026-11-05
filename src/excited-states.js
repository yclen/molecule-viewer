// Excited state of each molecule (site) in the scene.
//
// Every site carries one state E:
//   on     - is it excited right now (the only thing a time-based simulation needs to flip)
//   color  - glow colour, '#rrggbb'
//   gain   - brightness multiplier
//   sigma  - spread: Gaussian width of each glowing atom, in Å
//
// States are kept per object (an example or a file) and per site label, so they survive rebuilding the
// object (e.g. changing n on a wire: "DPA 2" keeps its state) and switching to another example and back.

/** Starting colour by site kind: warm for the sensitizer, blue for DPA mediators, cyan for the emitter. */
const KIND_COLOR = {
  sensitizer: '#ff8a3d',
  mediator: '#3d8bff',
  emitter: '#4de3ff',
};

export const defaultState = (kind) => ({
  on: false,
  color: KIND_COLOR[kind] ?? '#3d8bff',
  gain: 1,
  sigma: 1.1,
});

export class ExcitedStates {
  constructor() {
    this.objects = new Map(); // object key -> Map(label -> state)
    this.states = new Map(); // the current object's states
    this.sites = []; // the current object's sites: [{ label, kind, atoms }]
  }

  /** Switch to an object's sites; existing states are kept, new sites get defaults. */
  setSites(objectKey, sites) {
    if (!this.objects.has(objectKey)) this.objects.set(objectKey, new Map());
    this.states = this.objects.get(objectKey);
    this.sites = sites;
    for (const s of sites) if (!this.states.has(s.label)) this.states.set(s.label, defaultState(s.kind));
  }

  get(label) {
    return this.states.get(label);
  }

  /** Change part of a site's state, e.g. set('DPA 1', { on: true }). */
  set(label, patch) {
    const s = this.states.get(label);
    if (s) Object.assign(s, patch);
    return s;
  }

  /** Change every site of the current object, e.g. setAll({ on: false }). */
  setAll(patch) {
    for (const site of this.sites) this.set(site.label, patch);
  }

  /** Sites of the current object that are excited. */
  excited() {
    return this.sites.filter((s) => this.states.get(s.label)?.on);
  }
}
