// Excited states of each molecule (site) in the scene.
//
// Every site carries two excited states, the singlet S₁ ('S') and the triplet T₁ ('T'). Each one has:
//   on     - is it populated right now (the only thing a time-based simulation needs to flip)
//   color  - glow colour, '#rrggbb'
//   gain   - brightness multiplier
//   sigma  - spread: Gaussian width of each glowing atom, in Å
//   fade   - afterglow: time constant of the glow going off, in real seconds. 0 = the glow appears and vanishes
//            in the same frame, so a fast-hopping triplet leaves no trail
// The two are independent: either, both or neither can be on, and their glows add up.
//
// States are kept per object (an example or a file) and per site label, so they survive rebuilding the
// object (e.g. changing n on a wire: "DPA 2" keeps its states) and switching to another example and back.

export const STATES = ['S', 'T'];
export const STATE_NAME = { S: 'singlet', T: 'triplet' };
export const STATE_SYMBOL = { S: 'S₁', T: 'T₁' };

/**
 * Starting colours by site kind. Singlets: warm for the sensitizer, blue for DPA mediators, cyan for the emitter.
 * Triplets are redder: the sensitizer's a slightly redder orange, DPA's and perylene's a full red. A dye is JF525's
 * yellow-green (it emits at 549 nm), a fluorescent protein mNeonGreen's green, and retinal a purple to set it apart.
 */
const KIND_COLOR = {
  sensitizer: { S: '#ff8a3d', T: '#ff5f2a' },
  mediator: { S: '#3d8bff', T: '#ff3b2f' },
  emitter: { S: '#4de3ff', T: '#e8243c' },
  molecule: { S: '#3d8bff', T: '#ff3b2f' },
  dye: { S: '#b4ff3a', T: '#ff3b2f' },
  fp: { S: '#3dff8a', T: '#ff3b2f' },
  retinal: { S: '#c65dff', T: '#ff3b2f' },
};

/** Singlets keep a soft afterglow when they emit; triplets switch instantly so hops look fast. */
const DEFAULT_FADE = { S: 0.25, T: 0 };

const defaultOne = (color, fade) => ({ on: false, color, gain: 1, sigma: 1.1, fade });

export const defaultState = (kind) => {
  const c = KIND_COLOR[kind] ?? KIND_COLOR.molecule;
  return { S: defaultOne(c.S, DEFAULT_FADE.S), T: defaultOne(c.T, DEFAULT_FADE.T) };
};

/** 'S' | 'T' from 'S', 'T', 'singlet', 'triplet', 'S1', 'T1' (any case). Throws on anything else. */
export function stateKey(state = 'S') {
  const s = String(state).trim().toLowerCase();
  if (s === 's' || s === 's1' || s === 'singlet') return 'S';
  if (s === 't' || s === 't1' || s === 'triplet') return 'T';
  throw new Error(`unknown state '${state}' (use 'S' for the singlet or 'T' for the triplet)`);
}

export const copyState = (st) => ({ S: { ...st.S }, T: { ...st.T } });

export class ExcitedStates {
  constructor() {
    this.objects = new Map(); // object key -> Map(label -> { S, T })
    this.states = new Map(); // the current object's states
    this.sites = []; // the current object's sites: [{ label, kind, atoms }]
  }

  /** Switch to an object's sites; existing states are kept, new sites get defaults. */
  setSites(objectKey, sites) {
    if (!this.objects.has(objectKey)) this.objects.set(objectKey, new Map());
    this.states = this.objects.get(objectKey);
    this.sites = sites;
    for (const s of sites) {
      if (this.states.has(s.label)) continue;
      const st = defaultState(s.kind);
      st.S.on = Boolean(s.on); // an example can start lit, to show where its chromophore is
      this.states.set(s.label, st);
    }
  }

  /** get('DPA 1') → { S, T }; get('DPA 1', 'T') → the triplet alone. */
  get(label, state) {
    const st = this.states.get(label);
    return state === undefined ? st : st?.[state];
  }

  /** Change part of one state, e.g. set('DPA 1', 'T', { on: true }). */
  set(label, state, patch) {
    const s = this.get(label, state);
    if (s) Object.assign(s, patch);
    return s;
  }

  /** Is any state of this site on. */
  anyOn(label) {
    const st = this.states.get(label);
    return Boolean(st && (st.S.on || st.T.on));
  }

  /** Sites of the current object with at least one state on. */
  excited() {
    return this.sites.filter((s) => this.anyOn(s.label));
  }
}
