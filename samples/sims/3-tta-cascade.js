// TTA cascade
// Sensitizer absorbs → its triplet moves to DPA 1 → triplets hop along the DPA wire → when two triplets meet,
// they annihilate (TTA) and the perylene emits one upconverted photon. A triplet that reaches the perylene stays
// there (it is the lowest triplet) until it decays or a second triplet arrives.
// For TTA Flora 4: try n = 2–6 DPA mediators and watch how often two triplets manage to meet.

const kAbs = 1.5; // sensitizer absorptions per second (ground state only)
const tTET = 0.25; // sensitizer → DPA 1 triplet energy transfer time, s (mean)
const kHop = 3; // hops per second for a triplet on a DPA (to a random neighbour)
const tauT = 5; // triplet lifetime, s (mean)
const tauS = 0.3; // perylene singlet (emission) lifetime, s (mean)

const TRIPLET = { color: '#3d8bff', gain: 1 };
const SINGLET = { color: '#7df9ff', gain: 2.2 };
const SENS = { color: '#ff8a3d', gain: 1.2 };

const sens = sites('sensitizer')[0];
const emit = sites('emitter')[0];
const wire = sites('mediator');
if (!sens || !emit) {
  log('This script needs a sensitizer and an emitter: open the TTA Flora 4 example.');
  return;
}
const chain = [...wire, emit]; // where triplets live, in order from the sensitizer

set(sens, SENS);
for (const s of chain) off(s);
let singlet = false; // is the perylene emitting (rather than holding a triplet)

// 1. absorption, then triplet transfer into the chain
every(0.02, (dt) => {
  if (!isOn(sens) && chance(kAbs, dt)) on(sens);
});
onExcite(sens, () => {
  count('absorbed');
  after(exponential(tTET), () => {
    if (!isOn(sens)) return;
    off(sens);
    arrive(chain[0]);
  });
});

// 2. a triplet lands on a site: free → it sits there; taken → two triplets meet: TTA
function arrive(site) {
  if (!isOn(site)) {
    set(site, TRIPLET);
    on(site);
    decay(site, tauT);
    return;
  }
  if (site === emit && singlet) {
    off(site); // perylene is busy emitting: the arriving triplet just quenches here
    count('quenched');
    return;
  }
  off(site); // both triplets used up
  count('TTA events');
  emitPhoton();
}

function emitPhoton() {
  off(emit); // a triplet sitting on the perylene is replaced by the singlet
  singlet = true;
  set(emit, SINGLET);
  on(emit);
  count('photons out');
  decay(emit, tauS);
}
onDecay(emit, () => {
  singlet = false;
});

// 3. hopping along the DPA wire (the perylene is a trap: nothing hops out of it)
every(0.02, (dt) => {
  const moving = wire.filter((d) => isOn(d) && chance(kHop, dt));
  for (const d of moving) {
    if (!isOn(d)) continue; // already used up this step
    const k = chain.indexOf(d);
    const to = pick([chain[k - 1], chain[k + 1]].filter(Boolean));
    off(d);
    arrive(to);
  }
});

// upconversion yield = photons out / absorbed (2 absorbed per photon at best)
every(1, () => {
  const a = countOf('absorbed');
  if (a) show('photons per absorption', countOf('photons out') / a);
});
