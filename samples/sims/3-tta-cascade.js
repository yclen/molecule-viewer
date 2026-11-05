// TTA cascade
// Every site has a singlet S₁ and a triplet T₁ (their colours are set in the editor: triplets are the redder ones).
// Sensitizer absorbs (S₁) → intersystem crossing to its T₁ → triplet energy transfer to DPA 1 (T₁) → triplets hop
// along the DPA wire → when two triplets meet they annihilate (TTA): that site goes to S₁, the singlet passes to the
// perylene, and the perylene emits one upconverted photon. A triplet that reaches the perylene stays there (it is
// the lowest triplet) until it decays or a second triplet arrives.
// For TTA Flora 4: try n = 2–6 DPA mediators and watch how often two triplets manage to meet.

const kAbs = 1.5; // sensitizer absorptions per second (ground state only)
const tISC = 0.1; // sensitizer S₁ → T₁ intersystem crossing time, s (mean)
const tTET = 0.25; // sensitizer T₁ → DPA 1 T₁ triplet energy transfer time, s (mean)
const kHop = 3; // hops per second for a triplet on a DPA (to a random neighbour)
const tauT = 5; // triplet lifetime, s (mean)
const tFRET = 0.08; // DPA S₁ → perylene S₁ singlet energy transfer time, s (mean)
const tauS = 0.3; // perylene singlet (emission) lifetime, s (mean)

const sens = sites('sensitizer')[0];
const emit = sites('emitter')[0];
const wire = sites('mediator');
if (!sens || !emit) {
  log('This script needs a sensitizer and an emitter: open the TTA Flora 4 example.');
  return;
}
const chain = [...wire, emit]; // where triplets live, in order from the sensitizer

off(sens, '*');
for (const s of chain) off(s, '*');

// 1. absorption (singlet), intersystem crossing (→ triplet), then triplet transfer into the chain
every(0.02, (dt) => {
  if (!isOn(sens, '*') && chance(kAbs, dt)) on(sens, 'S');
});
onExcite(sens, () => {
  count('absorbed');
  after(exponential(tISC), () => {
    if (isOn(sens, 'S')) transfer(sens, sens, 'S', 'T');
  });
}, 'S');
onExcite(sens, () => {
  after(exponential(tTET), () => {
    if (!isOn(sens, 'T')) return;
    off(sens, 'T');
    arrive(chain[0]);
  });
}, 'T');

// 2. a triplet lands on a site: free → it sits there; taken → two triplets meet: TTA
function arrive(site) {
  if (site === emit && isOn(emit, 'S')) {
    count('quenched'); // perylene is busy emitting: the arriving triplet is lost
    return;
  }
  if (!isOn(site, 'T')) {
    on(site, 'T');
    decay(site, tauT, 'T');
    return;
  }
  off(site, 'T'); // both triplets used up ...
  count('TTA events');
  on(site, 'S'); // ... and one singlet made
  if (site === emit) emitPhoton();
  else {
    after(exponential(tFRET), () => {
      if (!isOn(site, 'S')) return;
      off(site, 'S');
      if (isOn(emit, 'T')) off(emit, 'T'); // the singlet replaces a triplet sitting on the perylene
      if (on(emit, 'S')) emitPhoton();
      else count('quenched');
    });
  }
}

function emitPhoton() {
  count('photons out');
  decay(emit, tauS, 'S');
}

// 3. hopping along the DPA wire (the perylene is a trap: nothing hops out of it)
every(0.02, (dt) => {
  const moving = wire.filter((d) => isOn(d, 'T') && chance(kHop, dt));
  for (const d of moving) {
    if (!isOn(d, 'T')) continue; // already used up this step
    const k = chain.indexOf(d);
    const to = pick([chain[k - 1], chain[k + 1]].filter(Boolean));
    off(d, 'T');
    arrive(to);
  }
});

// upconversion yield = photons out / absorbed (2 absorbed per photon at best)
every(1, () => {
  const a = countOf('absorbed');
  if (a) show('photons per absorption', countOf('photons out') / a);
});
