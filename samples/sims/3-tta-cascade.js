// TTA cascade
// The sensitizer gets excited at random and hands its triplet to the nearest DPA. Triplets random-walk along
// the DPA chain: each one hops to a neighbouring DPA (forward or back, never back onto the sensitizer). The
// perylene annihilator is a trap: a triplet that hops onto it stays there until a second triplet arrives.
//   - two triplets meet on a DPA (one hops onto the other): homoTTA. They annihilate and that DPA goes to its
//     singlet S₁. A higher hop rate (kHop) gets triplets to the perylene sooner, so there is less homoTTA.
//   - two triplets meet on the perylene: the perylene goes to its singlet S₁ and emits the upconverted photon
// For TTA Flora 4: try n = 2–6 DPA mediators and see where the triplets meet, then raise kHop to cut homoTTA.

const kExc = 1.5; // sensitizer excitations per second (only while it is in the ground state)
const tTET = 0.1; // sensitizer → nearest DPA triplet transfer time, s (mean): quick
const kHop = 3; // hops per second for a triplet on a DPA (higher = less homoTTA)
const tauT = 5; // triplet lifetime, s (mean); Infinity = triplets never decay on their own
const tauS = 0.3; // singlet lifetime (DPA or perylene emission), s (mean)

const sens = sites('sensitizer')[0];
const emit = sites('emitter')[0];
const wire = sites('mediator');
if (!sens || !emit) {
  log('This script needs a sensitizer and an annihilator: open the TTA Flora 4 example.');
  return;
}
const chain = [...wire, emit]; // where triplets can live, from the sensitizer's side; chain[0] is the nearest one

off(sens, '*');
for (const s of chain) off(s, '*');

// 1. the sensitizer gets excited (its triplet), then passes it to the nearest DPA (the perylene if n = 0)
every(0.02, (dt) => {
  if (!isOn(sens, '*') && chance(kExc, dt)) on(sens, 'T');
});
onExcite(sens, () => {
  count('sensitizer excitations');
  after(exponential(tTET), () => {
    if (!isOn(sens, 'T')) return;
    off(sens, 'T');
    arrive(chain[0]);
  });
}, 'T');

// 2. a triplet lands on a site: empty → it stays there; already holding a triplet → TTA, that site goes to S₁
function arrive(site) {
  if (!isOn(site, 'T')) {
    on(site, 'T');
    decay(site, tauT, 'T');
    return;
  }
  off(site, 'T'); // both triplets used up ...
  on(site, 'S'); // ... and that site is in its singlet
  decay(site, tauS, 'S');
  count(site === emit ? 'TTA on perylene (photon out)' : 'homoTTA (on a DPA)');
}

// 3. random walk: a triplet on a DPA hops to a random neighbour in the chain (the perylene never lets go)
every(0.02, (dt) => {
  const moving = wire.filter((d) => isOn(d, 'T') && chance(kHop, dt));
  for (const d of moving) {
    if (!isOn(d, 'T')) continue; // used up by an annihilation earlier in this step
    const k = chain.indexOf(d);
    const to = pick([chain[k - 1], chain[k + 1]].filter(Boolean)); // chain[-1] is undefined: no hop to the sensitizer
    off(d, 'T');
    arrive(to);
  }
});

// perylene photons per sensitizer excitation (two excitations per photon at best: 0.5)
every(1, () => {
  const a = countOf('sensitizer excitations');
  if (a) show('perylene photons per excitation', countOf('TTA on perylene (photon out)') / a);
});
