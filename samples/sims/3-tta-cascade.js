// TTA cascade
// The sensitizer gets excited at random into its singlet S₁ and stays there until intersystem crossing (ISC) turns
// it into its triplet T₁. That triplet hops to the nearest DPA at the same rate as any other hop. Triplets
// random-walk along the DPA chain: each one hops to a neighbouring DPA (forward or back, never back onto the
// sensitizer). The perylene annihilator is a trap: a triplet that hops onto it stays there until a second triplet arrives.
//   - two triplets meet on a DPA (one hops onto the other): homoTTA. They annihilate and that DPA goes to its
//     singlet S₁. A higher hop rate (kHop) gets triplets to the perylene sooner, so there is less homoTTA.
//   - two triplets meet on the perylene: the perylene goes to its singlet S₁ and emits the upconverted photon
// For TTA Flora 4: change n (the number of DPA mediators) to see where the triplets meet, and raise kHop to cut
// homoTTA.
// kExc, kISC, kHop and kDecay are the sliders under the controls: moving one changes a running simulation at once.
// Times and rates are in sim seconds (real seconds × Speed).

// sensitizer excitations per sim second (only while it is in the ground state)
const kExc = slider('kExc', 0.1, { min: 0.1, max: 20, log: true, unit: '/s' });
// sensitizer S₁ → T₁ conversions per sim second: 2.82 means it takes about 0.35 s, so you can watch the singlet
const kISC = slider('kISC', 2.82, { min: 0.1, max: 20, log: true, unit: '/s' });
// hops per sim second for a triplet, left and right together (higher = less homoTTA). Also the rate of the
// sensitizer's hop onto the nearest DPA
const kHop = slider('kHop', 54.2, { min: 1, max: 300, log: true, unit: '/s' });
// triplet decays per sim second: 0.001 means a triplet lasts about 1000 s, so few are lost before they meet
const kDecay = slider('kDecay', 0.001, { min: 0.001, max: 1, log: true, unit: '/s' });
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

// 1. the sensitizer gets excited into its singlet, and at each time step that singlet may convert to the triplet
every(0.002, (dt) => {
  if (!isOn(sens, '*') && chance(kExc(), dt)) on(sens, 'S');
  if (isOn(sens, 'S') && chance(kISC(), dt)) transfer(sens, sens, 'S', 'T');
});
onExcite(sens, () => count('sensitizer excitations'), 'S');

// 2. a triplet lands on a site: empty → it stays there; already holding a triplet → TTA, that site goes to S₁
function arrive(site) {
  if (!isOn(site, 'T')) {
    on(site, 'T');
    decay(site, 1 / kDecay(), 'T');
    return;
  }
  off(site, 'T'); // both triplets used up ...
  on(site, 'S'); // ... and that site is in its singlet
  decay(site, tauS, 'S');
  count(site === emit ? 'TTA on perylene (photon out)' : 'homoTTA (on a DPA)');
}

// 3. random walk: the moment a site gets a triplet, it gets a random time until its next hop (mean 1 / kHop). The
// sensitizer's triplet hops to the nearest DPA (the perylene if n = 0). A DPA's hops to a random neighbour in the
// chain. The perylene is not in this list, so it never lets go.
const nextHop = new Map(); // site -> timer of the triplet sitting on it
onExcite([sens, ...wire], (d) => {
  nextHop.get(d)?.stop();
  nextHop.set(d, after(exponential(1 / kHop()), () => {
    if (!isOn(d, 'T')) return; // it decayed or annihilated first
    const k = chain.indexOf(d);
    // chain[-1] is undefined, so a triplet on DPA 1 can't hop back to the sensitizer
    const to = d === sens ? chain[0] : pick([chain[k - 1], chain[k + 1]].filter(Boolean));
    off(d, 'T');
    arrive(to);
  }));
}, 'T');

// perylene photons per sensitizer excitation (two excitations per photon at best: 0.5)
every(1, () => {
  const a = countOf('sensitizer excitations');
  if (a) show('perylene photons per excitation', countOf('TTA on perylene (photon out)') / a);
});
