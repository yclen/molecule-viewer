// TTA cascade
// The same 10 rules as the standalone simulation (sim/README.md, sim/formulas.md), on one molecule under steady light.
// Every site is in its ground state, its triplet T₁ or its singlet S₁, never two at once.
//   1. excitation: the sensitizer goes to S₁ at kExc, only from its ground state
//   2. intersystem crossing: the sensitizer's S₁ turns into T₁ at kISC
//   3. triplet decay: a triplet on a DPA or the perylene decays at kDecay (the sensitizer's doesn't)
//   4. emission: a singlet on a DPA or the perylene emits a photon after tauS (on average)
//   5. hop: a triplet hops to a neighbour in its ground state. The sensitizer's goes to DPA 1, DPA 1's only forward,
//      the other DPAs' forward or back. kHop is the rate for both directions together. The perylene never lets go
//   6. blocked hop: a triplet can't hop onto a singlet. It stays put and tries again later
//   7. TTA: a triplet hops onto another triplet. With chance f the landing site goes to S₁ ...
//   8. ... otherwise one triplet survives, on the landing site
//   9. the sensitizer's S₁ can also just decay. tripletYield is the share that makes a triplet instead (rule 2)
//  10. back transfer (FRET): a singlet can go back to the sensitizer instead of emitting, with chance
//      E = 1 / (1 + (r / fretRadius)⁶), r = its distance from the sensitizer. Only while the sensitizer is in its ground
//      state: otherwise the singlet always emits
// For TTA Flora 4: change n (the number of DPA mediators) to see where the triplets meet, and raise kHop to cut
// homoTTA (TTA on a DPA instead of the perylene).
// All the sliders change a running simulation at once. Times and rates are in sim seconds (real seconds × Speed).

// sensitizer excitations per sim second (only while it is in the ground state)
const kExc = slider('kExc', 0.1, { min: 0.1, max: 20, log: true, unit: '/s' });
// sensitizer S₁ → T₁ conversions per sim second: 2.82 means it takes about 0.35 s, so you can watch the singlet
const kISC = slider('kISC', 2.82, { min: 0.1, max: 20, log: true, unit: '/s' });
// share of sensitizer singlets that become triplets: the rest decay (rule 9) at kISC × (1 − tripletYield) / tripletYield
const tripletYield = slider('tripletYield', 0.5, { min: 0.01, max: 1 });
// hops per sim second for a triplet, both directions together (higher = less homoTTA)
const kHop = slider('kHop', 54.2, { min: 1, max: 300, log: true, unit: '/s' });
// triplet decays per sim second: 0.001 means a triplet lasts about 1000 s, so few are lost before they meet
const kDecay = slider('kDecay', 0.001, { min: 0.001, max: 1, log: true, unit: '/s' });
// the spin factor: share of TTA events that make a singlet
const f = slider('f', 0.4, { min: 0, max: 1 });
// distance from the sensitizer where back transfer is 50% (0 = no back transfer)
const fretRadius = slider('fretRadius', 4, { min: 0, max: 10, unit: 'nm' });
const tauS = 0.3; // singlet lifetime (DPA or perylene emission), s (mean)

const sens = sites('sensitizer')[0];
const emit = sites('emitter')[0];
const wire = sites('mediator');
if (!sens || !emit) {
  log('This script needs a sensitizer and an annihilator: open the TTA Flora 4 example.');
  return;
}
const chain = [...wire, emit]; // where triplets can live, from the sensitizer's side; chain[0] is the nearest one
const n = wire.length;

// Distance in nm from the sensitizer to chain[k], between the same centres sim/rules.py uses (the sensitizer's core,
// each DPA's middle, the perylene's core). Matches its site_distances to 0.01 nm for every n
function distance(k) {
  if (k < n) return 1.403 + 1.3205 * k;
  return n ? 1.403 + 1.3205 * (n - 1) + 1.29 : 1.406;
}
function fret(k) {
  const R0 = fretRadius();
  return R0 > 0 ? 1 / (1 + (distance(k) / R0) ** 6) : 0;
}

const ground = (site) => !isOn(site, '*');

off(sens, '*');
for (const s of chain) off(s, '*');

// rules 1, 2 and 9, checked every time step
every(0.002, (dt) => {
  if (ground(sens)) {
    if (chance(kExc(), dt)) {
      on(sens, 'S');
      count('sensitizer excitations');
    }
  } else if (isOn(sens, 'S')) {
    const y = tripletYield();
    const kLoss = (kISC() * (1 - y)) / y; // rule 9
    if (chance(kISC() + kLoss, dt)) {
      if (random() * (kISC() + kLoss) < kISC()) {
        transfer(sens, sens, 'S', 'T'); // rule 2
        count('triplets made');
      } else {
        off(sens, 'S'); // rule 9
        count('sensitizer singlets lost');
      }
    }
  }
});

// rules 5, 7 and 8: the triplet on `from` hops onto `to`
function hop(from, to) {
  off(from, 'T');
  if (ground(to)) {
    on(to, 'T'); // rule 5
    decay(to, 1 / kDecay(), 'T'); // rule 3
    return;
  }
  const where = to === emit ? 'on the perylene' : 'on a DPA';
  if (random() < f()) {
    off(to, 'T'); // rule 7: both triplets used up, and the landing site is in its singlet
    on(to, 'S');
    count(`TTA singlet ${where}`);
  } else {
    count(`TTA failed ${where}`); // rule 8: the landing site keeps its triplet
  }
}

// The moment a site gets a triplet, it gets a random wait until its next hop try (mean 1 / kHop). A try picks a
// direction, and if that neighbour is in its singlet, nothing happens (rule 6) and the next wait starts. The perylene
// is not in this list, so it never lets go
const nextHop = new Map(); // site -> timer of the triplet sitting on it
function waitToHop(d) {
  nextHop.get(d)?.stop();
  nextHop.set(d, after(exponential(1 / kHop()), () => {
    if (!isOn(d, 'T')) return; // it decayed or hopped away first
    const k = chain.indexOf(d);
    // chain[-1] is undefined, so a triplet on DPA 1 can't hop back to the sensitizer
    const to = d === sens ? chain[0] : pick([chain[k - 1], chain[k + 1]].filter(Boolean));
    if (isOn(to, 'S')) waitToHop(d); // rule 6
    else hop(d, to);
  }));
}
onExcite([sens, ...wire], (d) => waitToHop(d), 'T');

// rules 4 and 10: a singlet on a DPA or the perylene leaves after tauS, as a photon or back to the sensitizer
const leave = new Map(); // site -> timer of its singlet
onExcite(chain, (site) => {
  leave.get(site)?.stop();
  leave.set(site, after(exponential(tauS), () => {
    if (!isOn(site, 'S')) return;
    off(site, 'S');
    if (ground(sens) && random() < fret(chain.indexOf(site))) {
      on(sens, 'S'); // rule 10
      count('back transfers');
    } else {
      count(site === emit ? 'photons from the perylene' : 'photons from DPAs'); // rule 4
    }
  }));
}, 'S');

// photons per sensitizer excitation, as the standalone simulation counts them (back transfers aren't new excitations)
every(1, () => {
  const a = countOf('sensitizer excitations');
  if (!a) return;
  const perylene = countOf('photons from the perylene');
  show('photons per excitation', (perylene + countOf('photons from DPAs')) / a);
  show('perylene photons per excitation', perylene / a);
});
