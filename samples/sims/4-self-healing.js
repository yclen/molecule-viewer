// Self-healing dye
// A dye with a triplet quencher linked to it (e.g. Cy5–COT), or a chlorophyll with its carotenoid (LHCII), under steady light. Run it on a dye–quencher molecule
// from the Self-healing dyes tab, then on the same dye alone (or set kTET to its minimum) to compare.
// Every site is in its ground state, its triplet T₁ or its singlet S₁, never two at once.
//   1. excitation: the dye goes to S₁ at kExc, only from its ground state
//   2. the dye's S₁ lasts tauS: it emits a photon, or with chance iscYield crosses to its T₁ (a dark state)
//   3. the dye's T₁ goes back to the ground state on its own at kT. While it sits there the dye is dark (blinking)
//   4. bleaching: from T₁ the dye reacts (with O₂, via singlet oxygen or radicals) at kBleach. A bleached dye is dead:
//      a fresh dye molecule replaces it, so the counters add up over many molecules
//   5. TTET: the dye's T₁ hands its energy to the quencher at kTET, if the quencher is in its ground state. Works
//      because the quencher's T₁ is lower than the dye's; the quencher's S₁ is higher, so singlets stay on the dye
//   6. the quencher's T₁ relaxes to the ground state at kQ (COT: it planarises and drops back fast)
// Photons per molecule (photons before bleaching) is the photon budget a self-healing dye is meant to raise.
// Rates are illustrative, in sim seconds; the real ones are ns (S₁), µs (TTET, COT) and µs–ms (dye T₁).

// dye excitations per sim second (only while it is in the ground state)
const kExc = slider('kExc', 2, { min: 0.1, max: 20, log: true, unit: '/s' });
// share of the dye's singlets that cross to the triplet (exaggerated so you can watch it; Cy5's is ~0.1-1%)
const iscYield = slider('iscYield', 0.05, { min: 0.001, max: 0.5, log: true });
// the dye's own triplet decay per sim second (slow: this is the dark time of a blink)
const kT = slider('kT', 0.2, { min: 0.01, max: 10, log: true, unit: '/s' });
// bleaching reactions per sim second while the dye is in its triplet
const kBleach = slider('kBleach', 0.05, { min: 0.001, max: 1, log: true, unit: '/s' });
// dye-to-quencher triplet transfers per sim second (the linker sets this: through space or through bonds)
const kTET = slider('kTET', 20, { min: 0.01, max: 300, log: true, unit: '/s' });
// the quencher's triplet decay per sim second
const kQ = slider('kQ', 30, { min: 0.1, max: 300, log: true, unit: '/s' });
const tauS = 0.3; // the dye's singlet lifetime, s (mean)

// In LHCII (14 chlorophylls, 4 carotenoids) it follows Chl a610, the lowest-energy chlorophyll, and Lutein 620 next to it
const chls = sites('chlorophyll');
const dye = sites('dye')[0] ?? (chls.includes('Chl a610') ? 'Chl a610' : chls[0]) ?? sites('emitter')[0] ?? sites()[0];
const quencher = sites('quencher')[0];
if (!quencher) log(`No quencher on this molecule: ${dye} blinks and bleaches on its own.`);

const ground = (site) => !isOn(site, '*');
for (const s of sites()) off(s, '*');

let darkSince = null; // when the dye's current triplet started
let dark = 0; // total time the dye spent in its triplet

every(0.002, (dt) => {
  if (ground(dye) && chance(kExc(), dt)) {
    on(dye, 'S'); // rule 1
    count('excitations');
  } else if (isOn(dye, 'T')) {
    const kOwn = kT();
    const kB = kBleach();
    const kQuench = quencher && ground(quencher) ? kTET() : 0;
    if (chance(kOwn + kB + kQuench, dt)) {
      const r = random() * (kOwn + kB + kQuench);
      off(dye, 'T');
      dark += now() - darkSince;
      if (r < kQuench) {
        on(quencher, 'T'); // rule 5
        decay(quencher, 1 / kQ(), 'T'); // rule 6
        count('triplets quenched (TTET)');
      } else if (r < kQuench + kB) {
        count('dyes bleached'); // rule 4: this molecule is gone, the next one starts in its ground state
      } else {
        count('triplets decayed on the dye'); // rule 3
      }
    }
  }
});

// the dark time starts whenever the dye gets a triplet (also when you switch it on by hand)
onExcite(dye, () => {
  darkSince = now();
}, 'T');

// rule 2: the singlet leaves after tauS, as a photon or into the triplet
onExcite(dye, () => {
  after(exponential(tauS), () => {
    if (!isOn(dye, 'S')) return;
    if (random() < iscYield()) {
      transfer(dye, dye, 'S', 'T');
      count('triplets made');
    } else {
      off(dye, 'S');
      count('photons');
    }
  });
}, 'S');

every(1, () => {
  const t = now();
  if (t <= 0) return;
  const darkNow = dark + (isOn(dye, 'T') ? t - darkSince : 0);
  show('time dark (in T₁)', `${((100 * darkNow) / t).toFixed(1)} %`);
  const bleached = countOf('dyes bleached');
  // photons per molecule: total photons over the molecules used up so far (the current one counts as one more)
  show('photons per molecule', Math.round(countOf('photons') / (bleached + 1)));
});
