// Walk along the molecule
// One excitation steps from site to site in molecule order, then starts over. Works on any molecule.

const stepTime = 0.6; // s per site

const all = sites();
let i = 0;
on(all[0]);

every(stepTime, () => {
  off(all[i]);
  i = (i + 1) % all.length;
  on(all[i]);
});
