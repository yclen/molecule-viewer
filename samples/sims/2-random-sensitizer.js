// Random sensitizer excitation
// The sensitizer absorbs at random (a fixed chance per second) and decays after a random lifetime.
// On a molecule without a sensitizer, the first site plays its part.

const kAbs = 0.8; // absorptions per second, while in the ground state
const tau = 1.2; // excited-state lifetime, s (mean)

const sens = sites('sensitizer')[0] ?? sites()[0];

every(0.02, (dt) => {
  if (!isOn(sens) && chance(kAbs, dt)) on(sens);
});

// also runs when you switch it on by hand
onExcite(sens, () => {
  count('absorptions');
  decay(sens, tau);
});
