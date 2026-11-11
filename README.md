# Molecule Viewer

A browser-based 3D viewer for `.xyz` molecule files, built with [Three.js](https://threejs.org) and [Vite](https://vite.dev).

## Open it (desktop shortcut)

Double-click **Molecule Viewer** on the Desktop. It opens `dist/index.html` in its own Chrome app window.
No server or terminal is needed: the whole app (code, styles, samples) is bundled into that one HTML file.

The shortcut shows the last *build*. After changing the code, rebuild:

```bash
npm run build
```

If you move the project folder, recreate the shortcut with `npm run shortcut`.

## Develop it

```bash
npm install      # first time only
npm run dev      # then open http://localhost:5173 (auto-reloads on save)
```

Drop an `.xyz` file onto the page, use **Open file…**, or pick one of the built-in samples.

## Controls

| Action | How |
| --- | --- |
| Rotate (free, any axis) | left-drag |
| Pan | right-drag |
| Zoom | scroll |
| Display style | buttons, or keys `1` ball & stick · `2` sticks · `3` space-filling · `4` wireframe |
| Show / hide hydrogens | checkbox or `H` |
| Reset view | button or `R` |
| Open file | button or `O` |
| Play / pause trajectory | button or `Space` |
| Atom info | hover an atom |
| Excite / un-excite | the **S** / **T** dots next to a molecule under **Excited states**, or `E` (singlet) / `T` (triplet) for the selected one |
| Run / pause simulation | **Simulate** mode: button or `Space` |
| Collapse a section | click its heading (remembered between visits) |

## How it works

- **`src/parser.js`**: reads XYZ (atom count line, comment line, then `Element x y z` lines in Å). Handles
  multi-frame trajectories, extended-XYZ extra columns, atom labels like `C12`, atomic numbers instead of symbols,
  and header-less files.
- **`src/bonds.js`**: XYZ has no bonds, so they're inferred: two atoms bond when
  `distance < r_cov(a) + r_cov(b) + 0.4 Å`. A spatial grid keeps it O(n).
- **`src/elements.js`**: covalent radii (Cordero 2008), van der Waals radii, Jmol/CPK colours for H–Rn.
- **`src/viewer.js`**: Three.js scene. Atoms and bonds are `InstancedMesh`es (one draw call each), so thousands
  of atoms stay smooth. Molecules are centred and rotated onto their principal axes, so flat molecules open face-on.
- **`src/main.js`**: UI wiring (panel, file loading, drag & drop, trajectory playback, tooltip, shortcuts).

## Excited state

The **Excited states** section puts a glowing, pulsing Gaussian cloud on excited chromophores. Parametric
objects list their chromophores (Sensitizer, DPA 1…n, Perylene); any other molecule has **Whole molecule**.

Each one has two excited states (`src/excited-states.js`), the **singlet S₁** and the **triplet T₁**, and each
state has its own **on/off**, **colour**, **brightness**, **spread** (Gaussian width in Å) and **afterglow**
(how long its glow takes to fade, in real seconds). The two are
independent: a molecule can have either, both or neither on, and their glows add up. Click a molecule in the
list (or one of its atoms in the 3D view) to select it; the editor below the list shows both of its states and
changes only that molecule, and every other molecule keeps what was set for it. The **S** and **T** dots beside
each name switch the two states on and off, `E` / `T` toggle the selected molecule's singlet / triplet, and
**All S₁ on / All T₁ on / All off** switch the whole object. States are remembered per example and per name, so
changing n keeps DPA 1's states, and switching to another example and back restores them.

Default colours: singlets are sensitizer orange, DPA blue and perylene cyan. Triplets are redder: the
sensitizer's a slightly redder orange, DPA's and perylene's red. A dye glows yellow-green with a red triplet, and a
triplet quencher (COT, a carotenoid) magenta when it holds a triplet, and a chlorophyll leaf green. Singlets have a 0.25 s afterglow. Triplets have
none: their glow appears and vanishes in the same frame, so a hopping triplet leaves no trail. Both states glow over the same atoms (S₁ and T₁
come from the same HOMO→LUMO pair). In dev builds `__setOn('DPA 2', true, 'T')` drives it from the console.

- **Where it glows**: each atom's share of the lowest π→π* excitation, (c²_HOMO + c²_LUMO) / 2, from Hückel
  theory on the 3D structure (`src/excitation.js`). Resonance integrals are scaled by p-orbital overlap
  (|n_i · n_j|), so twisted rings drop out by themselves: DPA's 90° phenyls stay dark and the glow sits on the
  anthracene (C9/C10 brightest, fused carbons zero), and fluorescein's benzoate stays dark while the xanthene
  glows. Heteroatoms use standard Hückel parameters. It is a qualitative picture, not a TD-DFT density.
- **How it's drawn** (`src/glow.js`): every glowing atom is a camera-facing sprite with a Gaussian falloff,
  added up in an HDR buffer (the line-of-sight integral of a 3D Gaussian cloud). Atoms and bonds in front hide
  it, a quarter-resolution blur adds bloom to the glow only, and a `1 − exp(−I)` curve keeps faint regions
  coloured and saturates dense ones toward white. Each chromophore's singlet and triplet glows have their own
  level in the shader. With an afterglow, a state fades in and out, and even a very short excitation shows as a
  brief flash. With afterglow 0, the level follows the state exactly.

## Simulation mode

**View / Simulate** at the top of **Excited states** switches between setting states by hand and letting a
script set them over time (`src/simulation.js`). A script is plain JavaScript, run once at t = 0. It starts
background loops that a clock then runs: **Run / Pause** (`Space`), **Step** (0.1 s), **Reset** (stops and
puts every state back as it was), and **Speed** (0.1×–20×). Time is in "sim seconds": the clock runs at real
seconds × Speed, so at 1× one sim second is one real second. They are arbitrary units, not real photophysical
times. Every rate and lifetime in a script is per sim second, and the glow's afterglow is in real seconds.

```js
const kAbs = 0.8; // absorptions per second
const tau = 1.2; // lifetime, s
every(0.02, (dt) => {
  if (!isOn('Sensitizer') && chance(kAbs, dt)) on('Sensitizer');
});
onExcite('Sensitizer', () => decay('Sensitizer', tau));
```

- **States**: `on(label, state?)`, `off`, `toggle`, `isOn`, `transfer(from, to, state?, toState?)`,
  `set(label, { color, gain, sigma, fade }, state?)`, `sites(kind?)`, `kindOf`, `neighbors`. `state` is `'S'`
  (singlet, the default, so scripts that never name a state drive the singlets) or `'T'` (triplet);
  `'singlet'` / `'triplet'` work too. `'*'` means either/both for `isOn`, `off` and `set`.
  `transfer(sens, sens, 'S', 'T')` is intersystem crossing, `transfer('DPA 1', 'DPA 2', 'T')` a triplet hop.
- **Time**: `now()`, `every(interval, fn(dt))`, `after(delay, fn)`, `at(time, fn)` (each returns `.stop()`),
  `onExcite(who, fn(label, state), state?)` / `onDecay(...)` (who = a label, a kind, an array or `'*'`;
  state defaults to `'S'`, `'*'` for either; clicks trigger them too), `decay(label, tau, state?)` (off after
  a random exponential lifetime, unless it was already switched off).
- **Chance**: `chance(rate, dt)` is true with probability 1 − e^(−rate·dt), so results don't depend on the step
  size. Also `random()`, `exponential(tau)`, `pick(array)`, `seed(n)` (repeatable runs).
- **Sliders**: `const kHop = slider('kHop', 30, { min: 1, max: 300, log: true, unit: '/s' })` puts a slider under
  the controls, and `kHop()` reads it. Moving it changes a running script at once. The app reads sliders from the
  script's text, so they show before Run, as long as the name and numbers are plain literals. Slider positions
  are saved per script, and **Revert** resets them too.
- **Output**: `count(name, n)`, `countOf(name)`, `show(name, value)`, `log(...)` show under the controls.
  `pause()` stops the clock.

The presets are in `samples/sims/` (the first comment line is the name): **Walk along the molecule**,
**Random sensitizer excitation**, **TTA cascade** and **Self-healing dye** (below). TTA cascade runs on TTA Flora 4 and follows the same 10
rules as the standalone simulation in `sim/` (listed in `sim/README.md`, as equations in `sim/formulas.md`), on one
molecule under steady light. Every site is in its ground state, T₁ or S₁, never two at once:

- The sensitizer gets excited at random into its S₁ (`kExc`). That singlet either turns into T₁ by intersystem
  crossing (`kISC`, default 2.82 per second, so about 0.35 s) or just decays. `tripletYield` (default 0.5) is the
  share that makes a triplet.
- The sensitizer's triplet hops to DPA 1, and triplets random-walk along the DPA chain, forward or back but never
  back onto the sensitizer. `kHop` (default 54.2 per second, both directions together) sets each triplet's random
  wait before its next hop. A triplet can't hop onto a singlet: it stays and tries again later. The perylene is a
  trap and never lets go.
- A triplet that hops onto another triplet annihilates with it (TTA). With chance `f` (default 0.4) the landing
  site goes to S₁; otherwise one triplet survives there. On a DPA this is homoTTA, which a higher `kHop` cuts by
  getting triplets to the perylene sooner.
- A singlet on a DPA or the perylene emits a photon after 0.3 s on average, unless it goes back to the sensitizer
  by FRET first: chance 1 / (1 + (r / `fretRadius`)⁶) at distance r, and only while the sensitizer is in its
  ground state.
- A triplet on a DPA or the perylene decays on its own at `kDecay` (default 0.001 per second, so about 1000 s).

`kExc`, `kISC`, `tripletYield`, `kHop`, `kDecay`, `f` and `fretRadius` are sliders. The counters show each kind of
event, and photons per excitation counts all photons (DPAs and perylene) per sensitizer excitation, as `sim/` does.
**Self-healing dye** runs on a dye–quencher molecule from the Self-healing dyes tab (or on any dye, without the
quencher). On LHCII it follows Chl a610 and Lutein 620, the pair in contact. The dye cycles ground → S₁ → photon; with chance `iscYield` its S₁ crosses to T₁ instead, and the dye is
dark until the triplet leaves. It leaves on its own (`kT`), by bleaching the dye (`kBleach`, after which a fresh
molecule takes over), or by triplet–triplet energy transfer to the quencher (`kTET`, only while the quencher is in its
ground state), whose triplet then relaxes at `kQ`. The counters show the share of time spent dark and the photons
per molecule before it bleaches. Set `kTET` to its minimum to see the dye without its quencher.

Edits are saved in the browser per script, and **Revert** goes back to the preset. Errors stop the run and
show the message (an unknown site, an `every(0)`, or handlers that keep switching each other). Changing the
molecule, or n, resets the run.

## Measuring

Click 2–4 atoms: the panel shows each distance, the angle(s) and the dihedral, and the 3D view labels each
segment. Click an atom again to unpick it, click empty space or press `Esc` to clear. Measurements update live
while a trajectory plays.

## Examples

Tabs (Basics, Triplet molecules, Self-healing dyes, Interesting photophysics, Xanthenes, Cyanines, Other dyes, Proteins, TTA upconversion, TADF, and two review tabs), all bundled into the app.
`↑` / `↓` steps through the current tab.

- **Basics** (`samples/basics/`): idealised geometries from `tools/make_samples.py` (`npm run samples`).
- **Everything else**: listed in `tools/catalog.py` and built by `tools/build_examples.py`
  (`npm run examples`, needs `pip install rdkit`):
  - *Small molecules*: structure (SMILES) from PubChem → RDKit ETKDGv3 conformer search → MMFF94s
    minimisation → lowest-energy conformer. Counter-ions are dropped. Bonds that MMFF gets wrong (C–I) or
    can't type (boron, so BODIPY uses UFF) are restrained to X-ray means. The PtOEP porphyrin core is
    placed from metalloporphyrin crystal averages, because MMFF puckers porphyrins.
  - *Proteins*: X-ray structures from the RCSB PDB (chain A, chromophore kept, solvent removed). HaloTag7–JF525 is a
    model, not a crystal structure of that dye: `tools/build_halotag.py` starts from HaloTag7 with TMR (6Y7A), keeps
    every crystal atom the two dyes share, and grows the difluoroazetidines onto them. Its glow marks the dye.
    Voltron525 and Ace2N-mNeon are composites too: `tools/build_voltron.py` joins rigid crystal domains (3AM6 with
    5LTR or the HaloTag7–JF525 model), regrows only the junction residues (sequence from US 11,598,792), and invents
    the pose of the two domains, so distances between them are illustrative. The rhodopsin is placed as OPM
    (Orientations of Proteins in Membranes) places 3AM6. "Voltron525 in a membrane" sets that model into a POPC
    bilayer (`tools/build_membrane.py`): the last frame of a published 270 ns MD run (Ollila, Ferreira and Topgaard,
    Zenodo 13279), tiled, cut to a disc, and stripped of lipids touching the protein. Nothing is relaxed afterwards
    and there is no water. Its index entry sets `surroundings` (the first n atoms, drawn faint) and `keepAxes`
    (file axes kept, so the membrane is seen side on).
- **Validation**: `docs/geometry-report.md` lists every bond type per molecule against typical
  crystallographic ranges, plus a bond-by-bond comparison with PubChem's own 3D conformer where one exists.

Downloads are cached in `tools/cache/`.

### Triplet molecules (parametric objects)

Objects you size in the app. Pick one and a control appears under the list (`−` / `+` buttons, slider, typed
number, or the `−` / `+` keys). The molecule rebuilds instantly and the panel shows its end-to-end length.

- **Molecular wire**: *n* 9,10-diphenylanthracene (DPA) units, each joined anthracene C2 → 1,4-phenylene → C6
  of the next one, with H on both ends (R = H). *n* = 1 to 20; *n* = 1 is plain DPA.
- **TTA Flora 4**: sensitizer – 1,4-phenylene – *n* DPA mediators (each followed by a 1,4-phenylene) –
  perylene (linked at C3). The sensitizer is 6,7-bis{4-[bis(4-methoxyphenyl)amino]phenyl}thieno[2,3-b]quinoxaline
  1,1-dioxide, linked at thiophene C2. *n* = 0 to 20; *n* = 0 is the direct sensitizer–linker–perylene dyad.
  The panel shows the sensitizer-to-perylene distance. The **Sensitizer** dropdown swaps in any sensitizer
  from the library (Review: TADF sensitizers, PtOEP, DMAC-TRZ; the list is `sensitizers.from` in
  `samples/triplet/index.json`). `src/attach.js` removes the original sensitizer from the end group and places
  the new one rigidly, using its own sample geometry, on the free phenylene carbon. It attaches through the
  outermost free aromatic C–H and picks the biaryl torsion that avoids clashes. The result isn't re-optimised.
  A few crowded sensitizers (PtOEP, MR-1, BN-2Cz-tBu, BBCz-SB-Br) keep one H contact of
  1.8–2.1 Å next to the new bond. The excited state of "Sensitizer" carries over when you swap.

How it's built: both anthracene-2,6-diyl and *p*-phenylene are centrosymmetric, so a straight wire is one
"DPA + linker" repeat unit copied along a fixed shift (13.2 Å). `tools/build_wire.py` (`npm run wire`)
builds an oligomer with RDKit + MMFF94s, takes its middle unit, and refines it in a relaxed pentamer until it is
self-consistent with its neighbours. The app (`src/parametric.js`) just stacks copies. `docs/triplet-report.md`
checks assembled wires: bond lengths in range, correct formula, and freely re-minimising an assembled wire moves
it by only ~0.01 Å RMSD. TTA Flora 4 reuses the same wire and adds two end groups from `tools/build_triplet.py` (`npm run triplet`): each gets its own conformer search, is attached where a wire linker would sit, and is relaxed against the fixed wire. To add a new object, give it a `param` entry in `tools/catalog.py` and a builder
in `src/parametric.js`.

### Self-healing dyes

A fluorescent dye with a triplet quencher bonded to it. The quencher's T₁ sits below the dye's T₁, so the dye's
triplet (a dark, reactive state: blinking, singlet oxygen, bleaching) drains into it by triplet–triplet energy transfer
(TTET), and its S₁ sits far above the dye's S₁, so the dye's fluorescence is left alone. Cyclooctatetraene (COT) is
the standard quencher: its T₁ relaxes to ~0.95 eV by flattening (Baird aromatic), its S₁ is ~4 eV and dark, and its
triplet drops back to the ground state. Cy5's T₁ is 1.60 eV and its S₁ 1.9 eV. Built with `blocks` in
`tools/catalog.py`, so the dye and the COT are separate sites. Each entry's About panel text is in
`tools/catalog_selfheal_info.py`.

- **Cy5–COT (Altman 2012)**: COT on a flexible tether from an indolenine nitrogen (Altman et al., *Nat. Methods* 9,
  68, 2012). TTET happens when the tether lets the two touch: Cy5's triplet lifetime falls from 63 µs to 1.1 µs with
  the shortest tether and 13 µs with a 13-atom one (Zheng et al., *J. Phys. Chem. Lett.* 3, 2200, 2012), and the
  photon budget rises up to ~70× without oxygen. The paper's 12-atom linker is in an SI figure we could not read, so
  the tether drawn here (hexanoyl amide to aminomethyl-COT) is illustrative.
- **Gentle Rhodamine GR555-mito**: TMR with COT esterified onto the pendant ring's carboxyl, COT–CH₂–O–C(=O)–Ar
  (Liu et al., *ACS Cent. Sci.* 10, 1933, 2024). Ten times less singlet oxygen and 4–5× less phototoxicity.
- **LHCII: chlorophylls + carotenoids (1RWT)**: nature's version, from the X-ray structure of spinach
  light-harvesting complex II (Liu et al., *Nature* 428, 287, 2004), chain A: 14 chlorophylls (8 a, 6 b) and 4
  carotenoids (2 lutein, violaxanthin, neoxanthin) in one protein, drawn faint (`tools/build_lhcii.py`;
  waters/lipids/detergent removed). X-ray maps at 2.72 Å can't see hydrogens (one electron each), so the deposited
  model has none; the pigments get theirs back from the PDB Chemical Component Dictionary (each atom's bonds, bond
  orders and H count, matched by atom name), placed by RDKit on the measured heavy atoms. Phytol atoms too disordered
  to be in the crystal stay out. The faint protein keeps no H, like the other proteins. Chlorophylls are chlorins (porphyrins with one ring reduced) with Mg in the middle
  and a phytol tail. Their triplet yield is high (~60% for Chl a in solution), and their triplets would make singlet oxygen; the carotenoids
  take them by TTET within nanoseconds (carotenoid T₁ ~0.9 eV is below chlorophyll's ~1.3 eV) and lose the energy as
  heat. There is no linker: the protein holds each lutein in van der Waals contact with the chlorophylls it guards
  (Lutein 620: Chl a610 3.6 Å, a612 3.7 Å; Lutein 621: Chl a602 3.4 Å, a603 3.7 Å, closest conjugated atoms), which
  is the Dexter overlap TTET needs. Each pigment is a site, but only its conjugated part glows: the macrocycle (Mg,
  4 N, 20 C) or the polyene chain. Chl a610, the lowest-energy chlorophyll and the main triplet donor, starts lit.
  Carotenoids also harvest light the other way (their S₂ absorbs blue-green and passes the energy to chlorophyll).
- **Chlorophyll + carotenoid pair (Chl a610 + Lutein 620)**: the same two pigments on their own, cut out of 1RWT at
  their crystal positions with the protein and everything else removed (whole molecules, phytol tail included, with their H: C55H72MgN4O5 and C40H56O2).
  It is the bare donor–quencher pair to compare with the self-healing dyes: the lutein's polyene lies across the edge
  of the chlorophyll's ring, 3.6 Å away, held there only by the protein that is no longer drawn.
- **Cy5–(phenylene)n–COT (our design)**: a parametric object. COT sits on Cy5's benzo ring at C5 (where sulfo-Cy5
  carries a sulfonate), at the end of a rigid rod of *n* 1,4-phenylenes, the linker of TTA Flora 4; *n* = 0 bonds COT
  straight to the dye. The other indolenine keeps its sulfonate and hexanoic-acid handle. The literature search found
  no published dye–phenylene–COT and no COT C–C bonded to a cyanine or rhodamine core (closest: COT on porphyrins,
  Ostapko et al., *Chem. Eur. J.* 26, 16666, 2020). `tools/build_selfheal.py` searches the dye's conformers once and
  keeps that pose for every *n*, so only the rod and the COT move with *n*; each *n* is a finished geometry
  (`frames` builder in `src/parametric.js`), checked in `docs/selfheal-report.md`.

- **Self-healing JF525 (JF525–phenylene–carotenoid)**: our design, with the triplet quencher photosynthesis uses.
  In light-harvesting complexes, carotenoids sit next to the chlorophylls and take their triplets before they make
  singlet oxygen. Here JF525's pendant ring carries a 1,4-phenylene para to the xanthene (between the two carboxyls'
  neighbours, so the lactone carboxyl and the HaloTag handle both stay), and the phenylene carries 7'-apo-β-carotene:
  β-carotene cut at C7' and bonded there, all-trans. That aryl-ended carotenoid is the one Gust and Moore joined to
  porphyrins in their carotenoporphyrin dyads. A normal catalog entry (RDKit + MMFF94s). The force field twists JF525's
  ring and the phenylene 64° apart (a carboxyl is next door) and the phenylene and the polyene 42°; the molecule is
  ~3 nm long.

  Levels (approximate literature values, not measured on this molecule): β-carotene's T₁ is ~0.9 eV, far below a
  rhodamine's (~1.8–2.0 eV; JF525's own is not reported that we know of), so the triplet condition holds easily. The
  singlet condition does not: a carotenoid's lowest singlet is the dark S₁ (2¹Ag), ~1.8 eV for this conjugation
  length, below JF525's S₁ (~2.3 eV, from 525 / 549 nm). Through the same bonds that carry TTET, JF525's singlet could
  drop into it, so this version would probably be dim. In photosynthesis the gap is small because chlorophyll's S₁ is
  low (~1.85 eV), and even there carotenoids can quench chlorophyll singlets; adding one conjugated C=C is enough to
  turn a carotenoid into a singlet quencher (Berera et al., *PNAS* 103, 5343, 2006). Shortening the polyene raises S₁
  (and T₁), so a shorter carotenoid is the version that could satisfy both conditions for JF525.

What the rigid rod changes, as estimates rather than results:

- **TTET goes through bonds, not by contact.** The rod holds the COT 1.5 Å (n = 0) to 18.7 Å (n = 4) from C5, so the
  two can't touch. Through-bond triplet transfer falls off about e^(−0.45 Å⁻¹ · r) along conjugated bridges (Eng et al.,
  *J. Phys. Chem. B* 110, 6483, 2006), so each phenylene (4.3 Å) costs ~7×: relative to n = 1, n = 2 is ~0.14× and
  n = 3 ~0.02×. The flexible Cy5–COT reaches ~10⁶ s⁻¹; a rod has to beat that to be worth it, so n = 0 or 1 is the
  range to test.
- **The twists decouple.** In the force-field geometry each biaryl twists ~55° and the COT stands ~80° to the ring it
  sits on, which keeps the dye's spectrum close to Cy5's but also weakens the coupling TTET needs.
- **COT is an electron acceptor.** Calculations on COT bonded straight to dyes put a dye⁺/COT⁻ charge-transfer state
  below the dye's S₁ (Chanmungkalakul et al., *Chin. Chem. Lett.* 35, 109227, 2024), which would quench fluorescence
  and even speed up triplet formation. A phenylene spacer weakens this as well as TTET, so n = 0 vs n = 1 is the real
  trade-off.
- **Conjugation shifts the levels.** Fusing or conjugating rings onto COT raises its T₁ (El Bakouri et al., *Chem.
  Sci.* 2025), and an aryl on C5 red-shifts Cy5 somewhat. Both narrow the gaps the design relies on; Cy5's T₁ (1.60 eV)
  still clears COT's (0.95–1.1 eV) with room to spare.

### Interesting photophysics

Systems that each teach one piece of photophysics. Every entry has an **About this molecule** panel (under the
examples list, collapsible like the other sections): what it is, the photophysics, why it matters for imaging, what
the model shows, and sources. The text is the entry's `info` in `tools/catalog_photophysics.py`, written in a small
Markdown subset (`### ` headings, `- ` lists with indented continuation lines, `**bold**`, `*italic*`, a small
"Sources:" paragraph). Any catalog entry can carry `info`.

- **Energy transfer**: *DNA FRET ruler: Cy3–Cy5* (parametric, 4–30 base pairs; the panel note gives the distance,
  the angle between the dyes' transition dipoles, κ², and FRET as drawn vs for tumbling dyes) and *Aequorin*
  (calcium-triggered chemiluminescence, BRET to GFP; 1EJ3).
- **Light harvesting**: *LH2* (the full 9-fold ring from the PDB's biological assembly: B800 and B850 rings of
  bacteriochlorophyll, 18 carotenoids; 1NKZ), *PCP* (peridinin–chlorophyll protein, carotenoids as the harvesters;
  1PPR) and *FMO* (an energy funnel of 7 BChl a at 1.3 Å; 3EOJ).
- **Triplets put to work**: *LOV2* (the flavin triplet bonds to a cysteine; 2V0U), *miniSOG* (no cysteine, so the
  triplet makes singlet oxygen; 6GPU) and *KillerRed* (a fluorescent protein with a channel to its chromophore that
  makes superoxide; 2WIQ).
- **Fluorogenicity and sensing**: *JF525 open* and *closed* (the zwitterion and its colourless spirolactone) and
  *GCaMP6m* (a calcium sensor that works by switching the chromophore's proton; 3WLD).

Protein entries are built by `tools/build_cofactors.py` from a catalog `cofactors` dict: the protein drawn faint,
each photoactive cofactor a site with its H restored from the PDB Chemical Component Dictionary (H the dictionary
marks as leaving when a residue joins the chain are skipped, and a chromophore's peptide bonds aim its H). Sites
can bring their own glow colour (`color`); new kinds `carotenoid` (orange, magenta triplet) and `flavin`. The
Self-healing dye simulation uses a `carotenoid` site as the quencher when there is no `quencher`.

The FRET ruler comes from `tools/build_fret.py`: the A5·T20 base pair of the Dickerson dodecamer (1BNA) and the
helical screw that carries it onto the next pair (37.4°, 3.13 Å; joined O3′–P 1.75–1.85 Å vs 1.60 real), stacked by
the `helix` builder in `src/parametric.js`. Cy3 and Cy5 (non-sulfonated, tethers not drawn) are stacked 3.4 Å
beyond the end pairs, each turned and slid (≤ 2 Å) to fit the DNA it meets (closest heavy atoms 2.5 and 3.1 Å),
then locked to its end pair. R₀ = 5.4 nm is assumed.

### Molecules from papers

The **Review: TADF sensitizers** and **Review: UC emitters** tabs collect the molecules from Uji & Yanai,
*J. Photochem. Photobiol. C* **66** (2026) 100741 (`tools/catalog_uji2026.py`). Each entry records where its
structure came from (`source`: PubChem, the original paper's SI, or the review's drawing), the review's
reference number and the original DOI. The geometry report shows each molecule's formula against the reference
formula. Only PubChem and SI (HRMS) formulas are independent checks; for drawings the formula is a hand count.

Workflow for a new paper: list the molecules → look up PubChem → for new compounds, read the SI (synthesis
section + HRMS) or the figures → add entries to a `catalog_<paper>.py` → `npm run examples` → check the report.

The icon in `assets/` comes from `tools/make_icon.py` (`npm run icon`).

## Roadmap

- v2: atom labels, screenshot export, background options; optional GFN2-xTB refinement of the dye geometries;
  animated energy transfer along the TTA chain
- v3: PDB / MOL / SDF formats (with explicit bond orders), double/triple bond rendering, ambient occlusion & outlines
