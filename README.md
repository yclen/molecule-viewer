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
| Excite / un-excite | the dot next to a molecule under **Excited states**, or `E` for the selected one |

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

Each one has its own excited state E (`src/excited-states.js`): **on/off**, **colour**, **brightness** and
**spread** (Gaussian width in Å). Click a molecule in the list (or one of its atoms in the 3D view) to select
it; the editor below the list changes only that one, and every other molecule keeps what was set for it. The
dot beside each name switches it on and off, `E` toggles the selected one, and **All on / All off** switch the
whole object. States are remembered per example and per name, so changing n keeps DPA 1's state, and
switching to another example and back restores them. Default colours: sensitizer orange, DPA blue, emitter
cyan. In dev builds `__excite.set('DPA 2', { on: true }); __updateGlow()` drives it from the console, the hook
for a later time-based simulation.

- **Where it glows**: each atom's share of the lowest π→π* excitation, (c²_HOMO + c²_LUMO) / 2, from Hückel
  theory on the 3D structure (`src/excitation.js`). Resonance integrals are scaled by p-orbital overlap
  (|n_i · n_j|), so twisted rings drop out by themselves: DPA's 90° phenyls stay dark and the glow sits on the
  anthracene (C9/C10 brightest, fused carbons zero), and fluorescein's benzoate stays dark while the xanthene
  glows. Heteroatoms use standard Hückel parameters. It is a qualitative picture, not a TD-DFT density.
- **How it's drawn** (`src/glow.js`): every glowing atom is a camera-facing sprite with a Gaussian falloff,
  added up in an HDR buffer (the line-of-sight integral of a 3D Gaussian cloud). Atoms and bonds in front hide
  it, a quarter-resolution blur adds bloom to the glow only, and a `1 − exp(−I)` curve keeps faint regions
  coloured and saturates dense ones toward white.

## Measuring

Click 2–4 atoms: the panel shows each distance, the angle(s) and the dihedral, and the 3D view labels each
segment. Click an atom again to unpick it, click empty space or press `Esc` to clear. Measurements update live
while a trajectory plays.

## Examples

Tabs (Basics, Triplet molecules, Xanthenes, Cyanines, Other dyes, Proteins, TTA upconversion, TADF, and two review tabs), all bundled into the app.
`↑` / `↓` steps through the current tab.

- **Basics** (`samples/basics/`): idealised geometries from `tools/make_samples.py` (`npm run samples`).
- **Everything else**: listed in `tools/catalog.py` and built by `tools/build_examples.py`
  (`npm run examples`, needs `pip install rdkit`):
  - *Small molecules*: structure (SMILES) from PubChem → RDKit ETKDGv3 conformer search → MMFF94s
    minimisation → lowest-energy conformer. Counter-ions are dropped. Bonds that MMFF gets wrong (C–I) or
    can't type (boron, so BODIPY uses UFF) are restrained to X-ray means. The PtOEP porphyrin core is
    placed from metalloporphyrin crystal averages, because MMFF puckers porphyrins.
  - *Proteins*: X-ray structures from the RCSB PDB (chain A, chromophore kept, solvent removed).
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
  per-state glow colours (singlet / triplet) and animated energy transfer along the TTA chain
- v3: PDB / MOL / SDF formats (with explicit bond orders), double/triple bond rendering, ambient occlusion & outlines
