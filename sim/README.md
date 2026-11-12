# TTA Flora 4 simulation

A standalone, data-only version of the TTA cascade. The app doesn't use anything in this folder.

The system is a string of `n + 2` letters, one per site: the sensitizer, DPA 1 to `n`, then the perylene. Each letter
is `G`, `T` or `S`. Rules rewrite one or two sites at a rate per ns, and there is no time step: the engine adds up
the rates of every rule that can fire, waits a random time, and fires one (the Gillespie method).

Needs Python with `numpy`, `scipy` and `numba`. The first run compiles for a few seconds. The easiest way to get
them, and to open the two windows, is the launcher in the project folder (`python launch.py`, or double-click
`Launch.bat` / `Launch.command`): see the main [README](../README.md#launch-it).

| File | What it is |
|---|---|
| `rules.py` | Default parameters, site distances, FRET and every rule. **The file to edit.** |
| `engine.py` | Monte Carlo, the exact solvers (long-run and pulse) and the metrics. |
| `plot_pulse.py` | The figures for pulse runs. |
| `plot_yield.py` | The quantum yield figure for sweeps. |
| `plot_intensity.py` | The figures for intensity sweeps: photons per molecule and quantum yield against `kExc`. |
| `tta_sim.py` | The commands and data saving. |
| `live.py` | The live histogram window (see "Live window"). Needs `PyQt6` and `pyqtgraph`. |
| `sweep_gui.py` | The window for running and plotting intensity sweeps (see "Sweep window"). Same needs, plus `matplotlib`. |
| `formulas.md` | The model as math: the state vector, the 10 rules, the master equation and how each method uses it. |

## Commands

Run these from the project folder:

```bash
python sim/tta_sim.py run --n 8 --excRatio 0.5
```

```bash
python sim/tta_sim.py sweep --excRatio 0.01:10:log7 --hopRatio 10,100,1000 --name ratios
```

```bash
python sim/tta_sim.py exact --n 4
```

```bash
python sim/tta_sim.py check
```

```bash
python sim/tta_sim.py pulse --excRatio 1 --hopRatio 100 --duration 5 --name first-pulse
```

```bash
python sim/tta_sim.py intensity --kExc 0.001:100:log31 --duration 5 --checks 3 --name qy-vs-kexc
```

```bash
python sim/tta_sim.py rebin 2026-11-06_first-pulse --bin 0.01
```

```bash
python sim/tta_sim.py plot 2026-11-06_first-pulse
```

- `run`: one setting, with a summary in the terminal.
- `sweep`: every combination of the values given.
- `exact`: the long-run answer with no random noise, up to `n = 8` (about 2 minutes there, seconds at `n` of 6 or less).
- `pulse`: flashes many molecules and records when every photon is emitted (see "Pulse runs").
- `intensity`: total photons and excitations per molecule over one pulse, against `kExc` (see "Intensity sweeps").
- `rebin`: builds a new histogram from a pulse run's saved photons, with no new simulation.
- `plot`: draws a run into a `plots/` folder. Needs `matplotlib`. A pulse run gets its emission histograms. A sweep of
  `excRatio` gets `quantum_yield.png`, with one line per value if it also varies one other parameter (up to 5 values).
  Points with under 30 photons are left out. An intensity run gets its four figures (see "Intensity sweeps").
- `check`: compares Monte Carlo with the exact solve and checks the bookkeeping.

Every parameter is an option, and any of them can be swept:

| Written as | Means |
|---|---|
| `5` | one value |
| `1,2,5` | a list |
| `0:12` | whole numbers from 0 to 12 |
| `1:300:6` | 6 evenly spaced values |
| `1:300:log6` | 6 values spaced by ratio |

Also: `--excitations` (per point, default 200,000), `--replicas` (default 64), `--seed` and `--name`. With 200,000
excitations, photons per excitation is good to about ±0.0005.

Pulse options: `--duration`, `--molecules` (default 100,000) and `--bin` (default `tauS / 10`). `pulse` also accepts sweeps,
for example `--excRatio 0.1,1,10 --duration 1:10:4`. `rebin` takes `--bin`, or `--logbins 60` for bins spaced by ratio.

## Parameters

Times are in ns and rates per ns, everywhere: the commands, the saved files and the live window. Only the wall-clock
figures (`wall_s`, `events_per_wall_s`, the `seconds` column of `index.csv`) are real seconds. On the command line,
hopping and excitation can be set relative to the singlet lifetime, or directly with `--kHop` and `--kExc`:

| Name | Default | Meaning |
|---|---|---|
| `n` | 8 | DPA mediators |
| `tauS` | 5 ns | singlet lifetime on a DPA or the perylene |
| `tauSens` | 5 ns | sensitizer singlet decay lifetime (rule 9), not counting ISC (`inf` = never decays) |
| `tripletYield` | 0.5 | share of sensitizer singlets that become triplets; sets `kISC` |
| `hopRatio` | 100 | how many times faster hopping is than the singlet lifetime: `kHop = hopRatio / tauS` |
| `excRatio` | 0.1 (1 in a pulse) | excitation rate relative to the sensitizer's lifetime: `kExc = excRatio / tauSens` |
| `kDecay` | 1e-6 /ns | triplet decay on a DPA or the perylene |
| `f` | 0.4 | the F factor: the share of TTA events that make a singlet |
| `fretRadius` | 4 nm | distance where back transfer is 50% (0 = no back transfer) |
| `spacing` | 0 | 0 = distances measured from the real molecule; above 0 = site `i` sits at `i × spacing` nm |
| `block` | 1 | 1 = a triplet can't hop onto a singlet; 0 = it lands and quenches it (for comparison only) |
| `duration` | 5 ns | pulse runs only: how long `kExc` stays on |

`kISC`, `kHop` and `kExc` are worked out from these and saved with the results. ISC and rule 9 compete, so
`tripletYield = kISC / (kISC + 1/tauSens)`, which gives `kISC = (1/tauSens) × tripletYield / (1 − tripletYield)`. At 50%
that makes `kISC` equal to `1/tauSens`, and the sensitizer's singlet actually lasts `tauSens / 2` on average.

Giving `--kISC`, `--kHop` or `--kExc` sets that rate directly, and the yield or ratio follows from it.

## Rules

Two-letter rules read as hopping site, then receiving site. `formulas.md` has the same rules as equations, and the
app's TTA cascade script (`samples/sims/3-tta-cascade.js`) follows the same rules on one molecule.

| # | Rule | Where | Change | Rate | Counted as |
|---|---|---|---|---|---|
| 1 | Excitation | sensitizer | `G → S` | `kExc` | excitation |
| 2 | Intersystem crossing | sensitizer | `S → T` | `kISC` | |
| 3 | Triplet decay | DPA, perylene | `T → G` | `kDecay` | triplet decayed |
| 4 | Singlet emission | DPA, perylene | `S → G` | `(1 − E)/tauS` with the sensitizer in `G`, otherwise `1/tauS` | photon |
| 5 | Hop | neighbouring sites | `T G → G T` | hop rate | |
| 6 | Blocked hop | neighbouring sites | `T S`: nothing happens | none | |
| 7 | TTA makes a singlet | neighbouring sites | `T T → G S` | `f` × hop rate | TTA singlet |
| 8 | TTA fails | neighbouring sites | `T T → G T` | `(1 − f)` × hop rate | TTA failed |
| 9 | Sensitizer singlet decay | sensitizer | `S → G` | `1/tauSens` | sensitizer singlet lost |
| 10 | Back transfer | sensitizer in `G` + DPA or perylene in `S` | sensitizer `G → S`, site `S → G` | `E/tauS` | back transfer |

**Hop rate:** the sensitizer hops to DPA 1 (or the perylene if `n = 0`) at `kHop`. DPA 1 only goes forward, at `kHop`.
Other DPAs go each way at `kHop/2`, and for the last DPA forward means the perylene. Nothing hops onto the sensitizer
or off the perylene.

**FRET:** `E = 1 / (1 + (r / fretRadius)⁶)`, worked out for each site from `r`, its distance to the sensitizer.

**Distances:** measured from the geometry in `samples/triplet/`, between the sensitizer's core atoms, the centre of each
DPA and the perylene's core atoms. These are the centres the app's "core centre to centre" note uses. DPA 1 is 1.40 nm
from the sensitizer and each further DPA adds 1.32 nm.

## Data

Each run gets its own folder in `sim/data/` (not tracked by git), and `sim/data/index.csv` lists them all.

Runs saved before the switch to ns (up to 2026-11-10) have the same numbers under `_s` names, such as `bin_start_s`.
`plot` and `rebin` still read them, and a rebinned file gets the new `_ns` names.

```text
sim/data/
  index.csv                   one line per run
  2026-11-06_wire-length/
    results.csv               one row per setting: parameters, then metrics
    sites.csv                 one row per site per setting
    params.json               the command, defaults, seed and what was swept
    rules.py                  the exact rules this run used
```

`point` links the rows of `results.csv` and `sites.csv`. Counts are **per excitation** unless the name says otherwise.

| `results.csv` column | Meaning |
|---|---|
| `photons_per_excitation`, `_se` | total photons, and its standard error |
| `perylene_photons`, `dpa_photons` | the same, split by where they came from |
| `back_transfers` | singlets taken back by the sensitizer |
| `isc_yield` | share of sensitizer singlets that became triplets |
| `sens_singlets_lost` | sensitizer singlets lost to rule 9 |
| `triplets_made`, `triplets_decayed` | triplets from ISC, and triplets lost to rule 3 |
| `tta_singlet_perylene`, `tta_singlet_dpa` | TTA events that made a singlet |
| `tta_failed_perylene`, `tta_failed_dpa` | TTA events that didn't |
| `excitations_per_ns`, `photons_per_ns` | rates per simulated ns |
| `mean_triplets_present`, `mean_singlets_present` | average number on the molecule at once |
| `time_to_photon_ns` | average time from excitation to photon. After TTA the clock follows the younger triplet |
| `triplet_balance`, `singlet_balance` | bookkeeping checks, should be 1 |
| `wall_s` | real seconds this point took to simulate |
| `events`, `events_per_wall_s` | rules fired after warmup, and how many per real second |
| `sim_ns` | simulated ns covered, added up across replicas |

The whole run's time is printed at the end and saved in the `seconds` column of `index.csv`.

Excitations mean rule 1 only: a back transfer doesn't count as a new one. `isc_yield` and `triplets_made` do include
singlets recycled through back transfer, so `triplets_made` can be above `isc_yield`. `exact` leaves `_se` and
`time_to_photon_ns` empty.

## Pulse runs

Every molecule starts in all `G`, and t = 0 is the start of the pulse. `kExc` is on until `duration`, then 0. Molecules
are independent, and each is followed until it's back in all `G`. Each molecule gets one pulse.

```text
sim/data/2026-11-06_first-pulse/
  histogram.csv               one row per time bin
  summary.csv                 one row per setting
  sites.csv                   one row per site per setting
  photons/point_0000.npz      every photon of setting 0
  params.json
  rules.py
```

`plot` adds `plots/point_0000.png` for each setting, and `plots/comparison.png` with shared axes when a run has 2 to 6
settings. Every figure shows `kHop` and the quantum yield, which is `photons_per_excitation`: photons emitted per photon
absorbed. Absorbed means rule 1 only, so a back transfer doesn't count as a new absorption.

The time axis ends once 99% of a setting's photons have arrived, and a comparison uses the middle setting's limit so one
long tail doesn't flatten every panel. A panel says so if about 1% or more of its photons land past the axis. `--xmax 20` sets
the end yourself.

Each setting has its own photon file, so a stopped sweep keeps the settings it finished. `rebin` writes a new file such
as `histogram_bin0.01.csv` and leaves `histogram.csv` as it is.

| `histogram.csv` column | Meaning |
|---|---|
| `point` | which setting |
| `bin_start_ns`, `bin_end_ns` | the bin, in ns after the pulse starts |
| `pulse_on` | 1 if the bin starts during the pulse |
| `photons`, `perylene_photons`, `dpa_photons` | counts in the bin |
| `photons_per_molecule_per_ns` | the count scaled so it doesn't depend on bin width or molecule count |
| `error` | the counting error on that scaled value, from the square root of the count |

| `summary.csv` column | Meaning |
|---|---|
| `excitations_per_molecule`, `photons_per_molecule` | averages over all molecules |
| `share_of_molecules_with_a_photon`, `share_of_photons_after_pulse` | fractions from 0 to 1 |
| `peak_bin_start_ns` | start of the busiest bin |
| `mean_emission_ns`, `median_emission_ns`, `last_photon_ns` | emission times after the pulse starts |
| `time_to_photon_ns` | average time from the younger triplet's excitation to the photon |
| `wall_s`, `events`, `events_per_wall_s` | real seconds, rules fired, and rules fired per real second |

`summary.csv` also has the parameters and the per-excitation columns from `results.csv`, such as
`photons_per_excitation`, `back_transfers` and `isc_yield`. A pulse ends with every molecule in all `G`, so here
`triplet_balance` and `singlet_balance` are exactly 1.

Each photon file holds `time`, `site`, `molecule` and `born` (when the younger triplet was excited), one entry per photon
sorted by time, plus `excitations`, one entry per molecule. Read one with:

```python
import numpy as np
photons = np.load('sim/data/2026-11-06_first-pulse/photons/point_0000.npz')
print(photons['time'][:10])
```

## Intensity sweeps

How many photons a molecule gives over one pulse, and how many times it absorbs, as the excitation rate goes up. This
is the data for photons against laser power and quantum yield against power: `kExc` is proportional to power, so the
conversion (an absorption cross-section) can be applied afterwards.

```bash
python sim/tta_sim.py intensity --kExc 0.001:100:log31 --duration 5,50 --checks 3 --name qy-vs-kexc
```

Every molecule starts in all `G`, gets one pulse of `duration` ns and is followed until it can't emit any more, the
same as `pulse`. Any parameter can be swept alongside `kExc`. Each extra value gives another curve.

- `--method exact` (the default): `engine.pulse_totals`, with no random noise at any `kExc`, so it holds up at very low
  power, where Monte Carlo would need enormous numbers of molecules (a photon needs two absorptions in one pulse, so
  photons fall as `kExc²`). Up to `n = 8`. At `n = 8` the first point takes about 20 s and each one after about 5 s;
  at `n` of 6 or less they take under a second. Long, bright pulses take longer. `formulas.md` has how it works.
- `--checks 3`: with the exact solve, also runs Monte Carlo at 3 points spread over the sweep, as a check.
- `--method mc`: Monte Carlo at every point, for `n` above 8. `--method both` runs both at every point.
- Monte Carlo adds batches of `--batch` molecules (default 100,000) until photons per molecule is known to
  `--precision` (default 0.01, a 1% standard error), or `--max-seconds` (default 600) runs out for that point.

It saves no photon times, only totals, so a run is a few kB. Each row is written as soon as its point is done, so a
stopped sweep keeps what it finished.

```text
sim/data/2026-11-10_qy-vs-kexc/
  results.csv                 one row per point per method
  params.json
  rules.py
  plots/                      after `plot`: photons_linear.png, photons_log.png, quantum_yield.png,
                              quantum_yield_log_kexc.png
```

Every count is **per molecule, per pulse**. Quantum yield is `photons_per_molecule / excitations_per_molecule`.

| `results.csv` column | Meaning |
|---|---|
| `point`, `method` | which setting, and `exact` or `mc` |
| the parameters | `n`, `kHop`, `kExc`, `duration` and the rest, as in `results.csv` of a sweep |
| `photons_per_molecule` | photons, the sum of the next two |
| `perylene_photons_per_molecule`, `dpa_photons_per_molecule` | where they came from |
| `excitations_per_molecule` | rule 1 events. Below `kExc × duration` once the sensitizer is often busy |
| `isc_per_molecule`, `sens_singlets_lost_per_molecule` | what sensitizer singlets did: ISC, or rule 9 |
| `back_transfers_per_molecule` | singlets taken back by the sensitizer |
| `tta_singlet_per_molecule`, `tta_failed_per_molecule` | TTA events that made a singlet, and ones that didn't |
| `triplets_decayed_per_molecule` | triplets lost to rule 3, including the one left over at the end, if any |
| `molecules` | Monte Carlo only: how many were flashed |
| `mean_photons_sq`, `mean_excitations_sq`, `mean_photons_x_excitations` | Monte Carlo only: averages over molecules, for error bars on any ratio |
| `photons_relative_error` | Monte Carlo only: the standard error of photons per molecule, as a share |
| `wall_s` | real seconds the point took |

The standard error of photons per molecule is `sqrt((mean_photons_sq − photons²) / molecules)`. `plot_intensity.py`
has the error on quantum yield, which also uses the excitation columns.

`plot` on an intensity run saves four figures. Several runs can be drawn together from the sweep window.

- **Photons per molecule, linear:** the total, plus perylene and DPA when there's one curve.
- **Photons per molecule, log, with slope:** log on both axes, and under it the local slope
  `d log(photons) / d log(kExc)`. It's 2 while a photon needs two absorptions in one pulse and falls to 0 at saturation,
  so the stages of saturation show as where it leaves 2 and where it reaches 0.
- **Quantum yield, linear,** in percent, with the limit `100 × tripletYield × f / (1 + f)`.
- **Quantum yield, log kExc:** the same, with `kExc` on a log scale and the yield still linear, so the rise at low
  power isn't squeezed against the axis.

Exact rows are lines, Monte Carlo rows open markers with error bars. When every curve has the same pulse length, a
second axis along the top shows `kExc × duration`, the absorptions per pulse if the sensitizer were always free.

**Laser power axis.** `kExc` is proportional to power, so any figure can have power in mW along the bottom instead,
with `kExc` along the top. `--kexc-per-mw 0.0389` does it for `plot`, saving files that end in `_power.png` beside
the `kExc` ones. The sweep window has a checkbox for it, and its default is the live window's `kExc per mW`, worked out
from the beam settings at the top of `live.py`. At 1 µm 1/e² radius, 635 nm and ε = 50,000 /M/cm (σ = 1.91 × 10⁻¹⁶ cm²)
that's 0.0389 /ns per mW, at the beam's peak intensity 2P/(πw²). Averaged over the 1/e² spot it would be half that.

## Sweep window

```bash
python sim/sweep_gui.py
```

**Run** has a box for every parameter, in the same units and with the same names as the live window. Every box takes
what the command line takes, so `5,50` in pulse length gives two curves. Set the `kExc` range, its number of points and
spacing (by ratio for the log plot, even for the linear ones), pick the method and press RUN SWEEP. It runs the
`intensity` command in its own process, shown underneath so you can copy it. Photons per molecule appear on a log plot
point by point, with progress and time left. STOP keeps the points finished so far. A box with a value out of range
turns red, and the exact solve isn't offered above `n = 8`.

**Plot** lists every intensity run in `sim/data`, newest first, with the `n` and pulse lengths it used. Select one, or
several with Ctrl or Shift, and open any of the four figures in its own matplotlib window, with zoom, pan and save.
Hovering shows the nearest point's values. "Log x axis on the linear plots" is there for sweeps spaced by ratio.
"Laser power on the x axis" switches every figure opened after it to power in mW, using the `kExc per mW` box
under it (see "Laser power axis" above), and "Save PNGs" then writes the `_power.png` files.
"Save PNGs" writes each selected run's `plots/` folder. A run that finishes opens here on its own.

## Live window

```bash
python sim/live.py
```

A pulse run you can watch. Fresh molecules are flashed continuously, and each photon lands in a histogram as it's
counted, so it builds up like a time-resolved measurement. It's the same simulation as `pulse`: every molecule starts in
all `G`, gets one pulse and is followed until it's back in all `G`, so nothing carries over from one pulse to the next.

**Everything in the window is in nanoseconds**, like the rest of the simulation. Only the wall-clock figures (molecules
per s, how long a solve takes) are real seconds. There are no ratios: `kHop` and `kExc` are both set directly.

The controls are laid out like an instrument, with the ones you change most at the top:

- **Laser:** START / STOP (or F5). Stop pauses and keeps the counts, and the count rate and molecules per s
  stay at their last values, marked "held". The slider sets **`kExc`** directly, in /ns, on a
  log scale from 0.0001 to 100 /ns. `kExc` is proportional to laser power, so this is the power knob. There are ÷10, ÷2,
  ×2 and ×10 buttons, or type a value. Dragging the slider only clears the histogram when you let go. **Pulse length**
  is here too.

  Under the slider is the **laser power** that gives this `kExc` at the centre of the beam, `P = kExc / (kExc per mW)`,
  with `kExc per mW = σ I0 / hν` for 1 mW: the absorption cross-section `σ = ln(10) · 1000 · ε / N_A` and the peak
  intensity of a Gaussian beam `I0 = 2P / (π w0²)`. The assumed setup is listed under it, and set at the top of
  `live.py`: beam radius `w0 = 1 µm` (1/e² of the intensity), `ε = 50,000 /M/cm` and 635 nm, which gives 0.0389 /ns per
  mW. So the default `kExc = 0.2 /ns` is about 5.1 mW. The power is only shown: the simulation runs on `kExc`.
- **Acquisition:** time window, bin width, throughput, log counts, the exact overlay, Clear and Save (or Ctrl+S).
- **Sample**, folded away with a one-line summary. **TTA | FLUORESCENCE** under its title picks what's in the beam, and
  stays in view when the panel is folded, so you can flip between the two traces with one click:
  - **TTA**: `n`, `kHop`, `f`, triplet yield, `τS`, `τSens`, `kDecay` and the FRET radius. Below them are the rates that
    follow: `kISC` and the hop time.
  - **Fluorescence**: a plain two-state fluorophore, `G` and `S`, on one site, to compare with (`fluorescence_rules` in
    `rules.py`). It's excited at `kExc` while the pulse is on and goes back to `G` at `1/τF`, giving a photon with chance
    `Φ_F` (default 1, so its quantum yield is `Φ_F`). Its trace is violet. Below the values is the saturation,
    `kExc / (kExc + 1/τF)`: the share a steady light would keep in `S`.

  Switching sample keeps the laser, pulse and time window, and starts a new histogram.
- **Readouts** along the top: the quantum yield with its counting error (and the exact solve's yield next to it when
  the overlay is on), the **count rate** in photons per real second, photons counted, molecules flashed with how many per
  second, and photons absorbed per molecule. The header has the acquisition clock (real time spent running since the
  last clear), the laser lamp and the exact solve's lamp.

The defaults are the old pulse defaults in flat units: `τS = τSens = 5 ns`, `kHop = 20 /ns`, `kExc = 0.2 /ns`, a 5 ns
pulse and a 50 ns window.

- **Parameters:** type a value and press Enter or Tab. Any change to the laser, pulse or sample clears the histogram, so
  it never mixes two experiments. Bin width, the log axis and throughput don't.
- **Time window:** photons arriving after it aren't in the histogram, but they do count towards the total and the quantum
  yield. Once 1% or more of the photons land after it, the Photons readout says so in amber. Make the window longer then.
- **Throughput:** molecules to aim for per real second, default 100,000. 0 means as fast as the machine can go. The rate
  actually reached is under "Molecules", averaged over 2 s.
- **Bin width, log axis:** these change only the view, so the counts stay. Bins are multiples of `τS / 100`.
- **Exact overlay:** works out the histogram the counts are heading for, with no random noise, while they pile up. It
  runs in its own process, with a progress bar and a time left, and when it's done it's drawn over the counts as dashed
  lines, scaled to the molecules counted so far. Its quantum yield, photons after the window included, appears under
  the counted one.
  Below the checkbox is how long it should take at the current values. Above `n = 9`, or when it would take over 10
  minutes, it isn't offered and says why. Changing a value starts a new solve; Clear, the bin width and the log axis
  don't. `formulas.md` has how it works. For the fluorophore it's instant, from the closed form
  (`FluorescenceExact` in `engine.py`).
- **Hover** over the histogram to read a bin, with what the exact solve expects there when it's shown. Drag and scroll to
  look along the time axis. Clicking a legend entry hides that curve, and clicking it again brings it back.
- **↕ AUTO Y** (on by default) keeps the counts axis fitted to every curve that's shown in the time range in view,
  the exact overlay included and hidden curves left out. Turn it off to zoom and pan the counts axis yourself (scroll
  over the axis to zoom only that axis). **⌂ HOME** (or the Home key) goes back to the whole time window with AUTO Y on.
- The exact overlay grows as the counts do, because it's drawn as the photons expected from the molecules counted so
  far. Its shape is fixed once the solve is done.
- **Save** writes `histogram.csv`, `summary.csv`, `params.json`, `rules.py` and `histogram.png` to a new
  `sim/data/YYYY-MM-DD_live/` folder. The columns are a pulse run's, in ns: `bin_start_ns`, `bin_end_ns` and
  `photons_per_molecule_per_ns`. It holds
  the histogram at the current bin width, not each photon, so `rebin` and `plot` don't work on it. With the exact solve
  shown, `histogram.csv` also has `exact_perylene_photons` and `exact_dpa_photons` (the photons it expects in each bin
  from the molecules counted) and `summary.csv` has `exact_photons_per_excitation`. A fluorescence run goes to a
  `YYYY-MM-DD_live-fluorescence/` folder, with `fluorescence_photons` (and `exact_fluorescence_photons`) in place of the
  perylene and DPA columns, and `sample` in its parameters says which it was.

The ± on the quantum yield is the counting error: the square root of the photons, over the excitations.

The exact solve has a row for every one of the `3^(n + 2)` strings, so `n` sets how long it takes. At the defaults on the
development machine, with counting running alongside: a few seconds up to `n = 7`, about 45 s at `n = 8` and several
minutes at `n = 9`. Hopping speed barely matters. A pulse longer than the time window gives no exact yield, because
the solve stops at the window.

Speed depends on how many excitations each molecule gets. With 5 ns pulses on the development machine, 100,000 molecules
per second held at `kExc` 0.2 and 2 /ns, and flat out it reached about 160,000 at 2 /ns. Long, bright pulses are far
slower: a 500 ns pulse at `kExc` 6 /ns runs about 1,200 molecules per second, though each gives about 25 photons.
