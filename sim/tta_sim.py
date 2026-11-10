"""Command line for the TTA Flora 4 simulation. See README.md for the commands and what's in the data folder.

Simulated times are in ns and rates per ns. Wall-clock figures (wall_s, seconds in index.csv) are real seconds.

  python sim/tta_sim.py run --n 8 --excRatio 0.5
  python sim/tta_sim.py sweep --excRatio 0.01:10:log7 --hopRatio 10,100,1000 --name ratios
  python sim/tta_sim.py exact --n 4
  python sim/tta_sim.py pulse --excRatio 1 --hopRatio 100 --duration 5 --name first-pulse
  python sim/tta_sim.py rebin 2026-11-06_first-pulse --bin 0.01
  python sim/tta_sim.py intensity --kExc 0.001:100:log31 --duration 5 --name qy-vs-kexc
  python sim/tta_sim.py plot 2026-11-06_first-pulse
  python sim/tta_sim.py check
"""

import argparse
import csv
import itertools
import json
import shutil
import sys
import time
from datetime import datetime
from pathlib import Path

import numba as nb
import numpy as np

from engine import EVENTS, exact, monte_carlo, per_site, pulse, pulse_totals, summarize
from rules import DEFAULTS, PULSE_DEFAULTS, build_rules, fret_efficiency, resolve, site_distances

HERE = Path(__file__).resolve().parent
DATA = HERE / 'data'
EXACT_MAX_N = 8  # 3^10 strings. Takes a couple of minutes
BINS_PER_LIFETIME = 10  # default histogram bin: a tenth of the singlet lifetime
UNITS = 'times in ns, rates per ns'  # saved in params.json


def parse_values(name, text):
    """'5' | '1,2,5' | '0:12' (whole numbers) | '1:300:6' (6 evenly spaced) | '1:300:log6' (6 spaced by ratio)."""
    if ',' in text:
        values = [float(v) for v in text.split(',')]
    elif ':' in text:
        parts = text.split(':')
        lo, hi = float(parts[0]), float(parts[1])
        if len(parts) == 2:
            values = list(np.arange(lo, hi + 1))
        elif parts[2].startswith('log'):
            values = list(np.geomspace(lo, hi, int(parts[2][3:])))
        else:
            values = list(np.linspace(lo, hi, int(parts[2])))
    else:
        values = [float(text)]
    return [int(v) if name in ('n', 'block') else float(v) for v in values]


def grid(args, defaults=DEFAULTS):
    """Every combination of the values given on the command line, as (parameters, names that vary)."""
    choices = {name: parse_values(name, getattr(args, name)) if getattr(args, name) is not None else [value]
               for name, value in defaults.items()}
    swept = [name for name, values in choices.items() if len(values) > 1]
    return [resolve(dict(zip(choices, combo))) for combo in itertools.product(*choices.values())], swept


def site_names(n):
    return ['sensitizer'] + [f'DPA {i}' for i in range(1, n + 1)] + ['perylene']


def new_run_folder(mode, name, swept):
    label = name or '_'.join([mode, *swept])
    stem = f'{datetime.now():%Y-%m-%d}_{label}'
    folder, k = DATA / stem, 1
    while folder.exists():
        k += 1
        folder = DATA / f'{stem}_{k}'
    folder.mkdir(parents=True)
    return folder


def add_to_index(folder, mode, points, swept, seconds):
    index = DATA / 'index.csv'
    fresh = not index.exists()
    with index.open('a', newline='') as out:
        writer = csv.writer(out)
        if fresh:
            writer.writerow(['folder', 'date', 'mode', 'points', 'swept', 'seconds', 'command'])
        writer.writerow([folder.name, f'{datetime.now():%Y-%m-%d %H:%M}', mode, points, ' '.join(swept),
                         round(seconds, 1), ' '.join(sys.argv[1:])])


def solve(mode, p, args, seed):
    rules = build_rules(p)
    sites = p['n'] + 2
    if mode == 'exact':
        return exact(rules, sites)
    return monte_carlo(rules, sites, args.excitations, args.replicas, seed)


def run_points(mode, args):
    points, swept = grid(args)
    if mode == 'exact' and max(p['n'] for p in points) > EXACT_MAX_N:
        sys.exit(f'exact only goes up to n = {EXACT_MAX_N}: it has 3^(n + 2) strings to solve')
    folder = new_run_folder(mode, args.name, swept)
    shutil.copy(HERE / 'rules.py', folder / 'rules.py')
    (folder / 'params.json').write_text(json.dumps({
        'mode': mode, 'command': ' '.join(sys.argv[1:]), 'date': f'{datetime.now():%Y-%m-%d %H:%M:%S}',
        'units': UNITS, 'defaults': DEFAULTS, 'swept': {name: sorted({p[name] for p in points}) for name in swept},
        'points': len(points), 'excitations_per_point': args.excitations, 'replicas': args.replicas, 'seed': args.seed,
        'threads': nb.get_num_threads(),
    }, indent=2))

    if mode != 'exact':  # load the compiled code now, so it isn't timed as part of the first point
        monte_carlo(build_rules({**DEFAULTS, 'n': 0}), 2, 2, 2, 1, 0)

    started = time.time()
    with (folder / 'results.csv').open('w', newline='') as results, (folder / 'sites.csv').open('w', newline='') as by_site:
        result_rows = site_rows = None
        for i, p in enumerate(points):
            began = time.perf_counter()
            tally = solve(mode, p, args, args.seed + i * args.replicas)
            wall = time.perf_counter() - began
            metrics = summarize(tally)
            events = tally.counts.sum() if mode != 'exact' else np.nan  # exact has rates, not counted events
            timing = {'wall_s': wall, 'events': events, 'events_per_wall_s': events / wall,
                      'sim_ns': tally.time if mode != 'exact' else np.nan}
            row = {'point': i, **p, **metrics, **timing}
            if result_rows is None:
                result_rows = csv.DictWriter(results, row.keys())
                result_rows.writeheader()
            result_rows.writerow(row)
            results.flush()  # a stopped sweep keeps what it finished

            r = np.r_[0.0, site_distances(p['n'], p['spacing'])]
            E = np.r_[0.0, fret_efficiency(r[1:], p['fretRadius'])]
            sites = per_site(tally)
            for s, name in enumerate(site_names(p['n'])):
                row = {'point': i, 'n': p['n'], 'site': s, 'name': name, 'r_nm': r[s], 'fret_efficiency': E[s], **sites[s]}
                if site_rows is None:
                    site_rows = csv.DictWriter(by_site, row.keys())
                    site_rows.writeheader()
                site_rows.writerow(row)
            if len(points) > 1 and (i + 1) % max(1, len(points) // 10) == 0:
                print(f'{i + 1}/{len(points)} points, {time.time() - started:.0f} s', flush=True)
    seconds = time.time() - started
    add_to_index(folder, mode, len(points), swept, seconds)
    if len(points) == 1:
        report(points[0], metrics, sites, r, E)
    print(f'\nSimulated {len(points)} point{"s" * (len(points) > 1)} in {seconds:.1f} s'
          f' on {nb.get_num_threads()} threads')
    print(f'Saved to {folder}')


def linear_edges(latest, width):
    return width * np.arange(int(latest / width) + 2)


def log_edges(times, bins):
    return np.geomspace(times[0], times[-1], bins + 1)


def histogram_rows(point, photons, perylene, molecules, duration, edges):
    """One dict per time bin. photons_per_molecule_per_ns doesn't depend on bin width or how many molecules ran."""
    total = np.histogram(photons['time'], edges)[0]
    from_perylene = np.histogram(photons['time'][photons['site'] == perylene], edges)[0]
    scale = 1 / (molecules * np.diff(edges))
    return [{
        'point': point, 'bin_start_ns': edges[k], 'bin_end_ns': edges[k + 1], 'pulse_on': int(edges[k] < duration),
        'photons': total[k], 'perylene_photons': from_perylene[k], 'dpa_photons': total[k] - from_perylene[k],
        'photons_per_molecule_per_ns': total[k] * scale[k], 'error': total[k] ** 0.5 * scale[k],
    } for k in range(len(total))]


def pulse_metrics(p, run, molecules, histogram):
    times = run.time
    m = summarize(run.tally)
    for name in ('excitations_per_ns', 'photons_per_ns', 'mean_triplets_present', 'mean_singlets_present',
                 'time_to_photon_ns'):
        del m[name]  # these describe a steady state
    has_photons = times.size > 0
    busiest = max(histogram, key=lambda row: row['photons'])
    return {
        'molecules': molecules,
        'excitations_per_molecule': run.excitations.sum() / molecules,
        'photons_per_molecule': times.size / molecules,
        **m,
        'share_of_molecules_with_a_photon': np.unique(run.molecule).size / molecules,
        'share_of_photons_after_pulse': (times >= p['duration']).mean() if has_photons else np.nan,
        'peak_bin_start_ns': busiest['bin_start_ns'] if has_photons else np.nan,
        'mean_emission_ns': times.mean() if has_photons else np.nan,
        'median_emission_ns': np.median(times) if has_photons else np.nan,
        'last_photon_ns': times[-1] if has_photons else np.nan,
        'time_to_photon_ns': (times - run.born).mean() if has_photons else np.nan,
    }


def run_pulses(args):
    points, swept = grid(args, PULSE_DEFAULTS)
    folder = new_run_folder('pulse', args.name, swept)
    (folder / 'photons').mkdir()
    shutil.copy(HERE / 'rules.py', folder / 'rules.py')
    (folder / 'params.json').write_text(json.dumps({
        'mode': 'pulse', 'command': ' '.join(sys.argv[1:]), 'date': f'{datetime.now():%Y-%m-%d %H:%M:%S}',
        'defaults': PULSE_DEFAULTS, 'swept': {name: sorted({p[name] for p in points}) for name in swept},
        'points': len(points), 'molecules': args.molecules, 'units': UNITS, 'bin_ns': args.bin or 'tauS / 10',
        'seed': args.seed,
        'threads': nb.get_num_threads(),
    }, indent=2))
    pulse(build_rules({**DEFAULTS, 'n': 0}), 2, 2, 0.1, 1)  # load the compiled code now, so it isn't timed

    started = time.time()
    outputs = {name: (folder / f'{name}.csv').open('w', newline='') for name in ('summary', 'histogram', 'sites')}
    writers = {}

    def write(name, rows):
        if name not in writers:
            writers[name] = csv.DictWriter(outputs[name], rows[0].keys())
            writers[name].writeheader()
        writers[name].writerows(rows)
        outputs[name].flush()  # a stopped sweep keeps what it finished

    for i, p in enumerate(points):
        began = time.perf_counter()
        run = pulse(build_rules(p), p['n'] + 2, args.molecules, p['duration'], args.seed + 256 * i)
        wall = time.perf_counter() - began

        np.savez_compressed(folder / 'photons' / f'point_{i:04d}.npz', time=run.time, site=run.site, molecule=run.molecule,
                            born=run.born, excitations=run.excitations)
        latest = max(run.time[-1] if run.time.size else 0.0, p['duration'])
        width = args.bin or p['tauS'] / BINS_PER_LIFETIME
        photons = {'time': run.time, 'site': run.site}
        histogram = histogram_rows(i, photons, p['n'] + 1, args.molecules, p['duration'], linear_edges(latest, width))
        metrics = pulse_metrics(p, run, args.molecules, histogram)
        events = run.tally.counts.sum()
        write('summary', [{'point': i, **p, **metrics, 'wall_s': wall, 'events': events, 'events_per_wall_s': events / wall}])
        write('histogram', histogram)

        r = np.r_[0.0, site_distances(p['n'], p['spacing'])]
        E = np.r_[0.0, fret_efficiency(r[1:], p['fretRadius'])]
        sites = per_site(run.tally)
        write('sites', [{'point': i, 'n': p['n'], 'site': s, 'name': name, 'r_nm': r[s], 'fret_efficiency': E[s],
                         **{k: v for k, v in sites[s].items() if not k.startswith('time_in')}}
                        for s, name in enumerate(site_names(p['n']))])
        if len(points) > 1 and (i + 1) % max(1, len(points) // 10) == 0:
            print(f'{i + 1}/{len(points)} points, {time.time() - started:.0f} s', flush=True)
    for output in outputs.values():
        output.close()
    seconds = time.time() - started
    add_to_index(folder, 'pulse', len(points), swept, seconds)
    if len(points) == 1:
        print(f"excitations per molecule  {metrics['excitations_per_molecule']:.3f}")
        print(f"photons per molecule      {metrics['photons_per_molecule']:.4f}")
        print(f"photons per excitation    {metrics['photons_per_excitation']:.4f}")
        print(f"photons after the pulse   {metrics['share_of_photons_after_pulse']:.1%}")
        print(f"peak bin starts at        {metrics['peak_bin_start_ns']:.2f} ns")
        print(f"mean / last photon        {metrics['mean_emission_ns']:.2f} ns / {metrics['last_photon_ns']:.2f} ns")
    print(f'\nSimulated {len(points)} point{"s" * (len(points) > 1)} of {args.molecules:,} molecules in {seconds:.1f} s'
          f' on {nb.get_num_threads()} threads')
    print(f'Saved to {folder}')


PER_MOLECULE = {  # intensity columns, each per molecule, and the event behind it
    'excitations_per_molecule': 'excitation', 'isc_per_molecule': 'isc',
    'sens_singlets_lost_per_molecule': 'sens_singlet_lost', 'back_transfers_per_molecule': 'back_transfer',
    'tta_singlet_per_molecule': 'tta_singlet', 'tta_failed_per_molecule': 'tta_failed',
    'triplets_decayed_per_molecule': 'triplet_lost',
}


def intensity_counts(counts, molecules):
    """The raw counts of an intensity point from a [event, site] table of events, all per molecule."""
    photons = counts[EVENTS.index('photon')] / molecules
    return {'photons_per_molecule': photons.sum(), 'perylene_photons_per_molecule': photons[-1],
            'dpa_photons_per_molecule': photons[:-1].sum(),
            **{name: counts[EVENTS.index(event)].sum() / molecules for name, event in PER_MOLECULE.items()}}


def intensity_mc(p, args, seed):
    """Monte Carlo for one intensity point: batches of molecules until photons per molecule is known to
    args.precision (relative standard error), or args.max_seconds runs out."""
    rules, sites = build_rules(p), p['n'] + 2
    started = time.perf_counter()
    molecules, counts = 0, np.zeros((len(EVENTS), sites))
    sums = np.zeros(3)  # photons², excitations² and photons × excitations, added up over molecules
    batch, error = args.batch, np.inf
    while True:
        run = pulse(rules, sites, batch, p['duration'], seed + molecules)
        photons = np.bincount(run.molecule, minlength=batch)
        excitations = run.excitations.astype(float)
        sums += [(photons ** 2).sum(), (excitations ** 2).sum(), (photons * excitations).sum()]
        counts += run.tally.counts
        molecules += batch
        mean = counts[EVENTS.index('photon')].sum() / molecules
        if mean > 0:
            error = np.sqrt(max(sums[0] / molecules - mean ** 2, 0) / molecules) / mean
        if error <= args.precision or time.perf_counter() - started > args.max_seconds:
            break
    return {**intensity_counts(counts, molecules), 'molecules': molecules, 'mean_photons_sq': sums[0] / molecules,
            'mean_excitations_sq': sums[1] / molecules, 'mean_photons_x_excitations': sums[2] / molecules,
            'photons_relative_error': error}


MC_ONLY = dict.fromkeys(['molecules', 'mean_photons_sq', 'mean_excitations_sq', 'mean_photons_x_excitations',
                         'photons_relative_error'], np.nan)  # empty on exact rows


def run_intensity(args):
    """Total photons and excitations per molecule over one pulse, against kExc (or anything else swept)."""
    points, swept = grid(args, PULSE_DEFAULTS)
    method = args.method
    if method in ('exact', 'both') and max(p['n'] for p in points) > EXACT_MAX_N:
        sys.exit(f'the exact solve only goes up to n = {EXACT_MAX_N}. Use --method mc above that')
    checks = set()
    if method == 'exact' and args.checks:
        checks = {int(round(k)) for k in np.linspace(0, len(points) - 1, args.checks + 2)[1:-1]}
    folder = new_run_folder('intensity', args.name, swept)
    shutil.copy(HERE / 'rules.py', folder / 'rules.py')
    (folder / 'params.json').write_text(json.dumps({
        'mode': 'intensity', 'command': ' '.join(sys.argv[1:]), 'date': f'{datetime.now():%Y-%m-%d %H:%M:%S}',
        'units': UNITS, 'defaults': PULSE_DEFAULTS, 'swept': {name: sorted({p[name] for p in points}) for name in swept},
        'points': len(points), 'method': method, 'checks': sorted(checks), 'precision': args.precision,
        'max_seconds': args.max_seconds, 'batch': args.batch, 'seed': args.seed, 'threads': nb.get_num_threads(),
    }, indent=2))
    print(f'Saving to {folder}', flush=True)
    if method != 'exact' or checks:
        pulse(build_rules({**DEFAULTS, 'n': 0}), 2, 2, 0.1, 1)  # load the compiled code now, so it isn't timed

    started = time.time()
    with (folder / 'results.csv').open('w', newline='') as out:
        writer = None
        for i, p in enumerate(points):
            runs = {'exact': ['exact'], 'mc': ['mc'], 'both': ['exact', 'mc']}[method] + ['mc'] * (i in checks)
            for how in runs:
                began = time.perf_counter()
                if how == 'exact':
                    found = {**intensity_counts(pulse_totals(build_rules(p), p).counts, 1), **MC_ONLY}
                else:
                    found = intensity_mc(p, args, args.seed + 1_000_000_000 * i)
                row = {'point': i, 'method': how, **p, **found, 'wall_s': time.perf_counter() - began}
                if writer is None:
                    writer = csv.DictWriter(out, row.keys())
                    writer.writeheader()
                writer.writerow(row)
                out.flush()  # a stopped sweep keeps what it finished
                print(f"point {i + 1}/{len(points)} {how} kExc {p['kExc']:.4g} photons {found['photons_per_molecule']:.4g}"
                      f" in {row['wall_s']:.1f} s, {time.time() - started:.0f} s so far", flush=True)
    seconds = time.time() - started
    add_to_index(folder, 'intensity', len(points), swept, seconds)
    print(f'\nDone: {len(points)} point{"s" * (len(points) > 1)} in {seconds:.1f} s')
    print(f'Saved to {folder}')


def rebin(args):
    """Build a new histogram from a pulse run's saved photons. The original histogram stays as it is."""
    folder = DATA / args.folder
    if not (folder / 'photons').is_dir():
        sys.exit(f'no saved photons in {folder}')
    rows = []
    with (folder / 'summary.csv').open(newline='') as summary:
        for entry in csv.DictReader(summary):
            point = int(entry['point'])
            with np.load(folder / 'photons' / f'point_{point:04d}.npz') as saved:
                photons = {'time': saved['time'], 'site': saved['site']}
            if photons['time'].size == 0:
                continue
            duration = float(entry['duration'])
            edges = (log_edges(photons['time'], args.logbins) if args.logbins else
                     linear_edges(max(photons['time'][-1], duration), args.bin))
            rows += histogram_rows(point, photons, int(entry['n']) + 1, int(entry['molecules']), duration, edges)
    name = f'histogram_log{args.logbins}.csv' if args.logbins else f'histogram_bin{args.bin:g}.csv'
    with (folder / name).open('w', newline='') as out:
        writer = csv.DictWriter(out, rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f'Saved {folder / name}')


def report(p, m, sites, r, E):
    se = m['photons_per_excitation_se']
    print(f"photons per excitation  {m['photons_per_excitation']:.4f}" + ('' if np.isnan(se) else f' +/- {se:.4f}'))
    print(f"  from perylene         {m['perylene_photons']:.4f}")
    print(f"  from DPA              {m['dpa_photons']:.4f}")
    print(f"ISC yield               {m['isc_yield']:.4f}")
    print(f"back transfers          {m['back_transfers']:.4f}")
    print(f"triplets decayed        {m['triplets_decayed']:.4f}")
    print(f"time to photon          {m['time_to_photon_ns']:.2f} ns")
    print(f"balance (should be 1)   triplets {m['triplet_balance']:.4f}, singlets {m['singlet_balance']:.4f}\n")
    print(f"{'site':<11}{'r nm':>7}{'FRET':>8}{'photons':>9}{'back':>8}{'TTA S':>8}{'TTA fail':>9}{'time T':>8}{'time S':>8}")
    for s, name in enumerate(site_names(p['n'])):
        v = sites[s]
        print(f"{name:<11}{r[s]:>7.2f}{E[s]:>8.3f}{v['photons']:>9.4f}{v['back_transfers']:>8.4f}{v['tta_singlet']:>8.4f}"
              f"{v['tta_failed']:>9.4f}{v['time_in_T']:>8.4f}{v['time_in_S']:>8.4f}")


def check(args):
    """Compare Monte Carlo with the exact solve, and check the bookkeeping, on a few settings."""
    failed = 0
    settings = ({'n': 0}, {'n': 2}, {'n': 4}, {'n': 4, 'kExc': 5.0}, {'n': 5, 'f': 1.0, 'fretRadius': 0.0})
    for i, changes in enumerate(settings):
        p = {**DEFAULTS, **changes}
        rules, sites = build_rules(p), p['n'] + 2
        truth = summarize(exact(rules, sites))
        tally = monte_carlo(rules, sites, 2_000_000, 64, args.seed + 1000 * i)  # own seed: shared ones correlate
        mc = summarize(tally)
        gap = abs(mc['photons_per_excitation'] - truth['photons_per_excitation'])
        sigmas = gap / mc['photons_per_excitation_se']
        books = max(abs(mc['triplet_balance'] - 1), abs(mc['singlet_balance'] - 1))
        ok = sigmas < 4 and books < 1e-3
        failed += not ok
        print(f"{'ok  ' if ok else 'FAIL'} {str(changes):<46} exact {truth['photons_per_excitation']:.5f}  "
              f"Monte Carlo {mc['photons_per_excitation']:.5f} ({sigmas:.1f} sigma)  books off by {books:.1e}")

    # a pulse run ends with every molecule back in all G, so its bookkeeping has to balance exactly
    p = {**PULSE_DEFAULTS, 'n': 4}
    run = pulse(build_rules(p), p['n'] + 2, 20_000, p['duration'], args.seed)
    count = {event: run.tally.counts[i].sum() for i, event in enumerate(EVENTS)}
    gaps = [count['isc'] - 2 * count['tta_singlet'] - count['tta_failed'] - count['triplet_lost'],
            count['excitation'] + count['back_transfer'] - count['isc'] - count['sens_singlet_lost'],
            count['tta_singlet'] - count['photon'] - count['back_transfer'],
            count['photon'] - run.time.size]
    ok = not any(gaps)
    failed += not ok
    print(f"{'ok  ' if ok else 'FAIL'} pulse bookkeeping: every triplet, singlet and photon accounted for {gaps}")
    sys.exit(1 if failed else 0)


def main():
    parser = argparse.ArgumentParser(description='TTA Flora 4 simulation')
    parser.add_argument('mode', choices=['run', 'sweep', 'exact', 'pulse', 'intensity', 'rebin', 'plot', 'check'])
    parser.add_argument('folder', nargs='?', help='rebin, plot: the run folder inside sim/data')
    for name, value in PULSE_DEFAULTS.items():
        parser.add_argument(f'--{name}', help=f'default {value}. One value, a,b,c, a:b, a:b:count or a:b:logcount')
    parser.add_argument('--excitations', type=int, default=200_000, help='excitations per point, across all replicas')
    parser.add_argument('--replicas', type=int, default=64)
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--name', help='run name, used in the folder name')
    parser.add_argument('--molecules', type=int, default=100_000, help='pulse: how many molecules to flash')
    parser.add_argument('--bin', type=float, help='pulse, rebin: histogram bin width in ns (pulse default: tauS / 10)')
    parser.add_argument('--logbins', type=int, help='rebin: this many bins spaced by ratio, instead of --bin')
    parser.add_argument('--xmax', type=float, help='plot: end of the time axis in ns')
    parser.add_argument('--kexc-per-mw', type=float,
                        help='plot, intensity run: put laser power on the x axis, at this kExc (/ns) per mW')
    parser.add_argument('--method', choices=['exact', 'mc', 'both'], default='exact',
                        help='intensity: the exact solve (n up to 8), Monte Carlo, or both at every point')
    parser.add_argument('--checks', type=int, default=0, help='intensity, exact: also run Monte Carlo at this many points')
    parser.add_argument('--precision', type=float, default=0.01,
                        help='intensity, Monte Carlo: stop at this relative standard error on photons per molecule')
    parser.add_argument('--max-seconds', type=float, default=600, help='intensity, Monte Carlo: most seconds per point')
    parser.add_argument('--batch', type=int, default=100_000, help='intensity, Monte Carlo: molecules per batch')
    args = parser.parse_args()
    if args.duration is not None and args.mode not in ('pulse', 'intensity'):
        parser.error('--duration only applies to pulse and intensity')
    if args.mode in ('rebin', 'plot') and not args.folder:
        parser.error(f'{args.mode} needs a run folder')
    if args.mode == 'rebin' and not (args.bin or args.logbins):
        parser.error('rebin needs --bin or --logbins')
    if args.mode == 'check':
        check(args)
    elif args.mode == 'pulse':
        run_pulses(args)
    elif args.mode == 'intensity':
        run_intensity(args)
    elif args.mode == 'rebin':
        rebin(args)
    elif args.mode == 'plot':
        from plot_pulse import plot_run  # only this command needs matplotlib
        from plot_yield import plot_sweep
        from plot_intensity import is_intensity, save_all
        folder = DATA / args.folder
        if is_intensity(folder):
            print(f'Saved to {save_all(folder, args.kexc_per_mw)}')
        else:
            print(f'Saved to {plot_run(folder, args.xmax) if (folder / "photons").is_dir() else plot_sweep(folder)}')
    else:
        run_points(args.mode, args)


if __name__ == '__main__':
    main()
