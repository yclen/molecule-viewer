"""Check the browser's simulation engine (src/sim/) against the Python one (sim/engine.py), which is the reference.

  python tools/check_js_engine.py

Needs Node and the sim/ packages (the launcher's .venv has them). Compares, at several parameter sets:
  1. the rules, exactly
  2. pulse_totals (every event per molecule over one pulse), to 1e-6 relative
  3. pulse_exact's histogram, to 1e-4 of its peak bin
  4. the JS Monte Carlo against the exact totals and histogram, within its statistical error
Exits non-zero if anything is off.
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'sim'))

from engine import EVENTS, pulse_exact, pulse_totals  # noqa: E402
from rules import PULSE_DEFAULTS, build_rules, resolve  # noqa: E402

failures = []


def js(job):
    done = subprocess.run(['node', str(ROOT / 'tools' / 'js_engine_cli.mjs')], input=json.dumps(job), capture_output=True,
                          text=True, cwd=ROOT)
    if done.returncode:
        raise RuntimeError(done.stderr)
    return json.loads(done.stdout)


def params(**changes):
    p = resolve({**PULSE_DEFAULTS, **changes})
    return {k: v for k, v in p.items()}


def check(name, ok, detail):
    print(f"{'ok  ' if ok else 'FAIL'} {name}: {detail}")
    if not ok:
        failures.append(name)


def compare_rules(p):
    py = [list(r) for r in build_rules(p)]
    mine = js({'kind': 'rules', 'p': p})['out']
    same = len(py) == len(mine) and all(
        a[:6] == b[:6] and a[7:] == b[7:] and abs(a[6] - b[6]) <= 1e-12 * max(abs(a[6]), 1) for a, b in zip(py, mine))
    check(f"rules n={p['n']} fret={p['fretRadius']}", same, f'{len(py)} rules')


def compare_totals(p):
    py = pulse_totals(build_rules(p), p).counts
    got = js({'kind': 'totals', 'p': p})
    mine = np.array([np.nan if v is None else v for v in got['out']]).reshape(py.shape)
    keep = np.isfinite(py)
    scale = np.maximum(np.abs(py[keep]), 1e-300)
    worst = np.max(np.abs(mine[keep] - py[keep]) / np.maximum(scale, 1e-12 * np.abs(py[keep]).max()))
    photons = py[EVENTS.index('photon')].sum()
    check(f"totals n={p['n']} kExc={p['kExc']:g} T={p['duration']:g}", worst < 1e-6,
          f'photons/molecule {photons:.6g}, worst relative gap {worst:.1e}, JS {got["seconds"]:.2f} s')


def compare_exact(p, window):
    edges = np.linspace(0, window, 101)
    py = pulse_exact(build_rules(p), p, window).binned(edges)
    got = js({'kind': 'exact', 'p': p, 'window': window, 'edges': edges.tolist()})
    mine = [np.array(a) for a in got['out']['binned']]
    peak = max(a.max() for a in py)
    gap = max(np.abs(a - b).max() for a, b in zip(py, mine)) / peak
    check(f"pulse_exact n={p['n']} kExc={p['kExc']:g}", gap < 1e-4,
          f'largest bin gap {gap:.1e} of the peak, JS {got["seconds"]:.2f} s')


def compare_mc(p, molecules=200_000):
    exact = pulse_totals(build_rules(p), p).counts
    got = js({'kind': 'mc', 'p': p, 'molecules': molecules, 'seed': 7})['out']
    photons = np.array(got['times']).size / molecules
    mean_sq = got['sums'][0] / molecules
    se = np.sqrt(max(mean_sq - photons ** 2, 0) / molecules)
    want = exact[EVENTS.index('photon')].sum()
    z = (photons - want) / se
    check(f"MC photons n={p['n']} kExc={p['kExc']:g}", abs(z) < 4, f'{photons:.5g} vs exact {want:.5g} ({z:+.1f} SE)')
    exc = got['excitations'] / molecules
    want_exc = exact[EVENTS.index('excitation')].sum()
    z = (exc - want_exc) / np.sqrt(got['sums'][1] / molecules - exc ** 2) * np.sqrt(molecules)
    check(f"MC excitations n={p['n']} kExc={p['kExc']:g}", abs(z) < 4, f'{exc:.5g} vs exact {want_exc:.5g} ({z:+.1f} SE)')
    # the histogram against the exact one: a chi-square over bins with enough counts
    window = 10 * p['tauS']
    edges = np.linspace(0, window, 51)
    expected = sum(pulse_exact(build_rules(p), p, window).binned(edges)) * molecules
    counts = np.histogram(got['times'], edges)[0]
    use = expected > 20
    chi2 = (((counts - expected) ** 2 / expected)[use]).sum()
    dof = use.sum()
    check(f"MC histogram n={p['n']} kExc={p['kExc']:g}", chi2 < dof + 5 * np.sqrt(2 * dof),
          f'chi² {chi2:.0f} over {dof} bins')


for p in (params(n=0), params(n=3), params(n=8), params(n=5, fretRadius=0), params(n=4, spacing=1.5, block=0)):
    compare_rules(p)
for p in (params(n=2, kExc=0.01), params(n=4, kExc=0.2), params(n=4, kExc=50.0), params(n=6, kExc=1.0, duration=25.0),
          params(n=3, kExc=1e-3)):
    compare_totals(p)
for p in (params(n=3, kExc=0.2, kHop=20.0), params(n=5, kExc=2.0, kHop=20.0)):
    compare_exact(p, 50.0)
for p in (params(n=3, kExc=0.5, kHop=20.0), params(n=6, kExc=0.2, kHop=20.0)):
    compare_mc(p)

print('\nAll checks passed.' if not failures else f'\n{len(failures)} failed: {", ".join(failures)}')
sys.exit(1 if failures else 0)
