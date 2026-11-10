"""The TTA Flora 4 model: default parameters, site distances and the rules. This is the file to edit.

The system is a string of n + 2 letters, one per site: the sensitizer (site 0), DPA 1 to n (sites 1 to n) and the
perylene (site n + 1). Each letter is G (ground), T (triplet) or S (singlet).

Times are in ns and rates per ns.

A rule rewrites one site, or two sites, at a rate per ns. There is no time step: the engine adds up the rates of
every rule that can fire, waits a random time, and fires one of them (see engine.py).
"""

import json
from pathlib import Path

import numpy as np

from engine import EVENTS

LETTERS = 'GTS'

DEFAULTS = {
    'n': 8,  # DPA mediators
    'tauS': 5.0,  # singlet lifetime on a DPA or the perylene, ns
    'tauSens': 5.0,  # sensitizer singlet decay lifetime, ns (rule 9), not counting ISC
    'tripletYield': 0.5,  # share of sensitizer singlets that become triplets. This sets kISC
    'hopRatio': 100.0,  # how many times faster hopping is than the singlet lifetime: kHop = hopRatio / tauS
    'excRatio': 0.1,  # excitation rate as a multiple of the sensitizer's lifetime: kExc = excRatio / tauSens
    'kDecay': 1e-6,  # triplet decay on a DPA or the perylene, /ns
    'f': 0.4,  # the F factor (spin-statistical factor): the share of TTA events that make a singlet
    'fretRadius': 4.0,  # nm: the distance from the sensitizer where back transfer is 50% (0 = no back transfer)
    'spacing': 0.0,  # nm: 0 = distances measured from the real molecule; above 0 = site i sits at i * spacing
    'block': 1,  # 1 = a triplet can't hop onto a singlet. 0 = it lands and quenches it (for comparison only)
    # normally worked out by resolve(). Giving one of these sets it directly, and the yield or ratio follows from it
    'kISC': None,  # sensitizer S -> T, /ns
    'kHop': None,  # triplet hops per ns, both directions together, /ns
    'kExc': None,  # sensitizer excitations per ns, /ns
}

# a pulse run: kExc is on from t = 0 to `duration`, then 0
PULSE_DEFAULTS = {**DEFAULTS, 'excRatio': 1.0, 'duration': 5.0}

# the fluorescence reference in live.py: a plain fluorophore with two states, G and S, on one site. No triplets, no hops
FLUORESCENCE_DEFAULTS = {
    'tauFl': 5.0,  # singlet lifetime, ns: S -> G at 1/tauFl in all
    'fluorYield': 1.0,  # share of those S -> G that give a photon. The rest are lost without one
}


def resolve(p):
    """Fill in kISC, kHop and kExc from the yield and ratios, or the other way round if a rate was given directly."""
    p = dict(p)
    decay = 1 / p['tauSens']
    if p['kISC'] is None:
        if p['tripletYield'] >= 1:
            raise SystemExit('tripletYield must be below 1: at 1 no singlet decays, so give --kISC instead')
        # ISC and decay compete, so the yield is kISC / (kISC + decay)
        p['kISC'] = decay * p['tripletYield'] / (1 - p['tripletYield'])
    else:
        p['tripletYield'] = p['kISC'] / (p['kISC'] + decay)
    if p['kHop'] is None:
        p['kHop'] = p['hopRatio'] / p['tauS']
    else:
        p['hopRatio'] = p['kHop'] * p['tauS']
    if p['kExc'] is None:
        p['kExc'] = p['excRatio'] * decay
    else:
        p['excRatio'] = p['kExc'] * p['tauSens']
    return p


SAMPLES = Path(__file__).resolve().parent.parent / 'samples' / 'triplet'


def _xyz(rows):
    return np.array([row[1:] for row in rows], float)


def site_distances(n, spacing=0.0):
    """Distance in nm from the sensitizer to sites 1 to n + 1 (index 0 of the result is DPA 1, or the perylene if n = 0).

    Measured from the same geometry the app builds, using the same centres as its "core centre to centre" note:
    the sensitizer's core atoms, the middle of each DPA (its anthracene centre) and the perylene's core atoms.
    """
    if spacing > 0:
        return spacing * np.arange(1, n + 2)
    flora = json.loads((SAMPLES / 'tta-flora-4.param.json').read_text())
    if n == 0:
        zero = _xyz(flora['zero'])
        gap = zero[flora['zeroRightCore']].mean(0) - zero[flora['zeroLeftCore']].mean(0)
        return np.array([np.linalg.norm(gap) / 10])
    wire = json.loads((SAMPLES / flora['wire']).read_text())
    shift = np.array(wire['t'])  # one DPA + linker repeat
    sens = _xyz(flora['left'])[flora['leftCore']].mean(0)
    dpa = _xyz(wire['dpa']).mean(0)
    perylene = _xyz(flora['right'])[flora['rightCore']].mean(0) + (n - 1) * shift
    centres = [dpa + k * shift for k in range(n)] + [perylene]
    return np.array([np.linalg.norm(c - sens) / 10 for c in centres])


def fret_efficiency(r, radius):
    """E = 1 / (1 + (r / R0)^6): the chance a singlet at distance r goes back to the sensitizer instead of emitting."""
    if radius <= 0:
        return np.zeros_like(r)
    return 1 / (1 + (r / radius) ** 6)


def hops(n, kHop):
    """Every allowed hop as (from site, to site, rate). Nothing hops onto the sensitizer or off the perylene."""
    moves = [(0, 1, kHop)]
    for site in range(1, n + 1):
        neighbours = [s for s in (site - 1, site + 1) if s >= 1]
        moves += [(site, to, kHop / len(neighbours)) for to in neighbours]
    return moves


def build_rules(p):
    """The rules for parameters p, as (site a, site b or -1, before a, before b, after a, after b, rate, event, where).

    `where` is the site an event is counted at. Rules with rate 0 are left out.
    """
    p = resolve(p)
    n = int(p['n'])
    E = fret_efficiency(site_distances(n, p['spacing']), p['fretRadius'])
    rules = []

    def one(site, before, after, rate, event):
        rules.append((site, -1, LETTERS.index(before), 0, LETTERS.index(after), 0, rate, EVENTS.index(event), site))

    def two(a, b, before, after, rate, event, where):
        (x, y), (u, v) = before.split(), after.split()
        rules.append((a, b, LETTERS.index(x), LETTERS.index(y), LETTERS.index(u), LETTERS.index(v), rate,
                      EVENTS.index(event), where))

    one(0, 'G', 'S', p['kExc'], 'excitation')  # rule 1
    one(0, 'S', 'T', p['kISC'], 'isc')  # rule 2
    one(0, 'S', 'G', 1 / p['tauSens'], 'sens_singlet_lost')  # rule 9

    emit = 1 / p['tauS']
    for site in range(1, n + 2):
        one(site, 'T', 'G', p['kDecay'], 'triplet_lost')  # rule 3
        back = E[site - 1]
        # rules 4 and 10: the sensitizer can only take the singlet back while it's in G. Otherwise FRET is off and
        # the site just emits. The first letter is the sensitizer, the second this site.
        two(0, site, 'G S', 'G G', emit * (1 - back), 'photon', site)
        two(0, site, 'G S', 'S G', emit * back, 'back_transfer', site)
        two(0, site, 'S S', 'S G', emit, 'photon', site)
        two(0, site, 'T S', 'T G', emit, 'photon', site)

    for a, b, rate in hops(n, p['kHop']):
        two(a, b, 'T G', 'G T', rate, 'hop', b)  # rule 5
        # rule 6: a hop onto a singlet is blocked, which means there's simply no 'T S' rule
        if not p['block']:
            two(a, b, 'T S', 'G T', rate, 'hop', b)
        two(a, b, 'T T', 'G S', p['f'] * rate, 'tta_singlet', b)  # rule 7
        two(a, b, 'T T', 'G T', (1 - p['f']) * rate, 'tta_failed', b)  # rule 8: one triplet survives, on site b

    return [rule for rule in rules if rule[6] > 0]


def fluorescence_rules(p):
    """A two-state fluorophore on site 0, in the same form as build_rules: excited at kExc, and back to G at 1/tauFl,
    as a photon with chance fluorYield. p needs kExc, tauFl and fluorYield.
    """
    G, S = LETTERS.index('G'), LETTERS.index('S')
    decay = 1 / p['tauFl']
    rules = [(0, -1, G, 0, S, 0, p['kExc'], EVENTS.index('excitation'), 0),
             (0, -1, S, 0, G, 0, decay * p['fluorYield'], EVENTS.index('photon'), 0),
             (0, -1, S, 0, G, 0, decay * (1 - p['fluorYield']), EVENTS.index('sens_singlet_lost'), 0)]
    return [rule for rule in rules if rule[6] > 0]
