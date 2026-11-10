"""Runs the rules from rules.py.

monte_carlo: a Gillespie simulation of the long-run behaviour. Add up the rates of every rule that can fire, wait a
random time, fire one.
exact: the same long-run answers by solving for the share of time in every possible string, with no random noise.
There are 3^(n + 2) strings, so this is for checking results at small n, not for sweeps.
pulse_exact: the expected pulse histogram the same way, by following the share of molecules in every string over time.
pulse: many molecules, each starting in all G, with a square pulse of light. Records every photon's emission time.

monte_carlo, exact and pulse give a Tally, so summarize() treats them alike. Simulated time is in ns; the solve
time estimates (exact_seconds, _setup_s, _step_s) are real seconds.
"""

from dataclasses import dataclass

import numba as nb
import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg
from scipy.interpolate import PchipInterpolator

EVENTS = ['excitation', 'isc', 'sens_singlet_lost', 'hop', 'triplet_lost', 'tta_singlet', 'tta_failed', 'photon',
          'back_transfer']
N_EVENTS = len(EVENTS)
EXCITATION, HOP, TTA_SINGLET, TTA_FAILED, PHOTON, BACK = (EVENTS.index(e) for e in (
    'excitation', 'hop', 'tta_singlet', 'tta_failed', 'photon', 'back_transfer'))


@dataclass
class Tally:
    counts: np.ndarray  # [event, site]
    time: float  # simulated ns these counts cover
    occupancy: np.ndarray  # [site, letter]: ns spent
    latency: float = np.nan  # mean ns from excitation to photon
    yield_se: float = np.nan  # standard error of photons per excitation


def rule_arrays(rules):
    table = np.array(rules, float)
    return tuple(table[:, k].astype(np.int64) for k in range(6)) + (table[:, 6], table[:, 7].astype(np.int64),
                                                                    table[:, 8].astype(np.int64))


@nb.njit(cache=True)
def _rates(a, b, fa, fb, rate, event, state, light_on, weights):
    """Fill in the rate of every rule that can fire now, and return their total."""
    total = 0.0
    for r in range(rate.size):
        can = state[a[r]] == fa[r] and (b[r] < 0 or state[b[r]] == fb[r]) and (light_on or event[r] != EXCITATION)
        weights[r] = rate[r] if can else 0.0
        total += weights[r]
    return total


@nb.njit(cache=True)
def _pick(weights, total):
    """Choose a rule, each with a chance in proportion to its rate."""
    u = np.random.random() * total
    r = 0
    acc = weights[0]
    while (acc < u or weights[r] == 0.0) and r < weights.size - 1:
        r += 1
        acc += weights[r]
    return r


@nb.njit(cache=True)
def _carry_birth(e, a, b, t, born):
    """Keep each excitation's birth time with it. After TTA the result carries the younger triplet's."""
    if e == EXCITATION:
        born[0] = t
    elif e == HOP:
        born[b] = born[a]
    elif e == TTA_SINGLET or e == TTA_FAILED:
        born[b] = max(born[a], born[b])
    elif e == BACK:
        born[0] = born[b]


@nb.njit(cache=True)
def _simulate(a, b, fa, fb, ta, tb, rate, event, where, sites, warmup, target, seed):
    """One replica from the all-G string. Throws away the first `warmup` excitations, then runs for `target` more."""
    np.random.seed(seed)
    state = np.zeros(sites, np.int64)
    born = np.zeros(sites)  # when the excitation on each site was made
    counts = np.zeros((N_EVENTS, sites), np.int64)
    occupancy = np.zeros((sites, 3))
    weights = np.zeros(rate.size)
    t = 0.0
    start = 0.0
    excited = 0
    measuring = warmup == 0
    latency = 0.0
    while True:
        total = _rates(a, b, fa, fb, rate, event, state, True, weights)
        dt = -np.log(1.0 - np.random.random()) / total
        if measuring:
            for s in range(sites):
                occupancy[s, state[s]] += dt
        t += dt

        r = _pick(weights, total)
        state[a[r]] = ta[r]
        if b[r] >= 0:
            state[b[r]] = tb[r]

        e = event[r]
        _carry_birth(e, a[r], b[r], t, born)
        if e == EXCITATION:
            excited += 1
            if not measuring and excited >= warmup:
                measuring = True
                start = t
                occupancy[:] = 0.0

        if measuring:
            counts[e, where[r]] += 1
            if e == PHOTON:
                latency += t - born[where[r]]
            if counts[EXCITATION, 0] >= target:
                return counts, occupancy, t - start, latency


@nb.njit(parallel=True, cache=True)
def _replicas(a, b, fa, fb, ta, tb, rate, event, where, sites, warmup, target, seed, replicas):
    counts = np.zeros((replicas, N_EVENTS, sites), np.int64)
    occupancy = np.zeros((replicas, sites, 3))
    time = np.zeros(replicas)
    latency = np.zeros(replicas)
    for k in nb.prange(replicas):
        c, o, t, lat = _simulate(a, b, fa, fb, ta, tb, rate, event, where, sites, warmup, target, seed + k)
        counts[k] = c
        occupancy[k] = o
        time[k] = t
        latency[k] = lat
    return counts, occupancy, time, latency


def monte_carlo(rules, sites, excitations=100_000, replicas=64, seed=1, warmup=200):
    """Run `replicas` independent copies in parallel, sharing out `excitations` between them."""
    target = max(1, excitations // replicas)
    counts, occupancy, time, latency = _replicas(*rule_arrays(rules), sites, warmup, target, seed, replicas)
    photons = counts[:, EVENTS.index('photon')].sum(1)
    yields = photons / counts[:, EXCITATION].sum(1)
    return Tally(counts.sum(0).astype(float), time.sum(), occupancy.sum(0),
                 latency.sum() / photons.sum() if photons.sum() else np.nan, yields.std(ddof=1) / replicas ** 0.5)


@dataclass
class Pulse:
    """One pulse run. The photon arrays hold one entry per photon, sorted by emission time."""
    time: np.ndarray  # emission time, s after the pulse starts
    site: np.ndarray  # site it came from
    molecule: np.ndarray  # which molecule
    born: np.ndarray  # when the younger triplet behind it was excited
    excitations: np.ndarray  # rule 1 events on each molecule
    tally: Tally


@nb.njit(cache=True, nogil=True)  # nogil: live.py runs this on several threads at once
def _pulse_chunk(a, b, fa, fb, ta, tb, rate, event, where, sites, duration, first, last, seed, times, from_site, from_mol,
                 from_born, counts, excitations):
    """Molecules first to last - 1, one after the other. Returns how many photons they gave, even past the arrays' end."""
    np.random.seed(seed)
    state = np.zeros(sites, np.int64)
    born = np.zeros(sites)
    weights = np.zeros(rate.size)
    found = 0
    for molecule in range(first, last):
        state[:] = 0
        t = 0.0
        while True:
            light_on = t < duration
            total = _rates(a, b, fa, fb, rate, event, state, light_on, weights)
            if total == 0.0:  # back to all G after the pulse, so nothing else can happen
                break
            dt = -np.log(1.0 - np.random.random()) / total
            if light_on and t + dt >= duration:
                t = duration  # the rates change here. Waits have no memory, so just draw a new one
                continue
            t += dt

            r = _pick(weights, total)
            state[a[r]] = ta[r]
            if b[r] >= 0:
                state[b[r]] = tb[r]
            e = event[r]
            _carry_birth(e, a[r], b[r], t, born)
            counts[e, where[r]] += 1
            if e == EXCITATION:
                excitations[molecule] += 1
            elif e == PHOTON:
                if found < times.size:
                    times[found] = t
                    from_site[found] = where[r]
                    from_mol[found] = molecule
                    from_born[found] = born[where[r]]
                found += 1
    return found


@nb.njit(parallel=True, cache=True)
def _pulse_chunks(a, b, fa, fb, ta, tb, rate, event, where, sites, duration, edges, seed, room):
    chunks = edges.size - 1
    times = np.zeros((chunks, room))
    from_site = np.zeros((chunks, room), np.int16)
    from_mol = np.zeros((chunks, room), np.int32)
    from_born = np.zeros((chunks, room))
    counts = np.zeros((chunks, N_EVENTS, sites), np.int64)
    excitations = np.zeros(edges[-1], np.int32)
    found = np.zeros(chunks, np.int64)
    for c in nb.prange(chunks):
        found[c] = _pulse_chunk(a, b, fa, fb, ta, tb, rate, event, where, sites, duration, edges[c], edges[c + 1], seed + c,
                                times[c], from_site[c], from_mol[c], from_born[c], counts[c], excitations)
    return times, from_site, from_mol, from_born, counts, excitations, found


def pulse(rules, sites, molecules, duration, seed=1, chunks=256):
    """Flash `molecules` molecules with the rules' kExc for `duration` ns, and follow each until nothing is left.

    Molecules are shared out in a fixed number of chunks, so the same seed gives the same photons on any machine.
    """
    edges = np.linspace(0, molecules, min(chunks, molecules) + 1).astype(np.int64)
    room = max(256, 2 * molecules // (edges.size - 1))
    while True:
        times, from_site, from_mol, from_born, counts, excitations, found = _pulse_chunks(
            *rule_arrays(rules), sites, duration, edges, seed, room)
        if found.max() <= room:
            break
        room = int(found.max())  # more photons than there was room for: same seed, so this run fits exactly

    keep = np.arange(room)[None, :] < found[:, None]
    order = np.argsort(times[keep], kind='stable')
    yields = counts[:, PHOTON].sum(1) / np.maximum(counts[:, EXCITATION].sum(1), 1)
    tally = Tally(counts.sum(0).astype(float), 1.0, np.zeros((sites, 3)), yield_se=yields.std(ddof=1) / yields.size ** 0.5)
    return Pulse(times[keep][order], from_site[keep][order], from_mol[keep][order], from_born[keep][order], excitations,
                 tally)


def pulse_alone(arrays, sites, molecules, duration, seed):
    """A few molecules on this thread only, for live.py: its photon times, their sites, and the excitations.

    `pulse` waits for its slowest thread every time, which costs a lot in small batches. Calling this from several
    threads at once doesn't, because none of them waits for another. `arrays` is rule_arrays(rules).
    """
    room = max(256, 4 * molecules)
    while True:
        times, from_site = np.zeros(room), np.zeros(room, np.int16)
        excitations = np.zeros(molecules, np.int32)
        found = _pulse_chunk(*arrays, sites, duration, 0, molecules, seed, times, from_site, np.zeros(room, np.int32),
                             np.zeros(room), np.zeros((N_EVENTS, sites), np.int64), excitations)
        if found <= room:
            return times[:found], from_site[:found], int(excitations.sum())
        room = found  # more photons than there was room for: same seed, so this run fits exactly


def _transitions(rules, sites):
    """Every rule as (strings it can fire from, the strings it leads to, rate, event, where), over all 3^sites strings.

    String number k has letter (k // 3^s) % 3 at site s, so all G is string 0.
    """
    index = np.arange(3 ** sites)
    digit = [(index // 3 ** s) % 3 for s in range(sites)]
    moves = []
    for a, b, fa, fb, ta, tb, rate, event, where in rules:
        can = digit[a] == fa
        if b >= 0:
            can = can & (digit[b] == fb)
        src = index[can]
        dst = src + (ta - fa) * 3 ** a + ((tb - fb) * 3 ** b if b >= 0 else 0)
        moves.append((src, dst, rate, event, where))
    return moves, digit


def exact(rules, sites):
    """Long-run behaviour by linear algebra: solve pi Q = 0 over all 3^sites strings."""
    size = 3 ** sites
    moves, digit = _transitions(rules, sites)
    rows, cols, vals, fires = [], [], [], []
    for src, dst, rate, event, where in moves:
        rows += [dst, src]
        cols += [src, src]
        vals += [np.full(src.size, rate), np.full(src.size, -rate)]
        fires.append((src, rate, event, where))
    Q = sparse.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(size, size))
    A = Q.tolil()
    A[0, :] = 1.0  # the shares add up to 1, in place of one redundant equation
    A = A.tocsc()
    rhs = np.zeros(size)
    rhs[0] = 1.0
    if sites <= 8:
        pi = linalg.spsolve(A, rhs)
    else:
        ilu = linalg.spilu(A, drop_tol=1e-4, fill_factor=30)
        pi, info = linalg.gmres(A, rhs, M=linalg.LinearOperator((size, size), ilu.solve), rtol=1e-10, restart=100,
                                maxiter=30)
        if info:
            raise RuntimeError(f'the exact solve did not converge (gmres info {info})')
    counts = np.zeros((len(EVENTS), sites))
    for src, rate, event, where in fires:
        counts[event, where] += rate * pi[src].sum()
    occupancy = np.array([[pi[digit[s] == letter].sum() for letter in range(3)] for s in range(sites)])
    return Tally(counts, 1.0, occupancy)


# pulse_exact: the share of molecules in every string, followed through time. It's a stiff system (hops are far
# faster than emission), so it steps with a 2-stage L-stable SDIRK, which takes steps much longer than a hop without
# going unstable. Each stage is a linear solve with the same matrix, done by GMRES with an incomplete LU of that matrix.
_GAMMA = 1 - 2 ** -0.5  # the SDIRK's diagonal
EXACT_ILU = {'drop_tol': 1e-3, 'fill_factor': 5}  # looser is quicker to set up, but needs more GMRES steps
EXACT_STEPS_PER_LIFETIME = 10  # of tauS or tauSens, whichever is shorter
EXACT_LEAST_STEPS = 10  # in the pulse, and after it
EXACT_MOST_STEPS = 2000  # in the time window. A longer window takes longer steps
EXACT_MOST_TAIL_STEPS = 2000  # after the window, until no molecule can still emit
QUIET = 1e-7  # the count after the window stops when the molecules that can still emit are this share of excitations


def exact_plan(p, window):
    """The time steps for resolved parameters p: [(light on, step length, number of steps)], up to the window.

    A step ends exactly where the pulse does.
    """
    step = max(min(p['tauS'], p['tauSens']) / EXACT_STEPS_PER_LIFETIME, window / EXACT_MOST_STEPS)
    plan = []
    for light, start, end in ((True, 0, min(p['duration'], window)), (False, p['duration'], window)):
        if end > start:
            count = max(EXACT_LEAST_STEPS, int(np.ceil((end - start) / step - 1e-9)))
            plan.append((light, (end - start) / count, count))
    return plan


SETUP_S = {6: 0.19, 7: 1.2, 8: 7.3, 9: 110}  # measured: the incomplete LU fills in faster and faster as n grows


def _setup_s(size):
    n = round(np.log(size) / np.log(3)) - 2
    return SETUP_S.get(n, 3e-8 * size ** 1.77 if n < 6 else SETUP_S[9] * 15 ** (n - 9))


def _step_s(size, stiffness):
    """One step. Stiffer steps (more hops per step) need more GMRES iterations."""
    return 1.5e-7 * size ** 1.2 * (1 + 0.35 * np.log10(max(stiffness / 5, 1)))


def exact_seconds(p, window):
    """Roughly how long pulse_exact takes, on the machine it was timed on. The count after the window is a guess."""
    size = 3 ** (p['n'] + 2)
    plan = exact_plan(p, window)
    seconds = sum(_setup_s(size) + _step_s(size, p['kHop'] * h) * count for _, h, count in plan)
    after = sum(count for light, _, count in plan if not light) / 4
    return seconds + after * _step_s(size, p['kHop'] * plan[-1][1])


@dataclass
class PulseExact:
    times: np.ndarray  # the end of each step, from 0 to the window, ns
    perylene: np.ndarray  # expected photons per molecule from the perylene, from t = 0 to each time
    dpa: np.ndarray  # the same from the DPAs
    excitations: float  # per molecule, over the whole pulse. nan if the pulse lasts past the window
    photons: float  # per molecule, counting the ones after the window. nan likewise
    still_active: float  # share of molecules that could still emit when the count after the window stopped

    def binned(self, edges):
        """Expected perylene and DPA photons per molecule in each bin between `edges`."""
        # each total only rises, and PCHIP keeps it rising, so no bin goes below 0
        edges = np.clip(edges, 0, self.times[-1])
        return [np.diff(PchipInterpolator(self.times, total)(edges)) for total in (self.perylene, self.dpa)]


def pulse_exact(rules, p, window, progress=None):
    """The expected pulse histogram, with no random noise: a pulse run with infinitely many molecules.

    `rules` are build_rules(p), with p resolved. Every molecule starts in all G, with kExc on until p['duration'].
    `progress(share done, what it's doing)` is called as it goes. exact_seconds says roughly how long it takes: at the
    defaults, a few seconds up to n = 7, under a minute at n = 8, and each further n is several times longer.
    """
    sites = p['n'] + 2
    size = 3 ** sites
    moves, _ = _transitions(rules, sites)
    counter = {PHOTON: lambda where: size if where == sites - 1 else size + 1, EXCITATION: lambda where: size + 2}

    def matrix(light):
        """Q, plus rows that add up perylene photons, DPA photons and excitations."""
        rows, cols, vals = [], [], []
        for src, dst, rate, event, where in moves:
            if event == EXCITATION and not light:
                continue
            rows += [dst, src]
            cols += [src, src]
            vals += [np.full(src.size, rate), np.full(src.size, -rate)]
            if event in counter:
                rows.append(np.full(src.size, counter[event](where)))
                cols.append(src)
                vals.append(np.full(src.size, rate))
        return sparse.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                                 shape=(size + 3, size + 3))

    plan = exact_plan(p, window)
    total = exact_seconds(p, window)
    spent = 0.0

    def report(seconds, doing):
        nonlocal spent
        spent += seconds
        if progress:
            progress(min(spent / total, 1.0), doing)

    # all G, or one triplet with nothing else excited: neither can make a photon without a new excitation
    quiet = [0] + [3 ** s for s in range(sites)]
    x = np.zeros(size + 3)
    x[0] = 1.0
    times, counts = [0.0], [x[size:].copy()]
    still_active = np.nan
    for light, h, count in plan:
        report(0, 'Preparing')
        M = matrix(light)
        A = (sparse.identity(size + 3, format='csr') - _GAMMA * h * M).tocsc()
        ilu = linalg.LinearOperator(A.shape, linalg.spilu(A, **EXACT_ILU).solve)
        report(_setup_s(size), 'Solving')
        each = _step_s(size, p['kHop'] * h)
        k1, k2 = np.zeros_like(x), np.zeros_like(x)

        def step():
            nonlocal x, k1, k2
            k1, info1 = linalg.gmres(A, M @ x, x0=k1, M=ilu, rtol=1e-8, atol=0)
            k2, info2 = linalg.gmres(A, M @ (x + (1 - _GAMMA) * h * k1), x0=k2, M=ilu, rtol=1e-8, atol=0)
            if info1 or info2:
                raise RuntimeError('the exact solve did not converge')
            x = x + h * ((1 - _GAMMA) * k1 + _GAMMA * k2)

        for _ in range(count):
            step()
            times.append(times[-1] + h)
            counts.append(x[size:].copy())
            report(each, 'Solving')
        if not light:
            for _ in range(EXACT_MOST_TAIL_STEPS):
                still_active = x[:size].sum() - x[quiet].sum()
                if still_active < QUIET * x[size + 2]:
                    break
                step()
                report(each, 'Counting the photons after the window')
    whole = p['duration'] <= window
    return PulseExact(np.array(times), np.array(counts)[:, 0], np.array(counts)[:, 1],
                      x[size + 2] if whole else np.nan, x[size] + x[size + 1] if whole else np.nan, still_active)


@dataclass
class FluorescenceExact:
    """pulse_exact's answer for a two-state fluorophore (rules.fluorescence_rules), which has a closed form.

    With a = kExc and k = 1/tauFl, the share in S during the pulse is P(t) = a/(a + k) (1 - exp(-(a + k) t)), and after
    it P decays at k. Photons come at fluorYield * k * P. Every excitation ends in S -> G, so over the whole pulse the
    photons are fluorYield times the excitations.
    """
    p: dict
    excitations: float = np.nan
    photons: float = np.nan
    still_active: float = 0.0  # the closed form counts every photon, however late

    def __post_init__(self):
        self.excitations = self.decays(np.inf)
        self.photons = self.p['fluorYield'] * self.excitations

    def decays(self, t):
        """S -> G from t = 0 to each time in t, per molecule."""
        a, k, pulse = self.p['kExc'], 1 / self.p['tauFl'], self.p['duration']
        s = a + k
        t = np.asarray(t, float)
        lit = np.minimum(t, pulse)
        during = k * a / s * (lit + np.expm1(-s * lit) / s)
        at_end = -a / s * np.expm1(-s * pulse)  # the share in S when the pulse ends
        dark = np.maximum(t - pulse, 0)
        return during + at_end * -np.expm1(-k * dark)

    def binned(self, edges):
        """Expected photons per molecule in each bin, in the same two-list form as PulseExact (the second is empty)."""
        photons = np.diff(self.p['fluorYield'] * self.decays(edges))
        return [photons, np.zeros_like(photons)]


def fluorescence_exact(p, window=None):
    """The expected pulse histogram for a two-state fluorophore. Instant, so there's no progress to report."""
    return FluorescenceExact(p)


_AFTER = {}  # the last after-pulse factorization. It doesn't depend on kExc or the pulse length, so a sweep of those
# factors once (about 15 s at n = 8)


def _solver(key, A):
    """A function solving A x = b: an LU when small, otherwise GMRES with an incomplete LU, as in exact()."""
    if _AFTER.get('key') != key:
        _AFTER.clear()
        if A.shape[0] <= 3 ** 8:
            _AFTER['solve'] = linalg.splu(A.tocsc()).solve
        else:
            ilu = linalg.LinearOperator(A.shape, linalg.spilu(A.tocsc(), drop_tol=1e-5, fill_factor=30).solve)

            def solve(b):
                x, info = linalg.gmres(A, b, M=ilu, rtol=1e-12, atol=0, restart=100, maxiter=50)
                if info:
                    raise RuntimeError(f'the after-pulse solve did not converge (gmres info {info})')
                return x
            _AFTER['solve'] = solve
        _AFTER['key'] = key
    return _AFTER['solve']


def pulse_totals(rules, p):
    """Every event one molecule gives over one whole pulse, on average, with no random noise and no cut-off.

    `rules` are build_rules(p), with p resolved. Every molecule starts in all G, with kExc on until p['duration'].
    Returns a Tally of events per molecule ([event, site]), for the intensity sweep.

    During the pulse: the matrix exponential of the master equation, with extra rows that add up each event as it
    fires (expm_multiply, which is accurate to rounding however small a count is, so it holds up at very low power).
    After the pulse: no more light, so one linear solve gives the expected time a molecule spends in each string from
    then on, and that times each rule's rate gives every later event. A molecule left with one triplet and nothing else
    can't make another photon, so it stops there and is counted as one triplet decayed, on the perylene (the triplet
    hops there long before it decays, since nothing hops off the perylene). Hops after the pulse aren't counted (nan).
    """
    sites = p['n'] + 2
    size = 3 ** sites
    moves, digit = _transitions(rules, sites)
    digit = np.array(digit)
    excited = (digit > 0).sum(0)
    lone_triplet = (excited == 1) & ((digit == 2).sum(0) == 0)
    active = (excited > 0) & ~lone_triplet  # can still emit with no more light

    counts = np.zeros((N_EVENTS, sites))
    x = np.zeros(size)
    x[0] = 1.0
    if p['duration'] > 0 and p['kExc'] > 0:
        rows, cols, vals = [], [], []
        for src, dst, rate, event, where in moves:
            rows += [dst, src, np.full(src.size, size + event * sites + where)]
            cols += [src, src, src]
            vals += [np.full(src.size, rate), np.full(src.size, -rate), np.full(src.size, rate)]
        extra = N_EVENTS * sites
        A = sparse.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                              shape=(size + extra, size + extra))
        start = np.zeros(size + extra)
        start[0] = 1.0
        end = linalg.expm_multiply(A * p['duration'], start)
        x, counts = end[:size], end[size:].reshape(N_EVENTS, sites).copy()

    after = [m for m in moves if m[3] != EXCITATION]
    key = (sites, tuple(rule for rule in rules if rule[7] != EXCITATION))
    keep = np.flatnonzero(active)
    if _AFTER.get('key') != key:
        rows = np.concatenate([np.r_[dst, src] for src, dst, *_ in after])
        cols = np.concatenate([np.r_[src, src] for src, *_ in after])
        vals = np.concatenate([np.r_[np.full(src.size, rate), np.full(src.size, -rate)] for src, _, rate, *_ in after])
        Q = sparse.csr_matrix((vals, (rows, cols)), shape=(size, size))
        _solver(key, -Q[keep][:, keep])
    time_in = np.zeros(size)  # ns spent in each string after the pulse, until all G or one triplet
    time_in[keep] = _AFTER['solve'](x[keep])
    ends_with_triplet = x[lone_triplet].sum()
    for src, dst, rate, event, where in after:
        counts[event, where] += rate * time_in[src].sum()
        ends_with_triplet += rate * time_in[src[lone_triplet[dst]]].sum()
    counts[EVENTS.index('triplet_lost'), -1] += ends_with_triplet
    counts[HOP] = np.nan
    return Tally(counts, 1.0, np.zeros((sites, 3)))


def summarize(tally):
    """The whole-system metrics, per excitation unless the name says otherwise."""
    per_event = {e: tally.counts[i] for i, e in enumerate(EVENTS)}
    excitations = per_event['excitation'].sum()
    isc = per_event['isc'].sum()
    lost_singlets = per_event['sens_singlet_lost'].sum()
    photons = per_event['photon']
    back = per_event['back_transfer'].sum()
    tta_s = per_event['tta_singlet']
    tta_f = per_event['tta_failed']
    decays = per_event['triplet_lost'].sum()
    x = 1 / excitations
    return {
        'photons_per_excitation': photons.sum() * x,
        'photons_per_excitation_se': tally.yield_se,
        'perylene_photons': photons[-1] * x,
        'dpa_photons': photons[:-1].sum() * x,
        'back_transfers': back * x,
        'isc_yield': isc / (isc + lost_singlets),
        'sens_singlets_lost': lost_singlets * x,
        'triplets_made': isc * x,
        'triplets_decayed': decays * x,
        'tta_singlet_perylene': tta_s[-1] * x,
        'tta_singlet_dpa': tta_s[:-1].sum() * x,
        'tta_failed_perylene': tta_f[-1] * x,
        'tta_failed_dpa': tta_f[:-1].sum() * x,
        'excitations_per_ns': excitations / tally.time,
        'photons_per_ns': photons.sum() / tally.time,
        'mean_triplets_present': tally.occupancy[:, 1].sum() / tally.time,
        'mean_singlets_present': tally.occupancy[:, 2].sum() / tally.time,
        'time_to_photon_ns': tally.latency,
        # bookkeeping: both should be 1 (a triplet made ends in TTA or decay; a singlet on the sensitizer ends in
        # ISC or decay). They only differ from 1 by whatever was still in flight when the run stopped
        'triplet_balance': (2 * tta_s.sum() + tta_f.sum() + decays) / isc if isc else np.nan,
        'singlet_balance': (isc + lost_singlets) / (excitations + back),
    }


def per_site(tally):
    """The same counts site by site: one dict per site, per excitation, plus the share of time in T and S."""
    x = 1 / tally.counts[EXCITATION].sum()
    return [{
        'photons': tally.counts[EVENTS.index('photon'), s] * x,
        'back_transfers': tally.counts[BACK, s] * x,
        'tta_singlet': tally.counts[TTA_SINGLET, s] * x,
        'tta_failed': tally.counts[TTA_FAILED, s] * x,
        'triplets_decayed': tally.counts[EVENTS.index('triplet_lost'), s] * x,
        'time_in_T': tally.occupancy[s, 1] / tally.time,
        'time_in_S': tally.occupancy[s, 2] / tally.time,
    } for s in range(tally.counts.shape[1])]
