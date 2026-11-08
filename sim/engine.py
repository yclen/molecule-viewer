"""Runs the rules from rules.py.

monte_carlo: a Gillespie simulation of the long-run behaviour. Add up the rates of every rule that can fire, wait a
random time, fire one.
exact: the same long-run answers by solving for the share of time in every possible string, with no random noise.
There are 3^(n + 2) strings, so this is for checking results at small n, not for sweeps.
pulse: many molecules, each starting in all G, with a square pulse of light. Records every photon's emission time.

All three give a Tally, so summarize() treats them alike.
"""

from dataclasses import dataclass

import numba as nb
import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg

EVENTS = ['excitation', 'isc', 'sens_singlet_lost', 'hop', 'triplet_lost', 'tta_singlet', 'tta_failed', 'photon',
          'back_transfer']
N_EVENTS = len(EVENTS)
EXCITATION, HOP, TTA_SINGLET, TTA_FAILED, PHOTON, BACK = (EVENTS.index(e) for e in (
    'excitation', 'hop', 'tta_singlet', 'tta_failed', 'photon', 'back_transfer'))


@dataclass
class Tally:
    counts: np.ndarray  # [event, site]
    time: float  # sim seconds these counts cover
    occupancy: np.ndarray  # [site, letter]: seconds spent
    latency: float = np.nan  # mean seconds from excitation to photon
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
    """Flash `molecules` molecules with the rules' kExc for `duration` seconds, and follow each until nothing is left.

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


def exact(rules, sites):
    """Long-run behaviour by linear algebra: solve pi Q = 0 over all 3^sites strings."""
    size = 3 ** sites
    index = np.arange(size)
    digit = [(index // 3 ** s) % 3 for s in range(sites)]
    rows, cols, vals, fires = [], [], [], []
    for a, b, fa, fb, ta, tb, rate, event, where in rules:
        can = digit[a] == fa
        if b >= 0:
            can = can & (digit[b] == fb)
        src = index[can]
        dst = src + (ta - fa) * 3 ** a + ((tb - fb) * 3 ** b if b >= 0 else 0)
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
        'excitations_per_s': excitations / tally.time,
        'photons_per_s': photons.sum() / tally.time,
        'mean_triplets_present': tally.occupancy[:, 1].sum() / tally.time,
        'mean_singlets_present': tally.occupancy[:, 2].sum() / tally.time,
        'time_to_photon_s': tally.latency,
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
