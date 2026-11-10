"""Figures for intensity runs: photons per molecule and quantum yield against kExc.

`draw(fig, kind, folders)` fills a matplotlib Figure, so sweep_gui.py can show it in a window and `plot` can save it.
Several runs can be drawn together. Each run gives one curve per value of whatever it swept besides kExc.

Exact rows are drawn as lines, Monte Carlo rows as points with their standard error. Quantum yield is
photons / excitations, worked out here from the saved counts.

With `per_mw` (kExc per mW of laser power), the x axis is laser power in mW instead of kExc, and the top axis
shows kExc. sweep_gui.py offers the live window's calibration, live.kexc_per_mw(), by default.
"""

import csv
import json
from pathlib import Path

import numpy as np

INK, MUTED, GRID, BASE = '#52514e', '#898781', '#e1e0d9', '#c3c2b7'  # the same as plot_pulse.py
DPA, PERYLENE = '#2a78d6', '#11d8f2'
SHADES = ['#85b7eb', '#5598e7', '#2a78d6', '#185fa5', '#042c53']  # one hue, light to dark: one ordered parameter
CATEGORICAL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4']  # anything else, always in this order
MARKERS = ['o', 's', '^', 'D', 'v']
MOST_CURVES = 5
KINDS = {
    'counts': 'Photons per molecule, linear',
    'counts_log': 'Photons per molecule, log, with slope',
    'yield': 'Quantum yield, linear',
    'yield_log': 'Quantum yield, log kExc',
}
FILES = {'counts': 'photons_linear.png', 'counts_log': 'photons_log.png', 'yield': 'quantum_yield.png',
         'yield_log': 'quantum_yield_log_kexc.png'}
LOG_KEXC = {'counts_log', 'yield_log'}  # the figures whose kExc axis is always on a log scale
CONTEXT = ['n', 'kHop', 'f', 'tripletYield', 'tauS', 'tauSens', 'kDecay', 'fretRadius', 'duration']
ALWAYS_TOGETHER = {'kExc', 'excRatio'}  # excRatio is kExc × tauSens, so it moves with kExc


def is_intensity(folder):
    try:
        return json.loads((Path(folder) / 'params.json').read_text())['mode'] == 'intensity'
    except (OSError, ValueError, KeyError):
        return False


def load(folder):
    """The run's rows, with numbers as floats, and the names it swept besides kExc."""
    folder = Path(folder)
    with (folder / 'results.csv').open(newline='') as source:
        rows = [{k: (v if k == 'method' else float(v) if v not in ('', 'nan') else np.nan) for k, v in row.items()}
                for row in csv.DictReader(source)]
    swept = [name for name in json.loads((folder / 'params.json').read_text())['swept'] if name not in ALWAYS_TOGETHER]
    return rows, swept


def curves(folders):
    """[{label, exact, mc}] with the rows of each sorted by kExc, plus the colours to draw them in."""
    found = []
    for folder in folders:
        rows, swept = load(folder)
        keys = sorted({tuple(row[name] for name in swept) for row in rows})
        for key in keys:
            mine = [row for row in rows if tuple(row[name] for name in swept) == key]
            parts = [f'{name} = {value:g}' for name, value in zip(swept, key)]
            if len(folders) > 1:
                parts.insert(0, Path(folder).name)
            found.append({
                'label': ', '.join(parts),  # empty for a run with one curve
                'exact': sorted((r for r in mine if r['method'] == 'exact'), key=lambda r: r['kExc']),
                'mc': sorted((r for r in mine if r['method'] == 'mc'), key=lambda r: r['kExc']),
                'rows': mine,
            })
        ordered = len(folders) == 1 and len(swept) == 1
    if not found:
        raise ValueError('no rows to plot yet')
    if len(found) > MOST_CURVES:
        raise ValueError(f'{len(found)} curves: at most {MOST_CURVES} can be told apart. Pick fewer runs')
    if len(found) == 1:
        colours = [CATEGORICAL[0]]
    elif ordered:
        colours = [SHADES[k] for k in np.linspace(0, len(SHADES) - 1, len(found)).round().astype(int)]
    else:
        colours = CATEGORICAL[:len(found)]
    for curve, colour, marker in zip(found, colours, MARKERS):
        curve['colour'], curve['marker'] = colour, marker
    return found


def column(rows, name):
    return np.array([row[name] for row in rows])


def yield_se(rows):
    """Standard error of photons / excitations for Monte Carlo rows, from the saved sums of squares."""
    N = column(rows, 'molecules')
    P, E = column(rows, 'photons_per_molecule'), column(rows, 'excitations_per_molecule')
    var_p = column(rows, 'mean_photons_sq') - P ** 2
    var_e = column(rows, 'mean_excitations_sq') - E ** 2
    cov = column(rows, 'mean_photons_x_excitations') - P * E
    Q = P / E
    return np.sqrt(np.maximum(var_p - 2 * Q * cov + Q ** 2 * var_e, 0) / N) / E


def photons_se(rows):
    N = column(rows, 'molecules')
    P = column(rows, 'photons_per_molecule')
    return np.sqrt(np.maximum(column(rows, 'mean_photons_sq') - P ** 2, 0) / N)


def local_slope(rows):
    """d log(photons) / d log(kExc): 2 while a photon needs two absorptions in one pulse, falling to 0 at saturation."""
    x, y = np.log(column(rows, 'kExc')), np.log(column(rows, 'photons_per_molecule'))
    ok = np.isfinite(y)
    if ok.sum() < 3:
        return None, None
    return np.exp(x[ok]), np.gradient(y[ok], x[ok])


def style(ax, logx, logy=False, xlabel='kExc, /ns'):
    ax.set_xscale('log' if logx else 'linear')
    ax.set_yscale('log' if logy else 'linear')
    if xlabel:
        ax.set_xlabel(f"{xlabel}{' (log scale)' if logx else ''}", color=MUTED)
    ax.grid(color=GRID, lw=1, which='major')
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('bottom', 'left'):
        ax.spines[side].set_color(BASE)
    ax.tick_params(colors=MUTED, which='both')


def shared(found, name):
    values = {row[name] for curve in found for row in curve['rows']}
    return values.pop() if len(values) == 1 else None


def plot_curve(ax, curve, y_of, se_of=None, split=False, scale=1.0):
    """Exact rows as a line with small markers, Monte Carlo rows as open markers with error bars. Returns the points.

    The x axis is kExc × scale: 1 for kExc itself, 1 / (kExc per mW) for laser power.
    """
    points = []
    for kind in ('exact', 'mc'):
        rows = curve[kind]
        if not rows:
            continue
        x, y = column(rows, 'kExc') * scale, y_of(rows)
        method = 'Exact solve' if kind == 'exact' else 'Monte Carlo'
        if not curve['label']:
            label = method
        else:
            label = curve['label'] + ('' if kind == 'exact' or not curve['exact'] else ', Monte Carlo')
        if kind == 'exact':
            ax.plot(x, y, color=curve['colour'], lw=2, marker=curve['marker'], ms=4, label=label, zorder=3)
        else:
            err = se_of(rows) if se_of else None
            ax.errorbar(x, y, yerr=err, color=curve['colour'], lw=0 if curve['exact'] else 2, elinewidth=1.2, capsize=3,
                        marker=curve['marker'], ms=8, mfc='white', mew=1.6, label=label, zorder=4)
        points += [(xi, yi, f"{curve['label'] or method}{' (Monte Carlo)' if curve['label'] and kind == 'mc' else ''}")
                   for xi, yi in zip(x, y)]
        if split and kind == 'exact':
            ax.plot(x, column(rows, 'perylene_photons_per_molecule'), color=PERYLENE, lw=1.5, ls=(0, (4, 2)),
                    label='from the perylene', zorder=2)
            ax.plot(x, column(rows, 'dpa_photons_per_molecule'), color=DPA, lw=1.5, ls=(0, (1, 2)),
                    label='from the DPAs', zorder=2)
    return points


def photons(rows):
    return column(rows, 'photons_per_molecule')


def quantum_yield(rows):
    return column(rows, 'photons_per_molecule') / column(rows, 'excitations_per_molecule')


def yield_percent(rows):
    return 100 * quantum_yield(rows)


def yield_percent_se(rows):
    return 100 * yield_se(rows)


def draw(fig, kind, folders, logx=None, per_mw=None):
    """Draw one figure of `kind` (a key of KINDS) for the runs in `folders`. Returns the points, for hover readouts.

    logx: the x axis on a log scale. The figures in LOG_KEXC always have it.
    per_mw: kExc per mW. Given, the x axis is laser power in mW = kExc / per_mw, and the top axis shows kExc.
    """
    found = curves(folders)
    logx = bool(logx) or kind in LOG_KEXC
    split = len(found) == 1
    fig.clear()
    fig.hover_unit = ''  # what the hover readout puts after a value
    scale = 1 / per_mw if per_mw else 1.0  # x axis value per unit of kExc
    xlabel = 'Laser power, mW' if per_mw else 'kExc, /ns'
    fig.hover_x = ('Power', 'mW') if per_mw else ('kExc', '/ns')
    fig.set_facecolor('white')
    points = {}
    if kind == 'counts_log':
        top, low = fig.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1.3]})
        points[top] = sum((plot_curve(top, c, photons, photons_se, split, scale) for c in found), [])
        style(top, True, True, xlabel=False)
        top.set_ylabel('Photons per molecule (log scale)', color=MUTED)
        for c in found:
            x, slope = local_slope(c['exact'] or c['mc'])
            if x is not None:
                x = x * scale
                low.plot(x, slope, color=c['colour'], lw=2, marker=c['marker'], ms=4)
                points.setdefault(low, []).extend((xi, si, f"{c['label'] or 'Local'} slope") for xi, si in zip(x, slope))
        for level, text in ((2, '2: needs two absorptions'), (1, '1'), (0, '0: saturated')):
            low.axhline(level, color=MUTED, lw=1, ls=(0, (4, 3)), zorder=1)
            low.text(1.005, level, text, color=MUTED, fontsize=8, va='center', transform=low.get_yaxis_transform())
        style(low, True, xlabel=xlabel)
        low.set_ylim(-0.15, 2.3)
        low.set_ylabel('Local slope', color=MUTED)
        legend_ax = top
    else:
        ax = fig.subplots()
        if kind == 'counts':
            points[ax] = sum((plot_curve(ax, c, photons, photons_se, split, scale) for c in found), [])
            ax.set_ylabel('Photons per molecule, per pulse', color=MUTED)
        else:
            points[ax] = sum((plot_curve(ax, c, yield_percent, yield_percent_se, scale=scale) for c in found), [])
            ax.set_ylabel('Quantum yield (%), photons per excitation', color=MUTED)
            fig.hover_unit = '%'
            triplet, f = shared(found, 'tripletYield'), shared(found, 'f')
            if triplet is not None and f is not None:
                limit = 100 * triplet * f / (1 + f)
                ax.axhline(limit, color=MUTED, lw=1, ls=(0, (4, 3)))
                ax.text(0.99, limit, f'Limit {limit:.3g}%: triplet yield × f/(1+f) ', color=MUTED, fontsize=9,
                        ha='right', va='bottom', transform=ax.get_yaxis_transform())
                ax.set_ylim(top=limit * 1.12)
        style(ax, logx, xlabel=xlabel)
        if not logx:
            ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        legend_ax = ax

    duration = shared(found, 'duration')
    first = fig.axes[0]
    top_axis = None
    if per_mw:  # power along the bottom, and the kExc it gives along the top
        top_axis = first.secondary_xaxis('top', functions=(lambda x: x * per_mw, lambda k: k / per_mw))
        top_axis.set_xlabel(f'kExc, /ns, at {per_mw:g} /ns per mW', color=MUTED)
    elif duration:  # the same pulse in every curve, so kExc also reads as absorptions per pulse
        top_axis = first.secondary_xaxis('top', functions=(lambda k: k * duration, lambda a: a / duration))
        top_axis.set_xlabel('kExc × pulse: absorptions per pulse if the sensitizer were always free', color=MUTED)
    if top_axis is not None:
        top_axis.tick_params(colors=MUTED, which='both')
        top_axis.spines['top'].set_color(BASE)
    handles, labels = legend_ax.get_legend_handles_labels()
    if len(handles) > 1:
        legend_ax.legend(frameon=False, labelcolor=INK, fontsize=9,
                         loc='lower right')
    context = ', '.join(f'{name} {value:g}' for name in CONTEXT if (value := shared(found, name)) is not None)
    fig.suptitle(f'{KINDS[kind]}\n{context}', x=0.01, ha='left', color=INK, fontsize=10.5, linespacing=1.5)
    fig.tight_layout()
    return points


def add_hover(fig, points):
    """Show the nearest point's values when the mouse is over a plot."""
    notes = {}
    for ax in points:
        note = ax.annotate('', (0, 0), xytext=(12, 12), textcoords='offset points', fontsize=9, color=INK,
                           bbox={'boxstyle': 'round,pad=0.4', 'fc': 'white', 'ec': BASE}, zorder=10)
        note.set_visible(False)
        notes[ax] = note

    def move(event):
        changed = False
        for ax, note in notes.items():
            show = False
            if event.inaxes is ax and points[ax]:
                where = ax.transData.transform([(x, y) for x, y, _ in points[ax]])
                gaps = np.hypot(where[:, 0] - event.x, where[:, 1] - event.y)
                k = int(np.nanargmin(gaps))
                if gaps[k] < 25:
                    x, y, label = points[ax][k]
                    note.xy = (x, y)
                    name, unit = getattr(fig, 'hover_x', ('kExc', '/ns'))
                    note.set_text(f"{label}\n{name} {x:.4g} {unit}\n{y:.5g}{getattr(fig, 'hover_unit', '')}")
                    show = True
            if note.get_visible() != show or show:
                note.set_visible(show)
                changed = True
        if changed:
            fig.canvas.draw_idle()

    fig.canvas.mpl_connect('motion_notify_event', move)


def save_all(folder, per_mw=None):
    """The figures in FILES, in plots/ of one run. Returns the plots folder.

    With per_mw, against laser power instead, in files ending _power.png, so the kExc ones stay.
    """
    from matplotlib.figure import Figure
    out = Path(folder) / 'plots'
    out.mkdir(exist_ok=True)
    for kind, name in FILES.items():
        fig = Figure(figsize=(10, 6.5 if kind == 'counts_log' else 5.6), dpi=150)
        draw(fig, kind, [folder], per_mw=per_mw)
        fig.savefig(out / (name.replace('.png', '_power.png') if per_mw else name), facecolor='white')
    return out
