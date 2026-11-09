"""Figures for a pulse run: one per setting in plots/, plus plots/comparison.png when a run has several settings.

Every figure shows kHop and the quantum yield (photons emitted per photon absorbed, where absorbed means rule 1
excitations). The comparison uses the same axes in every panel, so the panels can be read against each other.
"""

import csv
import json

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

DPA, PERYLENE = '#2a78d6', '#11d8f2'  # perylene emits cyan blue, DPA a deeper blue
INK, MUTED, GRID, BASE = '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
MOST_PANELS = 6


def load(folder):
    """The run's settings, as a list of (summary row, histogram arrays), and the names it swept."""
    with (folder / 'summary.csv').open(newline='') as source:
        summaries = list(csv.DictReader(source))
    with (folder / 'histogram.csv').open(newline='') as source:
        bins = list(csv.DictReader(source))
    points = []
    for summary in summaries:
        rows = [row for row in bins if row['point'] == summary['point']]
        start = np.array([float(row['bin_start_s']) for row in rows])
        end = np.array([float(row['bin_end_s']) for row in rows])
        scale = 1 / (int(summary['molecules']) * (end - start))
        points.append((summary, {
            'time': (start + end) / 2,
            'width': end[0] - start[0],
            'dpa': np.array([int(row['dpa_photons']) for row in rows]) * scale,
            'perylene': np.array([int(row['perylene_photons']) for row in rows]) * scale,
        }))
    return points, json.loads((folder / 'params.json').read_text())['swept']


def quantum_yield(summary):
    return f"{float(summary['photons_per_excitation']):.1%}"


def settings_text(summary, width, leave_out=()):
    parts = {
        'kHop': f"kHop = {float(summary['kHop']):.3g} /s",
        'n': f"n = {summary['n']}",
        'kExc': f"kExc = {float(summary['kExc']):.3g} /s",
        'duration': f"pulse = {float(summary['duration']):g} s",
    }
    text = [value for name, value in parts.items() if name not in leave_out]
    return ', '.join(text + [f"{int(summary['molecules']):,} molecules", f'{width:g} s bins'])


def time_limit(summary, h):
    """Long enough for 99% of a setting's photons, and for the pulse."""
    total = np.cumsum(h['dpa'] + h['perylene'])
    latest = h['time'][np.searchsorted(total, 0.99 * total[-1])] if total[-1] > 0 else 0.0
    return 1.05 * max(latest, float(summary['duration']))


def share_past(h, xmax):
    """The share of a setting's photons that land beyond the end of the time axis."""
    total = h['dpa'] + h['perylene']
    return total[h['time'] > xmax].sum() / total.sum() if total.sum() else 0.0


def note_share_past(ax, h, xmax):
    share = share_past(h, xmax)
    if share >= 0.005:
        ax.text(0.98, 0.5, f'{share:.0%} of photons arrive\nafter this axis ends', transform=ax.transAxes, ha='right',
                color=MUTED, fontsize=10)


def draw(ax, summary, h, xmax, ymax):
    ax.axvspan(0, float(summary['duration']), color=MUTED, alpha=0.14, lw=0, label='Pulse on')
    ax.plot(h['time'], h['dpa'], color=DPA, lw=2, solid_capstyle='round', label='DPA photons')
    ax.plot(h['time'], h['perylene'], color=PERYLENE, lw=2, ls=(0, (5, 2.5)), label='Perylene photons')
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    ax.set_xlabel('Time after pulse starts (s)', color=MUTED)
    ax.grid(axis='y', color=GRID, lw=1)
    ax.set_axisbelow(True)
    for side in ('top', 'right', 'left'):
        ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_color(BASE)
    ax.tick_params(colors=MUTED, length=0)


def peak(points):
    return 1.1 * max(max(h['dpa'].max(), h['perylene'].max()) for _, h in points)


def plot_run(folder, limit=None):
    """Draw every figure. `limit` sets the end of the time axis in seconds."""
    points, swept = load(folder)
    out = folder / 'plots'
    out.mkdir(exist_ok=True)

    for summary, h in points:
        fig, ax = plt.subplots(figsize=(9, 5.4), dpi=160)
        xmax = limit or time_limit(summary, h)
        draw(ax, summary, h, xmax, peak([(summary, h)]))
        note_share_past(ax, h, xmax)
        ax.set_ylabel('Photons per molecule per second', color=MUTED)
        ax.legend(frameon=False, labelcolor=INK, loc='upper right')
        fig.suptitle(f'Quantum yield {quantum_yield(summary)} (photons emitted / photons absorbed)', x=0.07, ha='left',
                     color=INK, fontsize=13)
        ax.set_title(settings_text(summary, h['width']), color=MUTED, fontsize=10, loc='left')
        fig.tight_layout()
        fig.savefig(out / f"point_{int(summary['point']):04d}.png", facecolor='white')
        plt.close(fig)

    if 1 < len(points) <= MOST_PANELS:
        fig, axes = plt.subplots(1, len(points), figsize=(4.4 * len(points), 5), dpi=160, sharey=True)
        # the middle setting's limit: one slow setting's long tail would flatten every panel
        xmax = limit or float(np.median([time_limit(summary, h) for summary, h in points]))
        ymax = peak(points)
        for ax, (summary, h) in zip(axes, points):
            draw(ax, summary, h, xmax, ymax)
            note_share_past(ax, h, xmax)
            shown = [f"{name} = {float(summary[name]):g}" for name in swept if name != 'kHop']
            ax.set_title('\n'.join([f"kHop = {float(summary['kHop']):.3g} /s", *shown,
                                    f'Quantum yield {quantum_yield(summary)}']), color=INK, fontsize=12, loc='left')
        axes[0].set_ylabel('Photons per molecule per second', color=MUTED)
        axes[-1].legend(frameon=False, labelcolor=INK, loc='upper right')
        summary, h = points[0]
        fig.suptitle(settings_text(summary, h['width'], ['kHop', *swept]), x=0.01, ha='left', color=MUTED, fontsize=10)
        fig.tight_layout()
        fig.savefig(out / 'comparison.png', facecolor='white')
        plt.close(fig)
    return out
