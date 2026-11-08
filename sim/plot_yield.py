"""Quantum yield against excitation rate for a sweep: plots/quantum_yield.png, on linear and log scales.

The sweep has to vary excRatio. It can vary one other parameter too, which gets one line per value.
"""

import csv
import json
import sys

import numpy as np

from plot_pulse import BASE, GRID, INK, MUTED, plt

SHADES = ['#b5d4f4', '#85b7eb', '#378add', '#185fa5', '#042c53']  # one hue, light to dark, for an ordered parameter
MARKERS = ['o', 's', '^', 'D', 'v']
FEWEST_PHOTONS = 30  # below this a point is mostly counting noise, so it's left out
CONTEXT = ['hopRatio', 'kDecay', 'tauS', 'tauSens', 'tripletYield', 'f', 'n']


def style(ax):
    ax.set_xscale('log')
    ax.set_xlabel('excRatio: excitation rate relative to the sensitizer lifetime (log scale)', color=MUTED)
    ax.grid(axis='y', color=GRID, lw=1)
    ax.set_axisbelow(True)
    for side in ('top', 'right', 'left'):
        ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_color(BASE)
    ax.tick_params(colors=MUTED, length=0, which='both')


def plot_sweep(folder):
    with (folder / 'results.csv').open(newline='') as source:
        rows = list(csv.DictReader(source))
    params = json.loads((folder / 'params.json').read_text())
    swept = list(params['swept'])
    others = [name for name in swept if name != 'excRatio']
    if 'excRatio' not in swept or len(others) > 1:
        sys.exit('plot needs a sweep of excRatio, plus at most one other parameter')
    group = others[0] if others else None
    values = sorted({float(row[group]) for row in rows}) if group else [None]
    if len(values) > len(SHADES):
        sys.exit(f'plot can show at most {len(SHADES)} values of {group}')
    picks = np.linspace(0, len(SHADES) - 1, len(values)).round().astype(int) if group else [2]

    fig, (linear, logged) = plt.subplots(1, 2, figsize=(13, 5.4), dpi=160)
    lowest = np.inf
    for value, pick in zip(values, picks):
        line = [row for row in rows if group is None or float(row[group]) == value]
        line.sort(key=lambda row: float(row['excRatio']))
        line = [row for row in line
                if float(row['photons_per_excitation']) * params['excitations_per_point'] >= FEWEST_PHOTONS]
        if not line:
            continue
        x = [float(row['excRatio']) for row in line]
        y = [100 * float(row['photons_per_excitation']) for row in line]
        lowest = min(lowest, min(y))
        for ax in (linear, logged):
            ax.plot(x, y, color=SHADES[pick], lw=2, marker=MARKERS[pick], ms=5,
                    label=f'{group} = {value:g}' if group else 'Quantum yield')

    first = rows[0]
    limit = 100 * float(first['tripletYield']) * float(first['f']) / (1 + float(first['f']))
    top = limit
    if not {'tripletYield', 'f'} & set(swept):
        for ax in (linear, logged):
            ax.axhline(limit, color=MUTED, lw=1, ls=(0, (4, 3)))
        linear.text(0.01, limit * 1.02, f'Limit {limit:.1f}%: triplet yield × f/(1+f)', color=MUTED, fontsize=10,
                    transform=linear.get_yaxis_transform())
    else:
        top = max(100 * float(row['photons_per_excitation']) for row in rows)
    for ax in (linear, logged):
        style(ax)
    linear.set_ylim(0, top * 1.15)
    linear.set_ylabel('Quantum yield (%)', color=MUTED)
    logged.set_yscale('log')
    logged.set_ylim(lowest / 2, top * 2)
    logged.set_ylabel('Quantum yield (%), log scale', color=MUTED)
    logged.legend(frameon=False, labelcolor=INK, loc='lower right')

    context = ', '.join(f'{name} = {float(first[name]):g}' for name in CONTEXT if name not in swept)
    fig.suptitle(f'Continuous light: {context}', x=0.01, ha='left', color=INK, fontsize=12)
    fig.tight_layout()
    out = folder / 'plots'
    out.mkdir(exist_ok=True)
    fig.savefig(out / 'quantum_yield.png', facecolor='white')
    plt.close(fig)
    return out
