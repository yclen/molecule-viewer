"""A live pulse experiment: fresh molecules are flashed all the time, and each photon lands in a histogram on screen.

  python sim/live.py

Like a time-resolved measurement, the histogram starts as a few scattered counts and smooths out as they build up. Each
molecule gets one pulse from the same code as `tta_sim.py pulse`, so a live run is that command in small batches.

Everything in the window is in nanoseconds, like the rest of the simulation: times in ns, rates per ns. Only the
wall-clock figures (molecules per s, how long a solve takes) are real seconds.

There are no ratios here. The hop rate is set directly, and so is the excitation rate, kExc, with the laser slider. kExc
is proportional to laser power, so the slider plays the part of the power knob.

The simulation runs in its own process, on threads that never wait for each other. Drawing the window would hold it up
otherwise, and so would the slower cores on a machine that has a mix of fast and slow ones.

"Exact solve" works out the histogram the counts are heading for, with engine.pulse_exact, in a second process. It's
drawn over the counts once it's done, scaled to the molecules counted so far.
"""

import csv
import json
import multiprocessing
import os
import shutil
import sys
import threading
import time
from collections import deque
from datetime import datetime

from PyQt6 import QtCore, QtGui, QtWidgets  # before pyqtgraph, so it picks this Qt and not any other that's installed

import numpy as np
import pyqtgraph as pg

from engine import QUIET, exact_seconds, fluorescence_exact, pulse_alone, pulse_exact, rule_arrays
from plot_pulse import DPA, PERYLENE
from rules import FLUORESCENCE_DEFAULTS, PULSE_DEFAULTS, build_rules, fluorescence_rules, resolve
from tta_sim import BINS_PER_LIFETIME, HERE, add_to_index, new_run_folder

# the values typed into the window: label, unit, lowest and highest allowed value
FIELDS = {
    'duration': ('Pulse length', 'ns', 1e-3, 1e6),
    'window': ('Time window', 'ns', 1e-2, 1e7),
    'n': ('DPA mediators, n', '', 0, 40),
    'kHop': ('Hop rate, kHop', '/ns', 1e-6, 1e7),
    'f': ('Spin factor, f', '', 0, 1),
    'tripletYield': ('Triplet yield', '', 0, 0.999),
    'tauS': ('Singlet lifetime, τS', 'ns', 1e-3, 1e6),
    'tauSens': ('Sensitizer lifetime, τSens', 'ns', 1e-3, 1e6),
    'kDecay': ('Triplet decay, kDecay', '/ns', 0, 1e3),
    'fretRadius': ('FRET radius, R0', 'nm', 0, 100),
    'tauFl': ('Fluorescence lifetime, τF', 'ns', 1e-3, 1e6),
    'fluorYield': ('Fluorescence yield, Φ_F', '', 0, 1),
}
SAMPLES = {  # the samples on offer, and the values each one uses. These are changed less often than the laser
    'tta': ('n', 'kHop', 'f', 'tripletYield', 'tauS', 'tauSens', 'kDecay', 'fretRadius'),
    'fluorescence': ('tauFl', 'fluorYield'),
}
KEXC_LOW, KEXC_HIGH = 1e-4, 1e2  # /ns, the ends of the laser slider, which is on a log scale
SLIDER_STEPS = 600  # 100 per decade
APPLY_MS = 300  # a slider change from the keyboard or mouse wheel waits this long for the next one before it clears

# the assumed setup, used only to show the laser power that gives the kExc set. Change these to match the experiment
BEAM_RADIUS_UM = 1.0  # Gaussian beam radius at the sample (1/e² of the intensity), µm
EXTINCTION = 50_000  # the sensitizer's molar extinction coefficient at the laser line, /M/cm
WAVELENGTH_NM = 635  # the laser line, nm
AVOGADRO, PLANCK, LIGHT = 6.02214076e23, 6.62607015e-34, 2.99792458e8


def kexc_per_mw():
    """The excitation rate (/ns) per mW of laser power, for a molecule at the centre of the beam.

    kExc = σ I0 / (h c / λ), with the absorption cross-section σ = ln(10) 1000 ε / N_A (cm²) and the peak intensity of a
    Gaussian beam I0 = 2 P / (π w0²). With the values above it's 0.0389 /ns per mW.
    """
    cross_section = np.log(10) * 1000 * EXTINCTION / AVOGADRO  # cm²
    peak = 2 * 1e-3 / (np.pi * (BEAM_RADIUS_UM * 1e-4) ** 2)  # W/cm² from 1 mW
    photon = PLANCK * LIGHT / (WAVELENGTH_NM * 1e-9)  # J
    return cross_section * peak / photon * 1e-9  # per ns


def power_text(mw):
    """A laser power in the unit that reads best: nW, µW, mW or W."""
    for unit, size in (('W', 1e3), ('mW', 1), ('µW', 1e-3)):
        if mw >= size:
            return f'{mw / size:.3g} {unit}'
    return f'{mw * 1e6:.3g} nW'
TARGET = 100_000  # molecules per second to aim for. 0 means as fast as the machine can go
WINDOW_LIFETIMES = 10  # default time window, in singlet lifetimes
CELLS_PER_LIFETIME = BINS_PER_LIFETIME * 10  # counts are kept this fine, so bins can change without losing anything
MOST_CELLS = 200_000
BIN_CELLS = (1, 2, 5, 10, 20, 50, 100)  # the bin widths on offer, in cells
DEFAULT_BIN_CELLS = 10  # a tenth of the singlet lifetime, the same default as a pulse run

THREADS = max(1, (os.cpu_count() or 2) // 2)  # more than this gains little, and the window gets jerky
BATCH_S = 0.02  # how long a thread's batch aims to take, so Stop and new parameters take effect quickly
MOST_MOLECULES = 2000  # in one batch
CATCH_UP_S = 0.25  # a thread that fell behind can make up this much of the schedule, then the rest is lost
FRAME_MS = 50
SEND_S = 0.05  # how often finished batches go to the window
FRAMES_PER_TEXT = 5  # the numbers change less often than the histogram, so they can be read
RATE_OVER_S = 2.0  # the rates on screen are averages over this long
EMPTY_BIN = 0.5  # where a bin with no counts sits on a log axis, which can't show 0
PAST_WARNING = 0.01  # the Photons readout warns once this share of photons lands after the time window
EXACT_MOST_N = 9  # 3^11 strings. Beyond this the exact solve doesn't fit in memory
EXACT_SLOW_S = 60  # an exact solve estimated to take longer than this is called slow
EXACT_MOST_S = 600  # and one estimated longer than this isn't started
EXACT_SLOWDOWN = 1.6  # the counting slows the exact solve this much, measured at n = 8
EXACT_PROGRESS_S = 0.1  # how often the exact solve reports how far it's got
HINT = 'Hover over the histogram to read a bin.'
STARTING = 'Starting the engine. The first run after a change to the code takes longer.'

# the instrument panel's colours. Perylene and DPA keep their colours from plot_pulse
BG, PANEL, EDGE, FIELD = '#0d1117', '#151b23', '#262e39', '#0a0e13'
TEXT, MUTED, DIM = '#d6dde6', '#7c8796', '#4a5361'
LASER = '#ff2100'  # 638 nm, the same red as the pulse band
GO, AMBER = '#2fd27a', '#f2a93b'
PERYLENE_EXACT, DPA_EXACT = '#b5f4fc', '#9cc3f5'  # lighter versions of the two colours, drawn dashed
FLUOR, FLUOR_EXACT = '#c792ff', '#e6d0ff'  # the fluorescence reference: violet, so it never passes for a TTA trace
PILES = {'tta': ('perylene', 'dpa'), 'fluorescence': ('fluorescence',)}  # each sample's photons, as Measurement keeps them
MONO = "'Cascadia Mono', Consolas, monospace"


class Job:
    """What the simulation threads share: the current settings, the schedule, and photons waiting to be sent."""

    def __init__(self):
        self.lock = threading.Lock()
        self.on = threading.Event()
        self.target = TARGET
        self.generation = 0
        self.seed = 1
        self.next_start = 0.0
        self.empty_outbox()

    def empty_outbox(self):
        self.molecules = self.excitations = 0
        self.times, self.sites = [], []

    def change(self, message):
        with self.lock:
            self.target = message.get('target', self.target)
            if 'p' in message:
                self.p = message['p']
                self.arrays = rule_arrays(rules_for(self.p))
                self.generation = message['generation']
                self.empty_outbox()
        if 'on' in message:
            self.on.set() if message['on'] else self.on.clear()

    def next_batch(self, speed, seen):
        """Claim a batch: its settings, its size, and when it's due to start.

        `speed` is what this thread managed on the settings `seen`. On new settings it starts from one molecule, because
        a batch sized for fast settings could take many seconds on slow ones.
        """
        with self.lock:
            size = int(np.clip(BATCH_S * speed if seen == self.generation else 1, 1, MOST_MOLECULES))
            if self.target:
                size = int(np.clip(size, 1, max(1, BATCH_S * self.target)))
            now = time.perf_counter()
            start = self.next_start
            if self.target:  # taken from the schedule, not from now, so a late wake-up doesn't lower the rate
                self.next_start = max(start, now - CATCH_UP_S) + size / self.target
            self.seed += 1
            return self.generation, self.arrays, sites_for(self.p), self.p['duration'], size, self.seed, start

    def run(self):
        speed, seen = 0.0, None  # molecules per second this thread managed on the settings `seen`
        while True:
            self.on.wait()
            generation, arrays, sites, duration, size, seed, start = self.next_batch(speed, seen)
            seen = generation
            time.sleep(max(0.0, start - time.perf_counter()))
            began = time.perf_counter()
            times, where, excitations = pulse_alone(arrays, sites, size, duration, seed)
            speed = size / (time.perf_counter() - began)
            with self.lock:
                if generation == self.generation:
                    self.molecules += size
                    self.excitations += excitations
                    self.times.append(times)
                    self.sites.append(where)

    def collect(self):
        with self.lock:
            batch = (self.generation, self.molecules, self.excitations, self.times, self.sites)
            self.empty_outbox()
        return batch


def acquire(link):
    """The simulation process. The window sends it dicts of changes: `on`, `target`, `close`, and `p` with `generation`."""
    job = Job()
    job.change(link.recv())  # the first parameters, which the threads need before they start
    for _ in range(THREADS):
        threading.Thread(target=job.run, daemon=True).start()
    while True:
        if link.poll(SEND_S):
            message = link.recv()
            if message.get('close'):
                return
            job.change(message)
        generation, molecules, excitations, times, sites = job.collect()
        if molecules:
            link.send((generation, molecules, excitations, np.concatenate(times), np.concatenate(sites)))


def solve_exactly(link):
    """An exact solve's process. Gets (key, p) to solve, sends ('progress', key, share, doing) as it goes, then
    ('done', key, result).

    It can't be interrupted, so to stop a solve the window ends the whole process and uses a spare one.
    """
    while True:
        key, p = link.recv()
        last = 0.0

        def progress(share, doing):
            nonlocal last
            if time.perf_counter() - last > EXACT_PROGRESS_S:
                last = time.perf_counter()
                link.send(('progress', key, share, doing))

        try:
            link.send(('done', key, exact_for(p, progress)))
        except Exception as error:  # shown in the window, which keeps running
            link.send(('failed', key, str(error)))


def exact_verdict(p):
    """Whether to offer the exact solve for p, and what to say about how long it takes."""
    if fluorescent(p):
        return True, 'Instant: a two-state molecule has a closed-form answer.'
    n = p['n']
    if n > EXACT_MOST_N:
        return False, f'n = {n} is too large to solve exactly: it has {3 ** (n + 2):,} strings.'
    seconds = exact_seconds(p, p['window']) * EXACT_SLOWDOWN
    if seconds > EXACT_MOST_S:
        return False, f'Too slow to solve exactly at these settings: about {how_long(seconds)}.'
    return True, f"Takes about {how_long(seconds)}{', which is slow' if seconds > EXACT_SLOW_S else ''}."


def how_long(seconds):
    if seconds < 1:
        return 'a second'
    if seconds < 90:
        return f'{seconds:.0f} s'
    if seconds < 90 * 60:
        return f'{seconds / 60:.0f} min'
    return f'{seconds / 3600:.0f} h'


def default_values():
    """What the window starts with: the pulse defaults, with the ratios turned into flat rates per ns."""
    p = resolve(PULSE_DEFAULTS)
    values = {name: p[name] for name in SAMPLES['tta'] + ('duration',)}
    values.update(FLUORESCENCE_DEFAULTS)
    values['window'] = WINDOW_LIFETIMES * p['tauS']
    values['kExc'] = p['kExc']
    values['sample'] = 'tta'
    return values


def parameters(values):
    """The engine's parameters for the window's values: the laser, the time window, and the chosen sample's values.

    kHop and kExc are given directly, so no ratio is used.
    """
    shared = {name: values[name] for name in ('sample', 'kExc', 'duration', 'window')}
    if values['sample'] == 'fluorescence':
        return {**shared, **{name: values[name] for name in SAMPLES['fluorescence']}}
    p = resolve({**PULSE_DEFAULTS, **shared, **{name: values[name] for name in SAMPLES['tta']}})  # works out kISC
    del p['excRatio'], p['hopRatio']
    return p


def fluorescent(p):
    return p['sample'] == 'fluorescence'


def rules_for(p):
    return fluorescence_rules(p) if fluorescent(p) else build_rules(p)


def sites_for(p):
    return 1 if fluorescent(p) else p['n'] + 2


def exact_for(p, progress=None):
    """The exact pulse histogram for p: in closed form for the fluorophore, solved step by step for TTA."""
    return fluorescence_exact(p) if fluorescent(p) else pulse_exact(build_rules(p), p, p['window'], progress)


def rate_to_step(kexc):
    return int(round(SLIDER_STEPS * np.log10(kexc / KEXC_LOW) / np.log10(KEXC_HIGH / KEXC_LOW)))


def step_to_rate(step):
    return float(f'{KEXC_LOW * (KEXC_HIGH / KEXC_LOW) ** (step / SLIDER_STEPS):.3g}')


def clock_text(seconds):
    seconds = int(seconds)
    return f'{seconds // 3600:02d}:{seconds // 60 % 60:02d}:{seconds % 60:02d}'


class Measurement:
    """Everything counted for one set of parameters.

    Photons are kept in two piles: `main`, from the emitter (the perylene, or the fluorophore), and `other`, from the
    DPAs (always empty for the fluorophore).
    """

    def __init__(self, p):
        self.p = p
        self.main_site = 0 if fluorescent(p) else p['n'] + 1
        lifetime = p['tauFl'] if fluorescent(p) else p['tauS']
        self.cell = max(lifetime / CELLS_PER_LIFETIME, p['window'] / MOST_CELLS)
        self.cells = int(np.ceil(p['window'] / self.cell))
        self.main = np.zeros(self.cells, np.int64)
        self.other = np.zeros(self.cells, np.int64)
        self.totals = {'molecules': 0, 'excitations': 0, 'photons': 0, 'past_window': 0}

    def add(self, molecules, excitations, times, sites):
        cell = (times / self.cell).astype(np.int64)
        inside = cell < self.cells
        from_main = sites == self.main_site
        self.main += np.bincount(cell[inside & from_main], minlength=self.cells)
        self.other += np.bincount(cell[inside & ~from_main], minlength=self.cells)
        self.totals['molecules'] += molecules
        self.totals['excitations'] += excitations
        self.totals['photons'] += times.size
        self.totals['past_window'] += int((~inside).sum())

    def binned(self, cells_per_bin):
        """The bin edges, then the main and other photons in each bin."""
        spare = -self.cells % cells_per_bin
        main = np.pad(self.main, (0, spare)).reshape(-1, cells_per_bin).sum(1)
        other = np.pad(self.other, (0, spare)).reshape(-1, cells_per_bin).sum(1)
        return self.cell * cells_per_bin * np.arange(main.size + 1), main, other


def label(text='', role=None):
    widget = QtWidgets.QLabel(text)
    if role:
        widget.setProperty('role', role)
    return widget


def panel(title=None, name='panel'):
    """A framed section of the control panel, with its title in small capitals."""
    frame = QtWidgets.QFrame()
    frame.setObjectName(name)
    layout = QtWidgets.QVBoxLayout(frame)
    layout.setContentsMargins(14, 12, 14, 14)
    layout.setSpacing(8)
    if title:
        layout.addWidget(label(title, 'section'))
    return frame, layout


class Lamp(QtWidgets.QLabel):
    """A status light with its name beside it."""

    def __init__(self, name):
        super().__init__()
        self.setProperty('role', 'lamp')
        self.set(None, name)

    def set(self, colour, text):
        self.setText(f"<span style='color:{colour or DIM}; font-size:16px'>●</span>&nbsp;"
                     f"<span style='color:{TEXT if colour else MUTED}'>{text}</span>")


def style(app):
    """A dark instrument-panel look: Fusion with a dark palette, then the accents from a style sheet."""
    app.setStyle('Fusion')
    palette = QtGui.QPalette()
    roles = QtGui.QPalette.ColorRole
    for role, colour in ((roles.Window, BG), (roles.WindowText, TEXT), (roles.Base, FIELD), (roles.AlternateBase, PANEL),
                         (roles.Text, TEXT), (roles.Button, PANEL), (roles.ButtonText, TEXT), (roles.Highlight, DPA),
                         (roles.HighlightedText, '#ffffff'), (roles.ToolTipBase, PANEL), (roles.ToolTipText, TEXT),
                         (roles.PlaceholderText, DIM)):
        palette.setColor(role, QtGui.QColor(colour))
    for role in (roles.Text, roles.WindowText, roles.ButtonText):
        palette.setColor(QtGui.QPalette.ColorGroup.Disabled, role, QtGui.QColor(DIM))
    app.setPalette(palette)
    app.setStyleSheet(f"""
        QWidget {{ color: {TEXT}; font-family: 'Segoe UI'; font-size: 13px; }}
        QMainWindow, QWidget#root, QWidget#column {{ background: {BG}; }}
        QToolTip {{ background: {PANEL}; color: {TEXT}; border: 1px solid {EDGE}; padding: 5px; }}
        QFrame#panel, QFrame#tile {{ background: {PANEL}; border: 1px solid {EDGE}; border-radius: 8px; }}
        QFrame#statusbar {{ background: {PANEL}; border: 1px solid {EDGE}; border-radius: 6px; }}
        QLabel {{ background: transparent; }}
        QLabel[role="title"] {{ font-size: 15px; font-weight: 600; letter-spacing: 3px; }}
        QLabel[role="subtitle"] {{ color: {MUTED}; font-size: 12px; letter-spacing: 1px; }}
        QLabel[role="section"] {{ color: {MUTED}; font-size: 11px; font-weight: 600; letter-spacing: 2px; }}
        QLabel[role="caption"] {{ color: {MUTED}; font-size: 12px; }}
        QLabel[role="note"] {{ color: {DIM}; font-size: 11px; }}
        QLabel[role="lamp"] {{ font-size: 11px; font-weight: 600; letter-spacing: 2px; }}
        QLabel[role="unit"] {{ color: {MUTED}; font-family: {MONO}; font-size: 12px; }}
        QLabel[role="value"] {{ font-family: {MONO}; font-size: 20px; }}
        QLabel[role="hero"] {{ font-family: {MONO}; font-size: 26px; color: {PERYLENE}; }}
        QLabel[role="big"] {{ font-family: {MONO}; font-size: 26px; }}
        QLabel[role="warn"] {{ font-family: {MONO}; color: {AMBER}; font-size: 12px; }}
        QLabel[role="sub"] {{ font-family: {MONO}; color: {MUTED}; font-size: 12px; }}
        QLabel[role="readout"] {{ font-family: {MONO}; color: {TEXT}; font-size: 12px; }}
        QLabel[role="clock"] {{ font-family: {MONO}; font-size: 15px; letter-spacing: 1px; }}
        QLabel[role="status"] {{ font-family: {MONO}; color: {MUTED}; font-size: 12px; }}
        QLineEdit, QComboBox {{ background: {FIELD}; border: 1px solid {EDGE}; border-radius: 4px; padding: 4px 6px;
                               font-family: {MONO}; selection-background-color: {DPA}; }}
        QLineEdit:focus, QComboBox:focus {{ border-color: {PERYLENE}; }}
        QLineEdit#power {{ font-size: 30px; color: #ff6a4a; background: transparent; border: none;
                           border-bottom: 1px solid {EDGE}; border-radius: 0; padding: 0 2px; }}
        QLineEdit#power:focus {{ border-bottom-color: {LASER}; }}
        QComboBox::drop-down {{ border: none; width: 18px; }}
        QComboBox QAbstractItemView {{ background: {FIELD}; border: 1px solid {EDGE}; selection-background-color: {EDGE}; }}
        QPushButton {{ background: #1c232d; border: 1px solid {EDGE}; border-radius: 5px; padding: 7px 12px;
                       font-weight: 600; letter-spacing: 1px; }}
        QPushButton:hover {{ background: #232c38; border-color: #3a4554; }}
        QPushButton:pressed {{ background: {FIELD}; }}
        QPushButton#step {{ padding: 4px 0; font-family: {MONO}; font-weight: normal; letter-spacing: 0; }}
        QPushButton#step:checked {{ color: {PERYLENE}; border-color: {PERYLENE}; background: #0f2a30; }}
        QPushButton#segment {{ padding: 6px 0; font-size: 11px; letter-spacing: 2px; color: {MUTED}; border-radius: 0; }}
        QPushButton#segment[side="left"] {{ border-top-left-radius: 5px; border-bottom-left-radius: 5px; }}
        QPushButton#segment[side="right"] {{ border-top-right-radius: 5px; border-bottom-right-radius: 5px; }}
        QPushButton#segment:checked {{ color: {TEXT}; background: #0f2a30; border-color: {PERYLENE}; }}
        QPushButton#fire {{ background: #10281b; border: 1px solid {GO}; color: {GO}; font-size: 15px; padding: 11px;
                            letter-spacing: 4px; }}
        QPushButton#fire:hover {{ background: #15361f; }}
        QPushButton#fire:checked {{ background: #360d07; border-color: {LASER}; color: #ff6a4a; }}
        QPushButton#fire:checked:hover {{ background: #45110a; }}
        QSlider {{ background: transparent; }}
        QSlider::groove:horizontal {{ height: 6px; background: {FIELD}; border: 1px solid {EDGE}; border-radius: 3px; }}
        QSlider::sub-page:horizontal {{ border-radius: 3px;
            background: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #3d0b04, stop:1 {LASER}); }}
        QSlider::handle:horizontal {{ background: #f1f4f8; border: 2px solid {LASER}; width: 14px; margin: -7px 0;
                                      border-radius: 9px; }}
        QCheckBox {{ spacing: 8px; background: transparent; }}
        QCheckBox::indicator {{ width: 14px; height: 14px; border: 1px solid #3a4554; border-radius: 3px;
                                background: {FIELD}; }}
        QCheckBox::indicator:checked {{ background: {PERYLENE}; border-color: {PERYLENE}; }}
        QCheckBox::indicator:disabled {{ border-color: {EDGE}; }}
        QProgressBar {{ background: {FIELD}; border: none; border-radius: 2px; max-height: 4px; }}
        QProgressBar::chunk {{ background: {AMBER}; border-radius: 2px; }}
        QToolButton#disclosure {{ background: transparent; border: none; color: {MUTED}; font-size: 11px;
                                  font-weight: 600; letter-spacing: 2px; padding: 0; }}
        QToolButton#disclosure:hover {{ color: {TEXT}; }}
        QScrollArea {{ background: transparent; border: none; }}
        QScrollBar:vertical {{ background: transparent; width: 8px; margin: 0; }}
        QScrollBar::handle:vertical {{ background: {EDGE}; border-radius: 4px; min-height: 30px; }}
        QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{ height: 0; }}
    """)


class Window(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Live pulse experiment')
        self.resize(1460, 900)
        self.values = default_values()
        self.shown_rate = self.values['kExc']  # what the kExc box shows, which a drag hasn't applied yet
        self.cells_per_bin = DEFAULT_BIN_CELLS
        self.recent = deque()  # (time, molecules, photons), for the rates on screen
        self.last_drawn = None
        self.hovered = None
        self.frame = 0
        self.running = False
        self.acquired_s = 0.0  # real time spent running since the histogram was last cleared
        self.last_tick = time.perf_counter()
        self.generation = 0  # batches still on their way from older parameters are ignored
        self.target = TARGET
        self.exact_key = None  # the parameters of the exact solve that's running, done or failed, as JSON
        self.exact = None  # its result, once it's done
        self.exact_doing = None  # what it last said it was doing, while it's running
        self.exact_began = None  # when it started work, from its first progress message

        self.link, far_end = multiprocessing.Pipe()
        self.process = multiprocessing.Process(target=acquire, args=(far_end,), daemon=True)
        self.process.start()
        # now, so they've loaded before they're needed. The spare takes over when a solve is stopped part way
        self.solver, self.spare = self.new_solver(), self.new_solver()
        self.build()
        self.restart()

        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.tick)
        self.timer.start(FRAME_MS)

    # ---- layout ----

    def build(self):
        self.boxes = {}
        column = QtWidgets.QWidget()
        column.setObjectName('column')
        left = QtWidgets.QVBoxLayout(column)
        left.setContentsMargins(0, 0, 6, 0)
        left.setSpacing(12)
        left.addWidget(self.build_laser())
        left.addWidget(self.build_acquisition())
        left.addWidget(self.build_sample())
        left.addStretch()
        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(column)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setFixedWidth(356)

        right = QtWidgets.QVBoxLayout()
        right.setSpacing(12)
        right.addLayout(self.build_tiles())
        right.addWidget(self.build_plot(), 1)
        status = QtWidgets.QFrame()
        status.setObjectName('statusbar')
        status_row = QtWidgets.QHBoxLayout(status)
        status_row.setContentsMargins(12, 6, 12, 6)
        self.status = label(HINT, 'status')
        status_row.addWidget(self.status)
        right.addWidget(status)

        body = QtWidgets.QHBoxLayout()
        body.setSpacing(14)
        body.addWidget(scroll)
        body.addLayout(right, 1)
        whole = QtWidgets.QVBoxLayout()
        whole.setContentsMargins(16, 12, 16, 14)
        whole.setSpacing(12)
        whole.addLayout(self.build_header())
        whole.addLayout(body, 1)
        root = QtWidgets.QWidget()
        root.setObjectName('root')
        root.setLayout(whole)
        self.setCentralWidget(root)

        QtGui.QShortcut(QtGui.QKeySequence('Ctrl+S'), self, self.save)
        QtGui.QShortcut(QtGui.QKeySequence('F5'), self, self.fire.click)
        QtGui.QShortcut(QtGui.QKeySequence('Home'), self, self.home)

    def build_header(self):
        row = QtWidgets.QHBoxLayout()
        row.addWidget(label('LIVE PULSE EXPERIMENT', 'title'))
        row.addSpacing(14)
        row.addWidget(label('time-resolved emission, Monte Carlo  ·  times in ns, rates per ns', 'subtitle'))
        row.addStretch()
        row.addWidget(label('ACQ', 'section'))
        self.clock = label(clock_text(0), 'clock')
        row.addWidget(self.clock)
        row.addSpacing(22)
        self.laser_lamp = Lamp('LASER OFF')
        self.exact_lamp = Lamp('EXACT')
        row.addWidget(self.laser_lamp)
        row.addSpacing(18)
        row.addWidget(self.exact_lamp)
        return row

    def build_laser(self):
        frame, box = panel(f'LASER  ·  {WAVELENGTH_NM} nm')
        self.fire = QtWidgets.QPushButton('▶   START')
        self.fire.setObjectName('fire')
        self.fire.setCheckable(True)
        self.fire.setToolTip('Start or stop flashing molecules (F5). Stopping keeps the counts.')
        self.fire.clicked.connect(self.toggle)
        box.addWidget(self.fire)
        box.addSpacing(4)

        box.addWidget(label('EXCITATION RATE, kExc', 'section'))
        reading = QtWidgets.QHBoxLayout()
        self.kexc_box = QtWidgets.QLineEdit()
        self.kexc_box.setObjectName('power')
        self.kexc_box.setToolTip(f'Excitation rate in /ns, from {KEXC_LOW:g} to {KEXC_HIGH:g}. It scales with laser '
                                 'power. Type a value and press Enter, or use the slider.')
        self.kexc_box.editingFinished.connect(self.rate_typed)
        reading.addWidget(self.kexc_box, 1)
        reading.addWidget(label('/ns', 'unit'), 0, QtCore.Qt.AlignmentFlag.AlignBottom)
        box.addLayout(reading)

        self.slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.slider.setRange(0, SLIDER_STEPS)
        self.slider.setPageStep(SLIDER_STEPS // 24)  # a quarter of a decade
        self.slider.setToolTip('kExc on a log scale. The histogram clears when you let go.')
        self.slider.valueChanged.connect(self.slider_moved)
        self.slider.sliderReleased.connect(self.apply_rate)
        box.addWidget(self.slider)
        ends = QtWidgets.QHBoxLayout()
        ends.addWidget(label(f'{KEXC_LOW:g} /ns', 'note'))
        ends.addStretch()
        ends.addWidget(label(f'{KEXC_HIGH:g} /ns', 'note'))
        box.addLayout(ends)
        self.apply_timer = QtCore.QTimer(self)
        self.apply_timer.setSingleShot(True)
        self.apply_timer.timeout.connect(self.apply_rate)

        steps = QtWidgets.QHBoxLayout()
        steps.setSpacing(6)
        for text, factor in (('÷10', 0.1), ('÷2', 0.5), ('×2', 2), ('×10', 10)):
            button = QtWidgets.QPushButton(text)
            button.setObjectName('step')
            button.setToolTip(f'Multiply kExc by {factor:g}')
            button.clicked.connect(lambda _, factor=factor: self.set_rate(self.values['kExc'] * factor))
            steps.addWidget(button)
        box.addLayout(steps)
        self.kexc_read = label('', 'readout')
        self.kexc_read.setToolTip('The laser power that gives this kExc at the centre of the beam, for the setup below. '
                                  'Change the setup at the top of live.py.')
        box.addWidget(self.kexc_read)
        setup = label(f'w0 {BEAM_RADIUS_UM:g} µm  ·  ε {EXTINCTION:,} /M/cm  ·  {WAVELENGTH_NM} nm  ·  '
                      f'{kexc_per_mw():.3g} /ns per mW', 'note')
        setup.setToolTip('Assumed: Gaussian beam radius w0 (1/e² of the intensity), the sensitizer\'s molar extinction '
                         'coefficient ε at the laser line, and the wavelength. kExc = σ I0 / hν at the beam centre.')
        setup.setWordWrap(True)
        box.addWidget(setup)

        box.addSpacing(4)
        grid = QtWidgets.QGridLayout()
        self.field(grid, 0, 'duration')
        box.addLayout(grid)
        return frame

    def build_acquisition(self):
        frame, box = panel('ACQUISITION')
        grid = QtWidgets.QGridLayout()
        self.field(grid, 0, 'window')
        self.bin_box = QtWidgets.QComboBox()
        self.bin_box.currentIndexChanged.connect(self.bin_changed)
        grid.addWidget(label('Bin width', 'caption'), 1, 0)
        grid.addWidget(self.bin_box, 1, 1, 1, 2)
        self.target_box = QtWidgets.QLineEdit(f'{TARGET}')
        self.target_box.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight)
        self.target_box.setToolTip('Molecules to flash per real second. 0 = as fast as the machine can go.')
        self.target_box.editingFinished.connect(self.target_changed)
        grid.addWidget(label('Throughput', 'caption'), 2, 0)
        grid.addWidget(self.target_box, 2, 1)
        grid.addWidget(label('mol/s', 'unit'), 2, 2)
        box.addLayout(grid)

        checks = QtWidgets.QHBoxLayout()
        self.log_box = QtWidgets.QCheckBox('Log counts')
        self.log_box.toggled.connect(self.log_changed)
        self.exact_box = QtWidgets.QCheckBox('Exact overlay')
        self.exact_box.setToolTip('Solve for the histogram the counts are heading for, and draw it over them.')
        self.exact_box.toggled.connect(self.update_exact)
        checks.addWidget(self.log_box)
        checks.addWidget(self.exact_box)
        checks.addStretch()
        box.addLayout(checks)
        self.exact_note = label('', 'caption')
        self.exact_note.setWordWrap(True)
        box.addWidget(self.exact_note)
        self.exact_bar = QtWidgets.QProgressBar()
        self.exact_bar.setRange(0, 1000)
        self.exact_bar.setTextVisible(False)
        self.exact_bar.hide()
        box.addWidget(self.exact_bar)

        buttons = QtWidgets.QHBoxLayout()
        clear = QtWidgets.QPushButton('CLEAR')
        clear.setToolTip('Empty the histogram. Keeps running if it was.')
        clear.clicked.connect(self.restart)
        save = QtWidgets.QPushButton('SAVE')
        save.setToolTip('Save the histogram, its settings and a picture to sim/data/ (Ctrl+S)')
        save.clicked.connect(self.save)
        buttons.addWidget(clear)
        buttons.addWidget(save)
        box.addLayout(buttons)
        return frame

    def build_sample(self):
        frame, box = panel()
        self.disclosure = QtWidgets.QToolButton()
        self.disclosure.setObjectName('disclosure')
        self.disclosure.setText('SAMPLE')
        self.disclosure.setCheckable(True)
        self.disclosure.setArrowType(QtCore.Qt.ArrowType.RightArrow)
        self.disclosure.setToolButtonStyle(QtCore.Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.disclosure.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.disclosure.toggled.connect(self.show_sample)
        box.addWidget(self.disclosure)

        # which sample is in the beam: always in view, so the two traces are one click apart
        choice = QtWidgets.QHBoxLayout()
        choice.setSpacing(0)
        self.sample_buttons = QtWidgets.QButtonGroup(self)
        for name, text, tip in (('tta', 'TTA', 'The TTA Flora 4 molecule: sensitizer, n DPAs and the perylene'),
                                ('fluorescence', 'FLUORESCENCE', 'A plain two-state fluorophore, G and S, to compare '
                                                                 'with: excited at kExc, emits at 1/τF')):
            button = QtWidgets.QPushButton(text)
            button.setObjectName('segment')
            button.setProperty('side', 'left' if name == 'tta' else 'right')
            button.setCheckable(True)
            button.setChecked(name == self.values['sample'])
            button.setToolTip(tip)
            button.clicked.connect(lambda _, name=name: self.choose_sample(name))
            self.sample_buttons.addButton(button)
            button.setProperty('sample', name)
            choice.addWidget(button)
        box.addLayout(choice)
        self.sample_summary = label('', 'readout')
        self.sample_summary.setWordWrap(True)
        box.addWidget(self.sample_summary)

        self.sample_body = QtWidgets.QWidget()
        inner = QtWidgets.QVBoxLayout(self.sample_body)
        inner.setContentsMargins(0, 4, 0, 0)
        self.sample_fields = {}  # each sample's values, shown only while it's chosen
        for sample, names in SAMPLES.items():
            self.sample_fields[sample] = QtWidgets.QWidget()
            grid = QtWidgets.QGridLayout(self.sample_fields[sample])
            grid.setContentsMargins(0, 0, 0, 0)
            for row, name in enumerate(names):
                self.field(grid, row, name)
            self.sample_fields[sample].setVisible(sample == self.values['sample'])
            inner.addWidget(self.sample_fields[sample])
        self.derived = label('', 'readout')
        inner.addWidget(self.derived)
        self.sample_body.hide()
        box.addWidget(self.sample_body)
        box.addWidget(label('Changing a value clears the histogram. Bin width,\nthe log axis and throughput do not.',
                            'note'))
        return frame

    def field(self, grid, row, name):
        """A row of a panel: the name, a box to type the value in, and its unit."""
        text, unit, low, high = FIELDS[name]
        box = QtWidgets.QLineEdit(f'{self.values[name]:g}')
        box.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight)
        box.setFixedWidth(104)
        box.setToolTip(f'From {low:g} to {high:g}. Press Enter or Tab to apply.')
        box.editingFinished.connect(lambda name=name: self.field_changed(name))
        self.boxes[name] = box
        grid.addWidget(label(text, 'caption'), row, 0)
        grid.addWidget(box, row, 1)
        unit_label = label(unit, 'unit')
        unit_label.setMinimumWidth(46)
        grid.addWidget(unit_label, row, 2)
        grid.setColumnStretch(0, 1)

    def build_tiles(self):
        self.tiles, self.subs = {}, {}
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(10)
        for name, caption, role, stretch, tip in (
                ('yield', 'QUANTUM YIELD', 'hero', 3, 'Photons emitted per photon absorbed, with its counting error, '
                                                     'and the exact solve\'s value when the overlay is on'),
                ('rate', 'COUNT RATE', 'big', 3, 'Photons detected per real second, averaged over 2 s'),
                ('photons', 'PHOTONS', 'value', 2, 'Photons counted since the last clear'),
                ('molecules', 'MOLECULES', 'value', 2, 'Molecules flashed since the last clear (one pulse each), and '
                                                       'how many per real second'),
                ('absorbed', 'ABSORBED / MOLECULE', 'value', 2, 'Photons absorbed (rule 1 excitations) per molecule')):
            frame, box = panel(caption, 'tile')
            frame.setToolTip(tip)
            box.setContentsMargins(12, 9, 12, 10)
            box.setSpacing(2)
            self.tiles[name] = label('–', role)
            box.addWidget(self.tiles[name])
            self.subs[name] = label(' ', 'sub')
            box.addWidget(self.subs[name])
            row.addWidget(frame, stretch)
        return row

    def build_plot(self):
        pg.setConfigOptions(antialias=True, background=PANEL, foreground=MUTED)
        frame, box = panel()
        box.setContentsMargins(6, 8, 10, 6)
        self.plot = pg.PlotWidget()
        self.plot.showGrid(x=True, y=True, alpha=0.12)
        self.plot.setLabel('left', 'Photons per bin')
        self.plot.setLabel('bottom', 'Time after the pulse starts (ns). Red: pulse on')
        self.plot.setMouseEnabled(x=True, y=False)
        for side in ('left', 'bottom'):
            self.plot.getAxis(side).setPen(pg.mkPen('#3a4554'))
            self.plot.getAxis(side).setTextPen(pg.mkPen(MUTED))
        self.plot.addLegend(offset=(-20, 12), labelTextColor=TEXT, labelTextSize='10pt', brush=pg.mkBrush(BG + 'cc'),
                            pen=pg.mkPen(EDGE))
        self.band = pg.LinearRegionItem((0, 1), movable=False, brush=pg.mkBrush(255, 33, 0, 46), pen=pg.mkPen(None))
        self.band.setZValue(-10)
        self.plot.addItem(self.band, ignoreBounds=True)
        # each sample's counted curves, main pile first, then its exact curves. A curve is on the plot only while its
        # sample is chosen (and an exact one only while the exact solve is shown), so the legend lists just those
        def counted(colour, name, alpha):
            return pg.PlotDataItem(stepMode='center', pen=pg.mkPen(colour, width=2), name=name, fillLevel=0,
                                   brush=pg.mkBrush(colour + alpha))

        def exact(colour, name):
            return pg.PlotDataItem(pen=pg.mkPen(colour, width=1.6, style=QtCore.Qt.PenStyle.DashLine), name=name,
                                   connect='finite')

        self.count_curves = {'tta': [counted(PERYLENE, 'Perylene photons', '24'), counted(DPA, 'DPA photons', '30')],
                             'fluorescence': [counted(FLUOR, 'Fluorescence photons', '28')]}
        self.exact_curves = {'tta': [exact(PERYLENE_EXACT, 'Perylene, exact'), exact(DPA_EXACT, 'DPA, exact')],
                             'fluorescence': [exact(FLUOR_EXACT, 'Fluorescence, exact')]}
        self.all_curves = [curve for curves in (*self.count_curves.values(), *self.exact_curves.values())
                           for curve in curves]
        self.cursor = pg.InfiniteLine(angle=90, pen=pg.mkPen(MUTED, style=QtCore.Qt.PenStyle.DotLine))
        self.cursor.hide()
        self.plot.addItem(self.cursor, ignoreBounds=True)
        self.plot.scene().sigMouseMoved.connect(self.hover)
        self.plot.hideButtons()  # HOME and AUTO Y do what pyqtgraph's small "A" did
        self.auto_y = False  # turned on below, once the button exists
        # the counts axis is fitted by fit_y, not pyqtgraph, so it can take in the exact curves and skip hidden ones
        self.plot.getViewBox().sigXRangeChanged.connect(lambda: self.auto_y and self.fit_y())
        for curve in self.all_curves:
            curve.visibleChanged.connect(lambda: self.auto_y and self.fit_y())  # a legend click hides or shows it

        bar = QtWidgets.QHBoxLayout()
        bar.setContentsMargins(8, 0, 0, 0)
        bar.addWidget(label('HISTOGRAM', 'section'))
        bar.addWidget(label('drag to pan  ·  scroll to zoom  ·  click a legend entry to hide it', 'note'))
        bar.addStretch()
        self.auto_button = QtWidgets.QPushButton('↕  AUTO Y')
        self.auto_button.setObjectName('step')
        self.auto_button.setCheckable(True)
        self.auto_button.setToolTip('On: the counts axis keeps fitting every curve in view, the exact overlay included.\n'
                                    'Off: zoom and pan the counts axis yourself (scroll over the axis to zoom only it).')
        self.auto_button.toggled.connect(self.set_auto_y)
        home = QtWidgets.QPushButton('⌂  HOME')
        home.setObjectName('step')
        home.setToolTip('Back to the whole time window, with AUTO Y on (Home key)')
        home.clicked.connect(self.home)
        for button in (self.auto_button, home):
            button.setMinimumWidth(104)
            bar.addWidget(button)
        box.addLayout(bar)
        box.addWidget(self.plot)
        self.auto_button.setChecked(True)
        return frame

    def home(self):
        """The original view: the whole time window, with the counts axis fitted to the data."""
        self.plot.setXRange(0, self.measurement.p['window'], padding=0)
        self.auto_button.setChecked(True)
        self.fit_y()

    def set_auto_y(self, on):
        self.auto_y = on
        self.plot.setMouseEnabled(x=True, y=not on)  # with AUTO Y off, the counts axis is yours to zoom
        if on:
            self.fit_y()

    def fit_y(self):
        """Fit the counts axis to every curve that's shown, within the time range in view."""
        box = self.plot.getViewBox()
        start, end = box.viewRange()[0]
        values = []
        for curve in self.all_curves:
            if not curve.isVisible() or curve not in self.plot.listDataItems():
                continue
            x, y = curve.getData()  # in log mode these are already log10 of the counts
            if x is None or y is None or not len(y):
                continue
            if len(x) == len(y) + 1:  # a step curve has bin edges: a bin is in view if any of it is
                inside = (x[1:] >= start) & (x[:-1] <= end)
            else:
                inside = (x >= start) & (x <= end)
            values.append(y[inside & np.isfinite(y)])
        values = np.concatenate(values) if values else np.array([])
        log = self.log_box.isChecked()
        bottom = np.log10(EMPTY_BIN) if log else 0.0
        top = values.max() if values.size else bottom
        if top <= bottom:
            top = bottom + 1
        box.setYRange(bottom, top + 0.06 * (top - bottom), padding=0)

    def show_sample(self, open_):
        self.disclosure.setArrowType(QtCore.Qt.ArrowType.DownArrow if open_ else QtCore.Qt.ArrowType.RightArrow)
        self.sample_body.setVisible(open_)
        self.sample_summary.setVisible(not open_)

    # ---- changing things ----

    def restart(self):
        """Start a new, empty histogram for the current values. If it's running, it stays running."""
        self.measurement = Measurement(parameters(self.values))
        p = self.measurement.p
        self.generation += 1
        self.link.send({'p': p, 'generation': self.generation})
        self.band.setRegion((0, p['duration']))
        self.home()
        self.bin_box.blockSignals(True)
        self.bin_box.clear()
        for cells in BIN_CELLS:
            self.bin_box.addItem(f'{cells * self.measurement.cell:.3g} ns', cells)
        self.bin_box.setCurrentIndex(BIN_CELLS.index(self.cells_per_bin))
        self.bin_box.blockSignals(False)
        self.recent.clear()
        self.acquired_s = 0.0
        self.last_drawn = None
        self.set_slider(self.values['kExc'])
        self.show_rate(self.values['kExc'])
        if fluorescent(p):
            self.sample_summary.setText(f"two-state fluorophore\nτF {p['tauFl']:g} ns  ·  Φ_F {p['fluorYield']:g}")
            # a steady light would leave this share in S: how close the pulse comes to saturating the fluorophore
            self.derived.setText(f"1/τF {1 / p['tauFl']:.3g} /ns  ·  saturation "
                                 f"{p['kExc'] / (p['kExc'] + 1 / p['tauFl']):.1%}")
        else:
            self.sample_summary.setText(f"n {p['n']}  ·  kHop {p['kHop']:g} /ns  ·  f {p['f']:g}\n"
                                        f"τS {p['tauS']:g} ns  ·  τSens {p['tauSens']:g} ns  ·  "
                                        f"Φ_T {p['tripletYield']:g}")
            self.derived.setText(f"kISC {p['kISC']:.3g} /ns  ·  hop time {1 / p['kHop']:.3g} ns")
        self.update_exact()

    def choose_sample(self, sample):
        """Put the other sample in the beam. The laser, pulse and time window stay as they are."""
        for name, fields in self.sample_fields.items():
            fields.setVisible(name == sample)
        for button in self.sample_buttons.buttons():
            button.setChecked(button.property('sample') == sample)
        if sample != self.values['sample']:
            self.values['sample'] = sample
            self.restart()

    def field_changed(self, name):
        box = self.boxes[name]
        text, unit, low, high = FIELDS[name]
        try:
            value = int(box.text()) if name == 'n' else float(box.text())
            allowed = low <= value <= high
        except ValueError:
            allowed = False
        if not allowed:
            self.status.setText(f'{text} has to be a number from {low:g} to {high:g} {unit}'.rstrip())
            box.setText(f'{self.values[name]:g}')
        elif value != self.values[name]:
            self.values[name] = value
            self.restart()

    def set_slider(self, kexc):
        self.slider.blockSignals(True)
        self.slider.setValue(rate_to_step(kexc))
        self.slider.blockSignals(False)

    def show_rate(self, kexc):
        """Show a kExc in the box, without applying it yet."""
        self.shown_rate = kexc
        self.kexc_box.setText(f'{kexc:#.3g}'.rstrip('.'))  # always three figures, like a meter: 0.200, 1.00
        self.kexc_read.setText(f'laser power  ≈ {power_text(kexc / kexc_per_mw())}')

    def slider_moved(self, step):
        self.show_rate(step_to_rate(step))
        if not self.slider.isSliderDown():  # the keyboard or the mouse wheel: wait in case another step follows
            self.apply_timer.start(APPLY_MS)

    def apply_rate(self):
        self.apply_timer.stop()
        if self.shown_rate != self.values['kExc']:
            self.values['kExc'] = self.shown_rate
            self.restart()

    def set_rate(self, kexc):
        """A kExc from the box or the step buttons, applied at once."""
        kexc = float(f'{np.clip(kexc, KEXC_LOW, KEXC_HIGH):.3g}')
        self.set_slider(kexc)
        self.show_rate(kexc)
        self.apply_rate()

    def rate_typed(self):
        try:
            kexc = float(self.kexc_box.text())
            if not kexc > 0:
                raise ValueError
        except ValueError:
            self.status.setText(f'kExc has to be a number from {KEXC_LOW:g} to {KEXC_HIGH:g} /ns')
            self.show_rate(self.values['kExc'])
            return
        if not KEXC_LOW <= kexc <= KEXC_HIGH:
            self.status.setText(f'kExc goes from {KEXC_LOW:g} to {KEXC_HIGH:g} /ns, so it was set to the nearest end')
        self.set_rate(kexc)

    def target_changed(self):
        try:
            self.target = max(0.0, float(self.target_box.text()))
            self.link.send({'target': self.target})
        except ValueError:
            self.status.setText('Throughput has to be a number of molecules per second')
        self.target_box.setText(f'{self.target:g}')
        self.recent.clear()

    def bin_changed(self):
        self.cells_per_bin = self.bin_box.currentData()
        self.last_drawn = None

    def log_changed(self, on):
        self.plot.setLogMode(y=on)
        floor = np.log10(EMPTY_BIN) if on else 0  # in log mode the fill level is in log units
        for curves in self.count_curves.values():
            for curve in curves:
                curve.setFillLevel(floor)
        self.plot.setXRange(0, self.measurement.p['window'], padding=0)  # setLogMode resets it
        self.draw()
        self.fit_y()  # the old counts range means nothing on the other kind of axis, so fit it even with AUTO Y off

    def toggle(self):
        self.running = self.fire.isChecked()
        self.link.send({'on': self.running})
        self.fire.setText('■   STOP' if self.running else '▶   START')
        if self.running:  # a fresh average. Stopping keeps the last one on screen
            self.recent.clear()
        self.update_lamps()

    # ---- the exact solve ----

    def update_exact(self):
        """Start the exact solve for the current values if it's wanted and not already running or done."""
        p = self.measurement.p
        allowed, verdict = exact_verdict(p)
        self.exact_box.setEnabled(allowed)
        key = json.dumps(p, sort_keys=True)
        if not (allowed and self.exact_box.isChecked()):
            if self.exact_doing:  # stop it, so it doesn't take time from the counting
                self.stop_solve()
                self.exact_key = None
            self.exact_note.setText(verdict)
            self.exact_bar.hide()
        elif key != self.exact_key:
            if self.exact_doing:
                self.stop_solve()
            self.solver[1].send((key, p))
            self.exact_key, self.exact, self.exact_doing, self.exact_began = key, None, 'Starting', None
            self.exact_note.setText(f'Starting. {verdict}')
            self.exact_bar.setValue(0)
            self.exact_bar.show()
        elif self.exact is not None:
            self.exact_note.setText('Exact solve done.')
        self.show_exact_curves()
        self.update_lamps()

    @staticmethod
    def new_solver():
        """A loaded, waiting exact solve process, as (process, link)."""
        link, far_end = multiprocessing.Pipe()
        process = multiprocessing.Process(target=solve_exactly, args=(far_end,), daemon=True)
        # one BLAS thread: more only compete with the simulation's threads, and the solve was slower for it
        threads = os.environ.get('OPENBLAS_NUM_THREADS')
        os.environ['OPENBLAS_NUM_THREADS'] = '1'
        process.start()
        if threads is None:
            del os.environ['OPENBLAS_NUM_THREADS']
        else:
            os.environ['OPENBLAS_NUM_THREADS'] = threads
        return process, link

    def stop_solve(self):
        """End the running solve's process at once, and switch to the spare, which has already loaded."""
        self.solver[0].terminate()
        self.solver, self.spare = self.spare, self.new_solver()
        self.exact_doing = None

    def exact_shown(self):
        return self.exact is not None and self.exact_box.isChecked() and self.exact_key == json.dumps(
            self.measurement.p, sort_keys=True)

    def show_exact_curves(self):
        """Put the chosen sample's curves on the plot, and its exact ones while they're shown. Take the rest off."""
        sample = self.measurement.p['sample']
        wanted = self.count_curves[sample] + (self.exact_curves[sample] if self.exact_shown() else [])
        there = [curve for curve in self.all_curves if curve in self.plot.listDataItems()]
        if there != wanted:  # all off and back on in order, so the legend lists them in order
            for curve in there:
                self.plot.removeItem(curve)
            for curve in wanted:
                self.plot.addItem(curve, ignoreBounds=curve not in self.count_curves[sample])
        self.last_drawn = None

    def read_exact(self):
        """Take in what the exact solve has sent: how far it's got, or its result."""
        while self.solver[1].poll():
            kind, key, *rest = self.solver[1].recv()
            if key != self.exact_key:  # from a solve that was dropped
                continue
            now = time.perf_counter()
            self.exact_began = self.exact_began or now
            spent = now - self.exact_began
            if kind == 'progress':
                share, self.exact_doing = rest
                left = spent / share * (1 - share) if share > 0.05 else None
                self.exact_bar.setValue(int(share * 1000))
                self.exact_note.setText(self.exact_doing + (f', about {how_long(left)} left.' if left else '.'))
                continue
            self.exact_doing = None
            self.exact_bar.hide()
            if kind == 'done':
                self.exact = rest[0]
                self.exact_note.setText(f'Exact solve done in {how_long(spent)}.')
                self.show_exact_curves()
                self.show_numbers()
            else:
                self.exact_note.setText(f'The exact solve failed: {rest[0]}')
            self.update_lamps()

    # ---- every frame ----

    def tick(self):
        try:
            while self.link.poll():
                generation, *batch = self.link.recv()
                if generation == self.generation:
                    self.measurement.add(*batch)
            self.read_exact()
        except (BrokenPipeError, EOFError):  # a process has gone, as when the window is closing. Nothing left to read
            self.timer.stop()
            return
        totals = self.measurement.totals
        molecules, photons = totals['molecules'], totals['photons']
        now = time.perf_counter()
        if self.running:
            self.acquired_s += now - self.last_tick
            self.recent.append((now, molecules, photons))
            while now - self.recent[0][0] > RATE_OVER_S:
                self.recent.popleft()
            if not molecules:
                self.status.setText(STARTING)
            elif self.status.text() == STARTING:
                self.status.setText(HINT)
        self.last_tick = now

        self.frame += 1
        if self.last_drawn != molecules:
            self.last_drawn = molecules
            self.draw()
        if self.frame % FRAMES_PER_TEXT == 0:
            self.show_numbers()

    def draw(self):
        edges, main, other = self.measurement.binned(self.cells_per_bin)
        sample = self.measurement.p['sample']
        floor = EMPTY_BIN if self.log_box.isChecked() else 0
        for curve, counts in zip(self.count_curves[sample], (main, other)):  # the fluorophore has no `other`
            curve.setData(edges, np.maximum(counts, floor))
        expected = None
        if self.exact_shown():
            # photons the exact solve expects in each bin from the molecules counted so far, drawn at the bin centres
            expected = [share * self.measurement.totals['molecules'] for share in self.exact.binned(edges)]
            centres = (edges[:-1] + edges[1:]) / 2
            for curve, counts in zip(self.exact_curves[sample], expected):
                curve.setData(centres, np.where(counts > 0, counts, np.nan) if self.log_box.isChecked() else counts)
        self.drawn = edges, main, other, expected
        if self.auto_y:
            self.fit_y()
        if self.hovered is not None:
            self.show_bin(self.hovered)

    def update_lamps(self):
        if self.running:
            self.laser_lamp.set(LASER, f'LASER ON  ·  kExc {self.values["kExc"]:.3g} /ns')
        else:
            self.laser_lamp.set(None, 'LASER OFF')
        if self.exact_doing:
            self.exact_lamp.set(AMBER, f'EXACT  ·  SOLVING {self.exact_bar.value() / 10:.0f}%')
        elif self.exact_shown():
            self.exact_lamp.set(PERYLENE, 'EXACT  ·  READY')
        else:
            self.exact_lamp.set(None, 'EXACT')

    def show_numbers(self):
        totals = self.measurement.totals
        molecules, photons, excitations = totals['molecules'], totals['photons'], totals['excitations']
        spent = self.recent[-1][0] - self.recent[0][0] if len(self.recent) > 1 else 0
        per_s = [(self.recent[-1][k] - self.recent[0][k]) / spent if spent else 0 for k in (1, 2)]
        past = totals['past_window'] / photons if photons else 0
        tiles = {'yield': '–', 'rate': f'{per_s[1]:,.0f}', 'photons': f'{photons:,}', 'molecules': f'{molecules:,}',
                 'absorbed': f'{excitations / molecules:.3g}' if molecules else '–'}
        held = '' if self.running or not per_s[0] else '  ·  held'  # the rates from just before Stop
        subs = {'yield': ' ', 'rate': f'photons / s{held}', 'molecules': f'{per_s[0]:,.0f} /s{held}',
                'absorbed': 'photons',
                'photons': f'{past:.0%} past the window' if past >= PAST_WARNING else 'in the histogram'}
        if excitations:
            tiles['yield'] = f'{photons / excitations:.2%}'
            subs['yield'] = f'± {photons ** 0.5 / excitations:.2%}'
        if self.exact_shown() and np.isfinite(self.exact.excitations):
            # at least, if it stopped counting after the window while some molecules could still emit
            unfinished = self.exact.still_active >= QUIET * self.exact.excitations
            subs['yield'] += f"  ·  exact {'≥ ' if unfinished else ''}{self.exact.photons / self.exact.excitations:.2%}"
        for name, text in tiles.items():
            self.tiles[name].setText(text)
        for name, text in subs.items():
            self.subs[name].setText(text)
        warn = 'warn' if past >= PAST_WARNING else 'sub'
        if self.subs['photons'].property('role') != warn:  # restyle only when it changes
            self.subs['photons'].setProperty('role', warn)
            self.subs['photons'].style().polish(self.subs['photons'])
        self.clock.setText(clock_text(self.acquired_s))
        self.update_lamps()

        p = self.measurement.p
        width = self.measurement.cell * self.cells_per_bin
        sample = (f"fluorophore, τF = {p['tauFl']:g} ns, Φ_F = {p['fluorYield']:g}" if fluorescent(p) else
                  f"TTA, kHop = {p['kHop']:.3g} /ns, n = {p['n']}")
        self.plot.setTitle(f"<span style='color:{TEXT}'>Quantum yield {tiles['yield']}</span>"
                           f"<span style='color:{MUTED}'>&nbsp;&nbsp; {sample}, kExc = {p['kExc']:.3g} /ns, "
                           f"pulse = {p['duration']:g} ns, {molecules:,} molecules, {width:.3g} ns bins</span>",
                           size='11pt')

    def hover(self, position):
        box = self.plot.getViewBox()
        if box.sceneBoundingRect().contains(position):
            self.hovered = box.mapSceneToView(position).x()
            self.cursor.setPos(self.hovered)
            self.cursor.show()
            self.show_bin(self.hovered)
        else:
            self.hovered = None
            self.cursor.hide()

    def show_bin(self, t):
        if not hasattr(self, 'drawn'):
            return
        edges, main, other, expected = self.drawn
        k = np.searchsorted(edges, t, side='right') - 1
        if 0 <= k < main.size:
            exact = [f' (exact {counts[k]:,.1f})' for counts in expected] if expected else ['', '']
            piles = PILES[self.measurement.p['sample']]
            self.status.setText('   ·   '.join([f'{edges[k]:.4g} to {edges[k + 1]:.4g} ns'] + [
                f'{name} {counts[k]:,}{note}' for name, counts, note in zip(piles, (main, other), exact)]))

    # ---- saving and closing ----

    def save(self):
        """Write the histogram as it stands, plus a picture of it. Times are in ns."""
        m = self.measurement
        edges, main, other = m.binned(self.cells_per_bin)
        piles = PILES[m.p['sample']]
        counted = (main, other)[:len(piles)]
        totals = m.totals
        if not totals['molecules']:
            self.status.setText('Nothing to save yet')
            return
        folder = new_run_folder('live', 'live-fluorescence' if fluorescent(m.p) else None, [])
        shutil.copy(HERE / 'rules.py', folder / 'rules.py')
        (folder / 'params.json').write_text(json.dumps({
            'mode': 'live', 'date': f'{datetime.now():%Y-%m-%d %H:%M:%S}', 'units': 'times in ns, rates per ns',
            'parameters': m.p, 'bin_ns': float(edges[1]),
        }, indent=2))
        width = np.diff(edges)
        scale = 1 / (totals['molecules'] * width)
        total = main + other
        # the exact solve's columns, if it's shown: the photons it expects in each bin from this many molecules
        exact = self.exact.binned(edges)[:len(piles)] if self.exact_shown() else None
        with (folder / 'histogram.csv').open('w', newline='') as out:
            writer = csv.writer(out)
            # TTA has perylene_photons and dpa_photons; the fluorophore has fluorescence_photons
            writer.writerow(['point', 'bin_start_ns', 'bin_end_ns', 'pulse_on', 'photons']
                            + [f'{name}_photons' for name in piles] + ['photons_per_molecule_per_ns', 'error']
                            + ([f'exact_{name}_photons' for name in piles] if exact else []))
            for k in range(total.size):
                writer.writerow([0, edges[k], edges[k + 1], int(edges[k] < m.p['duration']), total[k]]
                                + [counts[k] for counts in counted]
                                + [total[k] * scale[k], total[k] ** 0.5 * scale[k]]
                                + ([share[k] * totals['molecules'] for share in exact] if exact else []))
        with (folder / 'summary.csv').open('w', newline='') as out:
            row = {**m.p, **totals, 'photons_per_excitation': totals['photons'] / max(totals['excitations'], 1)}
            if exact:
                row['exact_photons_per_excitation'] = self.exact.photons / self.exact.excitations
            writer = csv.DictWriter(out, row.keys())
            writer.writeheader()
            writer.writerow(row)
        self.plot.grab().save(str(folder / 'histogram.png'))
        add_to_index(folder, 'live', 1, [], 0)
        self.status.setText(f'Saved to {folder}')

    def closeEvent(self, event):
        # no more frames: one that ran after this would poll pipes whose other ends are being closed
        self.timer.stop()
        self.apply_timer.stop()
        self.link.send({'close': True})
        self.process.join(2)
        self.process.terminate()  # only matters if it was stuck sending a batch nobody will read
        self.solver[0].terminate()
        self.spare[0].terminate()
        event.accept()


def main():
    app = QtWidgets.QApplication(sys.argv)
    style(app)
    window = Window()
    window.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
