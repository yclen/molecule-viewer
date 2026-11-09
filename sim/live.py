"""A live pulse experiment: fresh molecules are flashed all the time, and each photon lands in a histogram on screen.

  python sim/live.py

Like a time-resolved measurement, the histogram starts as a few scattered counts and smooths out as they build up. Each
molecule gets one pulse from the same code as `tta_sim.py pulse`, so a live run is that command in small batches.

The simulation runs in its own process, on threads that never wait for each other. Drawing the window would hold it up
otherwise, and so would the slower cores on a machine that has a mix of fast and slow ones.
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

from PyQt6 import QtCore, QtWidgets  # before pyqtgraph, so it picks this Qt and not any other that's installed

import numpy as np
import pyqtgraph as pg

from engine import pulse_alone, rule_arrays
from plot_pulse import BASE, DPA, GRID, INK, MUTED, PERYLENE
from rules import PULSE_DEFAULTS, build_rules, resolve
from tta_sim import BINS_PER_LIFETIME, HERE, add_to_index, new_run_folder

# what the left panel can change: label, lowest and highest allowed value
FIELDS = {
    'excRatio': ('Excitation (excRatio)', 1e-6, 1e6),
    'duration': ('Pulse length (s)', 1e-3, 1e6),
    'hopRatio': ('Hopping (hopRatio)', 1e-3, 1e7),
    'kDecay': ('Triplet decay (kDecay, /s)', 0, 1e3),
    'f': ('Spin factor (f)', 0, 1),
    'tripletYield': ('Triplet yield', 0, 0.999),
    'n': ('DPA mediators (n)', 0, 40),
    'tauS': ('Singlet lifetime (tauS, s)', 1e-3, 1e6),
    'tauSens': ('Sensitizer lifetime (tauSens, s)', 1e-3, 1e6),
    'fretRadius': ('FRET radius (nm)', 0, 100),
    'window': ('Time window (s)', 1e-2, 1e7),
}
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
HINT = 'Hover over the histogram to read a bin.'
STARTING = 'Starting the engine. The first run after a change to the code takes longer.'


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
                self.arrays = rule_arrays(build_rules(self.p))
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
            return self.generation, self.arrays, self.p['n'] + 2, self.p['duration'], size, self.seed, start

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


class Measurement:
    """Everything counted for one set of parameters."""

    def __init__(self, p):
        self.p = resolve(p)
        self.perylene_site = self.p['n'] + 1
        self.cell = max(self.p['tauS'] / CELLS_PER_LIFETIME, self.p['window'] / MOST_CELLS)
        self.cells = int(np.ceil(self.p['window'] / self.cell))
        self.perylene = np.zeros(self.cells, np.int64)
        self.dpa = np.zeros(self.cells, np.int64)
        self.totals = {'molecules': 0, 'excitations': 0, 'photons': 0, 'past_window': 0}

    def add(self, molecules, excitations, times, sites):
        cell = (times / self.cell).astype(np.int64)
        inside = cell < self.cells
        on_perylene = sites == self.perylene_site
        self.perylene += np.bincount(cell[inside & on_perylene], minlength=self.cells)
        self.dpa += np.bincount(cell[inside & ~on_perylene], minlength=self.cells)
        self.totals['molecules'] += molecules
        self.totals['excitations'] += excitations
        self.totals['photons'] += times.size
        self.totals['past_window'] += int((~inside).sum())

    def binned(self, cells_per_bin):
        """The bin edges, then the perylene and DPA photons in each bin."""
        spare = -self.cells % cells_per_bin
        perylene = np.pad(self.perylene, (0, spare)).reshape(-1, cells_per_bin).sum(1)
        dpa = np.pad(self.dpa, (0, spare)).reshape(-1, cells_per_bin).sum(1)
        return self.cell * cells_per_bin * np.arange(perylene.size + 1), perylene, dpa


class Window(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('TTA Flora 4: live pulse experiment')
        self.resize(1320, 780)
        self.values = {name: PULSE_DEFAULTS[name] for name in FIELDS if name != 'window'}
        self.values['window'] = WINDOW_LIFETIMES * self.values['tauS']
        self.cells_per_bin = DEFAULT_BIN_CELLS
        self.recent = deque()  # (time, molecules, photons), for the rates on screen
        self.last_drawn = None
        self.hovered = None
        self.frame = 0
        self.running = False
        self.generation = 0  # batches still on their way from older parameters are ignored
        self.target = TARGET

        self.link, far_end = multiprocessing.Pipe()
        self.process = multiprocessing.Process(target=acquire, args=(far_end,), daemon=True)
        self.process.start()
        self.build()
        self.restart()

        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.tick)
        self.timer.start(FRAME_MS)

    def build(self):
        self.setStyleSheet(f"""
            QWidget {{ background: white; color: {INK}; font-size: 13px; }}
            QLineEdit, QComboBox {{ border: 1px solid {BASE}; border-radius: 4px; padding: 4px; }}
            QPushButton {{ border: 1px solid {BASE}; border-radius: 4px; padding: 7px 12px; }}
            QPushButton:hover {{ background: {GRID}; }}
            QLabel[role="caption"] {{ color: {MUTED}; font-size: 12px; }}
            QLabel[role="tile"] {{ font-size: 22px; }}
            QCheckBox::indicator {{ width: 14px; height: 14px; border: 1px solid {BASE}; border-radius: 3px; }}
            QCheckBox::indicator:checked {{ background: {INK}; }}
        """)
        panel = QtWidgets.QVBoxLayout()
        self.start_button = QtWidgets.QPushButton('Start')
        self.start_button.clicked.connect(self.toggle)
        clear = QtWidgets.QPushButton('Clear')
        clear.clicked.connect(self.restart)
        save = QtWidgets.QPushButton('Save')
        save.clicked.connect(self.save)
        buttons = QtWidgets.QHBoxLayout()
        for button in (self.start_button, clear, save):
            buttons.addWidget(button)
        panel.addLayout(buttons)

        form = QtWidgets.QFormLayout()
        self.boxes = {}
        for name, (label, _, _) in FIELDS.items():
            self.boxes[name] = QtWidgets.QLineEdit(f'{self.values[name]:g}')
            self.boxes[name].editingFinished.connect(lambda name=name: self.field_changed(name))
            form.addRow(label, self.boxes[name])
        self.target_box = QtWidgets.QLineEdit(f'{TARGET}')
        self.target_box.editingFinished.connect(self.target_changed)
        form.addRow('Molecules per s (0 = flat out)', self.target_box)
        self.bin_box = QtWidgets.QComboBox()
        self.bin_box.currentIndexChanged.connect(self.bin_changed)
        form.addRow('Bin width', self.bin_box)
        self.log_box = QtWidgets.QCheckBox('Log counts axis')
        self.log_box.toggled.connect(self.log_changed)
        form.addRow('', self.log_box)
        panel.addLayout(form)

        self.rates = QtWidgets.QLabel()
        self.rates.setProperty('role', 'caption')
        panel.addWidget(self.rates)
        note = QtWidgets.QLabel('Changing a value above clears the histogram.\nBin width and the log axis do not.')
        note.setProperty('role', 'caption')
        panel.addWidget(note)
        panel.addStretch()

        self.tiles = {}
        tiles = QtWidgets.QHBoxLayout()
        for caption in ('Quantum yield', 'Photons', 'Photons per s', 'Molecules', 'Molecules per s', 'Past the window'):
            tile = QtWidgets.QVBoxLayout()
            name = QtWidgets.QLabel(caption)
            name.setProperty('role', 'caption')
            self.tiles[caption] = QtWidgets.QLabel('0')
            self.tiles[caption].setProperty('role', 'tile')
            tile.addWidget(name)
            tile.addWidget(self.tiles[caption])
            tiles.addLayout(tile)

        pg.setConfigOptions(antialias=True, background='w', foreground=MUTED)
        self.plot = pg.PlotWidget()
        self.plot.showGrid(y=True, alpha=0.25)
        self.plot.setLabel('left', 'Photons per bin')
        self.plot.setLabel('bottom', 'Time after the pulse starts (s). Shaded: pulse on')
        self.plot.setMouseEnabled(x=True, y=False)
        for side in ('left', 'bottom'):
            self.plot.getAxis(side).setPen(pg.mkPen(BASE))
            self.plot.getAxis(side).setTextPen(pg.mkPen(MUTED))
        self.plot.addLegend(offset=(-20, 12), labelTextColor=INK, labelTextSize='10pt')
        self.band = pg.LinearRegionItem((0, 1), movable=False, brush=pg.mkBrush(255, 33, 0, 40), pen=pg.mkPen(None))
        self.band.setZValue(-10)
        self.plot.addItem(self.band, ignoreBounds=True)
        self.perylene_curve = self.plot.plot(stepMode='center', pen=pg.mkPen(PERYLENE, width=2), name='Perylene photons')
        self.dpa_curve = self.plot.plot(stepMode='center', pen=pg.mkPen(DPA, width=2), name='DPA photons')
        self.cursor = pg.InfiniteLine(angle=90, pen=pg.mkPen(BASE))
        self.cursor.hide()
        self.plot.addItem(self.cursor, ignoreBounds=True)
        self.plot.scene().sigMouseMoved.connect(self.hover)

        self.status = QtWidgets.QLabel(HINT)
        self.status.setProperty('role', 'caption')

        right = QtWidgets.QVBoxLayout()
        right.addLayout(tiles)
        right.addWidget(self.plot, 1)
        right.addWidget(self.status)
        whole = QtWidgets.QHBoxLayout()
        whole.addLayout(panel)
        whole.addLayout(right, 1)
        root = QtWidgets.QWidget()
        root.setLayout(whole)
        self.setCentralWidget(root)

    def restart(self):
        """Start a new, empty histogram for the current values. If it's running, it stays running."""
        self.measurement = Measurement({**PULSE_DEFAULTS, **self.values})
        p = self.measurement.p
        self.generation += 1
        self.link.send({'p': p, 'generation': self.generation})
        self.band.setRegion((0, p['duration']))
        self.plot.setXRange(0, p['window'], padding=0)
        self.rates.setText(f"kExc = {p['kExc']:.3g} /s\nkHop = {p['kHop']:.3g} /s\nkISC = {p['kISC']:.3g} /s")
        self.bin_box.blockSignals(True)
        self.bin_box.clear()
        for cells in BIN_CELLS:
            self.bin_box.addItem(f'{cells * self.measurement.cell:.3g} s', cells)
        self.bin_box.setCurrentIndex(BIN_CELLS.index(self.cells_per_bin))
        self.bin_box.blockSignals(False)
        self.recent.clear()
        self.last_drawn = None

    def field_changed(self, name):
        box = self.boxes[name]
        _, low, high = FIELDS[name]
        try:
            value = int(box.text()) if name == 'n' else float(box.text())
            allowed = low <= value <= high
        except ValueError:
            allowed = False
        if not allowed:
            self.status.setText(f'{FIELDS[name][0]} has to be a number from {low:g} to {high:g}')
            box.setText(f'{self.values[name]:g}')
        elif value != self.values[name]:
            self.values[name] = value
            self.restart()

    def target_changed(self):
        try:
            self.target = max(0.0, float(self.target_box.text()))
            self.link.send({'target': self.target})
        except ValueError:
            self.status.setText('Molecules per s has to be a number')
        self.target_box.setText(f'{self.target:g}')
        self.recent.clear()

    def bin_changed(self):
        self.cells_per_bin = self.bin_box.currentData()
        self.last_drawn = None

    def log_changed(self, on):
        self.plot.setLogMode(y=on)
        self.plot.setXRange(0, self.measurement.p['window'], padding=0)  # setLogMode resets it
        self.last_drawn = None

    def toggle(self):
        self.running = not self.running
        self.link.send({'on': self.running})
        self.start_button.setText('Stop' if self.running else 'Start')
        self.recent.clear()

    def tick(self):
        while self.link.poll():
            generation, *batch = self.link.recv()
            if generation == self.generation:
                self.measurement.add(*batch)
        totals = self.measurement.totals
        molecules, photons = totals['molecules'], totals['photons']
        if self.running:
            now = time.perf_counter()
            self.recent.append((now, molecules, photons))
            while now - self.recent[0][0] > RATE_OVER_S:
                self.recent.popleft()
            if not molecules:
                self.status.setText(STARTING)
            elif self.status.text() == STARTING:
                self.status.setText(HINT)

        self.frame += 1
        if self.last_drawn != molecules:
            self.last_drawn = molecules
            self.draw()
        if self.frame % FRAMES_PER_TEXT == 0:
            self.show_numbers()

    def draw(self):
        edges, perylene, dpa = self.measurement.binned(self.cells_per_bin)
        floor = EMPTY_BIN if self.log_box.isChecked() else 0
        self.perylene_curve.setData(edges, np.maximum(perylene, floor))
        self.dpa_curve.setData(edges, np.maximum(dpa, floor))
        self.drawn = edges, perylene, dpa
        if self.hovered is not None:
            self.show_bin(self.hovered)

    def show_numbers(self):
        totals = self.measurement.totals
        molecules, photons = totals['molecules'], totals['photons']
        spent = self.recent[-1][0] - self.recent[0][0] if len(self.recent) > 1 else 0
        tiles = {'Photons': f'{photons:,}', 'Molecules': f'{molecules:,}', 'Quantum yield': '0'}
        tiles['Molecules per s'] = f'{(self.recent[-1][1] - self.recent[0][1]) / spent:,.0f}' if spent else '0'
        tiles['Photons per s'] = f'{(self.recent[-1][2] - self.recent[0][2]) / spent:,.0f}' if spent else '0'
        tiles['Past the window'] = f"{totals['past_window'] / photons:.1%}" if photons else '0'
        if totals['excitations']:
            tiles['Quantum yield'] = f"{photons / totals['excitations']:.2%} ± {photons ** 0.5 / totals['excitations']:.2%}"
        for caption, text in tiles.items():
            self.tiles[caption].setText(text)

        p = self.measurement.p
        width = self.measurement.cell * self.cells_per_bin
        self.plot.setTitle(f"<span style='color:{INK}'>Quantum yield {tiles['Quantum yield']}</span>"
                           f"<span style='color:{MUTED}'>&nbsp;&nbsp; kHop = {p['kHop']:.3g} /s, kExc = {p['kExc']:.3g} /s, "
                           f"pulse = {p['duration']:g} s, n = {p['n']}, {molecules:,} molecules, "
                           f"{width:.3g} s bins</span>", size='11pt')

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
        edges, perylene, dpa = self.drawn
        k = np.searchsorted(edges, t, side='right') - 1
        if 0 <= k < perylene.size:
            self.status.setText(f'{edges[k]:.4g} to {edges[k + 1]:.4g} s: {perylene[k]:,} perylene photons, '
                                f'{dpa[k]:,} DPA photons')

    def save(self):
        """Write the histogram as it stands, in the same columns as a pulse run, plus a picture of it."""
        m = self.measurement
        edges, perylene, dpa = m.binned(self.cells_per_bin)
        totals = m.totals
        if not totals['molecules']:
            self.status.setText('Nothing to save yet')
            return
        folder = new_run_folder('live', None, [])
        shutil.copy(HERE / 'rules.py', folder / 'rules.py')
        (folder / 'params.json').write_text(json.dumps({
            'mode': 'live', 'date': f'{datetime.now():%Y-%m-%d %H:%M:%S}', 'parameters': m.p, 'bin_s': float(edges[1]),
        }, indent=2))
        width = np.diff(edges)
        scale = 1 / (totals['molecules'] * width)
        total = perylene + dpa
        with (folder / 'histogram.csv').open('w', newline='') as out:
            writer = csv.writer(out)
            writer.writerow(['point', 'bin_start_s', 'bin_end_s', 'pulse_on', 'photons', 'perylene_photons', 'dpa_photons',
                             'photons_per_molecule_per_s', 'error'])
            for k in range(total.size):
                writer.writerow([0, edges[k], edges[k + 1], int(edges[k] < m.p['duration']), total[k], perylene[k], dpa[k],
                                 total[k] * scale[k], total[k] ** 0.5 * scale[k]])
        with (folder / 'summary.csv').open('w', newline='') as out:
            row = {**m.p, **totals, 'photons_per_excitation': totals['photons'] / max(totals['excitations'], 1)}
            writer = csv.DictWriter(out, row.keys())
            writer.writeheader()
            writer.writerow(row)
        self.plot.grab().save(str(folder / 'histogram.png'))
        add_to_index(folder, 'live', 1, [], 0)
        self.status.setText(f'Saved to {folder}')

    def closeEvent(self, event):
        self.link.send({'close': True})
        self.process.join(2)
        self.process.terminate()  # only matters if it was stuck sending a batch nobody will read
        event.accept()


def main():
    app = QtWidgets.QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
