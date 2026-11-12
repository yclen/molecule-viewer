# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy>=2.0", "scipy>=1.13", "numba>=0.60", "PyQt6>=6.6", "pyqtgraph>=0.13.7", "matplotlib>=3.8"]
# ///
"""A window for intensity sweeps: set the parameters, run the sweep in the background, then plot the results.

  python sim/sweep_gui.py

The Run tab builds a `tta_sim.py intensity` command from the boxes and runs it as its own process, so the window stays
free and a stopped sweep keeps every point it finished. The Plot tab lists every intensity run in sim/data and opens
matplotlib windows of photons per molecule and quantum yield against kExc, for one run or several together.

Every box takes what the command line takes: one value, a list (5,50), a:b, a:b:count or a:b:logcount. Giving a
second parameter several values gives one curve per value. Times are in ns and rates per ns.
"""

import csv
import json
import re
import sys
from pathlib import Path

from PyQt6 import QtCore, QtGui, QtWidgets  # before pyqtgraph and matplotlib, so both use this Qt

import numpy as np
import pyqtgraph as pg
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

from live import (AMBER, DIM, EDGE, FIELD, FIELDS, GO, LASER, MUTED, PANEL, TEXT, Lamp, big_axes, kexc_per_mw, label,
                  panel, style)
from plot_intensity import KINDS, add_hover, draw, is_intensity, save_all
from plot_pulse import PERYLENE
from rules import PULSE_DEFAULTS, resolve
from tta_sim import DATA, EXACT_MAX_N, HERE, parse_values

SAMPLE = ('n', 'kHop', 'f', 'tripletYield', 'tauS', 'tauSens', 'kDecay', 'fretRadius')  # the TTA model's boxes
METHODS = {'Exact solve': 'exact', 'Monte Carlo': 'mc', 'Both at every point': 'both'}
DEFAULT_RANGE = ('0.001', '100', '31', 'log')  # kExc from, to, points, spacing
LIVE_COLOURS = [PERYLENE, '#eb6834', '#1baf7a', '#eda100', '#e87ba4']  # the progress plot, one per curve, in order
POINT_LINE =re.compile(r'point (\d+)/(\d+) (\w+) kExc (\S+) photons (\S+)')


def note(text):
    widget = label(text, 'note')
    widget.setWordWrap(True)
    return widget


def default_values():
    p = resolve(PULSE_DEFAULTS)
    return {name: f'{p[name]:g}' for name in ('duration', *SAMPLE)}


def intensity_runs():
    """Every intensity run in sim/data, newest first, as (folder, params, its rows)."""
    runs = []
    for folder in DATA.iterdir() if DATA.is_dir() else []:
        if folder.is_dir() and is_intensity(folder):
            params = json.loads((folder / 'params.json').read_text())
            try:
                with (folder / 'results.csv').open(newline='') as source:
                    rows = list(csv.DictReader(source))
            except OSError:
                rows = []
            runs.append((folder, params, rows))
    return sorted(runs, key=lambda run: run[1].get('date', ''), reverse=True)


class PlotWindow(QtWidgets.QMainWindow):
    """One matplotlib figure in its own window, with the usual zoom and pan toolbar, and a readout on hover."""

    def __init__(self, kind, folders, logx, per_mw=None):
        super().__init__()
        self.setPalette(QtWidgets.QApplication.style().standardPalette())  # light, so the toolbar picks dark icons
        names = ', '.join(folder.name for folder in folders)
        self.setWindowTitle(f'{KINDS[kind]}: {names}')
        self.resize(1100, 760 if kind == 'counts_log' else 640)
        figure = Figure(dpi=100)
        canvas = FigureCanvasQTAgg(figure)
        add_hover(figure, draw(figure, kind, folders, logx, per_mw))
        holder = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(holder)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(NavigationToolbar2QT(canvas, self))
        layout.addWidget(canvas, 1)
        self.setCentralWidget(holder)
        self.setStyleSheet('QWidget { background: white; color: #52514e; }')


class Window(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Intensity sweep')
        self.resize(1300, 860)
        self.process = None
        self.folder = None  # the run being written
        self.expected = 0  # rows the running sweep will write
        self.written = 0
        self.started = None
        self.plots = []  # open plot windows, kept so they aren't closed straight away

        tabs = QtWidgets.QTabWidget()
        tabs.addTab(self.build_run(), 'Run')
        tabs.addTab(self.build_plot_tab(), 'Plot')
        self.tabs = tabs
        root = QtWidgets.QWidget()
        root.setObjectName('root')
        whole = QtWidgets.QVBoxLayout(root)
        whole.setContentsMargins(16, 12, 16, 14)
        header = QtWidgets.QHBoxLayout()
        title = QtWidgets.QVBoxLayout()
        title.setSpacing(0)
        title.addWidget(label('INTENSITY SWEEP', 'title'))
        title.addWidget(label('Photons and excitations per molecule, per pulse, against kExc', 'subtitle'))
        header.addLayout(title)
        header.addStretch()
        self.lamp = Lamp('IDLE')
        header.addWidget(self.lamp)
        whole.addLayout(header)
        whole.addWidget(tabs, 1)
        self.setCentralWidget(root)

        self.clock = QtCore.QTimer(self)
        self.clock.timeout.connect(self.show_progress)
        self.clock.start(500)
        self.check_boxes()

    # ---- the Run tab ----

    def build_run(self):
        self.boxes = {}
        column = QtWidgets.QWidget()
        column.setObjectName('column')
        left = QtWidgets.QVBoxLayout(column)
        left.setContentsMargins(0, 8, 6, 0)
        left.setSpacing(12)
        left.addWidget(self.build_sweep())
        left.addWidget(self.build_method())
        left.addWidget(self.build_sample())
        left.addStretch()
        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(column)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setFixedWidth(380)

        right = QtWidgets.QVBoxLayout()
        right.setContentsMargins(0, 8, 0, 0)
        right.setSpacing(10)
        frame, box = panel('PROGRESS')
        self.progress_text = label('Not running.', 'readout')
        box.addWidget(self.progress_text)
        self.bar = QtWidgets.QProgressBar()
        self.bar.setTextVisible(False)
        box.addWidget(self.bar)
        self.command_text = note('')
        self.command_text.setWordWrap(True)
        self.command_text.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        box.addWidget(self.command_text)
        right.addWidget(frame)

        pg.setConfigOptions(antialias=True)
        self.live_plot = pg.PlotWidget()
        self.live_plot.setBackground(PANEL)
        self.live_plot.setLogMode(True, True)
        self.live_plot.showGrid(x=True, y=True, alpha=0.15)
        big_axes(self.live_plot, 'Photons per molecule', 'kExc (/ns)')
        self.live_plot.getAxis('bottom').enableAutoSIPrefix(False)
        self.live_plot.getAxis('left').enableAutoSIPrefix(False)
        right.addWidget(self.live_plot, 1)

        self.log = QtWidgets.QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(2000)
        self.log.setFixedHeight(170)
        self.log.setStyleSheet(f"QPlainTextEdit {{ background: {FIELD}; border: 1px solid {EDGE}; border-radius: 6px;"
                               f" font-family: 'Cascadia Mono', Consolas, monospace; font-size: 12px; color: {MUTED}; }}")
        right.addWidget(self.log)

        tab = QtWidgets.QWidget()
        body = QtWidgets.QHBoxLayout(tab)
        body.setContentsMargins(0, 0, 0, 0)
        body.setSpacing(14)
        body.addWidget(scroll)
        body.addLayout(right, 1)
        return tab

    def box(self, grid, row, name, text, unit, value, tip):
        edit = QtWidgets.QLineEdit(value)
        edit.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight)
        edit.setFixedWidth(110)
        edit.setToolTip(tip)
        edit.textChanged.connect(self.check_boxes)
        self.boxes[name] = edit
        grid.addWidget(label(text, 'caption'), row, 0)
        grid.addWidget(edit, row, 1)
        grid.addWidget(label(unit, 'unit'), row, 2)

    def build_sweep(self):
        frame, box = panel('SWEEP OF kExc')
        grid = QtWidgets.QGridLayout()
        low, high, count, spacing = DEFAULT_RANGE
        self.box(grid, 0, 'from', 'From', '/ns', low, 'The lowest kExc')
        self.box(grid, 1, 'to', 'To', '/ns', high, 'The highest kExc')
        self.box(grid, 2, 'points', 'Points', '', count, 'How many values of kExc')
        self.spacing = QtWidgets.QComboBox()
        self.spacing.addItems(['spaced by ratio (log)', 'evenly spaced (linear)'])
        self.spacing.currentIndexChanged.connect(self.check_boxes)
        grid.addWidget(label('Spacing', 'caption'), 3, 0)
        grid.addWidget(self.spacing, 3, 1, 1, 2)
        values = default_values()
        self.box(grid, 4, 'duration', *FIELDS['duration'][:2], values['duration'],
                 'One value, or several (5,50) for one curve each')
        box.addLayout(grid)
        self.range_text = label('', 'readout')
        self.range_text.setWordWrap(True)
        box.addWidget(self.range_text)
        box.addWidget(note('Log spacing suits the log plot and the slope. For the linear plots, an even spacing '
                            'fills the axis better.'))
        return frame

    def build_method(self):
        frame, box = panel('METHOD')
        grid = QtWidgets.QGridLayout()
        self.method = QtWidgets.QComboBox()
        self.method.addItems(list(METHODS))
        self.method.currentIndexChanged.connect(self.check_boxes)
        grid.addWidget(label('Method', 'caption'), 0, 0)
        grid.addWidget(self.method, 0, 1, 1, 2)
        self.box(grid, 1, 'checks', 'Monte Carlo checks', 'points', '3',
                 'With the exact solve: also run Monte Carlo at this many points, spread over the sweep, as a check')
        self.box(grid, 2, 'precision', 'Monte Carlo precision', '%', '1',
                 'Monte Carlo adds molecules until photons per molecule is known to this relative standard error')
        self.box(grid, 3, 'max_seconds', 'Most time per point', 's', '600',
                 'Monte Carlo stops a point here even if it has not reached the precision')
        self.box(grid, 4, 'name', 'Run name', '', 'kexc-sweep',
                 'Used in the folder name: sim/data/YYYY-MM-DD_<name>')
        self.boxes['name'].setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft)
        box.addLayout(grid)
        self.method_text = label('', 'readout')
        self.method_text.setWordWrap(True)
        box.addWidget(self.method_text)

        buttons = QtWidgets.QHBoxLayout()
        self.run_button = QtWidgets.QPushButton('RUN SWEEP')
        self.run_button.setObjectName('fire')
        self.run_button.clicked.connect(self.run)
        self.stop_button = QtWidgets.QPushButton('STOP')
        self.stop_button.clicked.connect(self.stop)
        self.stop_button.setEnabled(False)
        buttons.addWidget(self.run_button, 1)
        buttons.addWidget(self.stop_button)
        box.addLayout(buttons)
        box.addWidget(note('A stopped sweep keeps every point it finished.'))
        return frame

    def build_sample(self):
        frame, box = panel('SAMPLE')
        grid = QtWidgets.QGridLayout()
        values = default_values()
        for row, name in enumerate(SAMPLE):
            text, unit, low, high = FIELDS[name]
            self.box(grid, row, name, text, unit, values[name],
                     f'From {low:g} to {high:g}. Several values (a,b or a:b:count) give one curve each')
        box.addLayout(grid)
        box.addWidget(note('Every box takes one value, a list (a,b), a:b:count or a:b:logcount.'))
        return frame

    def read(self):
        """The command-line arguments from the boxes, or raise ValueError naming the box that's wrong."""
        def values(name, text):
            try:
                found = parse_values(name, text.strip())
            except (ValueError, IndexError):
                raise ValueError(name)
            if name in FIELDS:
                low, high = FIELDS[name][2:]
                if not all(low <= v <= high for v in found):
                    raise ValueError(name)
            return found

        b = {name: edit.text().strip() for name, edit in self.boxes.items()}
        try:
            low, high, count = float(b['from']), float(b['to']), int(b['points'])
        except ValueError:
            raise ValueError('from')
        log = self.spacing.currentIndex() == 0
        if not (0 < low < high and count >= 2) and not (low == high and count == 1):
            raise ValueError('to')
        if not log and low < 0:
            raise ValueError('from')
        kexc = f"{low:g}:{high:g}:{'log' if log else ''}{count}"
        args = ['--kExc', kexc]
        counts = {'kExc': count}
        for name in ('duration', *SAMPLE):
            found = values(name, b[name])
            counts[name] = len(found)
            args += [f'--{name}', b[name]]
        method = METHODS[self.method.currentText()]
        args += ['--method', method]
        checks = 0
        if method == 'exact':
            try:
                checks = int(b['checks'])
            except ValueError:
                raise ValueError('checks')
            args += ['--checks', str(checks)]
        if method != 'exact' or checks:
            try:
                precision, seconds = float(b['precision']) / 100, float(b['max_seconds'])
            except ValueError:
                raise ValueError('precision')
            if not 0 < precision < 1:
                raise ValueError('precision')
            args += ['--precision', f'{precision:g}', '--max-seconds', f'{seconds:g}']
        if b['name']:
            args += ['--name', re.sub(r'[^\w.-]+', '-', b['name'])]
        points = int(np.prod(list(counts.values())))
        most_n = max(values('n', b['n']))
        return args, points, method, checks, most_n, [name for name, k in counts.items() if k > 1 and name != 'kExc']

    def check_boxes(self):
        for edit in self.boxes.values():
            edit.setStyleSheet('')
        try:
            args, points, method, checks, most_n, curves = self.read()
        except ValueError as wrong:
            name = str(wrong)
            if name in self.boxes:
                self.boxes[name].setStyleSheet(f'border-color: {LASER};')
            self.range_text.setText('Check the box marked in red.')
            self.run_button.setEnabled(False)
            return None
        too_big = method != 'mc' and most_n > EXACT_MAX_N
        self.boxes['checks'].setEnabled(method == 'exact')
        uses_mc = method != 'exact' or checks > 0
        self.boxes['precision'].setEnabled(uses_mc)
        self.boxes['max_seconds'].setEnabled(uses_mc)
        each = f"one curve per value of {', '.join(curves)}" if curves else 'one curve'
        self.range_text.setText(f'{points} point{"s" * (points > 1)}, {each}.')
        if too_big:
            self.method_text.setText(f'The exact solve only goes up to n = {EXACT_MAX_N}. Pick Monte Carlo.')
        elif method == 'mc':
            self.method_text.setText('Monte Carlo: each point runs until the precision or the time limit. Low kExc '
                                     'needs many molecules, since a photon needs two absorptions in one pulse.')
        else:
            self.method_text.setText('Exact: no noise at any kExc. About 5 s a point at n = 8, after a 15 s start; '
                                     'under a second at n of 6 or less. Long, bright pulses take longer.')
        self.method_text.setProperty('role', 'warn' if too_big else 'readout')
        self.method_text.style().polish(self.method_text)
        self.run_button.setEnabled(self.process is None and not too_big)
        self.command_text.setText('python sim/tta_sim.py intensity ' + ' '.join(args) if self.process is None else
                                  self.command_text.text())
        return args, points, method, checks

    def run(self):
        found = self.check_boxes()
        if not found or self.process is not None:
            return
        args, points, method, checks = found
        self.expected = points * (2 if method == 'both' else 1) + (checks if method == 'exact' else 0)
        self.written = 0
        self.folder = None
        self.started = QtCore.QElapsedTimer()
        self.started.start()
        self.log.clear()
        self.live_plot.clear()
        self.bar.setRange(0, self.expected)
        self.bar.setValue(0)
        self.command_text.setText('python sim/tta_sim.py intensity ' + ' '.join(args))

        process = QtCore.QProcess(self)
        process.setWorkingDirectory(str(HERE.parent))
        process.setProcessChannelMode(QtCore.QProcess.ProcessChannelMode.MergedChannels)
        process.readyReadStandardOutput.connect(self.read_output)
        process.finished.connect(self.finished)
        process.start(sys.executable, ['-u', str(HERE / 'tta_sim.py'), 'intensity', *args])
        self.process = process
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.lamp.set(AMBER, 'RUNNING')

    def stop(self):
        if self.process is not None:
            self.process.kill()
            self.log.appendPlainText('Stopped. The points finished so far are saved.')

    def read_output(self):
        text = bytes(self.process.readAllStandardOutput()).decode(errors='replace')
        for line in text.splitlines():
            if not line.strip():
                continue
            self.log.appendPlainText(line)
            if line.startswith('Saving to '):
                self.folder = Path(line[len('Saving to '):].strip())
            if POINT_LINE.match(line):
                self.written += 1
                self.bar.setValue(self.written)
                self.redraw_live()

    def redraw_live(self):
        """Photons per molecule against kExc for the run so far, straight from its results.csv."""
        if self.folder is None:
            return
        try:
            with (self.folder / 'results.csv').open(newline='') as source:
                rows = list(csv.DictReader(source))
        except OSError:
            return
        self.live_plot.clear()
        groups = {}
        for row in rows:
            curve = tuple(row[name] for name in ('duration', 'n', 'kHop', 'f', 'tripletYield', 'tauS', 'tauSens',
                                                 'kDecay', 'fretRadius'))
            groups.setdefault((curve, row['method']), []).append((float(row['kExc']), float(row['photons_per_molecule'])))
        curves = sorted({curve for curve, _ in groups})
        for (curve, method), points in sorted(groups.items()):
            points = [(x, y) for x, y in sorted(points) if y > 0]
            if not points:
                continue
            x, y = np.array(points).T
            exact = method == 'exact'
            colour = LIVE_COLOURS[curves.index(curve) % len(LIVE_COLOURS)]  # Monte Carlo shares its curve's colour
            self.live_plot.plot(x, y, pen=pg.mkPen(colour, width=2) if exact else None, symbol='o' if exact else 't',
                                symbolSize=6 if exact else 10, symbolBrush=colour if exact else None,
                                symbolPen=pg.mkPen(colour, width=1.5))

    def show_progress(self):
        if self.process is None or self.started is None:
            return
        seconds = self.started.elapsed() / 1000
        text = f'{self.written} of {self.expected} rows, {seconds:.0f} s'
        if self.written:
            left = seconds / self.written * (self.expected - self.written)
            text += f', about {left:.0f} s left' if left < 120 else f', about {left / 60:.0f} min left'
        self.progress_text.setText(text)

    def finished(self, code, status):
        stopped = status == QtCore.QProcess.ExitStatus.CrashExit
        self.process = None
        self.stop_button.setEnabled(False)
        self.check_boxes()
        seconds = self.started.elapsed() / 1000 if self.started else 0
        if stopped:
            self.lamp.set(AMBER, 'STOPPED')
            self.progress_text.setText(f'Stopped after {self.written} of {self.expected} rows, {seconds:.0f} s.')
        elif code == 0:
            self.lamp.set(GO, 'DONE')
            self.progress_text.setText(f'Done: {self.written} rows in {seconds:.0f} s. Saved to {self.folder.name}.')
        else:
            self.lamp.set(LASER, 'FAILED')
            self.progress_text.setText('The sweep failed. The log below says why.')
        self.refresh_runs(select=self.folder)

    # ---- the Plot tab ----

    def build_plot_tab(self):
        tab = QtWidgets.QWidget()
        layout = QtWidgets.QHBoxLayout(tab)
        layout.setContentsMargins(0, 8, 0, 0)
        layout.setSpacing(14)

        frame, box = panel('INTENSITY RUNS')
        self.runs = QtWidgets.QTreeWidget()
        self.runs.setHeaderLabels(['Run', 'Date', 'Rows', 'Swept', 'Method', 'n', 'Pulse, ns'])
        self.runs.setRootIsDecorated(False)
        self.runs.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        self.runs.itemDoubleClicked.connect(lambda *_: self.open_plot('counts_log'))
        self.runs.setStyleSheet(f"QTreeWidget {{ background: {FIELD}; border: 1px solid {EDGE}; border-radius: 6px; }}"
                                f"QHeaderView::section {{ background: {PANEL}; color: {MUTED}; border: none;"
                                f" padding: 4px 6px; }}"
                                f"QTreeWidget::item {{ padding: 3px; }}"
                                f"QTreeWidget::item:selected {{ background: #1f3346; color: {TEXT}; }}")
        box.addWidget(self.runs, 1)
        box.addWidget(note('Select one run, or several (Ctrl or Shift) to draw them together. Double-click opens '
                            'the log plot.'))
        layout.addWidget(frame, 1)

        side = QtWidgets.QWidget()
        side.setFixedWidth(340)
        column = QtWidgets.QVBoxLayout(side)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(12)
        frame, box = panel('PLOTS')
        for kind, text in KINDS.items():
            button = QtWidgets.QPushButton(text)
            button.clicked.connect(lambda _, kind=kind: self.open_plot(kind))
            box.addWidget(button)
        every = QtWidgets.QPushButton(f'All {len(KINDS)}')
        every.clicked.connect(lambda: [self.open_plot(kind) for kind in KINDS])
        box.addWidget(every)
        self.logx = QtWidgets.QCheckBox('Log x axis on the linear plots')
        box.addWidget(self.logx)
        self.power = QtWidgets.QCheckBox('Laser power on the x axis, not kExc')
        self.power.setToolTip('Power in mW = kExc / (kExc per mW). The top axis then shows kExc')
        self.power.toggled.connect(self.check_per_mw)
        box.addWidget(self.power)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(label('kExc per mW', 'caption'))
        self.per_mw = QtWidgets.QLineEdit(f'{kexc_per_mw():.4g}')
        self.per_mw.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight)
        self.per_mw.setFixedWidth(90)
        self.per_mw.setToolTip('The default is the live window\'s, worked out from the beam settings at the top of '
                               'live.py: a Gaussian beam at its peak intensity 2P/(πw²). Averaged over the 1/e² spot '
                               'it would be half this.')
        self.per_mw.textChanged.connect(self.check_per_mw)
        row.addWidget(self.per_mw)
        row.addWidget(label('/ns', 'unit'))
        row.addStretch()
        box.addLayout(row)
        box.addWidget(note('The default is the live window\'s, from the beam settings at the top of live.py. Hover '
                           'over the box for more.'))
        column.addWidget(frame)

        frame, box = panel('FILES')
        save = QtWidgets.QPushButton('Save PNGs to each run')
        save.clicked.connect(self.save_pngs)
        box.addWidget(save)
        reveal = QtWidgets.QPushButton('Open the run folder')
        reveal.clicked.connect(self.open_folder)
        box.addWidget(reveal)
        refresh = QtWidgets.QPushButton('Refresh the list')
        refresh.clicked.connect(lambda: self.refresh_runs())
        box.addWidget(refresh)
        box.addWidget(note('Each run is in sim/data. results.csv has the raw counts per molecule.'))
        column.addWidget(frame)
        self.plot_status = label('', 'readout')
        self.plot_status.setWordWrap(True)
        column.addWidget(self.plot_status)
        column.addStretch()
        layout.addWidget(side)
        self.refresh_runs()
        return tab

    def refresh_runs(self, select=None):
        self.runs.clear()
        for folder, params, rows in intensity_runs():
            swept = params.get('swept', {})

            def used(name):  # the values the run actually used, from its rows
                return ','.join(f'{v:g}' for v in sorted({float(row[name]) for row in rows}))
            item = QtWidgets.QTreeWidgetItem([
                folder.name, params.get('date', '')[:16], str(len(rows)),
                ', '.join(name for name in swept if name != 'excRatio'), params.get('method', ''),
                used('n'), used('duration')])
            item.setData(0, QtCore.Qt.ItemDataRole.UserRole, str(folder))
            self.runs.addTopLevelItem(item)
            if select is not None and folder == select:
                item.setSelected(True)
                self.runs.setCurrentItem(item)
        for k in range(self.runs.columnCount()):
            self.runs.resizeColumnToContents(k)
        if select is not None and self.runs.selectedItems():
            self.tabs.setCurrentIndex(1)

    def selected(self):
        return [Path(item.data(0, QtCore.Qt.ItemDataRole.UserRole)) for item in self.runs.selectedItems()]

    def open_plot(self, kind):
        folders = self.selected()
        if not folders:
            self.plot_status.setText('Select a run first.')
            return
        try:
            window = PlotWindow(kind, folders, self.logx.isChecked(), self.chosen_per_mw())
        except (ValueError, OSError, KeyError) as problem:
            self.plot_status.setText(f"Can't draw that: {problem}")
            return
        self.plot_status.setText('')
        self.plots = [w for w in self.plots if w.isVisible()] + [window]
        window.show()

    def chosen_per_mw(self):
        """kExc per mW when the power axis is on and the box holds a positive number, otherwise None."""
        if not self.power.isChecked():
            return None
        try:
            value = float(self.per_mw.text())
        except ValueError:
            return None
        return value if value > 0 else None

    def check_per_mw(self):
        wrong = self.power.isChecked() and self.chosen_per_mw() is None
        self.per_mw.setStyleSheet(f'border-color: {LASER};' if wrong else '')
        self.plot_status.setText('kExc per mW has to be a number above 0.' if wrong else '')

    def save_pngs(self):
        saved = []
        for folder in self.selected():
            try:
                saved.append(save_all(folder, self.chosen_per_mw()).parent.name)
            except (ValueError, OSError) as problem:
                self.plot_status.setText(f'{folder.name}: {problem}')
                return
        self.plot_status.setText(f"Saved plots/ in {', '.join(saved)}" if saved else 'Select a run first.')

    def open_folder(self):
        folders = self.selected()
        if folders:
            QtGui.QDesktopServices.openUrl(QtCore.QUrl.fromLocalFile(str(folders[0])))

    def closeEvent(self, event):
        if self.process is not None:
            answer = QtWidgets.QMessageBox.question(
                self, 'A sweep is running', 'Stop the sweep and close? The points finished so far are saved.')
            if answer != QtWidgets.QMessageBox.StandardButton.Yes:
                event.ignore()
                return
            self.process.kill()
            self.process.waitForFinished(2000)
        for window in self.plots:
            window.close()
        event.accept()


def main():
    app = QtWidgets.QApplication(sys.argv)
    style(app)
    app.setStyleSheet(app.styleSheet() + f"""
        QTabWidget::pane {{ border: none; }}
        QTabBar::tab {{ background: transparent; color: {MUTED}; padding: 8px 18px; font-weight: 600;
                        letter-spacing: 2px; border-bottom: 2px solid transparent; }}
        QTabBar::tab:selected {{ color: {TEXT}; border-bottom-color: {PERYLENE}; }}
        QPushButton:disabled {{ color: {DIM}; border-color: {EDGE}; }}
        QPushButton#fire:disabled {{ background: {PANEL}; border-color: {EDGE}; color: {DIM}; }}
        QLineEdit:disabled {{ color: {DIM}; }}
    """)
    window = Window()
    window.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
