"""One launcher for every part of the project. Needs only Python 3.10 or newer.

  python launch.py            a menu
  python launch.py viewer     the 3D molecule viewer, in your browser
  python launch.py live       the live pulse experiment window
  python launch.py sweep      the intensity sweep window

The first time a pulse window is opened, this makes a private Python environment in .venv/ next to this file and
installs requirements.txt into it (a few minutes, once). After that it starts right away. Delete .venv/ to start over.
The viewer needs nothing installed: it opens the online version, or your own build in dist/ if there is one.
"""

import hashlib
import os
import subprocess
import sys
import venv
import webbrowser
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VENV = ROOT / '.venv'
REQUIREMENTS = ROOT / 'requirements.txt'
STAMP = VENV / 'requirements.sha256'  # what was installed, so an edited requirements.txt is picked up
VIEWER_URL = 'https://yclen.github.io/molecule-viewer/'

APPS = {
    'viewer': ('Molecule Viewer', '3D molecules, excited states and on-molecule simulations (browser)'),
    'live': ('Live pulse experiment', 'flash molecules and watch the photon histogram build up'),
    'sweep': ('Intensity sweep', 'run photons and quantum yield against laser power, then plot them'),
}
SCRIPTS = {'live': ROOT / 'sim' / 'live.py', 'sweep': ROOT / 'sim' / 'sweep_gui.py'}


def venv_python() -> Path:
    return VENV / ('Scripts/python.exe' if os.name == 'nt' else 'bin/python')


def ensure_env() -> Path:
    """Make .venv and install requirements.txt into it, unless that's already done."""
    if sys.version_info < (3, 10):
        sys.exit(f'Python 3.10 or newer is needed (this is {sys.version.split()[0]}). Get it from python.org.')
    wanted = hashlib.sha256(REQUIREMENTS.read_bytes()).hexdigest()
    python = venv_python()
    if python.exists() and STAMP.exists() and STAMP.read_text().strip() == wanted:
        return python
    if not python.exists():
        print(f'Setting up a Python environment in {VENV} (first run only)...')
        venv.create(VENV, with_pip=True)
    print('Installing the simulation packages (numpy, scipy, numba, PyQt6, pyqtgraph, matplotlib)...')
    subprocess.run([python, '-m', 'pip', 'install', '--upgrade', 'pip'], check=True)
    subprocess.run([python, '-m', 'pip', 'install', '-r', REQUIREMENTS], check=True)
    STAMP.write_text(wanted)
    print('Done. The first run also compiles the simulation for a few seconds.\n')
    return python


def open_viewer() -> None:
    local = ROOT / 'dist' / 'index.html'
    target = local.as_uri() if local.exists() else VIEWER_URL
    print(f'Opening {target}')
    webbrowser.open(target)


def run(name: str) -> int:
    if name == 'viewer':
        open_viewer()
        return 0
    python = ensure_env()
    print(f'Starting the {APPS[name][0].lower()} window...')
    return subprocess.run([python, SCRIPTS[name]], cwd=ROOT).returncode


def menu() -> str:
    print('Molecule Viewer and TTA pulse experiments\n')
    names = list(APPS)
    for i, name in enumerate(names, 1):
        title, blurb = APPS[name]
        print(f'  {i}. {title:<24} {blurb}')
    while True:
        choice = input('\nPick 1-3 (or q to quit): ').strip().lower()
        if choice in ('q', ''):
            sys.exit(0)
        if choice in names:
            return choice
        if choice.isdigit() and 1 <= int(choice) <= len(names):
            return names[int(choice) - 1]


def main() -> None:
    args = sys.argv[1:]
    if args and args[0] not in APPS:
        sys.exit(f'Unknown part "{args[0]}". Choose one of: {", ".join(APPS)}')
    try:
        sys.exit(run(args[0] if args else menu()))
    except subprocess.CalledProcessError as err:
        sys.exit(f'\nSetup failed ({err}). Check your internet connection, delete the .venv folder and try again.')
    except (KeyboardInterrupt, EOFError):
        sys.exit(130)


if __name__ == '__main__':
    main()
