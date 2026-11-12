#!/bin/sh
# Double-click to open the launcher menu (macOS). On Linux, run: sh Launch.command
cd "$(dirname "$0")" || exit 1
exec python3 launch.py "$@"
