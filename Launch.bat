@echo off
rem Double-click to open the launcher menu (Windows). Needs Python 3.10+ from python.org.
cd /d "%~dp0"
where py >nul 2>nul
if %errorlevel%==0 (
  py -3 launch.py %*
) else (
  python launch.py %*
)
if errorlevel 1 pause
