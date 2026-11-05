# Creates a "Molecule Viewer" shortcut on the Desktop that opens the built app
# (dist/index.html) in its own Chrome/Edge app window, with no tabs or address bar, and no server.
#
# Run with:  npm run shortcut
# Re-run it if you move the project folder.

$ErrorActionPreference = 'Stop'

$root = Split-Path -Parent $PSScriptRoot
$html = Join-Path $root 'dist\index.html'
$icon = Join-Path $root 'assets\icon.ico'

if (-not (Test-Path $html)) {
    throw "dist\index.html not found. Run 'npm run build' first."
}

$browsers = @(
    "$env:ProgramFiles\Google\Chrome\Application\chrome.exe",
    "${env:ProgramFiles(x86)}\Google\Chrome\Application\chrome.exe",
    "$env:LOCALAPPDATA\Google\Chrome\Application\chrome.exe",
    "${env:ProgramFiles(x86)}\Microsoft\Edge\Application\msedge.exe",
    "$env:ProgramFiles\Microsoft\Edge\Application\msedge.exe"
)
$browser = $browsers | Where-Object { Test-Path $_ } | Select-Object -First 1
if (-not $browser) { throw 'Could not find Chrome or Edge.' }

$url = ([System.Uri]$html).AbsoluteUri   # file:///D:/Yaakov%20Vault/...
$desktop = [Environment]::GetFolderPath('Desktop')
$lnkPath = Join-Path $desktop 'Molecule Viewer.lnk'

$shell = New-Object -ComObject WScript.Shell
$lnk = $shell.CreateShortcut($lnkPath)
$lnk.TargetPath = $browser
$lnk.Arguments = "--app=`"$url`" --window-size=1400,900"
$lnk.WorkingDirectory = $root
$lnk.IconLocation = "$icon,0"
$lnk.Description = 'Molecule Viewer: 3D viewer for .xyz files'
$lnk.Save()

Write-Host "Created $lnkPath"
Write-Host "  -> $browser"
Write-Host "  -> $url"
