@echo off
TITLE Metadoon Launcher (WSLg Mode)
setlocal enabledelayedexpansion

:: --- CONFIGURATION ---
set IMAGE_NAME=engbio/metadoon:v1.0

echo ==========================================
echo      Metadoon (Docker via WSLg)
echo ==========================================

:: 1. Check if Docker is running
docker info >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Docker is not running. Please start Docker Desktop.
    pause
    exit /b
)

:: 2. Navigate to project root
cd /d "%~dp0.."

:: 3. Convert current Windows path (C:\...) to WSL path (/mnt/c/...)
:: This is required because we are running the command "inside" WSL context
for /f "usebackq tokens=*" %%a in (`wsl wslpath -a "%cd%"`) do set WSL_CURRENT_DIR=%%a

:: 4. Pull/Update the image
echo.
echo [INFO] Checking for updates...
docker pull %IMAGE_NAME%

:: 5. Run Container using WSLg (Native Windows Graphics)
echo.
echo [INFO] Launching Metadoon...
echo [INFO] Using Windows Native Graphics (WSLg)...

:: THE TRICK:
:: We use 'wsl -e' to run the docker command FROM WITHIN Linux (WSL).
:: We map /tmp/.X11-unix and /mnt/wslg to utilize Windows 11 native graphics.

wsl -e docker run --rm -it ^
  -v /tmp/.X11-unix:/tmp/.X11-unix ^
  -v /mnt/wslg:/mnt/wslg ^
  -e DISPLAY=:0 ^
  -e WAYLAND_DISPLAY=wayland-0 ^
  -e XDG_RUNTIME_DIR=/mnt/wslg/runtime-dir ^
  -v "%WSL_CURRENT_DIR%":/app ^
  -w /app ^
  %IMAGE_NAME%

if %errorlevel% neq 0 (
    echo.
    echo [ERROR] Execution failed.
    echo.
    echo Troubleshooting:
    echo 1. Ensure you have Windows 10 Build 19044+ or Windows 11.
    echo 2. Ensure 'wsl --update' has been run in PowerShell.
    echo 3. Docker Desktop must have 'Use WSL 2 based engine' enabled.
    pause
)