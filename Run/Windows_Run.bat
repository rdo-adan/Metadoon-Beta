@echo off
TITLE Metadoon Launcher (WSLg Mode)
setlocal enabledelayedexpansion

:: --- CONFIGURATION ---
set IMAGE_NAME=engbio/metadoon:v1.0

echo ==========================================
echo      Metadoon (Docker via WSLg)
echo ==========================================

:: 1. Check Docker Status
docker info >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Docker is not running. Please start Docker Desktop.
    pause
    exit /b
)

:: 2. Navigate to script directory
cd /d "%~dp0"

:: 3. Convert Windows paths to WSL format
:: Crucial for mounting volumes correctly in Docker via WSL
for /f "usebackq tokens=*" %%a in (`wsl wslpath -a "%cd%"`) do set WSL_CURRENT_DIR=%%a
for /f "usebackq tokens=*" %%a in (`wsl wslpath -a "%UserProfile%"`) do set WSL_USER_PROFILE=%%a

:: 4. Update Image (Optional - comment if offline)
echo.
echo [INFO] Checking for updates...
docker pull %IMAGE_NAME%

:: 5. Execute Container
echo.
echo [INFO] Launching Metadoon...
echo [INFO] Workspace mapped to: %WSL_CURRENT_DIR%

:: --- DOCKER RUN COMMAND ---
:: -v /tmp/.X11-unix...: Graphic drivers mapping for WSLg
:: -v "%WSL_CURRENT_DIR%":/workspace: Maps current folder to container workspace
:: -w /workspace: Sets working directory
:: python /app/metadoon.py: Executes the specific entry command

wsl -e docker run --rm -it ^
  --net=host ^
  -v /tmp/.X11-unix:/tmp/.X11-unix ^
  -v /mnt/wslg:/mnt/wslg ^
  -e DISPLAY=:0 ^
  -e WAYLAND_DISPLAY=wayland-0 ^
  -e XDG_RUNTIME_DIR=/mnt/wslg/runtime-dir ^
  -v "%WSL_CURRENT_DIR%":/workspace:rw ^
  -v "%WSL_USER_PROFILE%":/app/YOUR_DATA ^
  --workdir /workspace ^
  %IMAGE_NAME% ^
  python /app/metadoon.py

if %errorlevel% neq 0 (
    echo.
    echo [ERROR] Execution failed.
    pause
)