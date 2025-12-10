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
for /f "usebackq tokens=*" %%a in (`wsl wslpath -a "%cd%"`) do set WSL_CURRENT_DIR=%%a
for /f "usebackq tokens=*" %%a in (`wsl wslpath -a "%UserProfile%"`) do set WSL_USER_PROFILE=%%a

:: 4. Update Image (Optional)
echo.
echo [INFO] Checking for updates...
docker pull %IMAGE_NAME%

:: 5. Execute Container
echo.
echo [INFO] Launching Metadoon...
echo [TIP] Inside the app:
echo        - Go to "/app/YOUR_DATA" to see your User folders (Downloads, Documents).
echo        - Go to "/app/C_Drive" to see your entire C: disk.

:: --- DOCKER RUN COMMAND ---

wsl -e docker run --rm -it ^
  --net=host ^
  -v /tmp/.X11-unix:/tmp/.X11-unix ^
  -v /mnt/wslg:/mnt/wslg ^
  -e DISPLAY=:0 ^
  -e WAYLAND_DISPLAY=wayland-0 ^
  -e XDG_RUNTIME_DIR=/mnt/wslg/runtime-dir ^
  -v "%WSL_CURRENT_DIR%":/workspace:rw ^
  -v "%WSL_USER_PROFILE%":/app/YOUR_DATA ^
  -v /mnt/c:/app/C_Drive ^
  -v /mnt/d:/app/D_Drive ^
  --workdir /workspace ^
  %IMAGE_NAME% ^
  python /app/metadoon.py

if %errorlevel% neq 0 (
    echo.
    echo [ERROR] Execution failed.
    pause
)