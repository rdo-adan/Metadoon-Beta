@echo off
TITLE Metadoon Launcher (WSLg Mode)
setlocal enabledelayedexpansion

:: --- CONFIGURATION ---
set IMAGE_NAME=engbio/metadoon:v1.0

echo ==========================================
echo      Metadoon (Docker via WSLg)
echo ==========================================

:: 1. Check Docker Status
:: Verifies if Docker Desktop is running before proceeding
docker info >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Docker is not running. Please start Docker Desktop.
    pause
    exit /b
)

:: 2. Navigate to root directory
:: Sets the working directory to the script's location
cd /d "%~dp0"

:: 3. Convert Windows paths to WSL (Linux) format
:: This conversion is necessary for Docker to correctly map Windows volumes
for /f "usebackq tokens=*" %%a in (`wsl wslpath -a "%cd%"`) do set WSL_CURRENT_DIR=%%a
for /f "usebackq tokens=*" %%a in (`wsl wslpath -a "%UserProfile%"`) do set WSL_USER_PROFILE=%%a

:: 4. Update Image (Optional)
:: Pulls the latest version from Docker Hub. Comment out for offline mode.
echo.
echo [INFO] Checking for updates...
docker pull %IMAGE_NAME%

:: 5. Execute Container
echo.
echo [INFO] Launching Metadoon...
echo [INFO] Current directory mapped to: /workspace
echo [INFO] User Profile mapped to: /app/YOUR_DATA

:: --- DOCKER RUN COMMAND ---
:: --net=host: Uses the host network stack for better connectivity
:: -v /tmp/.X11-unix... & /mnt/wslg...: Graphic drivers mapping for WSLg (GUI support)
:: -v "%WSL_CURRENT_DIR%": Maps the current Windows folder to /workspace inside the container
:: -v "%WSL_USER_PROFILE%": Maps the Windows User folder to /app/YOUR_DATA for easy access
:: bash -c ...: Changes directory to /app and executes the python script

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
  bash -c "cd /app && python metadoon.py"

:: Check for execution errors
if %errorlevel% neq 0 (
    echo.
    echo [ERROR] Execution failed.
    pause
)