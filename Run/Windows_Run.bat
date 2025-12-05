@echo off
TITLE Metadoon Launcher (WSLg Mode)
setlocal enabledelayedexpansion

:: --- CONFIGURAÇÃO ---
set IMAGE_NAME=engbio/metadoon:v1.0

echo ==========================================
echo      Metadoon (Docker via WSLg)
echo ==========================================

:: 1. Verifica Docker
docker info >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Docker is not running. Please start Docker Desktop.
    pause
    exit /b
)

:: 2. Navega para raiz
cd /d "%~dp0.."

:: 3. Converte caminhos para WSL
for /f "usebackq tokens=*" %%a in (`wsl wslpath -a "%cd%"`) do set WSL_CURRENT_DIR=%%a
for /f "usebackq tokens=*" %%a in (`wsl wslpath -a "%UserProfile%"`) do set WSL_USER_PROFILE=%%a

:: 4. Atualiza Imagem
echo.
echo [INFO] Checking for updates...
docker pull %IMAGE_NAME%

:: 5. Executa
echo.
echo [INFO] Launching Metadoon...
echo [TIP] Your Windows files are inside the folder "YOUR_DATA" in the file selection window.

:: --- MUDANÇA AQUI: ---
:: Adicionamos -v "%WSL_USER_PROFILE%":/app/YOUR_DATA
:: Isso cria um atalho direto para seus documentos dentro da janela do programa.

wsl -e docker run --rm -it ^
  -v /tmp/.X11-unix:/tmp/.X11-unix ^
  -v /mnt/wslg:/mnt/wslg ^
  -e DISPLAY=:0 ^
  -e WAYLAND_DISPLAY=wayland-0 ^
  -e XDG_RUNTIME_DIR=/mnt/wslg/runtime-dir ^
  -v "%WSL_CURRENT_DIR%":/app ^
  -v "%WSL_USER_PROFILE%":/app/YOUR_DATA ^
  -w /app ^
  %IMAGE_NAME%

if %errorlevel% neq 0 (
    echo.
    echo [ERROR] Execution failed.
    pause
)