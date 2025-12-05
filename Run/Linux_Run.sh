@echo off
TITLE Metadoon Launcher (WSLg Mode)
setlocal enabledelayedexpansion

:: --- CONFIGURAÇÃO ---
set IMAGE_NAME=engbio/metadoon:v1.0

echo ==========================================
echo      Metadoon (Docker via WSLg)
echo ==========================================

:: 1. Verifica se o Docker está rodando
docker info >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Docker is not running. Please start Docker Desktop.
    pause
    exit /b
)

:: 2. Navega para a raiz do projeto
cd /d "%~dp0.."

:: 3. Converte o caminho atual do Windows (C:\...) para WSL (/mnt/c/...)
:: Isso é necessário porque vamos rodar o comando "dentro" do Linux
for /f "usebackq tokens=*" %%a in (`wsl wslpath -a "%cd%"`) do set WSL_CURRENT_DIR=%%a

:: 4. Baixa/Atualiza a imagem
echo.
echo [INFO] Checking for updates...
docker pull %IMAGE_NAME%

:: 5. Executa o Container usando WSLg (Sem VcXsrv)
echo.
echo [INFO] Launching Metadoon...
echo [INFO] Using Windows Native Graphics (WSLg)...

:: O TRUQUE ESTÁ AQUI:
:: Usamos 'wsl -e' para rodar o docker DE DENTRO do Linux.
:: Mapeamos /tmp/.X11-unix e /mnt/wslg para usar o gráfico nativo do Windows.

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