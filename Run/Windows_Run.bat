@echo off
TITLE Metadoon Launcher
setlocal enabledelayedexpansion

:: --- CONFIGURAÇÃO ---
:: Coloque aqui o nome exato da sua imagem no Docker Hub
set IMAGE_NAME=engbio/metadoon:v1.0

echo ==========================================
echo      Metadoon (Docker Edition)
echo ==========================================

:: 1. Verifica se o Docker está rodando
docker info >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Docker is not running. Please start Docker Desktop.
    pause
    exit /b
)

:: 2. Navega para a raiz do projeto (onde estão os dados do usuário)
cd /d "%~dp0.."

:: 3. Baixa/Atualiza a imagem (Pull)
echo.
echo [INFO] Checking for updates (%IMAGE_NAME%)...
docker pull %IMAGE_NAME%
if %errorlevel% neq 0 (
    echo [WARNING] Could not pull latest image. Trying to run existing version...
)

:: 4. Executa o Container
echo.
echo [INFO] Launching Metadoon...
echo [NOTE] Ensure X Server (VcXsrv) is running for the GUI to appear.

docker run --rm -it ^
  -e DISPLAY=host.docker.internal:0 ^
  -v "%cd%":/app ^
  -v "%UserProfile%":/host_home:ro ^
  -w /app ^
  %IMAGE_NAME%

if %errorlevel% neq 0 (
    echo.
    echo [ERROR] Execution failed.
    pause
)