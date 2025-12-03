@echo off
TITLE Metadoon Installer
echo ==========================================
echo      Installing Metadoon for Windows
echo ==========================================

:: 1. Verifica se Conda existe
call conda --version >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Conda not found. Please install Anaconda or Miniconda first.
    pause
    exit /b
)

:: 2. Volta para a pasta raiz (onde está o yaml)
cd ..\..

:: 3. Cria o ambiente
echo [INFO] Creating Conda environment...
call conda env create -f metadoon_env.yaml -n metadoon

:: 4. Instala pacotes R extras via Rscript do ambiente
echo [INFO] Installing extra R packages...
call conda run -n metadoon Rscript -e "if (!requireNamespace('devtools', quietly = TRUE)) install.packages('devtools', repos='https://cloud.r-project.org')"
call conda run -n metadoon Rscript -e "devtools::install_github('pmartinezarbizu/pairwiseAdonis/pairwiseAdonis')"

echo.
echo [SUCCESS] Installation Complete! You can now use the launcher in the 'Run' folder.
pause