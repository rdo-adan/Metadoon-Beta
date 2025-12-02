#!/bin/bash

# Setup.sh
# Script to create and configure the Metadoon Conda environment

echo "=== Step 1: Creating Conda environment from metadoon_env.yaml ==="

# Verifica se o conda está instalado
if ! command -v conda &> /dev/null; then
    echo "Conda is not installed or not in PATH. Please install Miniconda or Anaconda first."
    exit 1
fi

# Dá preferência ao mamba se estiver disponível (é mais rápido)
if command -v mamba &> /dev/null; then
    echo "[INFO] mamba detected. Creating environment with mamba..."
    mamba env create -f metadoon_env.yaml -n metadoon
else
    echo "[INFO] mamba not found. Using conda to create environment..."
    conda env create -f metadoon_env.yaml -n metadoon
fi

# Inicializa o conda dentro do script e ativa o ambiente
# IMPORTANTE: Isso só ativa o ambiente DURANTE a execução deste script.
eval "$(conda shell.bash hook)"
conda activate metadoon

echo "=== Step 2: Installing additional R packages from GitHub ==="

# Instala devtools se não existir e depois os pacotes do GitHub
Rscript -e "if (!requireNamespace('devtools', quietly = TRUE)) install.packages('devtools', repos='https://cloud.r-project.org')"
Rscript -e "devtools::install_github('pmartinezarbizu/pairwiseAdonis/pairwiseAdonis')"
Rscript -e "devtools::install_github('microbiome/microbiome')"

echo "✅ Metadoon environment setup complete."