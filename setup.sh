#!/bin/bash

# Setup.sh
# Script to create and configure the Metadoon Conda environment

echo "=== Step 1: Creating Conda environment from metadoon_env.yaml ==="

if ! command -v conda &> /dev/null; then
    echo "Conda is not installed or not in PATH. Please install Miniconda or Anaconda first."
    exit 1
fi

# Prefer mamba if available
if command -v mamba &> /dev/null; then
    echo "[INFO] mamba detected. Creating environment with mamba..."
    mamba env create -f metadoon_env.yaml
else
    echo "[INFO] mamba not found. Using conda to create environment..."
    conda env create -f metadoon_env.yaml
fi

# Activate environment (requires conda hook for shell)
eval "$(conda shell.bash hook)"
conda activate metadoon

echo "=== Step 2: Installing additional R packages from GitHub ==="

# Install devtools if needed, then GitHub packages
Rscript -e "if (!requireNamespace('devtools', quietly = TRUE)) install.packages('devtools', repos='http://cran.us.r-project.org')"
Rscript -e "devtools::install_github('vlubitch/pairwiseAdonis')"
Rscript -e "devtools::install_github('microbiome/microbiome')"

echo "✅ Metadoon environment setup complete."
