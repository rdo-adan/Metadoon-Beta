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
    mamba env create -f metadoon_env.yaml -n metadoon
else
    echo "[INFO] mamba not found. Using conda to create environment..."
    conda env create -f metadoon_env.yaml -n metadoon
fi

# Activate environment (use correct shell hook for bash/zsh)
eval "$(conda shell.$(basename $SHELL) hook)"
conda activate metadoon

echo "=== Step 2: Installing additional R packages from GitHub ==="

# Install devtools if not present, then install GitHub packages
Rscript -e "if (!requireNamespace('devtools', quietly = TRUE)) install.packages('devtools', repos='https://cloud.r-project.org')"
Rscript -e "devtools::install_github('pmartinezarbizu/pairwiseAdonis/pairwiseAdonis')"
Rscript -e "devtools::install_github('microbiome/microbiome')"

echo "✅ Metadoon environment setup complete."
