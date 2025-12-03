#!/bin/bash

# Setup.sh
# Script to create and configure the Metadoon Conda environment

# Exit immediately if a command exits with a non-zero status
set -e

echo "=== Step 1: Creating Conda environment from metadoon_env.yaml ==="

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    echo "Error: Conda is not installed or not in PATH. Please install Miniconda or Anaconda first."
    exit 1
fi

# Prefer Mamba if available (it is faster)
if command -v mamba &> /dev/null; then
    echo "[INFO] Mamba detected. Creating environment with mamba..."
    mamba env create -f metadoon_env.yaml -n metadoon
else
    echo "[INFO] Mamba not found. Using conda to create environment..."
    conda env create -f metadoon_env.yaml -n metadoon
fi

# Initialize Conda within the script and activate the environment
# IMPORTANT: This only activates the environment DURING the execution of this script.
# The user will still need to run 'conda activate metadoon' manually after this script finishes.
echo "[INFO] Activating 'metadoon' environment for configuration..."
eval "$(conda shell.bash hook)"
conda activate metadoon

echo "=== Step 2: Installing additional R packages from GitHub ==="

# Note: Most packages (phyloseq, deseq2, microbiome, etc.) are already installed via Conda.
# We only use Rscript here to install 'pairwiseAdonis' because it is not available in Conda channels.

echo "[INFO] Installing pairwiseAdonis..."
Rscript -e "if (!requireNamespace('devtools', quietly = TRUE)) install.packages('devtools', repos='https://cloud.r-project.org')"
Rscript -e "devtools::install_github('pmartinezarbizu/pairwiseAdonis/pairwiseAdonis')"

echo "✅ Metadoon environment setup complete."
echo "👉 To start using the tool, please run: conda activate metadoon"
