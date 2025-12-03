#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR/../../"

echo "Starting Metadoon..."

# Tenta carregar o conda profile para ativar ambientes
source ~/opt/anaconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null

conda activate metadoon
python metadoon.py