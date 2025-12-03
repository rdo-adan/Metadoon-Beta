#!/bin/bash
cd "$(dirname "$0")/../../"

# Tenta inicializar o conda para o script
eval "$(conda shell.bash hook)"
conda activate metadoon

python metadoon.py