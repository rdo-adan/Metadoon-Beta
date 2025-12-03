#!/bin/bash
cd "$(dirname "$0")/../"

echo "=== Launching Metadoon ==="

if [ ! -f "metadoon.py" ]; then
    echo "Error: metadoon.py not found in $(pwd)"
    read -p "Press Enter to exit..."
    exit 1
fi


eval "$(conda shell.bash hook)"
conda activate metadoon

python metadoon.py