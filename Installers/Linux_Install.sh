#!/bin/bash
# Garante que roda a partir da raiz
cd "$(dirname "$0")/../../"

echo "=== Metadoon Installer (Linux) ==="
bash setup.sh
read -p "Press Enter to exit..."