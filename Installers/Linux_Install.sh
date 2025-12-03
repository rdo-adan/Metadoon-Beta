#!/bin/bash
cd "$(dirname "$0")/../"

echo "=== Metadoon Installer (Linux) ==="

if [ -f "setup.sh" ]; then
    bash setup.sh
else
    echo "Error: setup.sh not found in $(pwd)"
fi

read -p "Press Enter to exit..."