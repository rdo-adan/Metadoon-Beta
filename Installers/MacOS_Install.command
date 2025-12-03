#!/bin/bash
# Pega o diretório onde o arquivo está e volta 2 níveis para a raiz
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR/../../"

echo "=== Metadoon Installer (macOS) ==="

# Verifica Conda
if ! command -v conda &> /dev/null; then
    # Tenta achar o conda nos locais padrões se não estiver no PATH
    if [ -f "$HOME/opt/anaconda3/bin/conda" ]; then
        source "$HOME/opt/anaconda3/bin/activate"
    elif [ -f "$HOME/miniconda3/bin/conda" ]; then
        source "$HOME/miniconda3/bin/activate"
    else
        echo "Error: Conda not found. Please install Anaconda first."
        exit 1
    fi
fi

# Roda o script principal
bash setup.sh

echo "✅ Installation Finished. You can close this window."