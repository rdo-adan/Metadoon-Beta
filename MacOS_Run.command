#!/bin/bash

# ==========================================
#      Metadoon Launcher (Native Conda)
# ==========================================

# 1. Navigate to the script's directory (Essential for relative paths)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"

echo "=========================================="
echo "      Metadoon Launcher"
echo "=========================================="
echo "[INFO] Working directory: $DIR"

# 2. Locate and Initialize Conda
# Conda is not automatically available in shell scripts. We must find 'conda.sh'.
# Common paths for Anaconda/Miniconda on macOS (Intel and Apple Silicon)
CONDA_PATHS=(
    "$HOME/opt/anaconda3/etc/profile.d/conda.sh"
    "$HOME/opt/miniconda3/etc/profile.d/conda.sh"
    "$HOME/anaconda3/etc/profile.d/conda.sh"
    "$HOME/miniconda3/etc/profile.d/conda.sh"
    "/opt/anaconda3/etc/profile.d/conda.sh"
    "/opt/miniconda3/etc/profile.d/conda.sh"
    "/usr/local/anaconda3/etc/profile.d/conda.sh"
    "/usr/local/miniconda3/etc/profile.d/conda.sh"
    "/opt/homebrew/anaconda3/etc/profile.d/conda.sh"
    "/opt/homebrew/miniconda3/etc/profile.d/conda.sh"
    "/opt/homebrew/Caskroom/miniconda/base/etc/profile.d/conda.sh"
)

CONDA_FOUND=false

for path in "${CONDA_PATHS[@]}"; do
    if [ -f "$path" ]; then
        echo "[INFO] Found Conda at: $path"
        source "$path"
        CONDA_FOUND=true
        break
    fi
done

# If Conda was not found in standard paths, try to assume it's in the PATH
if [ "$CONDA_FOUND" = false ]; then
    if command -v conda &> /dev/null; then
        echo "[INFO] Conda found in system PATH."
        eval "$(conda shell.bash hook)"
    else
        echo ""
        echo "[ERROR] Conda not found."
        echo "Please install Anaconda or Miniconda first."
        echo "Download: https://docs.conda.io/en/latest/miniconda.html"
        read -p "Press Enter to exit..."
        exit 1
    fi
fi

# 3. Check for 'metadoon' environment
echo "[INFO] Checking Conda environment 'metadoon'..."

if conda env list | grep -q "^metadoon "; then
    echo "[INFO] Environment 'metadoon' exists."
else
    echo "[WARN] Environment 'metadoon' NOT found."
    
    if [ -f "./setup.sh" ]; then
        echo "[INFO] Running setup.sh to create the environment..."
        chmod +x ./setup.sh
        ./setup.sh
        
        # Verify if creation was successful
        if conda env list | grep -q "^metadoon "; then
            echo "[SUCCESS] Environment created successfully."
        else
            echo "[ERROR] setup.sh failed to create the environment."
            read -p "Press Enter to exit..."
            exit 1
        fi
    else
        echo "[ERROR] 'setup.sh' not found in this directory."
        echo "Cannot create the environment automatically."
        read -p "Press Enter to exit..."
        exit 1
    fi
fi

# 4. Activate Environment and Run
echo "[INFO] Activating 'metadoon'..."
conda activate metadoon

echo "[INFO] Launching Metadoon.py..."
echo "------------------------------------------"

if [ -f "Metadoon.py" ]; then
    python Metadoon.py
else
    echo "[ERROR] Metadoon.py not found in current directory."
    read -p "Press Enter to exit..."
    exit 1
fi

echo "------------------------------------------"
echo "[INFO] Process finished."
# Wait for user input if it crashes immediately, so they can read the error
read -p "Press Enter to close..."