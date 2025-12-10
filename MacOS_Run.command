#!/bin/bash

# ==========================================
#      Metadoon Launcher (macOS)
# ==========================================

IMAGE_NAME="engbio/metadoon:v1.0"

# 1. Navigate to script directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"

echo "=========================================="
echo "      Metadoon (Docker via XQuartz)"
echo "=========================================="

# 2. Check Docker
if ! command -v docker &> /dev/null; then
    echo "[ERROR] Docker not found. Please install Docker Desktop for Mac."
    exit 1
fi

# 3. Check XQuartz
if ! command -v xhost &> /dev/null; then
    echo "[ERROR] XQuartz not found. Please install it (brew install --cask xquartz)."
    exit 1
fi

# 4. Configure X11 Permissions
echo "[INFO] Configuring X11 permissions..."
xhost + 127.0.0.1 > /dev/null 2>&1

# 5. Update Image
echo "[INFO] Checking for updates..."
docker pull $IMAGE_NAME

# 6. Execute Container
echo ""
echo "[INFO] Launching Metadoon..."
echo "[TIP] Workspace mapped to: $(pwd)"

# --- DOCKER RUN COMMAND ---
# -e DISPLAY=host.docker.internal:0: Connects to XQuartz on Mac
# -v "$(pwd)":/workspace:rw: Maps current folder
# python /app/metadoon.py: Explicit command

docker run --rm -it \
    -e DISPLAY=host.docker.internal:0 \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    -v "$(pwd)":/workspace:rw \
    -v "$HOME":/app/YOUR_DATA \
    --workdir /workspace \
    $IMAGE_NAME \
    python /app/metadoon.py

if [ $? -ne 0 ]; then
    echo ""
    echo "[ERROR] Execution failed."
    exit 1
fi