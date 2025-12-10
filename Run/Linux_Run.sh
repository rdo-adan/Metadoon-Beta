#!/bin/bash

# ==========================================
#      Metadoon Launcher (Linux)
# ==========================================

IMAGE_NAME="engbio/metadoon:v1.0"

# 1. Navigate to script directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"

echo "=========================================="
echo "      Metadoon (Linux Docker)"
echo "=========================================="

# 2. Check Docker
if ! command -v docker &> /dev/null; then
    echo "[ERROR] Docker not found. Please install Docker Engine."
    exit 1
fi

# 3. Configure X11 Permissions
echo "[INFO] Configuring X11 permissions..."
xhost +local:docker > /dev/null 2>&1

# 4. Update Image
echo "[INFO] Checking for updates..."
docker pull $IMAGE_NAME

# 5. Execute Container
echo ""
echo "[INFO] Launching Metadoon..."

# --- DOCKER RUN COMMAND ---
# --user $(id -u):$(id -g): Runs as current user to fix file permissions
# -v "$(pwd)":/workspace:rw: Maps current directory
# python /app/metadoon.py: Main execution command

docker run --rm -it \
    --net=host \
    -e DISPLAY=$DISPLAY \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    -v "$(pwd)":/workspace:rw \
    -v "$HOME":/app/YOUR_DATA \
    --user $(id -u):$(id -g) \
    --workdir /workspace \
    $IMAGE_NAME \
    python /app/metadoon.py

if [ $? -ne 0 ]; then
    echo ""
    echo "[ERROR] Execution failed."
    exit 1
fi