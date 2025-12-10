#!/bin/bash

# ==========================================
#      Metadoon Launcher (Linux)
# ==========================================

IMAGE_NAME="engbio/metadoon:v1.0"

# 1. Navigate to the script's directory
# Ensures the script runs from the project root regardless of where it's called
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"

echo "=========================================="
echo "      Metadoon (Docker for Linux)"
echo "=========================================="

# 2. Check Docker
if ! command -v docker &> /dev/null; then
    echo "[ERROR] Docker not found. Please install Docker Engine."
    exit 1
fi

# 3. Configure X11 Permissions (Required for GUI)
echo "[INFO] Configuring X11 permissions..."
# Allows the container to communicate with the local X11 server
xhost +local:docker > /dev/null 2>&1

# 4. Update Image
echo "[INFO] Checking for updates..."
docker pull $IMAGE_NAME

# 5. Execute Container
echo ""
echo "[INFO] Launching Metadoon..."
echo "[INFO] Current directory mapped to: /workspace"
echo "[INFO] Your Home folder mapped to: /app/YOUR_DATA"

# --- DOCKER RUN COMMAND ---
# --net=host: Shares network stack (often helps with X11 connectivity)
# --user $(id -u):$(id -g): Runs container as CURRENT USER to avoid 'root' file permission issues on outputs
# -e DISPLAY=$DISPLAY: Connects to the host display
# -v /tmp/.X11-unix...: Maps X11 socket
# -v "$(pwd)":/workspace: Maps current project folder for inputs/outputs
# -v "$HOME":/app/YOUR_DATA: Maps Linux Home folder for easy file selection
# bash -c ...: Enters /app and runs the tool

docker run --rm -it \
    --net=host \
    --user $(id -u):$(id -g) \
    -e DISPLAY=$DISPLAY \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    -v "$(pwd)":/workspace \
    -v "$HOME":/app/YOUR_DATA \
    --workdir /workspace \
    $IMAGE_NAME \
    bash -c "cd /app && python metadoon.py"

if [ $? -ne 0 ]; then
    echo ""
    echo "[ERROR] Execution failed."
    exit 1
fi