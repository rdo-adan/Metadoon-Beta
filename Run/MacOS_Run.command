#!/bin/bash

# Image Configuration
IMAGE_NAME="engbio/metadoon:v1.0"

# Navigate to project root
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR/.."

echo "=========================================="
echo "      Metadoon Launcher (macOS)           "
echo "=========================================="

# 1. Check Docker
if ! command -v docker &> /dev/null; then
    echo "[ERROR] Docker not found. Please install Docker Desktop."
    exit 1
fi

# 2. Pull/Update Image
echo "[INFO] Checking for updates ($IMAGE_NAME)..."
docker pull $IMAGE_NAME

if [[ "$(docker images -q $IMAGE_NAME 2> /dev/null)" == "" ]]; then
    echo "[ERROR] Image not found locally and download failed."
    exit 1
fi

# 3. Configure Display (XQuartz)
# Allow local network connections for XQuartz
echo "[INFO] Configuring XQuartz Display..."

# Optional: Verify if XQuartz is running
open -a XQuartz 2>/dev/null
xhost + 127.0.0.1 > /dev/null 2>&1

# 4. Run Container
# Mac uses 'host.docker.internal:0' to route video to XQuartz
echo "[INFO] Launching..."
echo "[NOTE] Ensure XQuartz is running and 'Allow connections from network clients' is enabled in settings."

docker run --rm -it \
    -e DISPLAY=host.docker.internal:0 \
    -v "$(pwd)":/app \
    -v "$HOME":/host_home:ro \
    -w /app \
    $IMAGE_NAME