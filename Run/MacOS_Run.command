#!/bin/bash
IMAGE_NAME="engbio/metadoon:v1.0"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR/.."

echo "=== Metadoon Launcher (macOS) ==="

if ! command -v docker &> /dev/null; then
    echo "[ERROR] Docker not found."
    exit 1
fi

echo "[INFO] Updating..."
docker pull $IMAGE_NAME

# Configura XQuartz
xhost + 127.0.0.1 > /dev/null 2>&1

echo "[INFO] Launching..."
echo "[TIP] Your Home files are inside the folder 'YOUR_DATA'"

docker run --rm -it \
    -e DISPLAY=host.docker.internal:0 \
    -v "$(pwd)":/app \
    -v "$HOME":/app/YOUR_DATA \
    -w /app \
    $IMAGE_NAME