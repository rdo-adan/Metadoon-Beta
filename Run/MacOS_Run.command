#!/bin/bash
IMAGE_NAME="engbio/metadoon:v1.0"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR/../"

echo "=== Metadoon Launcher (macOS) ==="

# 1. Verifica Docker
if ! command -v docker &> /dev/null; then
    echo "[ERROR] Docker not found. Please install Docker Desktop."
    exit 1
fi

# 2. Baixa Imagem
echo "[INFO] Updating image..."
docker pull $IMAGE_NAME

# 3. Configura Display (XQuartz)
xhost + 127.0.0.1 > /dev/null 2>&1

echo "[INFO] Launching..."
docker run --rm -it \
    -v "$(pwd)":/app \
    -e DISPLAY=host.docker.internal:0 \
    $IMAGE_NAME