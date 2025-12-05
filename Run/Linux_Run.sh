#!/bin/bash

IMAGE_NAME="engbio/metadoon:v1.0"
cd "$(dirname "$0")/.."

echo "=========================================="
echo "      Metadoon Launcher (Linux)           "
echo "=========================================="

# 1. Verifica Docker
if ! command -v docker &> /dev/null; then
    echo "[ERROR] Docker not found."
    exit 1
fi

# 2. Atualiza
echo "[INFO] Checking for updates..."
docker pull $IMAGE_NAME

if [[ "$(docker images -q $IMAGE_NAME 2> /dev/null)" == "" ]]; then
    echo "[ERROR] Image download failed."
    exit 1
fi

# 3. Configura X11
xhost +local:docker > /dev/null 2>&1

# 4. Executa
echo "[INFO] Launching..."
echo "[TIP] Your Home files are inside the folder 'YOUR_DATA'"

# --- MUDANÇA AQUI: ---
# -v "$HOME":/app/YOUR_DATA mapeia sua home para dentro da pasta de trabalho

docker run --rm -it \
  --user $(id -u):$(id -g) \
  -e DISPLAY=$DISPLAY \
  -v /tmp/.X11-unix:/tmp/.X11-unix \
  -v "$(pwd)":/app \
  -v "$HOME":/app/YOUR_DATA \
  -w /app \
  $IMAGE_NAME

if [ $? -ne 0 ]; then
    echo "[ERROR] Execution failed."
fi