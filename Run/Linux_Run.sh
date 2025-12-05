#!/bin/bash

# Define o nome da imagem
IMAGE_NAME="engbio/metadoon:v1.0"

# Garante que o script rode a partir da raiz do projeto
cd "$(dirname "$0")/.."

echo "=========================================="
echo "      Metadoon Launcher (Linux)           "
echo "=========================================="

# 1. Verifica se o Docker está instalado
if ! command -v docker &> /dev/null; then
    echo "[ERROR] Docker not found. Please install Docker first."
    read -p "Press Enter to exit..."
    exit 1
fi

# 2. Baixa/Atualiza a imagem
echo "[INFO] Pulling latest version ($IMAGE_NAME)..."
docker pull $IMAGE_NAME

# Verifica se a imagem existe (caso o pull falhe por falta de internet)
if [[ "$(docker images -q $IMAGE_NAME 2> /dev/null)" == "" ]]; then
    echo "[ERROR] Image not found locally and download failed."
    read -p "Press Enter to exit..."
    exit 1
fi

# 3. Configura permissão para o X11 (Interface Gráfica)
echo "[INFO] Configuring Display..."
xhost +local:docker > /dev/null 2>&1

# 4. Executa o Container
echo "[INFO] Launching..."

docker run --rm -it \
  --user $(id -u):$(id -g) \
  -e DISPLAY=$DISPLAY \
  -v /tmp/.X11-unix:/tmp/.X11-unix \
  -v "$(pwd)":/app \
  -v "$HOME":/host_home:ro \
  -w /app \
  $IMAGE_NAME

if [ $? -ne 0 ]; then
    echo "[ERROR] Execution failed."
    read -p "Press Enter to exit..."
fi