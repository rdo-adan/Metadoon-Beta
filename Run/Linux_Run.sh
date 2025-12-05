#!/bin/bash

# Image Configuration
IMAGE_NAME="engbio/metadoon:v1.0"

# Ensure the script runs from the project root
cd "$(dirname "$0")/.."

echo "=========================================="
echo "      Metadoon Launcher (Linux)           "
echo "=========================================="

# 1. Check if Docker is installed
if ! command -v docker &> /dev/null; then
    echo "[ERROR] Docker not found. Please install Docker first."
    echo "Try: sudo apt install docker.io"
    read -p "Press Enter to exit..."
    exit 1
fi

# 2. Pull/Update the image
echo "[INFO] Checking for updates ($IMAGE_NAME)..."
docker pull $IMAGE_NAME

# Check if image exists (fallback if pull fails due to offline status)
if [[ "$(docker images -q $IMAGE_NAME 2> /dev/null)" == "" ]]; then
    echo "[ERROR] Image not found locally and download failed."
    read -p "Press Enter to exit..."
    exit 1
fi

# 3. Configure X11 Display permissions
echo "[INFO] Configuring Display..."
xhost +local:docker > /dev/null 2>&1

# 4. Run Container
# --user: Runs with current user ID (prevents root-locked output files)
# -v /tmp/.X11-unix: Maps native X11 socket for GUI
# -w /app: Sets working directory inside container
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