#!/bin/bash

# Nome do ambiente Conda e arquivo YAML
ENV_NAME="metadoon"
YAML_FILE="metadoon.yaml"

# Verifique se o arquivo YAML existe
if [[ ! -f "$YAML_FILE" ]]; then
  echo "Arquivo $YAML_FILE não encontrado. Por favor, verifique o caminho do arquivo YAML."
  exit 1
fi

# Instalar pacotes adicionais do R, incluindo dependências do devtools
echo "Instalando pacotes adicionais do R (Bioconductor e dependências do devtools)..."
Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); BiocManager::install(c("scater", "phyloseq", "DESeq2")); install.packages(c("shiny", "miniUI", "curl", "httpuv", "devtools"), repos="https://cloud.r-project.org")'

# Instale as dependências do sistema necessárias para pacotes R e V8
echo "Instalando dependências do sistema para pacotes R..."
sudo apt-get update
sudo apt-get install -y libcurl4-openssl-dev libxml2-dev libssl-dev \
                        libharfbuzz-dev libfribidi-dev libfontconfig1-dev \
                        libgsl-dev liblapack-dev libmpfr-dev libv8-dev \
                        libcairo2-dev libuv-dev

# Criar o ambiente Conda usando o arquivo YAML
echo "Criando o ambiente Conda $ENV_NAME usando o arquivo YAML..."
conda env create -f "$YAML_FILE"

# Verificar se o ambiente foi criado com sucesso
if [[ $? -ne 0 ]]; then
  echo "Falha ao criar o ambiente Conda $ENV_NAME."
  exit 1
fi

# Ativar o ambiente Conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate "$ENV_NAME"

# Instalar pacotes adicionais do R, especialmente de Bioconductor
echo "Instalando pacotes adicionais do R (Bioconductor)..."
Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install(c("scater", "phyloseq", "DESeq2"))'

sudo apt-get update
sudo apt-get install libgmp-dev libmpfr-dev

#devtools
apt install r-cran-devtools
# Confirmação de finalização
echo "Ambiente $ENV_NAME configurado com sucesso. Para começar a usar, execute 'conda activate $ENV_NAME'."

