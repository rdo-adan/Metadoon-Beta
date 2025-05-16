#!/bin/bash
# Metadoon Linux Setup Script

# Atualizar pacotes
sudo apt-get update -y && sudo apt-get upgrade -y

# Dependências básicas
sudo apt-get install -y wget build-essential libssl-dev libcurl4-openssl-dev libxml2-dev python3 python3-pip python3-tk r-base r-base-dev git

# Instalar pacotes R necessários
cat <<EOF > install_r_packages.R
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
library(devtools)

# Tentar instalar pairwiseAdonis com fallback
tryCatch({
  devtools::install_github("vlubitch/pairwiseAdonis")
}, error = function(e) {
  message("Warning: Falha ao instalar 'pairwiseAdonis'. Verifique manualmente caso necessário.")
})

# Instalar outros pacotes
install.packages(c("tidyverse", "reshape2", "igraph", "foreach", "lme4"))
EOF

Rscript install_r_packages.R

# Ajustar permissões e configurar o diretório
chmod +x Setup.sh
chmod +x Analise.R
echo "Setup concluído! Execute a interface Metadoon com Python usando: python3 teste.py"

