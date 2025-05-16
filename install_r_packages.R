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
