################################### Carregamento de Pacotes ########################################

# Carrega todos os pacotes necessários
required_packages <- c("phyloseq", "vegan", "ggplot2", "ggpubr", "cowplot", "dplyr", 
                       "DESeq2", "scater", "rprojroot", "pairwiseAdonis", "pheatmap", "viridis", "ape", "microbiome", "rprojroot")

sapply(required_packages, require, character.only = TRUE) # Usando sapply para melhor tratamento de erros

########################################### Loading Files ###################################
script_dir <- find_root(has_file("Analise.R"))
print(script_dir)
########################################### Loading Files ###################################

# Full paths to files based on the script directory
otu_file <- file.path(script_dir, "OTUs/otutab.txt")
taxonomy_file <- file.path(script_dir, "Taxonomy/taxonomy.txt")
metadata_dir <- file.path(script_dir, "Metadata File")
tree_file_path <- file.path(script_dir, "Tree File/tree.nwk")

# Read the complete file as text
lines <- readLines(otu_file)
lines[1] <- gsub("^#", "", lines[1])
temp_file <- tempfile()
writeLines(lines, temp_file)
tabBIN <- as.matrix(read.table(temp_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))

# Load taxonomy table
tabTax <- as.matrix(read.csv(taxonomy_file, header = TRUE, row.names = 1, sep = "\t"))

# Load metadata file
metadata_files <- list.files(path = metadata_dir, pattern = "\\.tsv$|\\.csv$", full.names = TRUE)
if (length(metadata_files) > 0) {
  sep <- ifelse(grepl("\\.tsv$", metadata_files[1]), "\t", ",")
  metadataTab <- read.csv(metadata_files[1], header = TRUE, row.names = 1, sep = sep)
  message("Metadata file loaded: ", metadata_files[1])
} else {
  stop("No metadata file found in the specified folder.")
}

# Check if the tree file is present
if (file.exists(tree_file_path)) {
  Arv <- read.tree(tree_file_path)
  message("Tree file successfully loaded.")
} else {
  message("Tree file not found. Proceeding without it.")
  Arv <- NULL
}

########################################### Create `phyloseq` Object ###################################

if (is.null(Arv)) {
  ps <- phyloseq(
    otu_table(tabBIN, taxa_are_rows = TRUE),
    tax_table(tabTax),
    sample_data(metadataTab)
  )
} else {
  ps <- phyloseq(
    otu_table(tabBIN, taxa_are_rows = TRUE),
    tax_table(tabTax),
    sample_data(metadataTab),
    phy_tree(Arv)
  )
}

# Check the phyloseq object
print(ps)

########################################### Functions and Plots ###################################

# Output directory for saving plots
output_dir <- file.path(script_dir, "Output")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Function to save plots with a white background
save_plot <- function(plot, filename, width = 12, height = 8) {
  ggsave(filename = file.path(output_dir, filename), plot = plot, width = width, height = height, dpi = 300, bg = "white")
}


#Rarefaction
# Extract OTU table as matrix
otu <- as(otu_table(ps), "matrix")

# Ensure taxa are in rows
if (!taxa_are_rows(ps)) {
  otu <- t(otu)
}

# Plot rarefaction curve
vegan::rarecurve(otu, step = 100, cex = 0.6, label = FALSE)
title("Rarefaction Curve")




# Função para plotar abundância relativa
plot_abundance <- function(ps, taxrank, top_n = 15) {
  ps_tax <- tax_glom(ps, taxrank = taxrank)
  ps_rel <- transform_sample_counts(ps_tax, function(x) x / sum(x))
  taxa_sums <- taxa_sums(ps_rel)
  top_taxa <- names(sort(taxa_sums, decreasing = TRUE)[1:top_n])
  tax_table(ps_rel)[!(taxa_names(ps_rel) %in% top_taxa), taxrank] <- "Others"
  p <- ggplot(psmelt(ps_rel), aes(x = Sample, y = Abundance, fill = !!sym(taxrank))) +
    geom_bar(stat = "identity") +
    labs(title = paste("Relative Abundance by", taxrank), x = "Sample", y = "Relative Abundance") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  save_plot(p, paste0("abundance_", taxrank, ".png"), width = 14, height = 8)
  print(p)
}

# Generate abundance plots for Phylum, Class, Order, Family, and Genus
for (rank in c("Phylum", "Class", "Order", "Family", "Genus")) {
  plot_abundance(ps, rank, top_n = 15)
}

# Normalizar dados para torná-los composicionais
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Instalar e carregar o pacote microbiome e viridis
if (!requireNamespace("microbiome", quietly = TRUE)) {
  install.packages("devtools")
  devtools::install_github("microbiome/microbiome")
}
if (!requireNamespace("viridis", quietly = TRUE)) {
  install.packages("viridis")
}
# Define parameters for core microbiome analysis
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-4), log10(.2), length = 10)

# Função para criar heatmap do core microbioma para diferentes níveis taxonômicos
plot_core_heatmap <- function(ps, taxrank, top_n = 30) {
  # Aglomerar OTUs ao nível taxonômico desejado
  ps_tax <- tax_glom(ps, taxrank)
  
  # Normalizar dados para torná-los composicionais
  ps_tax_rel <- transform_sample_counts(ps_tax, function(x) x / sum(x))
  
  # Filtrar os táxons mais prevalentes no nível especificado
  prev_table <- microbiome::prevalence(ps_tax_rel) # Calcular prevalência
  top_taxa <- names(sort(prev_table, decreasing = TRUE)[1:top_n])
  ps_top_taxa <- prune_taxa(top_taxa, ps_tax_rel)  # Manter apenas os táxons selecionados
  
  # Parâmetros de prevalência e detecção
  prevalences <- seq(.05, 1, .05)
  detections <- 10^seq(log10(1e-4), log10(.2), length = 10)
  
  # Plotar o core microbioma
  p <- plot_core(ps_top_taxa, 
                 plot.type = "heatmap", 
                 prevalences = prevalences, 
                 detections = detections, 
                 min.prevalence = 0.5) +
    xlab("Detection Threshold (Relative Abundance (%))") +
    ylab(taxrank) +
    scale_fill_viridis() +
    theme_minimal() +
    ggtitle(paste("Core Microbiome at", taxrank, "Level"))
  
  # Salvar o gráfico com alta resolução e exibir
  save_plot(p, paste0("core_microbiome_heatmap_", taxrank, ".png"), width = 16, height = 10)
  print(p)
}

# Gerar heatmaps para os níveis Phylum, Class, Order, Family e Genus
for (rank in c("Phylum", "Class", "Order", "Family", "Genus")) {
  plot_core_heatmap(ps, rank, top_n = 30)
}


# Definir parâmetros para análise de core microbioma
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-4), log10(.2), length = 10)

# Função para criar heatmap do core microbioma para diferentes níveis taxonômicos
plot_core_heatmap <- function(ps, taxrank, top_n = 30) {
  # Aglomerar OTUs ao nível taxonômico desejado
  ps_tax <- tax_glom(ps, taxrank)
  
  # Normalizar dados para torná-los composicionais
  ps_tax_rel <- transform_sample_counts(ps_tax, function(x) x / sum(x))
  
  # Obter a prevalência dos táxons e selecionar os mais prevalentes
  prev_table <- microbiome::prevalence(ps_tax_rel)
  top_taxa <- names(sort(prev_table, decreasing = TRUE)[1:top_n])
  ps_top_taxa <- prune_taxa(top_taxa, ps_tax_rel)
  
  # Substituir identificadores OTU pelo nome do nível taxonômico desejado
  tax_table(ps_top_taxa)[, taxrank] <- as.character(tax_table(ps_top_taxa)[, taxrank])
  taxa_names(ps_top_taxa) <- tax_table(ps_top_taxa)[, taxrank]
  
  # Plotar o core microbioma
  p <- plot_core(ps_top_taxa, 
                 plot.type = "heatmap", 
                 prevalences = prevalences, 
                 detections = detections, 
                 min.prevalence = 0.5) +
    xlab("Detection Threshold (Relative Abundance (%))") +
    ylab(taxrank) +
    scale_fill_viridis() +
    theme_minimal() +
    ggtitle(paste("Core Microbiome at", taxrank, "Level"))
  
  # Salvar o gráfico com alta resolução e exibir
  save_plot(p, paste0("taxa_core_microbiome_heatmap_", taxrank, ".png"), width = 16, height = 10)
  print(p)
}

# Gerar heatmaps para os níveis Phylum, Class, Order, Family e Genus
for (rank in c("Phylum", "Class", "Order", "Family", "Genus")) {
  plot_core_heatmap(ps, rank, top_n = 30)
}

alpha_div <- estimate_richness(ps)
head(alpha_div)
alpha_merged <- merge(alpha_div, metadata, by = "SampleID")
head(alpha_merged)

any(is.na(alpha_div$SampleID))
any(is.na(meta$SampleID))


# Alpha Diversity Plot — Um gráfico por índice, contendo todas as amostras

plot_alpha_diversity <- function(ps) {
  
  # Calcula diversidade alfa por amostra
  alpha_div <- estimate_richness(ps, measures = c("Shannon", "Simpson", "Chao1", "ACE", "Observed"))
  alpha_div$SampleID <- rownames(alpha_div)
  
  # Verifica se SampleID está correto
  if (any(is.na(alpha_div$SampleID))) {
    stop("SampleID contains NA!")
  }
  
  indices <- c("Shannon", "Simpson", "Chao1", "ACE", "Observed")
  
  for (index in indices) {
    
    plot_df <- data.frame(
      Sample = alpha_div$SampleID,
      Diversity = alpha_div[[index]]
    )
    
    p <- ggplot(plot_df, aes(x = Sample, y = Diversity)) +
      geom_point(size = 3, color = "steelblue") +
      geom_line(group = 1, color = "gray")+
      labs(title = paste("Alpha Diversity -", index),
           x = "SampleID",
           y = index) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    ggsave(filename = paste0("Output/alpha_diversity_", index, ".png"),
           plot = p, width = 10, height = 6, dpi = 600, bg = "white")
    
    print(p)
  }
}
plot_alpha_diversity(ps)



# Function to plot NMDS with Bray-Curtis and sample labels
plot_beta_diversity <- function(ps, metadata_column) {
  ps_dist <- phyloseq::distance(ps, method = "bray")
  ord_nmds <- metaMDS(ps_dist, k = 2, trymax = 20)
  
  ord_data <- data.frame(ord_nmds$points)
  ord_data$SampleID <- rownames(ord_data)
  ord_data$Group <- sample_data(ps)[[metadata_column]]
  
  p <- ggplot(ord_data, aes(x = MDS1, y = MDS2, color = Group, label = SampleID)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.5, hjust = 0.5) +
    stat_ellipse() +
    labs(title = paste("Beta Diversity NMDS - Grouped by:", metadata_column),
         x = "NMDS1", y = "NMDS2") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right")
  
  save_plot(p, paste0("beta_diversity_NMDS_", metadata_column, ".png"), width = 10, height = 8)
  print(p)
}

# Loop to generate NMDS plots for each metadata column
for (metadata_column in colnames(sample_data(ps))) {
  plot_beta_diversity(ps, metadata_column)
}

# PERMANOVA Test
run_permanova <- function(ps, metadata_column) {
  metadata <- as(sample_data(ps), "data.frame")
  dist_bc <- phyloseq::distance(ps, method = "bray")
  if (metadata_column %in% colnames(metadata)) {
    formula <- as.formula(paste("dist_bc ~", metadata_column))
    permanova_result <- adonis2(formula, data = metadata, permutations = 999)
    result_file <- file.path(output_dir, paste0("permanova_result_", metadata_column, ".txt"))
    write.table(permanova_result, file = result_file, sep = "\t", col.names = NA, quote = FALSE)
    print(paste("PERMANOVA results for column:", metadata_column))
    print(permanova_result)
  } else {
    print(paste("Column", metadata_column, "not found in metadata"))
  }
}

for (metadata_column in colnames(sample_data(ps))) {
  run_permanova(ps, metadata_column)
}


########################################### Clean Environment ###################################

# Limpa o ambiente para liberar memória após a execução (opcional, mas recomendado)
rm(list = ls())
gc() # Coleta o lixo para liberar memória imediatamente
message("Script executado com sucesso. Todos os resultados foram salvos na pasta 'Output'.")
