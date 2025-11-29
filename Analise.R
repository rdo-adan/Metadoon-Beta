################################### Loading Packages ########################################

required_packages <- c("phyloseq", "vegan", "ggplot2", "ggpubr", "cowplot", "dplyr", 
                       "DESeq2", "scater", "rprojroot", "pairwiseAdonis", "pheatmap", "viridis", "ape", "microbiome", "rprojroot","wesanderson","RColorBrewer")

sapply(required_packages, require, character.only = TRUE) 


########################################### Loading Files ###################################
script_dir <- find_root(has_file("Analise.R"))
print(script_dir)
library(jsonlite)
`%||%` <- function(a, b) if (!is.null(a)) a else b 

params_file <- "pipeline_params.json"

if (file.exists(params_file)) {
  params <- fromJSON(params_file)
} else {
  stop("pipeline_params.json not found! Please configure parameters first.")
}

# Assign variables with fallbacks to defaults

stat_test <- params$stat_test %||% "anova"
dist_method <- params$dist_method %||% "bray"
color_palette <- params$color_palette %||% "viridis"
abundance_top_n <- params$abundance_top_n %||% 15
core_top_n <- params$core_top_n %||% 30
rarefaction_step <- params$rarefaction_step %||% 100
rarefaction_cex <- params$rarefaction_cex %||% 0.6
enable_rarefaction <- params$enable_rarefaction %||% FALSE
rarefaction_depth <- params$rarefaction_depth %||% 1000
get_palette <- function(palette_name) {
  switch(palette_name,
         viridis = viridis::viridis(10),
         plasma = viridis::plasma(10),
         inferno = viridis::inferno(10),
         magma = viridis::magma(10),
         cividis = viridis::cividis(10),
         wesanderson = wesanderson::wes_palette("Darjeeling1", 10),
         RColorBrewer::brewer.pal(10, "Set3")  # default fallback
  )
}
palette_colors <- get_palette(color_palette)

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
  ggsave(filename = file.path(output_dir, filename), plot = plot, width = width, height = height, dpi = 600, bg = "white")
}


#Rarefaction
# Extract OTU table as matrix
otu <- as(otu_table(ps), "matrix")

# Ensure taxa are in rows
if (!taxa_are_rows(ps)) {
  otu <- t(otu)
}

# Plot rarefaction curve
vegan::rarecurve(otu, step = rarefaction_step, cex = rarefaction_cex, label = rarefaction_label)
title("Rarefaction Curve")

#rarefaction for new data
rarefy_phyloseq <- function(ps, depth = rarefaction_depth) {
  if (enable_rarefaction) {
    message(paste("Applying rarefaction to depth:", depth))
    # Verifica se todos os samples têm reads suficientes
    min_reads <- min(sample_sums(ps))
    if (depth > min_reads) {
      warning(paste("Rarefaction depth", depth, "exceeds minimum sample reads (", min_reads, "). Adjusting..."))
      depth <- min_reads
    }
    ps_rarefied <- rarefy_even_depth(ps, sample.size = depth, replace = TRUE)
    message("Data rarefied successfully!")
    return(ps_rarefied)
  } else {
    message("Rarefaction disabled - using original data")
    return(ps)
  }
}

#check rarefied phyloseq object
ps <- rarefy_phyloseq(ps)
message("=== Rarefaction Settings ===")
message(paste("Enable rarefaction:", enable_rarefaction))
message(paste("Rarefaction depth:", rarefaction_depth))
message(paste("Samples before rarefaction:", nsamples(ps)))
# apply conditional rarefaction
ps <- rarefy_phyloseq(ps, rarefaction_depth)
message(paste("Samples after rarefaction:", nsamples(ps)))
message("============================\n")


# RElative abundance
plot_abundance <- function(ps, taxrank, top_n = abundance_top_n) {
  ps_tax <- tax_glom(ps, taxrank = taxrank)
  ps_rel <- transform_sample_counts(ps_tax, function(x) x / sum(x))
  taxa_sums <- taxa_sums(ps_rel)
  top_taxa <- names(sort(taxa_sums, decreasing = TRUE)[1:top_n])
  tax_table(ps_rel)[!(taxa_names(ps_rel) %in% top_taxa), taxrank] <- "Others"
  p <- ggplot(psmelt(ps_rel), aes(x = Sample, y = Abundance, fill = !!sym(taxrank))) +
    geom_bar(stat = "identity") +
    labs(title = paste("Relative Abundance by", taxrank), x = "Sample", y = "Relative Abundance") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    palette_colors <- get_palette(color_palette) +  # Define this function (see next step)
    p <- p + scale_fill_manual(values = palette_colors)
  save_plot(p, paste0("abundance_", taxrank, ".png"), width = 16, height = 10)
  print(p)
}

# Generate abundance plots for Phylum, Class, Order, Family, and Genus
for (rank in c("Phylum", "Class", "Order", "Family", "Genus")) {
  plot_abundance(ps, rank)
}

# Normalization
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))

# install_packages microbiome and viridis
if (!requireNamespace("microbiome", quietly = TRUE)) {
  install.packages("devtools")
  devtools::install_github("microbiome/microbiome")
}
if (!requireNamespace("viridis", quietly = TRUE)) {
  install.packages("viridis")
}
# Define parameters for core microbiome analysis
prevalences <- seq(0.05, 1, 0.05)
detections <- 10^seq(log10(1e-4), log10(.2), length = 10)

# Core_microbiome
plot_core_heatmap <- function(ps, taxrank, top_n = core_top_n) {
  ps_tax <- tax_glom(ps, taxrank)
  
  # Normalize data to make them compositional
  ps_tax_rel <- transform_sample_counts(ps_tax, function(x) x / sum(x))
  
  # Filtrar os táxons mais prevalentes no nível especificado
  prev_table <- microbiome::prevalence(ps_tax_rel) # Calcular prevalência
  top_taxa <- names(sort(prev_table, decreasing = TRUE)[1:top_n])
  ps_top_taxa <- prune_taxa(top_taxa, ps_tax_rel)  # Manter apenas os táxons selecionados
  
  # Detection parameters
  prevalences <- seq(0.05, 1, 0.05)
  detections <- 10^seq(log10(1e-4), log10(.2), length = 10)
  
  # Plot core
  p <- plot_core(ps_top_taxa, 
                 plot.type = "heatmap", 
                 prevalences = prevalences, 
                 detections = detections, 
                 min.prevalence = 0.5) +
    xlab("Detection Threshold (Relative Abundance (%))") +
    ylab(taxrank) +
    #scale_fill_viridis() +
    p <- p + scale_fill_manual(values = palette_colors)
    theme_minimal() +
    ggtitle(paste("Core Microbiome at", taxrank, "Level"))
  
  save_plot(p, paste0("core_microbiome_heatmap_", taxrank, ".png"), width = 16, height = 10)
  print(p)
}

# Gerar heatmaps para os níveis Phylum, Class, Order, Family e Genus
for (rank in c("Phylum", "Class", "Order", "Family", "Genus")) {
  plot_core_heatmap(ps, rank)
}


# Define parameters
prevalences <- seq(0.05, 1, 0.05)
detections <- 10^seq(log10(1e-4), log10(.2), length = 10)

# Function for core microbiome with taxa
plot_core_heatmap <- function(ps, taxrank, top_n = 30) {
  #aglomerate
  ps_tax <- tax_glom(ps, taxrank)
  
  # Normalize
  ps_tax_rel <- transform_sample_counts(ps_tax, function(x) x / sum(x))
  
  prev_table <- microbiome::prevalence(ps_tax_rel)
  top_taxa <- names(sort(prev_table, decreasing = TRUE)[1:top_n])
  ps_top_taxa <- prune_taxa(top_taxa, ps_tax_rel)
  
  tax_table(ps_top_taxa)[, taxrank] <- as.character(tax_table(ps_top_taxa)[, taxrank])
  taxa_names(ps_top_taxa) <- tax_table(ps_top_taxa)[, taxrank]
  
  # Plot core w tax
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
  
  save_plot(p, paste0("taxa_core_microbiome_heatmap_", taxrank, ".png"), width = 16, height = 10)
  print(p)
}

for (rank in c("Phylum", "Class", "Order", "Family", "Genus")) {
  plot_core_heatmap(ps, rank, top_n = 30)
}



# Alpha Diversity Plot grouped by any metadata column

plot_alpha_diversity <- function(ps) {
  
  if (!dir.exists("Output")) dir.create("Output")
  
  alpha_div <- estimate_richness(ps, measures = c("Shannon", "Simpson", "Chao1", "ACE", "Observed"))
  alpha_div$SampleID <- rownames(alpha_div)
  
  meta_file <- list.files(path = "Metadata File", pattern = "\\.tsv$|\\.csv$", full.names = TRUE)
  
  if (length(meta_file) == 0) stop("No metadata .tsv or .csv file found.")
  if (length(meta_file) > 1) warning("More than one metadata file found. Using the first: ", meta_file[1])
  
  meta_raw <- read.csv(meta_file[1], header = TRUE, sep = "\t", check.names = FALSE)
  if ("SampleID" %in% colnames(meta_raw)) {
    rownames(meta_raw) <- meta_raw$SampleID
    meta_raw$SampleID <- NULL
  }
  meta_data <- meta_raw
  meta_data$SampleID <- rownames(meta_data)
  
  alpha_merged <- merge(alpha_div, meta_data, by = "SampleID")
  
  indices <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson")
  meta_columns <- setdiff(colnames(meta_data), "SampleID")
  if(length(meta_columns) == 0) stop("No metadata columns found besides SampleID.")
  
  for(meta_var in meta_columns){
    message(paste0("Generating alpha diversity plots grouped by '", meta_var, "'"))
    for(index in indices){
      
      plot_df <- data.frame(
        Sample = alpha_merged$SampleID,
        Diversity = alpha_merged[[index]],
        Group = as.factor(alpha_merged[[meta_var]])
      )
      
      # Perform chosen statistical test - example: ANOVA here; adapt as desired
      stat_test_result <- compare_means(Diversity ~ Group, data = plot_df, method = tolower(stat_test))
      
      p <- ggboxplot(plot_df, x = "Group", y = "Diversity", color = "Group", add = "jitter") +
        labs(title = paste("Alpha Diversity -", index, "\nGrouped by:", meta_var),
             x = meta_var, y = index) +
        theme_minimal(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        stat_pvalue_manual(stat_test_result, label = "p.format") 
      
      ggsave(filename = paste0("Output/alpha_diversity_", index, "_", meta_var, ".png"),
             plot = p, width = 8, height = 6, dpi = 600, bg = "white")
      
      print(p)
    }
  }
}

# Run the function
plot_alpha_diversity(ps)


plot_beta_diversity <- function(ps, dist_method = "bray") {
  
  if (!dir.exists("Output")) dir.create("Output")
  
  meta_data <- as(sample_data(ps), "data.frame")
  meta_data$SampleID <- rownames(meta_data)
  
  meta_columns <- setdiff(colnames(meta_data), "SampleID")
  
  dist <- phyloseq::distance(ps, method = dist_method)
  
  for (meta_var in meta_columns) {
    
    message(paste0("Generating beta diversity NMDS plots grouped by ", meta_var, " using method ", dist_method))
    
    ord_nmds <- metaMDS(dist, k = 2, trymax = 50)
    nmds_df <- data.frame(ord_nmds$points)
    nmds_df$SampleID <- rownames(nmds_df)
    nmds_df$Group <- meta_data[match(nmds_df$SampleID, meta_data$SampleID), meta_var]
    
    # Calculate PERMANOVA for p-value
    formula <- as.formula(paste("dist ~", meta_var))
    permanova_result <- adonis2(dist ~ Group, data = nmds_df, permutations = 999)
    p_value <- permanova_result$`Pr(>F)`[1]
    
    p <- ggplot(nmds_df, aes(x = MDS1, y = MDS2, color = Group, label = SampleID)) +
      geom_point(size = 3) +
      geom_text(vjust = -0.5, hjust = 0.5) +
      stat_ellipse() +
      labs(title = paste("Beta Diversity - NMDS grouped by", meta_var),
           subtitle = paste("PERMANOVA p =", signif(p_value, digits = 3)),
           x = "NMDS1", y = "NMDS2") +
      theme_minimal(base_size = 14)
    
    ggsave(filename = paste0("Output/beta_diversity_NMDS_", meta_var, ".png"),
           plot = p, width = 10, height = 8, dpi = 600, bg = "white")
    
    print(p)
  }
}

plot_beta_diversity(ps)


# Beta Diversity Plot using PCoA — grouped by any metadata column

plot_pcoa_diversity <- function(ps, dist_method = "bray") {
  
  if (!dir.exists("Output")) dir.create("Output")
  
  meta_data <- as(sample_data(ps), "data.frame")
  meta_data$SampleID <- rownames(meta_data)
  
  meta_columns <- setdiff(colnames(meta_data), "SampleID")
  
  dist <- phyloseq::distance(ps, method = dist_method)
  ord_pcoa <- ape::pcoa(as.matrix(dist))
  pcoa_df <- data.frame(ord_pcoa$vectors[, 1:2])
  colnames(pcoa_df) <- c("PCoA1", "PCoA2")
  pcoa_df$SampleID <- rownames(pcoa_df)
  
  for (meta_var in meta_columns) {
    
    message(paste0("Generating PCoA plots grouped by '", meta_var, "' using method ", dist_method))
    
    pcoa_df$Group <- meta_data[match(pcoa_df$SampleID, meta_data$SampleID), meta_var]
    
    # Calculate PERMANOVA for p-value
    permanova_result <- adonis2(dist ~ Group, data = meta_data, permutations = 999)
    p_value <- permanova_result$`Pr(>F)`[1]
    
    p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Group, label = SampleID)) +
      geom_point(size = 3) +
      geom_text(vjust = -0.5, hjust = 0.5) +
      stat_ellipse() +
      labs(title = paste("Beta Diversity - PCoA grouped by", meta_var),
           subtitle = paste("PERMANOVA p =", signif(p_value, digits = 3)),
           x = "PCoA1", y = "PCoA2") +
      theme_minimal(base_size = 14)
    
    ggsave(filename = paste0("Output/beta_diversity_PCoA_", meta_var, ".png"),
           plot = p, width = 10, height = 8, dpi = 600, bg = "white")
    
    print(p)
  }
}

plot_pcoa_diversity(ps)

# PERMANOVA Test
run_permanova <- function(ps) {
  
  if (!dir.exists("Output")) dir.create("Output")
  
  meta_data <- as(sample_data(ps), "data.frame")
  
  dist_bc <- phyloseq::distance(ps, method = "bray")
  
  meta_columns <- colnames(meta_data)
  
  for (meta_var in meta_columns) {
    
    formula <- as.formula(paste("dist_bc ~", meta_var))
    permanova_result <- adonis2(formula, data = meta_data, permutations = 999)
    
    result_file <- paste0("Output/permanova_result_", meta_var, ".txt")
    
    write.table(permanova_result, file = result_file, sep = "\t", col.names = NA, quote = FALSE)
    
    message(paste0("✔️ PERMANOVA result saved for ", meta_var))
    print(permanova_result)
  }
}

run_permanova(ps)

run_deseq2 <- function(ps) {
  
  if (!dir.exists("Output")) dir.create("Output")
  
  meta_data <- as(sample_data(ps), "data.frame")
  meta_data$SampleID <- rownames(meta_data)
  
  meta_columns <- setdiff(colnames(meta_data), "SampleID")
  
  taxa_table <- as(tax_table(ps), "matrix")
  
  get_tax_label <- function(tax_mat) {
    apply(tax_mat, 1, function(x) {
      taxa_levels <- c("Genus", "Family", "Order", "Class", "Phylum", "Kingdom")
      name <- NA
      for (level in taxa_levels) {
        if (!is.na(x[level]) && x[level] != "") {
          name <- x[level]
          break
        }
      }
      if (is.na(name)) {
        name <- rownames(x)
      }
      return(name)
    })
  }
  
  tax_labels <- get_tax_label(taxa_table)
  
  for (meta_var in meta_columns) {
    
    message(paste0("🚀 Running DESeq2 for metadata variable: ", meta_var))
    
    formula_deseq <- as.formula(paste("~", meta_var))
    
    dds <- phyloseq_to_deseq2(ps, formula_deseq)
    dds <- DESeq(dds, test = "Wald", fitType = "parametric")
    
    res <- results(dds, cooksCutoff = FALSE)
    res <- res[order(res$padj, na.last = NA), ]
    res_sig <- subset(res, padj < 0.05)
    
    res_df <- as.data.frame(res)
    res_df$Taxa <- tax_labels[rownames(res_df)]
    
    res_sig_df <- as.data.frame(res_sig)
    res_sig_df$Taxa <- tax_labels[rownames(res_sig_df)]
    
    write.csv(res_df, file.path("Output", paste0("DESeq2_all_results_", meta_var, ".csv")))
    write.csv(res_sig_df, file.path("Output", paste0("DESeq2_significant_results_", meta_var, ".csv")))
    
    if (nrow(res_sig_df) > 0) {
      p <- ggplot(res_sig_df, aes(x = reorder(Taxa, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        labs(title = paste("Differential Abundance (", meta_var, ")"),
             x = "Taxa", y = "Log2 Fold Change") +
        theme_minimal() +
        scale_fill_manual(values = c("red", "blue"))
      
      ggsave(file.path("Output", paste0("deseq2_differential_taxa_", meta_var, ".png")),
             plot = p, width = 10, height = 8, dpi = 600)
      
      print(p)
    } else {
      message(paste0("⚠️ No significant taxa found for ", meta_var))
    }
  }
}
run_deseq2(ps)


########################################### Clean Environment ###################################

rm(list = ls())
gc()
message("Script executed successfully. All results have been saved in the 'Output' folder.")
