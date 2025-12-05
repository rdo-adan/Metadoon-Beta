################################### Loading Packages ########################################
install.packages("rprojroot")
required_packages <- c("phyloseq", "vegan", "ggplot2", "ggpubr", "cowplot", "dplyr", 
                       "DESeq2", "scater", "rprojroot", "pairwiseAdonis", "pheatmap", 
                       "viridis", "ape", "microbiome", "wesanderson", "RColorBrewer", 
                       "ANCOMBC")

# Instala/Carrega pacotes silenciosamente
suppressPackageStartupMessages({
  sapply(required_packages, require, character.only = TRUE)
})

########################################### Configuration ###################################
library(jsonlite)
library(rprojroot)

# Localiza a raiz do projeto
script_dir <- tryCatch(find_root(has_file("Analise.R")), error = function(e) getwd())
print(paste("Working directory:", script_dir))

`%||%` <- function(a, b) if (!is.null(a)) a else b 

params_file <- "pipeline_params.json"

if (file.exists(params_file)) {
  params <- fromJSON(params_file)
} else {
  # Defaults de segurança se não houver JSON
  warning("pipeline_params.json not found! Using defaults.")
  params <- list()
}

# Assign variables with fallbacks
stat_test <- "kruskal.test"  # Default non-parametric test for Alpha Diversity
dist_method <- "bray"        # Default distance for Beta Diversity
color_palette <- params$color_palette %||% "viridis"
abundance_top_n <- params$abundance_top_n %||% 15
core_top_n <- params$core_top_n %||% 30
rarefaction_step <- params$rarefaction_step %||% 100
rarefaction_cex <- params$rarefaction_cex %||% 0.6
enable_rarefaction <- params$enable_rarefaction %||% FALSE
rarefaction_depth <- params$rarefaction_depth %||% 1000

# --- FIX: PALETTE ERROR ---
get_palette <- function(palette_name) {
  switch(palette_name,
         viridis = viridis::viridis(10),
         plasma = viridis::plasma(10),
         inferno = viridis::inferno(10),
         magma = viridis::magma(10),
         cividis = viridis::cividis(10),
         # Wesanderson 'Darjeeling1' tem apenas 5 cores. Pedir 10 causava erro.
         wesanderson = wesanderson::wes_palette("Darjeeling1", 5, type = "discrete"),
         RColorBrewer::brewer.pal(10, "Set3")
  )
}
palette_base <- get_palette(color_palette)

get_dynamic_palette <- function(base_palette, n_colors) {
  if (n_colors > length(base_palette)) {
    return(colorRampPalette(base_palette)(n_colors))
  } else {
    return(base_palette[1:n_colors])
  }
}

########################################### Loading Files ###################################

# Definição dos caminhos (CRUCIAL PARA O SEU ERRO)
otu_file <- file.path(script_dir, "OTUs/otutab.txt")
taxonomy_file <- file.path(script_dir, "Taxonomy/taxonomy.txt")
metadata_dir <- file.path(script_dir, "Metadata File")
tree_file_path <- file.path(script_dir, "Tree File/tree.nwk")

# Output directory
output_dir <- file.path(script_dir, "Output")
if (!dir.exists(output_dir)) dir.create(output_dir)

save_plot <- function(plot, filename, width = 12, height = 8) {
  ggsave(filename = file.path(output_dir, filename), plot = plot, width = width, height = height, dpi = 600, bg = "white")
}

# Load Metadata
metadata_files <- list.files(path = metadata_dir, pattern = "\\.tsv$|\\.csv$", full.names = TRUE)
if (length(metadata_files) > 0) {
  sep <- ifelse(grepl("\\.tsv$", metadata_files[1]), "\t", ",")
  metadataTab <- read.csv(metadata_files[1], header = TRUE, row.names = 1, sep = sep, check.names = FALSE)
  message("Metadata loaded.")
} else {
  stop("No metadata file found.")
}

########################################### Create `phyloseq` Object ###################################

# 1. Load Tables (COM CORREÇÃO PARA CABEÇALHOS)
message("Loading OTU and Taxonomy tables...")

# 'comment.char=""' é crucial porque o VSEARCH coloca # no cabeçalho
otu_mat <- as.matrix(read.table(otu_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, comment.char = ""))
tax_mat <- as.matrix(read.table(taxonomy_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, fill = TRUE, comment.char = ""))

# 2. Sync IDs (Fixes "Component taxa/OTU names do not match")
common_taxa <- intersect(rownames(otu_mat), rownames(tax_mat))

if(length(common_taxa) == 0) {
  # Debug info se der erro
  print("Head OTU Table:")
  print(head(rownames(otu_mat)))
  print("Head Tax Table:")
  print(head(rownames(tax_mat)))
  stop("CRITICAL ERROR: No common OTUs found between OTU Table and Taxonomy.")
}

message(paste("Syncing: Keeping", length(common_taxa), "common OTUs."))
otu_mat_clean <- otu_mat[common_taxa, ]
tax_mat_clean <- tax_mat[common_taxa, ]

# 3. Sync Samples
if(exists("metadataTab")){
  common_samples <- intersect(colnames(otu_mat_clean), rownames(metadataTab))
  if(length(common_samples) == 0){
    warning("Warning: No common samples between OTU table and Metadata!")
    metadataTab_clean <- metadataTab
  } else {
    otu_mat_clean <- otu_mat_clean[, common_samples]
    metadataTab_clean <- metadataTab[common_samples, , drop=FALSE]
    message(paste("Syncing: Keeping", length(common_samples), "common samples."))
  }
} else {
  metadataTab_clean <- NULL
}

# 4. Load Tree (Optional)
Arv <- NULL
if (file.exists(tree_file_path)) {
  message("Tree file found. Loading...")
  tryCatch({
    Arv <- read.tree(tree_file_path)
    message("Tree loaded successfully.")
  }, error = function(e) {
    warning(paste("Error loading tree:", e$message))
    Arv <- NULL
  })
}

# 5. Create Phyloseq Object
if (is.null(Arv)) {
  ps <- phyloseq(
    otu_table(otu_mat_clean, taxa_are_rows = TRUE),
    tax_table(tax_mat_clean),
    sample_data(metadataTab_clean)
  )
} else {
  ps <- phyloseq(
    otu_table(otu_mat_clean, taxa_are_rows = TRUE),
    tax_table(tax_mat_clean),
    sample_data(metadataTab_clean),
    phy_tree(Arv)
  )
}

print(ps)
########################################### Functions and Plots ###################################

# --- RAREFACTION ---
rarefy_phyloseq <- function(ps, depth = rarefaction_depth) {
  if (enable_rarefaction) {
    message(paste("Applying rarefaction to depth:", depth))
    min_reads <- min(sample_sums(ps))
    if (depth > min_reads) {
      warning(paste("Adjusting rarefaction depth to minimum reads:", min_reads))
      depth <- min_reads
    }
    ps_rarefied <- rarefy_even_depth(ps, sample.size = depth, replace = TRUE, verbose = FALSE)
    return(ps_rarefied)
  } else {
    message("Rarefaction disabled.")
    return(ps)
  }
}

# Generate Curve (using raw data)
otu_curve <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu_curve <- t(otu_curve)
tryCatch({
  png(filename = file.path(output_dir, "rarefaction_curve.png"), width = 12, height = 8, units = "in", res = 300)
  vegan::rarecurve(t(otu_curve), step = rarefaction_step, cex = rarefaction_cex, label = TRUE)
  title("Rarefaction Curve")
  dev.off()
}, error = function(e) message("Error plotting rarefaction curve"))

# Apply rarefaction for downstream analysis
ps <- rarefy_phyloseq(ps)


# --- RELATIVE ABUNDANCE ---
plot_abundance <- function(ps, taxrank, top_n = abundance_top_n) {
  tryCatch({
    ps_tax <- tax_glom(ps, taxrank = taxrank)
    ps_rel <- transform_sample_counts(ps_tax, function(x) x / sum(x))
    taxa_sums <- taxa_sums(ps_rel)
    top_taxa <- names(sort(taxa_sums, decreasing = TRUE)[1:top_n])
    tax_table(ps_rel)[!(taxa_names(ps_rel) %in% top_taxa), taxrank] <- "Others"
    
    n_colors_needed <- length(unique(psmelt(ps_rel)[[taxrank]]))
    pal <- get_dynamic_palette(palette_base, n_colors_needed)
    
    p <- ggplot(psmelt(ps_rel), aes(x = Sample, y = Abundance, fill = !!sym(taxrank))) +
      geom_bar(stat = "identity") +
      labs(title = paste("Relative Abundance by", taxrank), x = "Sample", y = "Relative Abundance") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_manual(values = pal)
    
    save_plot(p, paste0("abundance_", taxrank, ".png"))
  }, error = function(e) message(paste("Error in abundance plot for", taxrank)))
}

for (rank in c("Phylum", "Class", "Order", "Family", "Genus")) {
  plot_abundance(ps, rank)
}


# --- CORE MICROBIOME ---
prevalences <- seq(0.05, 1, 0.05)
detections <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

plot_core_heatmap <- function(ps, taxrank, top_n = core_top_n) {
  tryCatch({
    ps_tax <- tax_glom(ps, taxrank)
    ps_tax_rel <- transform_sample_counts(ps_tax, function(x) x / sum(x))
    prev_table <- microbiome::prevalence(ps_tax_rel)
    
    # Check if we have taxa
    if(length(prev_table) == 0) return(NULL)
    
    top_taxa <- names(sort(prev_table, decreasing = TRUE)[1:min(length(prev_table), top_n)])
    ps_top_taxa <- prune_taxa(top_taxa, ps_tax_rel)
    
    # Replace names
    tax_mat <- as(tax_table(ps_top_taxa), "matrix")
    tax_names <- tax_mat[, taxrank]
    tax_names[is.na(tax_names)] <- taxa_names(ps_top_taxa)[is.na(tax_names)]
    taxa_names(ps_top_taxa) <- tax_names
    
    n_colors <- length(taxa_names(ps_top_taxa))
    pal <- get_dynamic_palette(palette_base, n_colors)
    
    p <- plot_core(ps_top_taxa, 
                   plot.type = "heatmap", 
                   prevalences = prevalences, 
                   detections = detections, 
                   min.prevalence = 0.2) + # Ajustado para ser menos estrito
      xlab("Detection Threshold (Relative Abundance)") + 
      ylab(taxrank) +
      scale_fill_gradientn(colours = pal) + 
      theme_minimal() +
      scale_x_discrete(breaks = as.character(detections)) +
      ggtitle(paste("Core Microbiome -", taxrank))
    
    save_plot(p, paste0("core_microbiome_heatmap_", taxrank, ".png"))
  }, error = function(e) message(paste("Skipping core heatmap for", taxrank, "(insufficient data)")))
}

for (rank in c("Phylum", "Class", "Order", "Family", "Genus")) {
  plot_core_heatmap(ps, rank, top_n = core_top_n)
}


# --- ALPHA DIVERSITY ---
plot_alpha_diversity <- function(ps) {
  alpha_div <- estimate_richness(ps, measures = c("Shannon", "Simpson", "Chao1", "ACE", "Observed"))
  alpha_div$SampleID <- rownames(alpha_div)
  
  meta_data <- as(sample_data(ps), "data.frame")
  meta_data$SampleID <- rownames(meta_data)
  
  alpha_merged <- merge(alpha_div, meta_data, by = "SampleID")
  meta_columns <- setdiff(colnames(meta_data), "SampleID")
  
  # Clean test name
  clean_test <- tolower(stat_test)
  if (clean_test %in% c("t-test", "ttest")) clean_test <- "t.test"
  if (clean_test %in% c("wilcox")) clean_test <- "wilcox.test"
  if (clean_test %in% c("kruskal")) clean_test <- "kruskal.test"
  if (clean_test %in% c("anova")) clean_test <- "anova"
  method_for_title <- ifelse(clean_test %in% c("t.test", "anova"), "anova", "kruskal.test")
  
  for (meta_var in meta_columns) {
    if (length(unique(meta_data[[meta_var]])) < 2) next
    
    pal <- get_dynamic_palette(palette_base, length(unique(meta_data[[meta_var]])))
    
    for (index in c("Observed", "Chao1", "ACE", "Shannon", "Simpson")) {
      plot_df <- data.frame(Sample = alpha_merged$SampleID, Diversity = alpha_merged[[index]], Group = as.factor(alpha_merged[[meta_var]]))
      
      p_label <- ""
      tryCatch({
        stat_res <- compare_means(Diversity ~ Group, data = plot_df, method = method_for_title)
        p_label <- paste0("(p = ", stat_res$p.format, ")")
      }, error = function(e) { p_label <<- "" })
      
      p <- ggboxplot(plot_df, x = "Group", y = "Diversity", color = "Group", add = "jitter") +
        labs(title = paste("Alpha Diversity -", index, p_label), x = meta_var, y = index) +
        theme_minimal() +
        scale_color_manual(values = pal)
      
      save_plot(p, paste0("alpha_diversity_", index, "_", meta_var, ".png"), width = 8, height = 6)
    }
  }
}
plot_alpha_diversity(ps)


# --- BETA DIVERSITY (NMDS & PCoA) ---
plot_beta <- function(ps, method_ord, filename_prefix) {
  meta_data <- as(sample_data(ps), "data.frame")
  meta_data$SampleID <- rownames(meta_data)
  meta_columns <- setdiff(colnames(meta_data), "SampleID")
  dist_mat <- phyloseq::distance(ps, method = dist_method)
  
  # Ordination
  ord <- tryCatch(ordinate(ps, method = method_ord, distance = dist_mat), error=function(e) NULL)
  if(is.null(ord)) return()
  
  for (meta_var in meta_columns) {
    if (length(unique(meta_data[[meta_var]])) < 2) next
    
    pal <- get_dynamic_palette(palette_base, length(unique(meta_data[[meta_var]])))
    
    # Permanova
    perm_p <- tryCatch({
      f <- as.formula(paste("dist_mat ~", meta_var))
      res <- vegan::adonis2(f, data = meta_data, permutations = 999)
      signif(res$`Pr(>F)`[1], 3)
    }, error = function(e) "NA")
    
    p <- plot_ordination(ps, ord, color = meta_var) +
      stat_ellipse(geom = "polygon", alpha = 0.1, show.legend = FALSE) +
      geom_point(size = 3) +
      labs(title = paste("Beta Diversity -", method_ord, "by", meta_var),
           subtitle = paste("PERMANOVA p =", perm_p)) +
      theme_minimal() +
      scale_color_manual(values = pal)
    
    save_plot(p, paste0("beta_diversity_", filename_prefix, "_", meta_var, ".png"), width = 10, height = 8)
  }
}

plot_beta(ps, "NMDS", "NMDS")
plot_beta(ps, "PCoA", "PCoA")


# --- DESEQ2 ---
run_deseq2 <- function(ps) {
  meta_data <- as(sample_data(ps), "data.frame")
  meta_columns <- setdiff(colnames(meta_data), "SampleID")
  
  tax_labels <- apply(as(tax_table(ps), "matrix"), 1, function(x) {
    paste(na.omit(x[c("Genus", "Family", "Order", "Class")]), collapse=";")
  })
  
  gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)) }
  
  for (meta_var in meta_columns) {
    if (length(unique(meta_data[[meta_var]])) < 2) next
    
    tryCatch({
      # Convert phyloseq to DESeq2
      dds <- phyloseq_to_deseq2(ps, as.formula(paste("~", meta_var)))
      
      # Handle zeros
      geoMeans = apply(counts(dds), 1, gm_mean)
      estimateSizeFactors(dds, geoMeans = geoMeans)
      
      # Run DESeq
      dds <- DESeq(dds, test = "Wald", fitType = "parametric", quiet = TRUE)
      res <- results(dds, cooksCutoff = FALSE)
      res_sig <- subset(res, padj < 0.05)
      
      if (nrow(res_sig) > 0) {
        res_sig_df <- as.data.frame(res_sig)
        # Add simpler Taxa names for plot
        res_sig_df$Taxa <- rownames(res_sig_df) 
        
        write.csv(res_sig_df, file.path(output_dir, paste0("DESeq2_significant_", meta_var, ".csv")))
        
        p <- ggplot(res_sig_df, aes(x = reorder(Taxa, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
          geom_bar(stat = "identity") +
          coord_flip() +
          labs(title = paste("Differential Abundance -", meta_var), x = "Taxa", y = "Log2 Fold Change") +
          theme_minimal() +
          scale_fill_manual(values = c("red", "blue"), labels = c("Decrease", "Increase"), name = "Change")
        
        save_plot(p, paste0("deseq2_differential_", meta_var, ".png"))
      }
    }, error = function(e) message(paste("DESeq2 failed for", meta_var, ":", e$message)))
  }
}

run_deseq2(ps)


# --- ANCOM-BC ANALYSIS ---
run_ancombc <- function(ps) {
  
  meta_data <- as(sample_data(ps), "data.frame")
  meta_columns <- setdiff(colnames(meta_data), "SampleID")
  
  for (meta_var in meta_columns) {
    
    # 1. First check: Do we have at least 2 groups defined?
    if (length(unique(na.omit(meta_data[[meta_var]]))) < 2) next
    
    cat(paste0("\n#### ANCOM-BC Analysis for: ", meta_var, "\n"))
    
    tryCatch({
      
      # 2. Filter Groups with n < 2 (CRITICAL FIX)
      # Calculate sample counts per group
      group_counts <- table(meta_data[[meta_var]])
      valid_groups <- names(group_counts[group_counts >= 2])
      
      # If less than 2 valid groups remain, skip analysis
      if(length(valid_groups) < 2) {
        cat(paste0("\n> **Skipped:** Not enough groups with replication (n >= 2). ",
                   "Groups detected with n=1: ", 
                   paste(names(group_counts[group_counts < 2]), collapse=", "), ".\n\n"))
        next # Skip to next variable
      }
      
      # 3. Subset data to keep only valid groups
      # We verify if the sample belongs to a valid group
      keep_samples <- sample_data(ps)[[meta_var]] %in% valid_groups
      ps_sub <- prune_samples(keep_samples, ps)
      
      # Remove empty taxa after subsetting
      ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
      
      # Inform user if groups were dropped
      dropped_groups <- setdiff(names(group_counts), valid_groups)
      if(length(dropped_groups) > 0) {
        cat(paste0("\n> *Note: Groups dropped due to low sample size (n<2): ", 
                   paste(dropped_groups, collapse=", "), "*\n\n"))
      }
      
      # 4. Run ANCOM-BC
      out <- ancombc(data = ps_sub, 
                     formula = meta_var, 
                     p_adj_method = "holm", 
                     lib_cut = 0, 
                     group = meta_var, 
                     struc_zero = TRUE, 
                     neg_lb = TRUE,
                     tol = 1e-5, 
                     max_iter = 100, 
                     conserve = TRUE, 
                     alpha = 0.05, 
                     global = TRUE)
      
      res <- out$res
      
      # 5. Process Results
      diff_df <- res$diff_abn
      sig_taxa_indices <- apply(diff_df, 1, any)
      
      if (sum(sig_taxa_indices) > 0) {
        
        beta_df <- res$lfc[sig_taxa_indices, , drop = FALSE]
        qval_df <- res$q_val[sig_taxa_indices, , drop = FALSE]
        
        df_plot <- data.frame(
          Taxa = rownames(beta_df),
          LogFoldChange = beta_df[,1], 
          QValue = qval_df[,1]
        )
        
        df_plot <- df_plot[order(df_plot$LogFoldChange), ]
        df_plot$Taxa <- factor(df_plot$Taxa, levels = df_plot$Taxa)
        
        write.csv(res$lfc, file.path(output_dir, paste0("ANCOMBC_lfc_", meta_var, ".csv")))
        
        p <- ggplot(df_plot, aes(x = Taxa, y = LogFoldChange, fill = LogFoldChange > 0)) +
          geom_bar(stat = "identity") +
          coord_flip() +
          labs(title = paste("ANCOM-BC -", meta_var),
               subtitle = "Significant Taxa (Bias Corrected)",
               x = "Taxa", y = "Log Fold Change") +
          theme_minimal() +
          scale_fill_manual(values = c("red", "blue"), labels = c("Decrease", "Increase"), name = "Direction")
        
        save_plot(p, paste0("ancombc_differential_", meta_var, ".png"))
        print(p)
        
        cat(paste0("\n\n> **Result:** ANCOM-BC identified ", nrow(df_plot), " significant taxa.\n\n"))
        
      } else {
        cat(paste0("\n\n> No significant taxa found by ANCOM-BC for **", meta_var, "**.\n\n"))
      }
      
    }, error = function(e) {
      cat(paste0("\n\n> **⚠️ ANCOM-BC Error for ", meta_var, ":** ", e$message, "\n\n"))
    })
  }
}

run_ancombc(ps)


message("Analysis Script Completed Successfully!")
rm(list = ls())
gc()
