################################### Loading Packages ########################################

required_packages <- c("phyloseq", "vegan", "ggplot2", "ggpubr", "cowplot", "dplyr", 
                       "DESeq2", "scater", "rprojroot", "pairwiseAdonis", "pheatmap", 
                       "viridis", "ape", "microbiome", "wesanderson", "RColorBrewer", 
                       "ANCOMBC", "ggrepel", "tibble", "tidyr")

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

metadata_files <- list.files(path = metadata_dir, pattern = "\\.(tsv|csv|txt)$", full.names = TRUE, ignore.case = TRUE)

if (length(metadata_files) > 0) {
  target_file <- metadata_files[1] 
  message(paste("Processing metadata file:", basename(target_file)))
  con_raw <- file(target_file, "rb")
  bom <- readBin(con_raw, "raw", 2)
  close(con_raw)
  is_utf16le <- (length(bom) >= 2 && bom[1] == as.raw(0xff) && bom[2] == as.raw(0xfe))
  if (is_utf16le) {
    message("Encoding detected: UTF-16LE (Windows format). Converting to UTF-8...")
    con <- file(target_file, open = "r", encoding = "UTF-16LE")
    file_lines <- readLines(con, warn = FALSE)
    close(con)
  } else {
    message("Encoding detected: Standard (UTF-8/ASCII).")
    file_lines <- readLines(target_file, warn = FALSE)
  }
  
  file_lines <- iconv(file_lines, to = "UTF-8")
  file_lines <- gsub("\r", "", file_lines) 
  header_line <- file_lines[1]
  raw_header <- utf8ToInt(header_line)
  n_tab <- sum(raw_header == 9)
  n_comma <- sum(raw_header == 44)
  n_semi <- sum(raw_header == 59)
  
  if (n_tab > n_comma && n_tab > n_semi) {
    detected_sep <- "\t"
    sep_name <- "TAB"
  } else if (n_semi > n_comma) {
    detected_sep <- ";"
    sep_name <- "SEMICOLON"
  } else {
    detected_sep <- ","
    sep_name <- "COMMA"
  }
  
  message(paste("Separator detected:", sep_name))
  
  
  metadataTab <- read.table(text = file_lines, 
                            header = TRUE, 
                            row.names = 1, 
                            sep = detected_sep, 
                            check.names = FALSE, 
                            fill = TRUE, 
                            quote = "",
                            comment.char = "")
  
  
  rownames(metadataTab) <- trimws(rownames(metadataTab))
  
  message("Metadata loaded successfully.")
  
} else {
  stop("No metadata file (csv, tsv, txt) found in the specified directory.")
}
########################################### Create `phyloseq` Object ###################################

# 1. Load Tables
message("Loading OTU and Taxonomy tables...")

otu_mat <- as.matrix(read.table(otu_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, comment.char = ""))
tax_mat <- as.matrix(read.table(taxonomy_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, fill = TRUE, comment.char = ""))

colnames(otu_mat) <- trimws(colnames(otu_mat))
rownames(tax_mat) <- trimws(rownames(tax_mat))
rownames(otu_mat) <- trimws(rownames(otu_mat))

message("--- Debugging Sample Names ---")
message("OTU Table Columns (First 5):")
print(head(colnames(otu_mat), 5))


common_taxa <- intersect(rownames(otu_mat), rownames(tax_mat))

if(length(common_taxa) == 0) {
  # Debug info se der erro
  print("Head OTU Table IDs:")
  print(head(rownames(otu_mat)))
  print("Head Tax Table IDs:")
  print(head(rownames(tax_mat)))
  stop("CRITICAL ERROR: No common OTUs found between OTU Table and Taxonomy.")
}

message(paste("Syncing: Keeping", length(common_taxa), "common OTUs."))
otu_mat_clean <- otu_mat[common_taxa, ]
tax_mat_clean <- tax_mat[common_taxa, ]


if(exists("metadataTab")){
  # Limpa nomes do metadado também
  rownames(metadataTab) <- trimws(rownames(metadataTab))
  
  message("Metadata Rows (First 5):")
  print(head(rownames(metadataTab), 5))
  
  common_samples <- intersect(colnames(otu_mat_clean), rownames(metadataTab))
  
  if(length(common_samples) == 0){
    message("\n!!! SAMPLE MISMATCH DETECTED !!!")
    message("Samples in OTU Table:")
    print(colnames(otu_mat_clean))
    message("Samples in Metadata:")
    print(rownames(metadataTab))
    stop("Error: No common samples between OTU table and Metadata! Check names above.")
  } else {
    otu_mat_clean <- otu_mat_clean[, common_samples]
    metadataTab_clean <- metadataTab[common_samples, , drop=FALSE]
    message(paste("Syncing: Keeping", length(common_samples), "common samples."))
  }
} else {
  metadataTab_clean <- NULL
}


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


# --- ALPHA DIVERSITY (Intelligent: Boxplot for Groups / Scatter + Gradient for Numeric) ---
plot_alpha_diversity <- function(ps) {
  
  if (!dir.exists("Output")) dir.create("Output")
  
  # Calculate indices
  alpha_div <- estimate_richness(ps, measures = c("Shannon", "Simpson", "Chao1", "ACE", "Observed"))
  alpha_div$SampleID <- rownames(alpha_div)
  
  meta_data <- as(sample_data(ps), "data.frame")
  meta_data$SampleID <- rownames(meta_data)
  
  alpha_merged <- merge(alpha_div, meta_data, by = "SampleID")
  meta_columns <- setdiff(colnames(meta_data), "SampleID")
  
  target_indices <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson")
  available_indices <- intersect(target_indices, colnames(alpha_merged))
  
  # Determine Statistical Method for Categorical Data
  clean_test <- tolower(stat_test)
  if (clean_test %in% c("t-test", "ttest", "t_test")) clean_test <- "t.test"
  if (clean_test %in% c("wilcox", "wilcoxon")) clean_test <- "wilcox.test"
  if (clean_test %in% c("kruskal", "kruskal.test")) clean_test <- "kruskal.test"
  if (clean_test %in% c("anova")) clean_test <- "anova"
  
  # Fallback for title label
  method_for_title <- ifelse(clean_test %in% c("t.test", "anova"), "anova", "kruskal.test")
  
  for (meta_var in meta_columns) {
    
    # Extract the column data to check type
    raw_vals <- meta_data[[meta_var]]
    # Remove NAs just to check properties
    clean_vals <- na.omit(raw_vals)
    
    # Skip if empty or only 1 value
    if (length(unique(clean_vals)) < 2) next
    
    # CHECK: Is it Numeric or Categorical?
    is_numeric_var <- is.numeric(clean_vals)
    
    cat(paste0("Generating Alpha Diversity for variable: ", meta_var, " (Type: ", ifelse(is_numeric_var, "Numeric", "Categorical"), ")\n"))
    
    # Dynamic Palette Setup
    if (!is_numeric_var) {
      pal <- get_dynamic_palette(palette_base, length(unique(clean_vals)))
    }
    
    for (index in available_indices) {
      
      # Create plotting dataframe
      # Note: We use the raw column for 'Group' to preserve numeric type if applicable
      plot_df <- data.frame(
        Sample = alpha_merged$SampleID, 
        Diversity = alpha_merged[[index]], 
        Variable = alpha_merged[[meta_var]] 
      )
      plot_df <- na.omit(plot_df)
      
      # --- PLOTTING LOGIC ---
      
      if (is_numeric_var) {
        # === NUMERIC VARIABLES: SCATTER PLOT + GRADIENT ===
        
        p <- ggplot(plot_df, aes(x = Variable, y = Diversity, color = Variable)) +
          # Scatter points with gradient color
          geom_point(size = 3, alpha = 0.8) +
          
          # Trend line (Linear Model) - shows correlation direction
          geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE, alpha = 0.2) +
          
          # Automatic Correlation Stats (R and p-value) on the plot
          ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
          
          # Gradient Scale using user's palette base
          scale_color_gradientn(colours = palette_base) +
          
          labs(
            title = paste("Alpha Diversity (Correlation) -", index), 
            subtitle = paste("Variable:", meta_var),
            x = meta_var,
            y = index,
            color = meta_var
          ) +
          theme_minimal(base_size = 14)
        
      } else {
        # === CATEGORICAL VARIABLES: BOXPLOT (Old Logic) ===
        
        # Ensure it is a factor for boxplot
        plot_df$Variable <- as.factor(plot_df$Variable)
        
        # Calc Stats
        p_label <- ""
        tryCatch({
          stat_res <- compare_means(Diversity ~ Variable, data = plot_df, method = method_for_title)
          p_label <- paste0("(p = ", stat_res$p.format, ")")
        }, error = function(e) { p_label <<- "" })
        
        p <- ggboxplot(plot_df, x = "Variable", y = "Diversity", color = "Variable", add = "jitter") +
          labs(
            title = paste("Alpha Diversity (Group Comparison) -", index, p_label), 
            x = meta_var,
            y = index,
            color = meta_var
          ) +
          theme_minimal(base_size = 14) +
          scale_color_manual(values = pal) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
      
      save_plot(p, paste0("alpha_diversity_", index, "_", meta_var, ".png"), width = 8, height = 6)
    }
  }
}

plot_alpha_diversity(ps)


# --- BETA DIVERSITY (NMDS & PCoA) ---
plot_beta <- function(ps, method_ord, filename_prefix) {
  
  if (!dir.exists("Output")) dir.create("Output")

  meta_data <- as(sample_data(ps), "data.frame")
  meta_data$SampleID <- rownames(meta_data)
  meta_columns <- setdiff(colnames(meta_data), "SampleID")
  

  dist_mat <- phyloseq::distance(ps, method = dist_method)

  ord <- tryCatch({
    ordinate(ps, method = method_ord, distance = dist_mat)
  }, error = function(e) NULL)
  
  if (is.null(ord)) {
    message(paste("Skipping", method_ord, "- Ordination failed."))
    return()
  }

  stress_info <- ""
  if (method_ord == "NMDS") {
    stress_val <- tryCatch(ord$stress, error = function(e) NA)
    if (!is.na(stress_val)) {
      stress_info <- paste("| Stress =", signif(stress_val, 3))
    }
  }
  
  for (meta_var in meta_columns) {

    valid_idx <- !is.na(meta_data[[meta_var]])
    if (sum(valid_idx) < 2) next
    

    val_vector <- meta_data[[meta_var]][valid_idx]
    is_numeric <- is.numeric(val_vector)

    if (!is_numeric && length(unique(val_vector)) < 2) next
    
    perm_p <- tryCatch({
      dist_sub <- dist_subset(dist_mat, valid_idx) 
      meta_sub <- meta_data[valid_idx, ]
      
      f <- as.formula(paste("dist_sub ~", meta_var))
      res <- vegan::adonis2(f, data = meta_sub, permutations = 999)
      signif(res$`Pr(>F)`[1], 3)
    }, error = function(e) "NA")

    p <- plot_ordination(ps, ord, type = "samples", color = meta_var) +
      geom_point(size = 3, alpha = 0.8) +
      labs(title = paste(method_ord, "-", meta_var),
           subtitle = paste("PERMANOVA p =", perm_p, stress_info),
           color = meta_var) +
      theme_minimal(base_size = 14)
    

    if (is_numeric) {

      p <- p + scale_color_gradientn(colours = palette_base)
    } else {

      n_groups <- length(unique(val_vector))
      pal <- get_dynamic_palette(palette_base, n_groups)
      
      p <- p + 
        stat_ellipse(geom = "polygon", alpha = 0.1, show.legend = FALSE, level = 0.95) +
        scale_color_manual(values = pal)
    }
    
    save_plot(p, paste0("beta_diversity_", filename_prefix, "_", meta_var, ".png"), width = 10, height = 8)
    print(p)
  }
}


dist_subset <- function(d, idx) {
  m <- as.matrix(d)
  as.dist(m[idx, idx])
}

plot_beta(ps, "NMDS", "NMDS")
plot_beta(ps, "PCoA", "PCoA")


# --- GLOBAL RDA (Numeric Arrows + Categorical Colors) ---
run_global_rda <- function(ps) {
  
  if (!dir.exists("Output")) dir.create("Output")
  if(!require("ggrepel", quietly=TRUE)) library(ggrepel)
  
  # 1. Prepare Data
  ps_hel <- transform_sample_counts(ps, function(x) sqrt(x / sum(x)))
  otu_mat <- t(as(otu_table(ps_hel), "matrix"))
  meta_df <- as(sample_data(ps_hel), "data.frame")
  
  # 2. Separate Numeric vs Categorical Variables
  numeric_cols <- c()
  factor_cols <- c()
  
  for(col in setdiff(colnames(meta_df), "SampleID")) {
    val <- meta_df[[col]]
    if(all(is.na(val))) next
    
    if(is.numeric(val) && length(unique(val)) > 1) {
      numeric_cols <- c(numeric_cols, col)
    } else if(length(unique(val)) > 1 && length(unique(val)) < nrow(meta_df)) {
      factor_cols <- c(factor_cols, col)
    }
  }
  
  # Need at least 1 numeric variable for arrows to make sense in this context
  if (length(numeric_cols) < 1) {
    message("Skipping Global RDA: No numeric environmental variables found for arrows.")
    return(NULL)
  }
  
  cat(paste0("Running RDA with Environmental Variables (Arrows): ", paste(numeric_cols, collapse=", "), "\n"))
  
  # Subset Data (Complete cases for numeric vars)
  meta_env <- meta_df[, numeric_cols, drop=FALSE]
  valid_idx <- complete.cases(meta_env)
  
  if(sum(valid_idx) < 3) {
    message("Skipping RDA: Too few complete samples.")
    return(NULL)
  }
  
  otu_use <- otu_mat[valid_idx, ]
  meta_env_use <- meta_env[valid_idx, , drop=FALSE]
  meta_cat_use <- meta_df[valid_idx, factor_cols, drop=FALSE]
  
  tryCatch({
    # 3. Run RDA constrained ONLY by Numeric Variables (Environmental)
    # This creates the arrows for pH, Temp, etc.
    rda_model <- vegan::rda(otu_use ~ ., data = meta_env_use, scale = TRUE)
    
    # Check
    if(is.null(rda_model$CCA) || rda_model$CCA$rank == 0) return()
    
    # 4. Stats
    anova_global <- vegan::anova.cca(rda_model, permutations = 999)
    p_val <- anova_global$`Pr(>F)`[1]
    r2_adj <- vegan::RsquareAdj(rda_model)$adj.r.squared
    r2_txt <- if(is.na(r2_adj) || r2_adj < 0) "< 1%" else paste0(round(r2_adj * 100, 1), "%")
    
    # 5. Extract Scores
    sc <- scores(rda_model, display = c("sites", "bp"), scaling = 2)
    
    df_sites <- as.data.frame(sc$sites)
    df_arrows <- as.data.frame(sc$biplot)
    df_arrows$Label <- rownames(df_arrows)
    
    # Scale arrows to fit plot
    mult <- min(
      (max(df_sites$RDA1) - min(df_sites$RDA1)) / (max(df_arrows$RDA1) - min(df_arrows$RDA1)),
      (max(df_sites$RDA2) - min(df_sites$RDA2)) / (max(df_arrows$RDA2) - min(df_arrows$RDA2))
    ) * 0.8
    
    df_arrows$RDA1 <- df_arrows$RDA1 * mult
    df_arrows$RDA2 <- df_arrows$RDA2 * mult
    
    # Axis Labels
    eig <- rda_model$CCA$eig
    var_pct <- round(eig / sum(rda_model$tot.chi) * 100, 1)
    
    # 6. Plot Loop (Color by each Categorical Variable)
    # We generate one plot for each categorical variable (Type, Biome, etc.)
    # showing the same RDA arrows (pH, Temp) but colored differently.
    
    # If no categorical vars, create a dummy group
    if(length(factor_cols) == 0) {
      factor_cols <- c("None")
      meta_cat_use$None <- "All Samples"
    }
    
    for (cat_var in factor_cols) {
      
      df_sites$Group <- meta_cat_use[[cat_var]]
      n_groups <- length(unique(df_sites$Group))
      pal <- get_dynamic_palette(palette_base, n_groups)
      
      p <- ggplot(df_sites, aes(x = RDA1, y = RDA2)) +
        geom_vline(xintercept = 0, linetype="dashed", color="grey80") +
        geom_hline(yintercept = 0, linetype="dashed", color="grey80") +
        
        # Arrows (Environmental Vectors)
        geom_segment(data = df_arrows, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                     arrow = arrow(length = unit(0.2, "cm")), color = "red", linewidth = 0.8) +
        
        geom_text_repel(data = df_arrows, aes(x = RDA1, y = RDA2, label = Label),
                        color = "red", fontface = "bold", bg.color = "white", bg.r = 0.15) +
        
        # Samples
        geom_point(aes(color = as.factor(Group), fill = as.factor(Group)), size = 3, shape = 21, alpha = 0.8) +
        
        # Ellipses
        stat_ellipse(aes(color = as.factor(Group)), geom="polygon", alpha=0.1, show.legend = FALSE, level=0.95) +
        
        scale_color_manual(values = pal, name = cat_var) +
        scale_fill_manual(values = pal, name = cat_var) +
        
        labs(title = paste("RDA - Environment vs", cat_var),
             subtitle = paste("Vectors: ", paste(numeric_cols, collapse=", "), "\nStats: p =", signif(p_val, 3), "| R2.adj =", r2_txt),
             x = paste0("RDA1 (", var_pct[1], "%)"),
             y = paste0("RDA2 (", var_pct[2], "%)")) +
        theme_minimal(base_size = 14)
      
      save_plot(p, paste0("beta_diversity_RDA_Env_vs_", cat_var, ".png"), width = 10, height = 8)
      print(p)
    }
    
  }, error = function(e) {
    message(paste("Global RDA Failed:", e$message))
  })
}

run_global_rda(ps)


# --- DESEQ2 (Smart Labels for Numeric Variables) ---
run_deseq2 <- function(ps) {
  
  ranks_to_analyze <- c("OTU", "Genus", "Family", "Order", "Class", "Phylum")
  meta_data <- as(sample_data(ps), "data.frame")
  meta_columns <- setdiff(colnames(meta_data), "SampleID")
  full_tax_table <- as.data.frame(tax_table(ps))
  
  get_pretty_name_vectorized <- function(ids, tax_df, current_rank) {
    sapply(ids, function(id) {
      if(!id %in% rownames(tax_df)) return(id)
      if (current_rank != "OTU") {
        val <- tax_df[id, current_rank]
        if (!is.na(val) && val != "" && val != "NA") return(val)
      }
      row <- tax_df[id, ]
      for (r in c("Species", "Genus", "Family", "Order", "Class", "Phylum")) {
        val <- row[[r]]
        if (!is.na(val) && val != "" && val != "NA") {
          if (current_rank == "OTU") return(paste0(val, " (", id, ")"))
          return(val)
        }
      }
      return(id)
    })
  }
  
  for (rank in ranks_to_analyze) {
    cat(paste0("\n[DESeq2] Processing Rank: ", rank, "\n"))
    
    if (rank == "OTU") { ps_run <- ps } else {
      ps_run <- tryCatch(tax_glom(ps, taxrank = rank, NArm = FALSE), error = function(e) NULL)
    }
    if (is.null(ps_run)) next
    
    local_tax_table <- as.data.frame(tax_table(ps_run))
    
    for (meta_var in meta_columns) {
      
    
      is_numeric <- is.numeric(meta_data[[meta_var]])
      
     
      if (!is_numeric && length(unique(na.omit(meta_data[[meta_var]]))) < 2) next
      
  
      if (is_numeric) {
        pairs <- list(c(meta_var, "Numeric")) # Dummy pair for loop
      } else {
        groups <- unique(as.character(na.omit(meta_data[[meta_var]])))
        pairs <- combn(groups, 2, simplify = FALSE)
      }
      
      tryCatch({
       
        if (!is_numeric) {
          ps_sub <- subset_samples(ps_run, !is.na(get(meta_var)))
          dds <- phyloseq_to_deseq2(ps_sub, as.formula(paste("~", meta_var)))
        } else {
          
          valid_samples <- rownames(meta_data)[!is.na(meta_data[[meta_var]])]
          ps_sub <- prune_samples(valid_samples, ps_run)
          dds <- phyloseq_to_deseq2(ps_sub, as.formula(paste("~", meta_var)))
        }
        
        dds <- estimateSizeFactors(dds, type = "poscounts")
        dds <- DESeq(dds, test="Wald", fitType="parametric", quiet=TRUE)
        
        for (pair in pairs) {
          
          if (is_numeric) {
            
            res <- results(dds, name = meta_var, cooksCutoff = FALSE)
            comp_name <- paste0(meta_var, "_correlation")
            title_text <- paste("DESeq2 Correlation:", meta_var)
            legend_labels <- c(paste("Decreases with", meta_var), paste("Increases with", meta_var))
            
          } else {
            g1 <- pair[1]; g2 <- pair[2]
            res <- results(dds, contrast = c(meta_var, g1, g2), cooksCutoff = FALSE)
            comp_name <- paste0(meta_var, "_", g1, "_vs_", g2)
            title_text <- paste("DESeq2:", g1, "vs", g2)
            legend_labels <- c(paste("Higher in", g2), paste("Higher in", g1))
          }
          
          res_df <- as.data.frame(res)
          
        
          ids <- rownames(res_df)
          if (rank == "OTU") {
            names_vec <- get_pretty_name_vectorized(ids, full_tax_table, "OTU")
          } else {
            names_vec <- get_pretty_name_vectorized(ids, local_tax_table, rank)
          }
          res_df$Taxon <- names_vec
          res_df$OTU_ID <- ids
          
          write.csv(res_df %>% select(Taxon, OTU_ID, everything()), 
                    file.path(output_dir, paste0("DESeq2_ALL_", comp_name, "_", rank, ".csv")), row.names=FALSE)
          
          # Plot Signif
          sig <- subset(res_df, padj < 0.05)
          
          if (nrow(sig) > 0) {
            if (nrow(sig) > 20) sig <- head(sig[order(abs(sig$log2FoldChange), decreasing=TRUE),], 20)
            sig$Taxon <- factor(sig$Taxon, levels = sig$Taxon[order(sig$log2FoldChange)])
            
            p <- ggplot(sig, aes(x = Taxon, y = log2FoldChange, fill = log2FoldChange > 0)) +
              geom_bar(stat = "identity") + coord_flip() +
              labs(title = title_text, subtitle = paste("Rank:", rank), x = "", y = "Log2 Fold Change") +
              scale_fill_manual(values = c("red", "blue"), labels = legend_labels, name = "Effect") +
              theme_minimal()
            
            save_plot(p, paste0("deseq2_diff_", comp_name, "_", rank, ".png"))
            print(p)
          }
        }
      }, error=function(e) message(paste("DESeq2 error:", e$message)))
    }
  }
}
run_deseq2(ps)

# --- ANCOM-BC (Smart Labels for Numeric) ---
run_ancombc <- function(ps) {
  require(dplyr); require(tidyr)
  
  ranks <- c("OTU", "Genus", "Family", "Order", "Class", "Phylum")
  meta_data <- as(sample_data(ps), "data.frame")
  meta_columns <- setdiff(colnames(meta_data), "SampleID")
  full_tax <- as.data.frame(tax_table(ps))
  
  get_pretty_name_vectorized <- function(ids, tax_df, current_rank) {
    sapply(ids, function(id) {
      if(!id %in% rownames(tax_df)) return(id)
      if (current_rank != "OTU") {
        val <- tax_df[id, current_rank]
        if (!is.na(val) && val != "" && val != "NA") return(val)
      }
      row <- tax_df[id, ]
      for (r in c("Species", "Genus", "Family", "Order", "Class", "Phylum")) {
        val <- row[[r]]
        if (!is.na(val) && val != "" && val != "NA") {
          if (current_rank == "OTU") return(paste0(val, " (", id, ")"))
          return(val)
        }
      }
      return(id)
    })
  }
  
  for (rank in ranks) {
    cat(paste0("\n[ANCOM-BC] Processing Rank: ", rank, "\n"))
    
    if (rank == "OTU") { ps_run <- ps } else {
      ps_run <- tryCatch(tax_glom(ps, taxrank = rank, NArm = FALSE), error=function(e) NULL)
    }
    if (is.null(ps_run)) next
    
    local_tax <- as.data.frame(tax_table(ps_run))
    
    for (var in meta_columns) {
      
      # Check if Numeric or Categorical
      is_numeric <- is.numeric(meta_data[[var]])
      
      # Pre-check
      if (!is_numeric && length(unique(na.omit(meta_data[[var]]))) < 2) next
      
      tryCatch({
        # Filter small groups only if categorical
        if (!is_numeric) {
          valid <- names(which(table(meta_data[[var]]) >= 2))
          if(length(valid) < 2) next
          keep <- sample_data(ps_run)[[var]] %in% valid
          ps_sub <- prune_samples(keep, ps_run)
        } else {
          # For numeric, just remove NAs
          valid_s <- rownames(meta_data)[!is.na(meta_data[[var]])]
          ps_sub <- prune_samples(valid_s, ps_run)
        }
        ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
        
        # Run
        out <- ancombc(data=ps_sub, formula=var, p_adj_method="holm", group=var, struc_zero=TRUE, neg_lb=TRUE, tol=1e-5, max_iter=100, conserve=TRUE, alpha=0.05, global=TRUE)
        res <- out$res
        
        # Translate Names
        lfc_df <- as.data.frame(res$lfc)
        ids <- if("taxon" %in% colnames(lfc_df)) lfc_df$taxon else rownames(lfc_df)
        
        if (rank == "OTU") names_vec <- get_pretty_name_vectorized(ids, full_tax, "OTU")
        else names_vec <- get_pretty_name_vectorized(ids, local_tax, rank)
        
        lfc_df$Taxon <- names_vec; lfc_df$ID <- ids
        write.csv(lfc_df %>% select(Taxon, ID, everything()), file.path(output_dir, paste0("ANCOMBC_LFC_", var, "_", rank, ".csv")), row.names=FALSE)
        
        # Plot
        sig_idx <- apply(res$diff_abn, 1, any, na.rm=TRUE)
        if (sum(sig_idx) > 0) {
          beta <- lfc_df[sig_idx, , drop=FALSE]
          
          # If Numeric: Simple Bar Plot with Correlation Logic
          if (is_numeric) {
            df_plot <- data.frame(Label = beta$Taxon, LFC = beta[[var]]) # Numeric vars usually have 1 column named after var
            
            if(nrow(df_plot) > 20) df_plot <- head(df_plot[order(abs(df_plot$LFC), decreasing=TRUE),], 20)
            df_plot$Label <- factor(df_plot$Label, levels = df_plot$Label[order(df_plot$LFC)])
            
            p <- ggplot(df_plot, aes(x=Label, y=LFC, fill=LFC>0)) +
              geom_bar(stat="identity", width=0.7) + coord_flip() +
              scale_fill_manual(values=c("red", "blue"), 
                                labels=c(paste("Decreases with", var), paste("Increases with", var)), 
                                name="Correlation") +
              labs(title=paste("ANCOM-BC Correlation:", var, "-", rank), x="", y="Log Fold Change (per unit)") +
              theme_minimal()
            
            save_plot(p, paste0("ancombc_diff_", var, "_", rank, ".png"))
            print(p)
            
          } else {
            # --- CATEGORICAL LOGIC (Groups) ---
            
            # 1. Prepare data columns (remove Taxon names and IDs for calculation)
            data_cols <- beta %>% select(-Taxon, -ID)
            
            # 2. Determine Plot Type based on number of contrasts
            # If there is only 1 column of results, it means 2 groups (Case vs Control).
            # If there are multiple columns, it means >2 groups (A vs Ref, B vs Ref...).
            
            if (ncol(data_cols) == 1) {
              # === CASE A: BAR PLOT (2 Groups) ===
              
              val_col_name <- colnames(data_cols)[1]
              
              df_plot <- data.frame(
                Label = beta$Taxon,
                LFC = data_cols[[1]]
              )
              
              # Filter Top 20 by absolute change
              if(nrow(df_plot) > 20) {
                df_plot <- df_plot[order(abs(df_plot$LFC), decreasing = TRUE), ]
                df_plot <- head(df_plot, 20)
              }
              
              # Order factors for plotting
              df_plot$Label <- factor(df_plot$Label, levels = df_plot$Label[order(df_plot$LFC)])
              
              p <- ggplot(df_plot, aes(x = Label, y = LFC, fill = LFC > 0)) +
                geom_bar(stat = "identity", width = 0.7) +
                coord_flip() +
                labs(title = paste("ANCOM-BC -", var, "(", rank, ")"),
                     subtitle = paste("Contrast:", val_col_name),
                     x = "", y = "Log Fold Change") +
                theme_minimal(base_size = 12) +
                scale_fill_manual(values = c("red", "blue"), labels = c("Decrease", "Increase"), name = "Direction") +
                theme(axis.text.y = element_text(size=10))
              
              save_plot(p, paste0("ancombc_bar_", var, "_", rank, ".png"))
              print(p)
              
            } else {
              # === CASE B: HEATMAP (>2 Groups) ===
              
              # Clean column names (remove variable prefix for cleaner heatmap)
              # e.g., "BiomeUrban" becomes "Urban"
              clean_cols <- data_cols
              colnames(clean_cols) <- gsub(var, "", colnames(clean_cols))
              
              # Bind taxa names back
              df_heatmap <- bind_cols(Taxon = beta$Taxon, clean_cols)
              
              # Filter Top 30 by variance (most variable taxa across groups)
              if(nrow(df_heatmap) > 30) {
                row_vars <- apply(clean_cols, 1, var)
                top_30_idx <- head(order(row_vars, decreasing = TRUE), 30)
                df_heatmap <- df_heatmap[top_30_idx, ]
              }
              
              # Pivot to long format for ggplot
              df_long <- df_heatmap %>%
                pivot_longer(cols = -Taxon, names_to = "Group", values_to = "LFC")
              
              # Order Taxa by mean LFC for hierarchical-like visualization
              tax_order <- df_long %>% 
                group_by(Taxon) %>% 
                summarise(mean = mean(LFC)) %>% 
                arrange(mean) %>% 
                pull(Taxon)
              
              df_long$Taxon <- factor(df_long$Taxon, levels = tax_order)
              
              # Define color limits centered at 0
              limit <- max(abs(df_long$LFC)) * c(-1, 1)
              
              p <- ggplot(df_long, aes(x = Group, y = Taxon, fill = LFC)) +
                geom_tile(color = "white") +
                scale_fill_gradient2(low = "firebrick", mid = "white", high = "royalblue", 
                                     midpoint = 0, limit = limit, name = "LFC") +
                geom_text(aes(label = round(LFC, 1)), color = "black", size = 3) +
                labs(title = paste("ANCOM-BC Heatmap:", var, "-", rank),
                     subtitle = "Log Fold Change vs Reference Group",
                     x = "Comparison", y = "") +
                theme_minimal(base_size = 12) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
              
              save_plot(p, paste0("ancombc_heatmap_", var, "_", rank, ".png"))
              print(p)
            }
          }
        }
      }, error = function(e) message(paste("ANCOM-BC Error:", e$message)))
    }
  }
}
run_ancombc(ps)



message("Analysis Script Completed Successfully!")
rm(list = ls())
gc()
