################################### Loading Packages ########################################

required_packages <- c("phyloseq", "vegan", "ggplot2", "ggpubr", "cowplot", "dplyr", 
                       "DESeq2", "scater", "rprojroot", "pheatmap", 
                       "viridis", "ape", "microbiome", "wesanderson", "RColorBrewer", 
                       "ANCOMBC", "ggrepel", "tibble", "tidyr")

suppressPackageStartupMessages({
  sapply(required_packages, require, character.only = TRUE)
})

library(jsonlite)
library(rprojroot)

# --- SETUP ---
script_dir <- tryCatch(find_root(has_file("Analise.R")), error = function(e) getwd())
output_dir <- file.path(script_dir, "Output")
if (!dir.exists(output_dir)) dir.create(output_dir)

# Load Params
params_file <- "pipeline_params.json"
params <- if (file.exists(params_file)) fromJSON(params_file) else list()

# Variables
`%||%` <- function(a, b) if (!is.null(a)) a else b 
stat_test <- "kruskal.test"
dist_method <- "bray"
color_palette <- params$color_palette %||% "viridis"
abundance_top_n <- params$abundance_top_n %||% 15
core_top_n <- params$core_top_n %||% 30
rarefaction_depth <- params$rarefaction_depth %||% 1000
enable_rarefaction <- params$enable_rarefaction %||% FALSE

# --- PALETTE FUNCTION ---
get_base_palette <- function(name) {
  switch(name,
         "viridis" = viridis::viridis(10),
         "plasma" = viridis::plasma(10),
         "magma" = viridis::magma(10),
         "wesanderson" = wesanderson::wes_palette("Darjeeling1", 5, type = "discrete"),
         "Set3" = RColorBrewer::brewer.pal(10, "Set3"),
         viridis::viridis(10)
  )
}
palette_base <- get_base_palette(color_palette)

get_dynamic_palette <- function(n_colors) {
  if (is.null(n_colors) || is.na(n_colors) || n_colors < 1) n_colors <- 1
  if (n_colors <= length(palette_base)) {
    return(palette_base[1:n_colors])
  } else {
    return(colorRampPalette(palette_base)(n_colors))
  }
}

################################### DATA LOADING ###################################

otu_file <- file.path(script_dir, "OTUs/otutab.txt")
tax_file <- file.path(script_dir, "Taxonomy/taxonomy.txt")
meta_dir <- file.path(script_dir, "Metadata File")
tree_file <- file.path(script_dir, "Tree File/tree.nwk")

# 1. Metadata Load (UTF-16 Safe)
meta_files <- list.files(meta_dir, pattern = "\\.(tsv|csv|txt)$", full.names = TRUE)
if (length(meta_files) == 0) stop("No metadata file found!")
target_file <- meta_files[1]

message(paste("Reading Metadata:", basename(target_file)))

con_raw <- file(target_file, "rb")
bom <- readBin(con_raw, "raw", 2)
close(con_raw)

is_utf16 <- (length(bom) >= 2 && bom[1] == as.raw(0xff) && bom[2] == as.raw(0xfe))

if(is_utf16) {
  message("Detected Windows UTF-16LE encoding. Converting...")
  con <- file(target_file, open="r", encoding="UTF-16LE")
  file_content <- readLines(con, warn=FALSE)
  close(con)
} else {
  file_content <- readLines(target_file, warn=FALSE)
}

file_content <- iconv(file_content, to="UTF-8")
file_content <- gsub("\r", "", file_content)

header <- file_content[1]
sep_char <- if(grepl("\t", header)) "\t" else if(grepl(";", header)) ";" else ","

meta_df <- read.table(text=file_content, header=TRUE, row.names=1, sep=sep_char, 
                      check.names=FALSE, comment.char="", quote="", stringsAsFactors=TRUE)
rownames(meta_df) <- trimws(rownames(meta_df))

message("Metadata Loaded Successfully.")

# 2. OTU & Tax Load
otu_mat <- as.matrix(read.table(otu_file, header=T, row.names=1, sep="\t", check.names=F, comment.char=""))
tax_mat <- as.matrix(read.table(tax_file, header=T, row.names=1, sep="\t", check.names=F, fill=T, comment.char=""))

# 3. Sync
common_samples <- intersect(colnames(otu_mat), rownames(meta_df))
if(length(common_samples) == 0) stop("CRITICAL: No matching samples between Metadata and OTU table.")

otu_mat <- otu_mat[, common_samples, drop=FALSE]
meta_df <- meta_df[common_samples, , drop=FALSE]

common_taxa <- intersect(rownames(otu_mat), rownames(tax_mat))
otu_mat <- otu_mat[common_taxa, , drop=FALSE]
tax_mat <- tax_mat[common_taxa, , drop=FALSE]

# 4. Tree
phy_tree_obj <- NULL
if(file.exists(tree_file)) try({ phy_tree_obj <- read.tree(tree_file) }, silent=TRUE)

# 5. Phyloseq
if(is.null(phy_tree_obj)) {
  ps <- phyloseq(otu_table(otu_mat, taxa_are_rows=T), tax_table(tax_mat), sample_data(meta_df))
} else {
  ps <- phyloseq(otu_table(otu_mat, taxa_are_rows=T), tax_table(tax_mat), sample_data(meta_df), phy_tree(phy_tree_obj))
}

################################### ANALYSIS ###################################

save_p <- function(p, name) ggsave(file.path(output_dir, name), p, width=10, height=8, bg="white")

# --- RAREFACTION ---
tryCatch({
  if(enable_rarefaction) {
    min_r <- min(sample_sums(ps))
    if(rarefaction_depth > min_r) rarefaction_depth <- min_r
    ps_rare <- rarefy_even_depth(ps, sample.size=rarefaction_depth, verbose=F, replace=T)
    ps <- ps_rare
  }
  otu_curve <- as(otu_table(ps), "matrix")
  if(taxa_are_rows(ps)) otu_curve <- t(otu_curve)
  png(filename = file.path(output_dir, "rarefaction_curve.png"), width = 12, height = 8, units = "in", res = 300)
  vegan::rarecurve(otu_curve, step = 100, label = TRUE)
  dev.off()
}, error=function(e) message("Rarefaction Plot Failed (Non-critical)"))

# --- ABUNDANCE ---
for(rank in c("Phylum", "Class", "Order", "Family", "Genus")) {
  tryCatch({
    ps_glom <- tax_glom(ps, rank)
    ps_rel <- transform_sample_counts(ps_glom, function(x) x/sum(x))
    top <- names(sort(taxa_sums(ps_rel), decreasing=T)[1:min(ntaxa(ps_rel), abundance_top_n)])
    tax_table(ps_rel)[!taxa_names(ps_rel) %in% top, rank] <- "Others"
    
    df <- psmelt(ps_rel)
    pal <- get_dynamic_palette(length(unique(df[[rank]])))
    
    p <- ggplot(df, aes(x=Sample, y=Abundance, fill=!!sym(rank))) +
      geom_bar(stat="identity") + scale_fill_manual(values=pal) +
      theme_minimal() + theme(axis.text.x=element_text(angle=90))
    save_p(p, paste0("abundance_", rank, ".png"))
  }, error=function(e) NULL)
}

# --- CORE MICROBIOME (FIXED) ---
plot_core_heatmap <- function(ps, taxrank, top_n = core_top_n) {
  tryCatch({
    ps_tax <- tax_glom(ps, taxrank)
    ps_rel <- transform_sample_counts(ps_tax, function(x) x / sum(x))
    
    # Relaxed Prevalence check
    prev_table <- microbiome::prevalence(ps_rel, detection = 0, sort = TRUE)
    if(length(prev_table) == 0) return(NULL)
    
    top_taxa <- names(prev_table)[1:min(length(prev_table), top_n)]
    ps_sub <- prune_taxa(top_taxa, ps_rel)
    
    tax_mat <- as(tax_table(ps_sub), "matrix")
    tn <- tax_mat[, taxrank]
    tn[is.na(tn)] <- taxa_names(ps_sub)[is.na(tn)]
    taxa_names(ps_sub) <- tn
    
    pal <- get_dynamic_palette(length(taxa_names(ps_sub)))
    
    p <- plot_core(ps_sub, plot.type = "heatmap", 
                   prevalences = seq(0.05, 1, 0.05), 
                   detections = c(0.01, 0.05, 0.1, 0.2), 
                   min.prevalence = 0.01) + # Very relaxed threshold
      xlab("Detection Threshold") + ylab(taxrank) +
      scale_fill_gradientn(colours = pal) + theme_minimal() +
      ggtitle(paste("Core Microbiome -", taxrank))
    
    save_p(p, paste0("core_microbiome_heatmap_", taxrank, ".png"))
  }, error = function(e) message(paste("Skipping core heatmap for", taxrank)))
}
for (rank in c("Phylum", "Class", "Order", "Family", "Genus")) plot_core_heatmap(ps, rank, top_n = core_top_n)

# --- ALPHA DIVERSITY ---
tryCatch({
  alpha <- estimate_richness(ps, measures = c("Shannon", "Simpson", "Chao1", "ACE", "Observed"))
  alpha$SampleID <- rownames(alpha)
  meta <- as(sample_data(ps), "data.frame")
  meta$SampleID <- rownames(meta)
  merged <- merge(alpha, meta, by = "SampleID")
  cols <- setdiff(colnames(meta), "SampleID")
  
  for (var in cols) {
    clean_vals <- na.omit(merged[[var]])
    if (length(unique(clean_vals)) < 2) next
    
    is_num <- is.numeric(clean_vals)
    pal <- if(!is_num) get_dynamic_palette(length(unique(clean_vals))) else NULL
    
    for (idx in c("Observed", "Chao1", "ACE", "Shannon", "Simpson")) {
      df <- data.frame(Sample=merged$SampleID, Div=merged[[idx]], Var=merged[[var]])
      df <- na.omit(df)
      
      if(is_num) {
        p <- ggplot(df, aes(x=Var, y=Div, color=Var)) + geom_point(size=3) + 
          geom_smooth(method="lm") + ggpubr::stat_cor() + 
          scale_color_gradientn(colours=palette_base) +
          labs(title=paste(idx, "by", var), x=var, color=var)
      } else {
        df$Var <- as.factor(df$Var)
        p <- ggboxplot(df, x="Var", y="Div", color="Var", add="jitter") + 
          scale_color_manual(values=pal) + 
          ggpubr::stat_compare_means(label="p.format") +
          labs(title=paste(idx, "by", var), x=var, color=var)
      }
      p <- p + theme_minimal()
      save_p(p, paste0("alpha_diversity_", idx, "_", var, ".png"))
    }
  }
}, error=function(e) message("Alpha Diversity Failed"))

# --- BETA DIVERSITY (SAFE DATAFRAME) ---
dist_obj <- tryCatch(phyloseq::distance(ps, method=dist_method), error=function(e) NULL)

if(!is.null(dist_obj)) {
  for(method in c("NMDS", "PCoA")) {
    ord <- tryCatch(ordinate(ps, method=method, distance=dist_obj), error=function(e) NULL)
    if(is.null(ord)) next
    
    meta_cols <- colnames(sample_data(ps))
    
    for(var in meta_cols) {
      raw_val <- sample_data(ps)[[var]]
      if(all(is.na(raw_val))) next 
      
      # 1. Manual Coordinate Extraction
      coord_df <- NULL
      if(method == "NMDS") {
        coord_df <- as.data.frame(ord$points)
        colnames(coord_df) <- c("Axis1", "Axis2")
      } else { 
        coord_df <- as.data.frame(ord$vectors[, 1:2])
        colnames(coord_df) <- c("Axis1", "Axis2")
      }
      
      # 2. Add Metadata SAFELY (Merge Logic)
      coord_df$SampleID <- rownames(coord_df)
      
      # Create a mini meta dataframe to merge
      meta_mini <- data.frame(SampleID = rownames(sample_data(ps)), Group = raw_val)
      
      # Merge to ensure alignment
      plot_df <- merge(coord_df, meta_mini, by="SampleID", sort=FALSE)
      
      # 3. Determine Type
      is_numeric_col <- is.numeric(plot_df$Group) && length(unique(plot_df$Group)) > 5
      if(!is_numeric_col) plot_df$Group <- as.factor(plot_df$Group)
      
      # 4. Plot using clean 'plot_df'
      p <- ggplot(plot_df, aes(x=Axis1, y=Axis2, color=Group)) + 
        geom_point(size=4, alpha=0.8) + theme_minimal()
      
      if(is_numeric_col) {
        p <- p + scale_color_gradientn(colours=viridis::viridis(10))
      } else {
        n_grps <- length(unique(plot_df$Group))
        pal <- get_dynamic_palette(n_grps)
        p <- p + scale_color_manual(values=pal)
        
        try({
          if(nrow(plot_df) >= 5 && n_grps > 1 && n_grps < nrow(plot_df)) {
            p <- p + stat_ellipse(level=0.95, linetype=2)
          }
        }, silent=TRUE)
      }
      
      # Stats
      perm_title <- ""
      try({
        if(length(unique(na.omit(raw_val))) > 1 && sum(!is.na(raw_val)) > 1) {
          # Use plot_df to ensure we match valid samples only
          d_mat <- as.matrix(dist_obj)
          valid_s <- plot_df$SampleID
          d_sub <- as.dist(d_mat[valid_s, valid_s])
          
          res_perm <- vegan::adonis2(d_sub ~ plot_df$Group)$`Pr(>F)`[1]
          p_val_str <- if(res_perm < 0.001) "< 0.001" else signif(res_perm, 3)
          perm_title <- paste("| PERMANOVA p =", p_val_str)
        }
      }, silent=TRUE)
      
      p <- p + geom_text_repel(aes(label = SampleID), size = 3.5, show.legend=FALSE)
      p <- p + labs(title=paste(method, "-", var), subtitle=perm_title, color=var, x="Axis 1", y="Axis 2")
      save_p(p, paste0("beta_diversity_", method, "_", var, ".png"))
    }
  }
}

# --- RDA ---
run_global_rda <- function(ps) {
  meta <- as(sample_data(ps), "data.frame")
  nums <- sapply(meta, is.numeric)
  num_cols <- names(nums)[nums]
  cat_cols <- names(nums)[!nums] 
  
  if(length(num_cols) == 0) return()
  meta_env <- meta[, num_cols, drop=FALSE]
  idx <- complete.cases(meta_env)
  if(sum(idx) < 3) return()
  
  tryCatch({
    ps_hel <- transform_sample_counts(ps, function(x) sqrt(x / sum(x)))
    otu <- t(as(otu_table(ps_hel), "matrix"))
    rda_mod <- vegan::rda(otu[idx,] ~ ., data=meta_env[idx,,drop=FALSE], scale=TRUE)
    sc <- scores(rda_mod, display=c("sites", "bp"), scaling=2)
    sites <- as.data.frame(sc$sites)
    arrows <- as.data.frame(sc$biplot)
    arrows$Label <- rownames(arrows)
    
    if(length(cat_cols) == 0) { sites$Group <- "All Samples"; cat_cols <- c("Group") }
    
    for(var in cat_cols) {
      if(!var %in% colnames(meta) && var != "Group") next
      sites$Group <- if(var == "Group") "All Samples" else meta[[var]][idx]
      sites$Group <- as.factor(sites$Group)
      pal <- get_dynamic_palette(length(unique(sites$Group)))
      
      p <- ggplot(sites, aes(x=RDA1, y=RDA2)) +
        geom_segment(data=arrows, aes(x=0, y=0, xend=RDA1, yend=RDA2), arrow=arrow(length=unit(0.2,"cm")), color="red") +
        ggrepel::geom_text_repel(data=arrows, aes(label=Label), color="red") +
        geom_point(aes(color=Group), size=3) +
        scale_color_manual(values=pal) + theme_minimal() + labs(title=paste("RDA - Env vs", var), color=var)
      save_p(p, paste0("beta_diversity_RDA_Env_vs_", var, ".png"))
    }
  }, error=function(e) NULL)
}
run_global_rda(ps)

# --- DESEQ2 (FIXED: EXPLICIT LABELS) ---
meta_cols <- colnames(sample_data(ps))
ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "OTU")

for(rank in ranks) {
  if(rank=="OTU") ps_run <- ps else ps_run <- tryCatch(tax_glom(ps, rank, NArm=F), error=function(e) NULL)
  if(is.null(ps_run)) next
  
  for(var in meta_cols) {
    # Skip Numeric
    raw_vals <- sample_data(ps_run)[[var]]
    if(is.numeric(raw_vals) && length(unique(raw_vals)) > 5) next
    
    vals <- as.factor(raw_vals[!is.na(raw_vals)])
    group_counts <- table(vals)
    valid_groups <- names(group_counts[group_counts >= 2])
    if(length(valid_groups) < 2) next 
    
    sample_data(ps_run)[[var]] <- as.factor(sample_data(ps_run)[[var]])
    ps_clean <- prune_samples(sample_data(ps_run)[[var]] %in% valid_groups, ps_run)
    
    tryCatch({
      dds <- phyloseq_to_deseq2(ps_clean, as.formula(paste("~", var)))
      dds <- estimateSizeFactors(dds, type="poscounts")
      dds <- DESeq(dds, test="Wald", fitType="parametric", quiet=TRUE)
      
      pairs <- combn(valid_groups, 2, simplify=FALSE)
      
      for(pair in pairs) {
        # Contrast: Numerator = pair[1], Denominator = pair[2]
        res <- results(dds, contrast=c(var, pair[1], pair[2]), cooksCutoff=FALSE)
        comp <- paste0(var, "_", pair[1], "_vs_", pair[2])
        
        res_df <- as.data.frame(res)
        write.csv(res_df, file.path(output_dir, paste0("DESeq2_ALL_", comp, "_", rank, ".csv")))
        
        sig <- subset(res_df, padj < 0.05)
        if(nrow(sig) > 0) {
          if(rank!="OTU") sig$Taxon <- tax_table(ps_clean)[rownames(sig), rank] else sig$Taxon <- rownames(sig)
          top <- head(sig[order(abs(sig$log2FoldChange), decreasing=T),], 20)
          
          # EXPLICIT LEGEND
          # LFC > 0: Higher in pair[1] (Blue)
          # LFC < 0: Higher in pair[2] (Red)
          labels_vec <- c(paste("Enriched in", pair[2]), paste("Enriched in", pair[1]))
          
          p <- ggplot(top, aes(x=reorder(Taxon, log2FoldChange), y=log2FoldChange, fill=log2FoldChange>0)) +
            geom_col() + coord_flip() + theme_minimal() + 
            scale_fill_manual(values=c("red", "blue"), labels=labels_vec, name="Direction") +
            labs(title=paste("DESeq2:", comp))
          save_p(p, paste0("deseq2_diff_", comp, "_", rank, ".png"))
        }
      }
    }, error=function(e) NULL)
  }
}

# --- ANCOM-BC (FIXED: EXPLICIT LABELS) ---
for(rank in ranks) {
  if(rank=="OTU") ps_run <- ps else ps_run <- tryCatch(tax_glom(ps, rank, NArm=F), error=function(e) NULL)
  if(is.null(ps_run)) next
  
  for(var in meta_cols) {
    raw_vals <- sample_data(ps_run)[[var]]
    if(is.numeric(raw_vals) && length(unique(raw_vals)) > 5) next
    
    vals <- as.factor(raw_vals[!is.na(raw_vals)])
    group_counts <- table(vals)
    valid_groups <- names(group_counts[group_counts >= 2])
    if(length(valid_groups) < 2) next 
    
    sample_data(ps_run)[[var]] <- as.factor(sample_data(ps_run)[[var]])
    ps_clean <- prune_samples(sample_data(ps_run)[[var]] %in% valid_groups, ps_run)
    
    tryCatch({
      out <- ancombc(data=ps_clean, formula=var, p_adj_method="holm", group=var, struc_zero=T, neg_lb=T, tol=1e-5)
      res <- out$res
      lfc <- as.data.frame(res$lfc)
      lfc$ID <- rownames(lfc)
      write.csv(lfc, file.path(output_dir, paste0("ANCOMBC_LFC_", var, "_", rank, ".csv")))
      
      if(sum(apply(res$diff_abn, 1, any)) > 0) {
        beta <- lfc[apply(res$diff_abn, 1, any), , drop=FALSE]
        df_long <- pivot_longer(beta %>% select(-ID), cols=everything(), names_to="Group", values_to="LFC")
        
        if(rank != "OTU") {
          tn <- as.data.frame(tax_table(ps_run))[beta$ID, rank]
          df_long$Taxon <- rep(tn, each=ncol(beta)-1)
        } else {
          df_long$Taxon <- rep(beta$ID, each=ncol(beta)-1)
        }
        
        # BARPLOT IF PAIRWISE
        if(length(valid_groups) == 2) {
          col_name <- colnames(beta %>% select(-ID))[1]
          
          # Determine Reference
          ref_group <- setdiff(valid_groups, gsub(var, "", col_name))
          comp_group <- gsub(var, "", col_name)
          
          # Legend: 
          # LFC > 0: Higher in Comparison (Blue)
          # LFC < 0: Higher in Reference (Red)
          labels_vec <- c(paste("Enriched in", ref_group), paste("Enriched in", comp_group))
          
          df_plot <- beta
          if(nrow(df_plot) > 20) df_plot <- head(df_plot[order(abs(df_plot[[col_name]]), decreasing=T),], 20)
          
          if(rank != "OTU") df_plot$Taxon <- as.data.frame(tax_table(ps_run))[df_plot$ID, rank] else df_plot$Taxon <- df_plot$ID
          
          p <- ggplot(df_plot, aes(x=reorder(Taxon, .data[[col_name]]), y=.data[[col_name]], fill=.data[[col_name]]>0)) +
            geom_col() + coord_flip() + theme_minimal() +
            scale_fill_manual(values=c("red", "blue"), labels=labels_vec, name="Direction") +
            labs(title=paste("ANCOM-BC:", var), y="Log Fold Change")
          save_p(p, paste0("ancombc_heatmap_", var, "_", rank, ".png"))
        } else {
          p <- ggplot(df_long, aes(x=Group, y=Taxon, fill=LFC)) + geom_tile() +
            scale_fill_gradient2() + theme_minimal() + labs(title=paste("ANCOM-BC:", var))
          save_p(p, paste0("ancombc_heatmap_", var, "_", rank, ".png"))
        }
      }
    }, error=function(e) NULL)
  }
}

message("Analysis Script Completed Successfully!")
rm(list = ls())
gc()