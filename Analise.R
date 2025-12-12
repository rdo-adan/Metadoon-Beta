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

# 1. Metadata Load
meta_files <- list.files(meta_dir, pattern = "\\.(tsv|csv|txt)$", full.names = TRUE)
if (length(meta_files) == 0) stop("No metadata file found!")

line1 <- readLines(meta_files[1], n=1)
sep_char <- if(grepl("\t", line1)) "\t" else if(grepl(";", line1)) ";" else ","

meta_df <- read.table(meta_files[1], header=TRUE, row.names=1, sep=sep_char, 
                      check.names=FALSE, comment.char="", quote="", stringsAsFactors=TRUE)
rownames(meta_df) <- trimws(rownames(meta_df))

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

# --- CORE MICROBIOME (RESTORED & FIXED) ---
prevalences <- seq(0.05, 1, 0.05)
detections <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

plot_core_heatmap <- function(ps, taxrank, top_n = core_top_n) {
  tryCatch({
    ps_tax <- tax_glom(ps, taxrank)
    ps_rel <- transform_sample_counts(ps_tax, function(x) x / sum(x))
    prev_table <- microbiome::prevalence(ps_rel)
    
    # Check if we have taxa
    if(length(prev_table) == 0) return(NULL)
    
    top_taxa <- names(sort(prev_table, decreasing = TRUE)[1:min(length(prev_table), top_n)])
    ps_sub <- prune_taxa(top_taxa, ps_rel)
    
    # Name Fix
    tax_mat <- as(tax_table(ps_sub), "matrix")
    tn <- tax_mat[, taxrank]
    tn[is.na(tn)] <- taxa_names(ps_sub)[is.na(tn)]
    taxa_names(ps_sub) <- tn
    
    # FIX: Palette call (Only n_colors)
    pal <- get_dynamic_palette(length(taxa_names(ps_sub)))
    
    p <- plot_core(ps_sub, plot.type = "heatmap", prevalences = prevalences, detections = detections, min.prevalence = 0.1) +
      xlab("Detection Threshold (Relative Abundance)") + ylab(taxrank) +
      scale_fill_gradientn(colours = pal) + theme_minimal() +
      scale_x_discrete(breaks = as.character(detections)) +
      ggtitle(paste("Core Microbiome -", taxrank))
    
    save_p(p, paste0("core_microbiome_heatmap_", taxrank, ".png"))
  }, error = function(e) message(paste("Skipping core heatmap for", taxrank)))
}

for (rank in c("Phylum", "Class", "Order", "Family", "Genus")) {
  plot_core_heatmap(ps, rank, top_n = core_top_n)
}

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
          scale_color_gradientn(colours=palette_base)
      } else {
        df$Var <- as.factor(df$Var)
        p <- ggboxplot(df, x="Var", y="Div", color="Var", add="jitter") + 
          scale_color_manual(values=pal) + 
          ggpubr::stat_compare_means(label="p.format")
      }
      p <- p + labs(title=paste(idx, "by", var)) + theme_minimal()
      save_p(p, paste0("alpha_diversity_", idx, "_", var, ".png"))
    }
  }
}, error=function(e) message("Alpha Diversity Failed"))

# --- BETA DIVERSITY (MANUAL CONSTRUCTION) ---
dist_obj <- tryCatch(phyloseq::distance(ps, method=dist_method), error=function(e) NULL)

if(!is.null(dist_obj)) {
  for(method in c("NMDS", "PCoA")) {
    ord <- tryCatch(ordinate(ps, method=method, distance=dist_obj), error=function(e) NULL)
    if(is.null(ord)) next
    
    meta_cols <- colnames(sample_data(ps))
    
    for(var in meta_cols) {
      raw_val <- sample_data(ps)[[var]]
      if(all(is.na(raw_val))) next 
      
      coord_df <- NULL
      if(method == "NMDS") {
        coord_df <- as.data.frame(ord$points)
        colnames(coord_df) <- c("Axis1", "Axis2")
      } else { 
        coord_df <- as.data.frame(ord$vectors[, 1:2])
        colnames(coord_df) <- c("Axis1", "Axis2")
      }
      
      coord_df$Group <- raw_val[match(rownames(coord_df), rownames(sample_data(ps)))]
      is_numeric_col <- is.numeric(coord_df$Group) && length(unique(coord_df$Group)) > 5
      
      p <- ggplot(coord_df, aes(x=Axis1, y=Axis2)) + theme_minimal()
      
      if(is_numeric_col) {
        p <- p + geom_point(aes(color = Group), size=4, alpha=0.8) +
          scale_color_gradientn(colours=viridis::viridis(10))
      } else {
        coord_df$Group <- as.factor(coord_df$Group)
        n_grps <- length(unique(coord_df$Group))
        pal <- get_dynamic_palette(n_grps)
        
        p <- p + geom_point(aes(color = Group), size=4, alpha=0.8) +
          scale_color_manual(values=pal)
        
        try({
          if(nrow(coord_df) >= 5 && n_grps > 1 && n_grps < nrow(coord_df)) {
            p <- p + stat_ellipse(aes(color = Group), level=0.95, linetype=2)
          }
        }, silent=TRUE)
      }
      
      perm_title <- ""
      try({
        if(length(unique(na.omit(raw_val))) > 1) {
          idx <- !is.na(raw_val)
          if(sum(idx) > 1) {
            valid_s <- rownames(sample_data(ps))[idx]
            d_mat <- as.matrix(dist_obj)
            d_sub <- as.dist(d_mat[valid_s, valid_s])
            groups_sub <- raw_val[idx]
            res_perm <- vegan::adonis2(d_sub ~ groups_sub)$`Pr(>F)`[1]
            perm_title <- paste("| PERMANOVA p =", signif(res_perm, 3))
          }
        }
      }, silent=TRUE)
      
      p <- p + labs(title=paste(method, "-", var), subtitle=perm_title, color=var)
      save_p(p, paste0("beta_diversity_", method, "_", var, ".png"))
    }
  }
}

# --- RDA (NUMERIC ONLY) ---
run_global_rda <- function(ps) {
  meta <- as(sample_data(ps), "data.frame")
  nums <- sapply(meta, is.numeric)
  num_cols <- names(nums)[nums]
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
    
    sites$Group <- "All Samples"
    pal <- get_dynamic_palette(1)
    
    p <- ggplot(sites, aes(x=RDA1, y=RDA2)) +
      geom_segment(data=arrows, aes(x=0, y=0, xend=RDA1, yend=RDA2), arrow=arrow(length=unit(0.2,"cm")), color="red") +
      ggrepel::geom_text_repel(data=arrows, aes(label=Label), color="red") +
      geom_point(aes(color=as.factor(Group)), size=3) +
      scale_color_manual(values=pal) + theme_minimal() + labs(title="RDA")
    save_p(p, "beta_diversity_RDA_Env_Vectors.png")
  }, error=function(e) NULL)
}
run_global_rda(ps)

# --- DESEQ2 (FIXED SCOPE) ---
meta_cols <- colnames(sample_data(ps))
ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "OTU")

for(rank in ranks) {
  if(rank=="OTU") ps_run <- ps else ps_run <- tryCatch(tax_glom(ps, rank, NArm=F), error=function(e) NULL)
  if(is.null(ps_run)) next
  
  for(var in meta_cols) {
    valid_samples <- rownames(sample_data(ps_run))[!is.na(sample_data(ps_run)[[var]])]
    vals <- sample_data(ps_run)[[var]][!is.na(sample_data(ps_run)[[var]])]
    group_counts <- table(vals)
    valid_groups <- names(group_counts[group_counts >= 2])
    
    if(length(valid_groups) < 2) next 
    
    ps_clean <- prune_samples(sample_data(ps_run)[[var]] %in% valid_groups, ps_run)
    
    tryCatch({
      dds <- phyloseq_to_deseq2(ps_clean, as.formula(paste("~", var)))
      dds <- estimateSizeFactors(dds, type="poscounts")
      dds <- DESeq(dds, test="Wald", fitType="parametric", quiet=TRUE)
      
      pairs <- combn(valid_groups, 2, simplify=FALSE)
      
      for(pair in pairs) {
        res <- results(dds, contrast=c(var, pair[1], pair[2]), cooksCutoff=FALSE)
        comp <- paste0(var, "_", pair[1], "_vs_", pair[2])
        
        res_df <- as.data.frame(res)
        write.csv(res_df, file.path(output_dir, paste0("DESeq2_ALL_", comp, "_", rank, ".csv")))
        
        sig <- subset(res_df, padj < 0.05)
        if(nrow(sig) > 0) {
          if(rank!="OTU") sig$Taxon <- tax_table(ps_clean)[rownames(sig), rank] else sig$Taxon <- rownames(sig)
          top <- head(sig[order(abs(sig$log2FoldChange), decreasing=T),], 20)
          p <- ggplot(top, aes(x=reorder(Taxon, log2FoldChange), y=log2FoldChange, fill=log2FoldChange>0)) +
            geom_col() + coord_flip() + theme_minimal() + scale_fill_manual(values=c("firebrick", "royalblue")) +
            labs(title=paste("DESeq2:", comp))
          save_p(p, paste0("deseq2_diff_", comp, "_", rank, ".png"))
        }
      }
    }, error=function(e) NULL)
  }
}

# --- ANCOM-BC (FIXED) ---
for(rank in ranks) {
  if(rank=="OTU") ps_run <- ps else ps_run <- tryCatch(tax_glom(ps, rank, NArm=F), error=function(e) NULL)
  if(is.null(ps_run)) next
  
  for(var in meta_cols) {
    vals <- sample_data(ps_run)[[var]]
    group_counts <- table(vals)
    valid_groups <- names(group_counts[group_counts >= 2])
    
    if(length(valid_groups) < 2) next 
    
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
        
        p <- ggplot(df_long, aes(x=Group, y=Taxon, fill=LFC)) + geom_tile() +
          scale_fill_gradient2() + theme_minimal() + labs(title=paste("ANCOM-BC:", var))
        save_p(p, paste0("ancombc_heatmap_", var, "_", rank, ".png"))
      }
    }, error=function(e) NULL)
  }
}

message("Analysis Script Completed Successfully!")
rm(list = ls())
gc()