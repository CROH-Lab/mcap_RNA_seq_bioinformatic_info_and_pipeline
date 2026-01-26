#!/usr/bin/env Rscript
# ==============================================================================
# Publication-Ready Figures: M. capitata and D. trenchii OA Response
# ==============================================================================

# Set working directory
setwd("/home/darmstrong4/mc_rework/12_publication_figures")

# ==============================================================================
# LOAD LIBRARIES
# ==============================================================================

cat("Loading required packages...\n")

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
    library(ggridges)
    library(VennDiagram)
    library(pheatmap)
    library(circlize)
    library(DESeq2)
    library(RColorBrewer)
    library(viridis)
    library(scales)
    library(cowplot)
    library(grid)
    library(gridExtra)
    library(ggrepel)
})

# Set publication theme
theme_pub <- theme_bw() +
    theme(
        text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_blank()
    )

# Color palettes
host_color <- "#E69F00"
symbiont_color <- "#56B4E9"
summer_color <- "#D55E00"
winter_color <- "#0072B2"
up_color <- "#D73027"
down_color <- "#4575B4"

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("\n==============================================================================\n")
cat("Loading DESeq2 Results and Annotations\n")
cat("==============================================================================\n\n")

# DESeq2 results
host_summer <- read.csv("/home/darmstrong4/mc_rework/07_deseq2/results/Mcap_Summer_BvsD.csv", row.names = 1)
host_winter <- read.csv("/home/darmstrong4/mc_rework/07_deseq2/results/Mcap_Winter_BvsD.csv", row.names = 1)
sym_summer <- read.csv("/home/darmstrong4/mc_rework/07_deseq2/results/Dtre_Summer_BvsD.csv", row.names = 1)
sym_winter <- read.csv("/home/darmstrong4/mc_rework/07_deseq2/results/Dtre_Winter_BvsD.csv", row.names = 1)

# Add gene IDs as column
host_summer$gene_id <- rownames(host_summer)
host_winter$gene_id <- rownames(host_winter)
sym_summer$gene_id <- rownames(sym_summer)
sym_winter$gene_id <- rownames(sym_winter)

# Annotations
host_annot <- read.delim("/home/darmstrong4/mc_rework/08_host_deg_annotation/results/all_annotations_full.tsv", 
                          stringsAsFactors = FALSE)
sym_annot <- read.delim("/home/darmstrong4/mc_rework/09_symbiont_deg_annotation/results/all_annotations_full.tsv",
                         stringsAsFactors = FALSE)

# Sample metadata
sample_info <- read.delim("/home/darmstrong4/mc_rework/sample_info.txt", stringsAsFactors = FALSE)

# VSD matrices (stored as matrices, not DESeqTransform objects)
host_summer_vsd <- readRDS("/home/darmstrong4/mc_rework/07_deseq2/objects/Mcapitata_Summer_vsd.rds")
host_winter_vsd <- readRDS("/home/darmstrong4/mc_rework/07_deseq2/objects/Mcapitata_Winter_vsd.rds")
sym_summer_vsd <- readRDS("/home/darmstrong4/mc_rework/07_deseq2/objects/Dtrenchii_Summer_vsd.rds")
sym_winter_vsd <- readRDS("/home/darmstrong4/mc_rework/07_deseq2/objects/Dtrenchii_Winter_vsd.rds")

# Define DEGs (padj < 0.05)
host_summer_degs <- host_summer %>% filter(!is.na(padj) & padj < 0.05)
host_winter_degs <- host_winter %>% filter(!is.na(padj) & padj < 0.05)
sym_summer_degs <- sym_summer %>% filter(!is.na(padj) & padj < 0.05)
sym_winter_degs <- sym_winter %>% filter(!is.na(padj) & padj < 0.05)

cat("DEG counts:\n")
cat("  Host Summer:", nrow(host_summer_degs), "\n")
cat("  Host Winter:", nrow(host_winter_degs), "\n")
cat("  Symbiont Summer:", nrow(sym_summer_degs), "\n")
cat("  Symbiont Winter:", nrow(sym_winter_degs), "\n")

# ==============================================================================
# FIGURE 1: Density Ridgeline Plot of log2FoldChange
# ==============================================================================

cat("\n==============================================================================\n")
cat("Figure 1: Density Ridgeline Plot\n")
cat("==============================================================================\n\n")

ridgeline_data <- bind_rows(
    host_summer %>% 
        filter(!is.na(log2FoldChange)) %>%
        mutate(Organism = "Host (M. capitata)", Season = "Summer") %>%
        select(log2FoldChange, Organism, Season),
    host_winter %>% 
        filter(!is.na(log2FoldChange)) %>%
        mutate(Organism = "Host (M. capitata)", Season = "Winter") %>%
        select(log2FoldChange, Organism, Season),
    sym_summer %>% 
        filter(!is.na(log2FoldChange)) %>%
        mutate(Organism = "Symbiont (D. trenchii)", Season = "Summer") %>%
        select(log2FoldChange, Organism, Season),
    sym_winter %>% 
        filter(!is.na(log2FoldChange)) %>%
        mutate(Organism = "Symbiont (D. trenchii)", Season = "Winter") %>%
        select(log2FoldChange, Organism, Season)
)

ridgeline_data <- ridgeline_data %>%
    mutate(Group = paste(Organism, "-", Season),
           Group = factor(Group, levels = c(
               "Symbiont (D. trenchii) - Winter",
               "Symbiont (D. trenchii) - Summer",
               "Host (M. capitata) - Winter",
               "Host (M. capitata) - Summer"
           )))

fig1_ridgeline <- ggplot(ridgeline_data, aes(x = log2FoldChange, y = Group, fill = interaction(Organism, Season))) +
    geom_density_ridges(alpha = 0.8, scale = 1.5, rel_min_height = 0.01) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray60", linewidth = 0.5) +
    scale_fill_manual(values = c(
        "Host (M. capitata).Summer" = "#E69F00",
        "Host (M. capitata).Winter" = "#F0E442",
        "Symbiont (D. trenchii).Summer" = "#56B4E9",
        "Symbiont (D. trenchii).Winter" = "#009E73"
    )) +
    scale_x_continuous(limits = c(-10, 10), breaks = seq(-10, 10, 2)) +
    labs(
        x = expression(Log[2]~Fold~Change~(OA~vs~Ambient)),
        y = "",
        title = "Distribution of Gene Expression Changes Under Ocean Acidification"
    ) +
    theme_pub +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 11, face = "bold"))

ggsave("figures/Fig1_ridgeline_log2FC.pdf", fig1_ridgeline, width = 10, height = 6, dpi = 300)
ggsave("figures/Fig1_ridgeline_log2FC.png", fig1_ridgeline, width = 10, height = 6, dpi = 300)
cat("Saved: Fig1_ridgeline_log2FC.pdf/png\n")

# ==============================================================================
# FIGURE 2: Venn Diagrams - Seasonal DEG Overlap
# ==============================================================================

cat("\n==============================================================================\n")
cat("Figure 2: Venn Diagrams\n")
cat("==============================================================================\n\n")

# Host Venn Diagram
host_summer_ids <- host_summer_degs$gene_id
host_winter_ids <- host_winter_degs$gene_id
host_overlap <- length(intersect(host_summer_ids, host_winter_ids))

pdf("figures/Fig2A_host_venn.pdf", width = 6, height = 6)
grid.newpage()
venn_host <- draw.pairwise.venn(
    area1 = length(host_summer_ids),
    area2 = length(host_winter_ids),
    cross.area = host_overlap,
    category = c("Summer", "Winter"),
    fill = c(summer_color, winter_color),
    alpha = 0.6,
    col = c(summer_color, winter_color),
    lwd = 2,
    cex = 2,
    fontface = "bold",
    cat.cex = 1.5,
    cat.fontface = "bold",
    cat.pos = c(-30, 30),
    cat.dist = c(0.05, 0.05),
    margin = 0.1
)
grid.text("Host (M. capitata) DEGs", y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()

# Symbiont Venn Diagram
sym_summer_ids <- sym_summer_degs$gene_id
sym_winter_ids <- sym_winter_degs$gene_id
sym_overlap <- length(intersect(sym_summer_ids, sym_winter_ids))

pdf("figures/Fig2B_symbiont_venn.pdf", width = 6, height = 6)
grid.newpage()
venn_sym <- draw.pairwise.venn(
    area1 = length(sym_summer_ids),
    area2 = length(sym_winter_ids),
    cross.area = sym_overlap,
    category = c("Summer", "Winter"),
    fill = c(summer_color, winter_color),
    alpha = 0.6,
    col = c(summer_color, winter_color),
    lwd = 2,
    cex = 2,
    fontface = "bold",
    cat.cex = 1.5,
    cat.fontface = "bold",
    cat.pos = c(-30, 30),
    cat.dist = c(0.05, 0.05),
    margin = 0.1
)
grid.text("Symbiont (D. trenchii) DEGs", y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()

cat("Saved: Fig2A_host_venn.pdf, Fig2B_symbiont_venn.pdf\n")
cat("  Host overlap:", host_overlap, "DEGs\n")
cat("  Symbiont overlap:", sym_overlap, "DEGs\n")

# ==============================================================================
# FIGURE 3: Heatmaps of Top DEGs
# ==============================================================================

cat("\n==============================================================================\n")
cat("Figure 3: Heatmaps\n")
cat("==============================================================================\n\n")

# Function to create heatmap from matrix
create_deg_heatmap <- function(vsd_mat, deseq_results, title_text, n_genes = 50) {
    
    # Get top DEGs by absolute log2FC
    top_degs <- deseq_results %>%
        filter(!is.na(padj) & padj < 0.05) %>%
        arrange(desc(abs(log2FoldChange))) %>%
        head(n_genes)
    
    if (nrow(top_degs) < 5) {
        cat("  Warning: Fewer than 5 DEGs for", title_text, "\n")
        return(NULL)
    }
    
    # Filter to top DEGs
    common_genes <- intersect(top_degs$gene_id, rownames(vsd_mat))
    if (length(common_genes) < 5) {
        cat("  Warning: Fewer than 5 matching genes for", title_text, "\n")
        return(NULL)
    }
    
    heatmap_mat <- vsd_mat[common_genes, , drop = FALSE]
    
    # Scale by row (z-score)
    heatmap_mat_scaled <- t(scale(t(heatmap_mat)))
    
    # Get treatment from column names (B = OA, D = Ambient)
    sample_names <- colnames(heatmap_mat)
    treatments <- ifelse(grepl("B", sample_names), "OA", "Ambient")
    
    anno_col <- data.frame(
        Treatment = treatments,
        row.names = sample_names
    )
    
    anno_colors <- list(Treatment = c("OA" = "#E69F00", "Ambient" = "#56B4E9"))
    
    return(list(
        matrix = heatmap_mat_scaled,
        annotation = anno_col,
        colors = anno_colors,
        title = title_text
    ))
}

cat("Creating heatmaps...\n")

# Host Summer
hm_host_summer <- create_deg_heatmap(host_summer_vsd, host_summer, "Host Summer", n_genes = 50)
if (!is.null(hm_host_summer)) {
    pdf("figures/Fig3A_heatmap_host_summer.pdf", width = 8, height = 10)
    pheatmap::pheatmap(hm_host_summer$matrix,
             main = "Host (M. capitata) - Summer\nTop 50 DEGs",
             color = colorRampPalette(c(down_color, "white", up_color))(100),
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             show_rownames = FALSE,
             show_colnames = TRUE,
             annotation_col = hm_host_summer$annotation,
             annotation_colors = hm_host_summer$colors,
             fontsize = 10,
             border_color = NA)
    dev.off()
    cat("  Saved: Fig3A_heatmap_host_summer.pdf\n")
}

# Host Winter
hm_host_winter <- create_deg_heatmap(host_winter_vsd, host_winter, "Host Winter", n_genes = 50)
if (!is.null(hm_host_winter)) {
    pdf("figures/Fig3B_heatmap_host_winter.pdf", width = 8, height = 10)
    pheatmap::pheatmap(hm_host_winter$matrix,
             main = "Host (M. capitata) - Winter\nTop 50 DEGs",
             color = colorRampPalette(c(down_color, "white", up_color))(100),
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             show_rownames = FALSE,
             show_colnames = TRUE,
             annotation_col = hm_host_winter$annotation,
             annotation_colors = hm_host_winter$colors,
             fontsize = 10,
             border_color = NA)
    dev.off()
    cat("  Saved: Fig3B_heatmap_host_winter.pdf\n")
}

# Symbiont Summer (may have few DEGs)
hm_sym_summer <- create_deg_heatmap(sym_summer_vsd, sym_summer, "Symbiont Summer", n_genes = 50)
if (!is.null(hm_sym_summer)) {
    pdf("figures/Fig3C_heatmap_symbiont_summer.pdf", width = 8, height = 10)
    pheatmap::pheatmap(hm_sym_summer$matrix,
             main = "Symbiont (D. trenchii) - Summer\nDEGs",
             color = colorRampPalette(c(down_color, "white", up_color))(100),
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             show_rownames = FALSE,
             show_colnames = TRUE,
             annotation_col = hm_sym_summer$annotation,
             annotation_colors = hm_sym_summer$colors,
             fontsize = 10,
             border_color = NA)
    dev.off()
    cat("  Saved: Fig3C_heatmap_symbiont_summer.pdf\n")
}

# Symbiont Winter
hm_sym_winter <- create_deg_heatmap(sym_winter_vsd, sym_winter, "Symbiont Winter", n_genes = 50)
if (!is.null(hm_sym_winter)) {
    pdf("figures/Fig3D_heatmap_symbiont_winter.pdf", width = 8, height = 10)
    pheatmap::pheatmap(hm_sym_winter$matrix,
             main = "Symbiont (D. trenchii) - Winter\nTop 50 DEGs",
             color = colorRampPalette(c(down_color, "white", up_color))(100),
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             show_rownames = FALSE,
             show_colnames = TRUE,
             annotation_col = hm_sym_winter$annotation,
             annotation_colors = hm_sym_winter$colors,
             fontsize = 10,
             border_color = NA)
    dev.off()
    cat("  Saved: Fig3D_heatmap_symbiont_winter.pdf\n")
}

# ==============================================================================
# FIGURE 4: Chord Diagram - Calcification/Ion Transport Genes
# ==============================================================================

cat("\n==============================================================================\n")
cat("Figure 4: Chord Diagram - Calcification/Ion Transport\n")
cat("==============================================================================\n\n")

calc_keywords <- c("calcif", "calcium", "carbonate", "bicarbon", "ion.transport", 
                   "homeostasis", "channel", "pump", "exchanger", "ATPase", 
                   "voltage.gated", "sodium", "potassium", "chloride", "proton",
                   "pH", "acid", "carbonic.anhydrase", "SLC", "transporter")

keyword_pattern <- paste(calc_keywords, collapse = "|")

# Find host DEGs with calcification-related annotations
host_all_degs <- bind_rows(
    host_summer_degs %>% mutate(season = "Summer"),
    host_winter_degs %>% mutate(season = "Winter")
) %>% distinct(gene_id, .keep_all = TRUE)

host_calc_degs <- host_all_degs %>%
    left_join(host_annot, by = "gene_id") %>%
    filter(grepl(keyword_pattern, description, ignore.case = TRUE)) %>%
    mutate(organism = "Host",
           short_desc = str_extract(description, "^[^,]+"))

# Find symbiont DEGs with calcification-related annotations
sym_all_degs <- bind_rows(
    sym_summer_degs %>% mutate(season = "Summer"),
    sym_winter_degs %>% mutate(season = "Winter")
) %>% distinct(gene_id, .keep_all = TRUE)

sym_calc_degs <- sym_all_degs %>%
    left_join(sym_annot, by = "gene_id") %>%
    filter(grepl(keyword_pattern, description, ignore.case = TRUE)) %>%
    mutate(organism = "Symbiont",
           short_desc = str_extract(description, "^[^,]+"))

cat("Calcification-related DEGs found:\n")
cat("  Host:", nrow(host_calc_degs), "\n")
cat("  Symbiont:", nrow(sym_calc_degs), "\n")

# Save calcification gene lists
if (nrow(host_calc_degs) > 0) {
    write.csv(host_calc_degs, "data/host_calcification_degs.csv", row.names = FALSE)
}
if (nrow(sym_calc_degs) > 0) {
    write.csv(sym_calc_degs, "data/symbiont_calcification_degs.csv", row.names = FALSE)
}

# Create chord diagram if we have genes
if (nrow(host_calc_degs) > 0 || nrow(sym_calc_degs) > 0) {
    
    # Extract key functional categories
    extract_category <- function(desc) {
        if (is.na(desc)) return("Other")
        desc_lower <- tolower(desc)
        if (grepl("calcium.channel|voltage.gated.*calcium|voltage.dependent.*calcium", desc_lower)) return("Ca2+ Channel")
        if (grepl("calcium.pump|calcium.transport.*atpase|plasma membrane calcium", desc_lower)) return("Ca2+ ATPase")
        if (grepl("calcium.bind|calmodulin|ef.hand|calcium.dependent.protein", desc_lower)) return("Ca2+ Binding")
        if (grepl("sodium|potassium|na\\+|k\\+", desc_lower)) return("Na+/K+ Transport")
        if (grepl("chloride|anion", desc_lower)) return("Anion Transport")
        if (grepl("carbonic.anhydrase", desc_lower)) return("Carbonic Anhydrase")
        if (grepl("exchanger|antiporter", desc_lower)) return("Ion Exchanger")
        if (grepl("proton|h\\+.*atpase|v.type", desc_lower)) return("H+ ATPase")
        if (grepl("transporter|solute.carrier", desc_lower)) return("Solute Carrier")
        return("Other Ion/Transport")
    }
    
    if (nrow(host_calc_degs) > 0) {
        host_calc_degs$category <- sapply(host_calc_degs$description, extract_category)
    }
    if (nrow(sym_calc_degs) > 0) {
        sym_calc_degs$category <- sapply(sym_calc_degs$description, extract_category)
    }
    
    # Create summary table for chord diagram
    host_cat_summary <- if(nrow(host_calc_degs) > 0) {
        host_calc_degs %>% 
            group_by(category) %>% 
            summarise(count = n(), avg_lfc = mean(log2FoldChange, na.rm = TRUE)) %>%
            mutate(organism = "Host")
    } else { data.frame() }
    
    sym_cat_summary <- if(nrow(sym_calc_degs) > 0) {
        sym_calc_degs %>% 
            group_by(category) %>% 
            summarise(count = n(), avg_lfc = mean(log2FoldChange, na.rm = TRUE)) %>%
            mutate(organism = "Symbiont")
    } else { data.frame() }
    
    all_cat_summary <- bind_rows(host_cat_summary, sym_cat_summary)
    
    if (nrow(all_cat_summary) > 0) {
        write.csv(all_cat_summary, "data/calcification_category_summary.csv", row.names = FALSE)
        
        # Create bar plot instead of chord if limited data
        fig4_bar <- ggplot(all_cat_summary, aes(x = reorder(category, count), y = count, fill = organism)) +
            geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
            coord_flip() +
            scale_fill_manual(values = c("Host" = host_color, "Symbiont" = symbiont_color)) +
            labs(
                x = "",
                y = "Number of DEGs",
                title = "Calcification & Ion Transport DEGs",
                fill = "Organism"
            ) +
            theme_pub
        
        ggsave("figures/Fig4_calcification_degs.pdf", fig4_bar, width = 10, height = 6)
        ggsave("figures/Fig4_calcification_degs.png", fig4_bar, width = 10, height = 6, dpi = 300)
        cat("Saved: Fig4_calcification_degs.pdf/png\n")
    }
} else {
    cat("No calcification-related DEGs found\n")
}

# ==============================================================================
# FIGURE 5: Volcano Plots
# ==============================================================================

cat("\n==============================================================================\n")
cat("Figure 5: Volcano Plots\n")
cat("==============================================================================\n\n")

create_volcano <- function(deseq_df, title_text, padj_cutoff = 0.05, lfc_cutoff = 1) {
    
    df <- deseq_df %>%
        filter(!is.na(pvalue) & !is.na(log2FoldChange)) %>%
        mutate(
            neg_log10_pval = -log10(pvalue),
            significance = case_when(
                padj < padj_cutoff & log2FoldChange > lfc_cutoff ~ "Up",
                padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Down",
                padj < padj_cutoff ~ "Sig (|LFC|<1)",
                TRUE ~ "NS"
            ),
            significance = factor(significance, levels = c("Up", "Down", "Sig (|LFC|<1)", "NS"))
        )
    
    # Cap extreme values for visualization
    df$neg_log10_pval <- pmin(df$neg_log10_pval, 50)
    df$log2FoldChange <- pmax(pmin(df$log2FoldChange, 10), -10)
    
    n_up <- sum(df$significance == "Up", na.rm = TRUE)
    n_down <- sum(df$significance == "Down", na.rm = TRUE)
    
    p <- ggplot(df, aes(x = log2FoldChange, y = neg_log10_pval, color = significance)) +
        geom_point(alpha = 0.5, size = 1) +
        geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "gray50") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
        scale_color_manual(
            values = c("Up" = up_color, "Down" = down_color, 
                       "Sig (|LFC|<1)" = "gray50", "NS" = "gray80"),
            name = "Expression"
        ) +
        labs(
            x = expression(Log[2]~Fold~Change),
            y = expression(-Log[10]~P-value),
            title = title_text,
            subtitle = paste0("Up: ", n_up, " | Down: ", n_down)
        ) +
        theme_pub +
        theme(legend.position = "right")
    
    return(p)
}

vol_host_summer <- create_volcano(host_summer, "Host (M. capitata) - Summer")
vol_host_winter <- create_volcano(host_winter, "Host (M. capitata) - Winter")
vol_sym_summer <- create_volcano(sym_summer, "Symbiont (D. trenchii) - Summer")
vol_sym_winter <- create_volcano(sym_winter, "Symbiont (D. trenchii) - Winter")

ggsave("figures/Fig5A_volcano_host_summer.pdf", vol_host_summer, width = 8, height = 6)
ggsave("figures/Fig5B_volcano_host_winter.pdf", vol_host_winter, width = 8, height = 6)
ggsave("figures/Fig5C_volcano_symbiont_summer.pdf", vol_sym_summer, width = 8, height = 6)
ggsave("figures/Fig5D_volcano_symbiont_winter.pdf", vol_sym_winter, width = 8, height = 6)

# Combined volcano plot
vol_combined <- plot_grid(
    vol_host_summer + theme(legend.position = "none"),
    vol_host_winter + theme(legend.position = "none"),
    vol_sym_summer + theme(legend.position = "none"),
    vol_sym_winter + theme(legend.position = "none"),
    ncol = 2, nrow = 2,
    labels = c("A", "B", "C", "D"),
    label_size = 16
)

legend <- get_legend(vol_host_summer + theme(legend.position = "bottom"))
vol_combined_legend <- plot_grid(vol_combined, legend, ncol = 1, rel_heights = c(1, 0.1))

ggsave("figures/Fig5_volcano_combined.pdf", vol_combined_legend, width = 12, height = 10)
ggsave("figures/Fig5_volcano_combined.png", vol_combined_legend, width = 12, height = 10, dpi = 300)
cat("Saved: Fig5_volcano_*.pdf/png\n")

# ==============================================================================
# FIGURE 6: PCA Plots
# ==============================================================================

cat("\n==============================================================================\n")
cat("Figure 6: PCA Plots\n")
cat("==============================================================================\n\n")

create_pca_from_matrix <- function(vsd_mat, title_text) {
    
    # Perform PCA on transposed matrix (samples as rows)
    pca_result <- prcomp(t(vsd_mat), scale. = FALSE)
    
    # Calculate percent variance
    percent_var <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)
    
    # Create data frame for plotting
    pca_data <- data.frame(
        PC1 = pca_result$x[, 1],
        PC2 = pca_result$x[, 2],
        sample = rownames(pca_result$x)
    )
    
    # Extract treatment from sample names (B = OA, D = Ambient)
    pca_data$Treatment <- ifelse(grepl("B", pca_data$sample), "OA", "Ambient")
    
    p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Treatment, shape = Treatment)) +
        geom_point(size = 4, alpha = 0.8) +
        stat_ellipse(level = 0.95, linetype = "dashed") +
        scale_color_manual(values = c("OA" = "#E69F00", "Ambient" = "#56B4E9"),
                           name = "Treatment") +
        scale_shape_manual(values = c("OA" = 16, "Ambient" = 17),
                           name = "Treatment") +
        labs(
            x = paste0("PC1: ", percent_var[1], "% variance"),
            y = paste0("PC2: ", percent_var[2], "% variance"),
            title = title_text
        ) +
        theme_pub +
        theme(legend.position = "right")
    
    return(p)
}

pca_host_summer <- create_pca_from_matrix(host_summer_vsd, "Host (M. capitata) - Summer")
pca_host_winter <- create_pca_from_matrix(host_winter_vsd, "Host (M. capitata) - Winter")
pca_sym_summer <- create_pca_from_matrix(sym_summer_vsd, "Symbiont (D. trenchii) - Summer")
pca_sym_winter <- create_pca_from_matrix(sym_winter_vsd, "Symbiont (D. trenchii) - Winter")

ggsave("figures/Fig6A_pca_host_summer.pdf", pca_host_summer, width = 7, height = 6)
ggsave("figures/Fig6B_pca_host_winter.pdf", pca_host_winter, width = 7, height = 6)
ggsave("figures/Fig6C_pca_symbiont_summer.pdf", pca_sym_summer, width = 7, height = 6)
ggsave("figures/Fig6D_pca_symbiont_winter.pdf", pca_sym_winter, width = 7, height = 6)

# Combined PCA plot
pca_combined <- plot_grid(
    pca_host_summer + theme(legend.position = "none"),
    pca_host_winter + theme(legend.position = "none"),
    pca_sym_summer + theme(legend.position = "none"),
    pca_sym_winter + theme(legend.position = "none"),
    ncol = 2, nrow = 2,
    labels = c("A", "B", "C", "D"),
    label_size = 16
)

pca_legend <- get_legend(pca_host_summer + theme(legend.position = "bottom"))
pca_combined_legend <- plot_grid(pca_combined, pca_legend, ncol = 1, rel_heights = c(1, 0.1))

ggsave("figures/Fig6_pca_combined.pdf", pca_combined_legend, width = 12, height = 10)
ggsave("figures/Fig6_pca_combined.png", pca_combined_legend, width = 12, height = 10, dpi = 300)
cat("Saved: Fig6_pca_*.pdf/png\n")

# ==============================================================================
# SUMMARY TABLE
# ==============================================================================

cat("\n==============================================================================\n")
cat("Creating Summary Statistics\n")
cat("==============================================================================\n\n")

summary_stats <- data.frame(
    Organism = c("Host", "Host", "Symbiont", "Symbiont"),
    Season = c("Summer", "Winter", "Summer", "Winter"),
    Total_Genes = c(nrow(host_summer), nrow(host_winter), nrow(sym_summer), nrow(sym_winter)),
    DEGs_total = c(nrow(host_summer_degs), nrow(host_winter_degs), 
                   nrow(sym_summer_degs), nrow(sym_winter_degs)),
    DEGs_up = c(sum(host_summer_degs$log2FoldChange > 0),
                sum(host_winter_degs$log2FoldChange > 0),
                sum(sym_summer_degs$log2FoldChange > 0),
                sum(sym_winter_degs$log2FoldChange > 0)),
    DEGs_down = c(sum(host_summer_degs$log2FoldChange < 0),
                  sum(host_winter_degs$log2FoldChange < 0),
                  sum(sym_summer_degs$log2FoldChange < 0),
                  sum(sym_winter_degs$log2FoldChange < 0))
)

write.csv(summary_stats, "data/DEG_summary_statistics.csv", row.names = FALSE)
print(summary_stats)

# ==============================================================================
# SESSION INFO
# ==============================================================================

cat("\n==============================================================================\n")
cat("Analysis Complete!\n")
cat("==============================================================================\n\n")

cat("Output files saved to:\n")
cat("  figures/ - All PDF and PNG figures\n")
cat("  data/ - Summary tables and gene lists\n")

cat("\n=== Session Info ===\n")
sessionInfo()
