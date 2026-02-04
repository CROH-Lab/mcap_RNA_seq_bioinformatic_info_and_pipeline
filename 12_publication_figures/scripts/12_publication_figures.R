#!/usr/bin/env Rscript
# ==============================================================================
# Publication-Ready Figures: M. capitata and D. trenchii OA Response
# Version 5 - Updated dot plots for calcification pathways
# ==============================================================================

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
    library(ComplexHeatmap)
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

host_summer <- read.csv("/home/darmstrong4/mc_rework/07_deseq2/results/Mcap_Summer_BvsD.csv", row.names = 1)
host_winter <- read.csv("/home/darmstrong4/mc_rework/07_deseq2/results/Mcap_Winter_BvsD.csv", row.names = 1)
sym_summer <- read.csv("/home/darmstrong4/mc_rework/07_deseq2/results/Dtre_Summer_BvsD.csv", row.names = 1)
sym_winter <- read.csv("/home/darmstrong4/mc_rework/07_deseq2/results/Dtre_Winter_BvsD.csv", row.names = 1)

host_summer$gene_id <- rownames(host_summer)
host_winter$gene_id <- rownames(host_winter)
sym_summer$gene_id <- rownames(sym_summer)
sym_winter$gene_id <- rownames(sym_winter)

host_annot <- read.delim("/home/darmstrong4/mc_rework/08_host_deg_annotation/results/all_annotations_full.tsv",
                          stringsAsFactors = FALSE)
sym_annot <- read.delim("/home/darmstrong4/mc_rework/09_symbiont_deg_annotation/results/all_annotations_full.tsv",
                         stringsAsFactors = FALSE)

host_summer_vsd <- readRDS("/home/darmstrong4/mc_rework/07_deseq2/objects/Mcapitata_Summer_vsd.rds")
host_winter_vsd <- readRDS("/home/darmstrong4/mc_rework/07_deseq2/objects/Mcapitata_Winter_vsd.rds")
sym_summer_vsd <- readRDS("/home/darmstrong4/mc_rework/07_deseq2/objects/Dtrenchii_Summer_vsd.rds")
sym_winter_vsd <- readRDS("/home/darmstrong4/mc_rework/07_deseq2/objects/Dtrenchii_Winter_vsd.rds")

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
# UPDATED FIGURE 1-2: Venn Diagrams with Legend + Ridgeline Plots
# Changes: Consistent font, removed category text, added color legend
# ==============================================================================

cat("\n==============================================================================\n")
cat("Figure 1-2 Combined: Venn Diagrams + Sample L2FC Ridgeline Plots\n")
cat("==============================================================================\n\n")

# ------------------------------------------------------------------------------
# Color Palettes
# ------------------------------------------------------------------------------

# Host Summer: Warm palette (for ridgeline and Venn summer)
host_summer_colors <- c(
    "1BS" = "#FDBB84",
    "2BS" = "#EF6548",
    "3BS" = "#D7301F"
)

# Host Winter: Cool palette
host_winter_colors <- c(
    "1BW" = "#9ECAE1",
    "2BW" = "#4292C6",
    "3BW" = "#08519C"
)

# Symbiont Winter: Green palette
sym_winter_colors <- c(
    "1BW" = "#A1D99B",
    "2BW" = "#41AB5D",
    "3BW" = "#006D2C"
)

# Venn colors - match the ridgeline mid-tones
venn_summer_color <- "#EF6548"  # Warm (matches host summer)
venn_winter_color <- "#4292C6"  # Cool (matches host winter)

# ------------------------------------------------------------------------------
# Prepare DEG IDs
# ------------------------------------------------------------------------------

host_summer_ids <- host_summer_degs$gene_id
host_winter_ids <- host_winter_degs$gene_id
host_overlap <- length(intersect(host_summer_ids, host_winter_ids))

sym_summer_ids <- sym_summer_degs$gene_id
sym_winter_ids <- sym_winter_degs$gene_id
sym_overlap <- length(intersect(sym_summer_ids, sym_winter_ids))

cat("DEG counts:\n")
cat("  Host Summer:", length(host_summer_ids), "\n")
cat("  Host Winter:", length(host_winter_ids), "\n")
cat("  Symbiont Winter:", length(sym_winter_ids), "\n")

# ------------------------------------------------------------------------------
# Prepare Ridgeline Data - Per-sample L2FC relative to Ambient mean
# ------------------------------------------------------------------------------

prepare_ridgeline_l2fc <- function(vsd_mat, deg_ids, season = "S") {
    # Filter to DEGs present in VST matrix
    common_genes <- intersect(deg_ids, rownames(vsd_mat))
    cat("  Using", length(common_genes), "DEGs\n")

    # Define sample patterns based on season
    if (season == "S") {
        oa_pattern <- "BS$"
        ambient_pattern <- "DS$"
    } else {
        oa_pattern <- "BW$"
        ambient_pattern <- "DW$"
    }

    # Get OA (B) and Ambient (D) samples
    oa_samples <- grep(oa_pattern, colnames(vsd_mat), value = TRUE)
    ambient_samples <- grep(ambient_pattern, colnames(vsd_mat), value = TRUE)

    cat("    OA samples:", paste(oa_samples, collapse = ", "), "\n")
    cat("    Ambient samples:", paste(ambient_samples, collapse = ", "), "\n")

    # Calculate mean ambient expression per gene
    ambient_mean <- rowMeans(vsd_mat[common_genes, ambient_samples, drop = FALSE])

    # Calculate L2FC for each OA sample relative to ambient mean
    # Since VST is already log2-scale, subtraction gives log2 fold change
    expr_long <- lapply(oa_samples, function(samp) {
        l2fc <- vsd_mat[common_genes, samp] - ambient_mean
        data.frame(
            gene_id = common_genes,
            sample = samp,
            log2FC = l2fc,
            sample_num = gsub("[^0-9]", "", samp),
            sample_label = paste0(gsub("[^0-9]", "", samp), "B", season)
        )
    }) %>% bind_rows()

    return(expr_long)
}

cat("\nPreparing Host Summer L2FC data...\n")
ridge_host_summer <- prepare_ridgeline_l2fc(host_summer_vsd, host_summer_ids, "S")
ridge_host_summer$organism <- "M. capitata"
ridge_host_summer$group <- "Host Summer"

cat("Preparing Host Winter L2FC data...\n")
ridge_host_winter <- prepare_ridgeline_l2fc(host_winter_vsd, host_winter_ids, "W")
ridge_host_winter$organism <- "M. capitata"
ridge_host_winter$group <- "Host Winter"

cat("Preparing Symbiont Winter L2FC data...\n")
ridge_sym_winter <- prepare_ridgeline_l2fc(sym_winter_vsd, sym_winter_ids, "W")
ridge_sym_winter$organism <- "D. trenchii"
ridge_sym_winter$group <- "Symbiont Winter"

# ------------------------------------------------------------------------------
# Create Ridgeline Plots
# ------------------------------------------------------------------------------

# Panel C: Host Summer ridgeline (L2FC)
ridge_host_summer <- ridge_host_summer %>%
    mutate(sample_label = factor(sample_label, levels = c("3BS", "2BS", "1BS")))

p_ridge_summer <- ggplot(ridge_host_summer,
                          aes(x = log2FC, y = sample_label, fill = sample_label)) +
    geom_density_ridges(alpha = 0.8, scale = 1.2, rel_min_height = 0.01) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray60", linewidth = 0.4) +
    scale_fill_manual(values = host_summer_colors) +
    scale_x_continuous(limits = c(-4, 4), breaks = seq(-4, 4, 1)) +
    labs(x = expression(Log[2]~Fold~Change), y = "",
         title = expression(italic("M. capitata")~"- Summer")) +
    theme_bw(base_size = 11) +
    theme(
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, color = "#D7301F"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.title.x = element_text(size = 11),
        panel.grid.minor = element_blank()
    )

# Panel D: Winter combined ridgeline (Host + Symbiont)
ridge_winter_combined <- bind_rows(
    ridge_host_winter %>% mutate(group_label = paste0(sample_label, "\n(Host)")),
    ridge_sym_winter %>% mutate(group_label = paste0(sample_label, "\n(Sym)"))
)

ridge_winter_combined <- ridge_winter_combined %>%
    mutate(
        plot_order = case_when(
            organism == "D. trenchii" & sample_label == "3BW" ~ 6,
            organism == "D. trenchii" & sample_label == "2BW" ~ 5,
            organism == "D. trenchii" & sample_label == "1BW" ~ 4,
            organism == "M. capitata" & sample_label == "3BW" ~ 3,
            organism == "M. capitata" & sample_label == "2BW" ~ 2,
            organism == "M. capitata" & sample_label == "1BW" ~ 1
        ),
        display_label = case_when(
            organism == "D. trenchii" ~ paste0(sample_label, " (Sym)"),
            organism == "M. capitata" ~ paste0(sample_label, " (Host)")
        ),
        color_key = paste(organism, sample_label, sep = "_")
    )

ridge_winter_combined <- ridge_winter_combined %>%
    mutate(display_label = factor(display_label,
                                   levels = c("3BW (Sym)", "2BW (Sym)", "1BW (Sym)",
                                             "3BW (Host)", "2BW (Host)", "1BW (Host)")))

winter_combined_colors <- c(
    "1BW (Host)" = "#9ECAE1",
    "2BW (Host)" = "#4292C6",
    "3BW (Host)" = "#08519C",
    "1BW (Sym)" = "#A1D99B",
    "2BW (Sym)" = "#41AB5D",
    "3BW (Sym)" = "#006D2C"
)

p_ridge_winter <- ggplot(ridge_winter_combined,
                          aes(x = log2FC, y = display_label, fill = display_label)) +
    geom_density_ridges(alpha = 0.8, scale = 1.2, rel_min_height = 0.01) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray60", linewidth = 0.4) +
    scale_fill_manual(values = winter_combined_colors) +
    scale_x_continuous(limits = c(-4, 4), breaks = seq(-4, 4, 1)) +
    labs(x = expression(Log[2]~Fold~Change), y = "",
         title = "Winter") +
    theme_bw(base_size = 11) +
    theme(
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 11),
        panel.grid.minor = element_blank()
    )

# ------------------------------------------------------------------------------
# Create Combined Figure with Venn Legend
# Layout:
#   Row 1: [A] Host Venn    [B] Symbiont Venn    [Legend]
#   Row 2: [C] Host Summer Ridge  [D] Winter Combined Ridge
# ------------------------------------------------------------------------------

cat("\nCreating combined figure with Venn legend...\n")

pdf("figures/Fig1_combined_venn_ridgeline.pdf", width = 12, height = 10)

grid.newpage()

# Layout: 2 rows, 3 columns (extra column for legend)
pushViewport(viewport(layout = grid.layout(
    nrow = 2, ncol = 3,
    heights = unit(c(0.45, 0.55), "npc"),
    widths = unit(c(0.42, 0.42, 0.16), "npc")
)))

# ------------------------------------------------------------------------------
# Panel A: Host Venn (row 1, col 1) - NO category labels
# ------------------------------------------------------------------------------
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
pushViewport(viewport(width = 0.9, height = 0.9))

grid.rect(gp = gpar(col = NA, fill = "white"))
venn_host <- draw.pairwise.venn(
    area1 = length(host_summer_ids), area2 = length(host_winter_ids),
    cross.area = host_overlap, 
    category = c("", ""),  # Remove category labels
    fill = c(venn_summer_color, venn_winter_color), alpha = 0.6,
    col = c(venn_summer_color, venn_winter_color), lwd = 2, cex = 1.8,
    fontface = "bold", fontfamily = "sans",  # Match ggplot font
    cat.cex = 0, cat.default.pos = "text",  # Hide category text
    margin = 0.05, ind = FALSE
)
grid.draw(venn_host)

upViewport()
# Panel label and title
grid.text("A", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "sans"))
grid.text(expression(italic("M. capitata")), x = unit(0.5, "npc"), y = unit(0.97, "npc"),
          just = "center", gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "sans"))
upViewport()

# ------------------------------------------------------------------------------
# Panel B: Symbiont Venn (row 1, col 2) - NO category labels
# ------------------------------------------------------------------------------
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
pushViewport(viewport(width = 0.9, height = 0.9))

grid.rect(gp = gpar(col = NA, fill = "white"))
venn_sym <- draw.pairwise.venn(
    area1 = length(sym_summer_ids), area2 = length(sym_winter_ids),
    cross.area = sym_overlap, 
    category = c("", ""),  # Remove category labels
    fill = c(venn_summer_color, venn_winter_color), alpha = 0.6,
    col = c(venn_summer_color, venn_winter_color), lwd = 2, cex = 1.8,
    fontface = "bold", fontfamily = "sans",  # Match ggplot font
    cat.cex = 0, cat.default.pos = "text",  # Hide category text
    margin = 0.05, ind = FALSE
)
grid.draw(venn_sym)

upViewport()
grid.text("B", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "sans"))
grid.text(expression(italic("D. trenchii")), x = unit(0.5, "npc"), y = unit(0.97, "npc"),
          just = "center", gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "sans"))
upViewport()

# ------------------------------------------------------------------------------
# Legend Panel (row 1, col 3)
# ------------------------------------------------------------------------------
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
pushViewport(viewport(width = 0.8, height = 0.5, y = 0.5))

# Draw legend box
grid.rect(gp = gpar(col = "gray70", fill = "white", lwd = 1))

# Legend title
grid.text("Season", x = unit(0.5, "npc"), y = unit(0.85, "npc"),
          gp = gpar(fontsize = 12, fontface = "bold", fontfamily = "sans"))

# Summer circle and label
grid.circle(x = unit(0.25, "npc"), y = unit(0.55, "npc"), r = unit(0.08, "npc"),
            gp = gpar(fill = venn_summer_color, col = venn_summer_color, alpha = 0.6))
grid.text("Summer", x = unit(0.55, "npc"), y = unit(0.55, "npc"),
          just = "left", gp = gpar(fontsize = 11, fontfamily = "sans"))

# Winter circle and label
grid.circle(x = unit(0.25, "npc"), y = unit(0.25, "npc"), r = unit(0.08, "npc"),
            gp = gpar(fill = venn_winter_color, col = venn_winter_color, alpha = 0.6))
grid.text("Winter", x = unit(0.55, "npc"), y = unit(0.25, "npc"),
          just = "left", gp = gpar(fontsize = 11, fontfamily = "sans"))

upViewport(2)

# ------------------------------------------------------------------------------
# Panel C: Host Summer Ridgeline (row 2, col 1)
# ------------------------------------------------------------------------------
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
pushViewport(viewport(width = 0.95, height = 0.92))

p_ridge_summer_grob <- ggplotGrob(p_ridge_summer)
grid.draw(p_ridge_summer_grob)

upViewport()
grid.text("C", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "sans"))
upViewport()

# ------------------------------------------------------------------------------
# Panel D: Winter Combined Ridgeline (row 2, col 2-3)
# ------------------------------------------------------------------------------
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2:3))
pushViewport(viewport(width = 0.95, height = 0.92))

p_ridge_winter_grob <- ggplotGrob(p_ridge_winter)
grid.draw(p_ridge_winter_grob)

upViewport()
grid.text("D", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "sans"))
upViewport()

dev.off()

# ------------------------------------------------------------------------------
# PNG version
# ------------------------------------------------------------------------------

png("figures/Fig1_combined_venn_ridgeline.png", width = 12, height = 10, units = "in", res = 300)

grid.newpage()

pushViewport(viewport(layout = grid.layout(
    nrow = 2, ncol = 3,
    heights = unit(c(0.45, 0.55), "npc"),
    widths = unit(c(0.42, 0.42, 0.16), "npc")
)))

# Panel A
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
pushViewport(viewport(width = 0.9, height = 0.9))
grid.rect(gp = gpar(col = NA, fill = "white"))
venn_host <- draw.pairwise.venn(
    area1 = length(host_summer_ids), area2 = length(host_winter_ids),
    cross.area = host_overlap, category = c("", ""),
    fill = c(venn_summer_color, venn_winter_color), alpha = 0.6,
    col = c(venn_summer_color, venn_winter_color), lwd = 2, cex = 1.8,
    fontface = "bold", fontfamily = "sans", cat.cex = 0,
    margin = 0.05, ind = FALSE
)
grid.draw(venn_host)
upViewport()
grid.text("A", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "sans"))
grid.text(expression(italic("M. capitata")), x = unit(0.5, "npc"), y = unit(0.97, "npc"),
          just = "center", gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "sans"))
upViewport()

# Panel B
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
pushViewport(viewport(width = 0.9, height = 0.9))
grid.rect(gp = gpar(col = NA, fill = "white"))
venn_sym <- draw.pairwise.venn(
    area1 = length(sym_summer_ids), area2 = length(sym_winter_ids),
    cross.area = sym_overlap, category = c("", ""),
    fill = c(venn_summer_color, venn_winter_color), alpha = 0.6,
    col = c(venn_summer_color, venn_winter_color), lwd = 2, cex = 1.8,
    fontface = "bold", fontfamily = "sans", cat.cex = 0,
    margin = 0.05, ind = FALSE
)
grid.draw(venn_sym)
upViewport()
grid.text("B", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "sans"))
grid.text(expression(italic("D. trenchii")), x = unit(0.5, "npc"), y = unit(0.97, "npc"),
          just = "center", gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "sans"))
upViewport()

# Legend Panel
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
pushViewport(viewport(width = 0.8, height = 0.5, y = 0.5))
grid.rect(gp = gpar(col = "gray70", fill = "white", lwd = 1))
grid.text("Season", x = unit(0.5, "npc"), y = unit(0.85, "npc"),
          gp = gpar(fontsize = 12, fontface = "bold", fontfamily = "sans"))
grid.circle(x = unit(0.25, "npc"), y = unit(0.55, "npc"), r = unit(0.08, "npc"),
            gp = gpar(fill = venn_summer_color, col = venn_summer_color, alpha = 0.6))
grid.text("Summer", x = unit(0.55, "npc"), y = unit(0.55, "npc"),
          just = "left", gp = gpar(fontsize = 11, fontfamily = "sans"))
grid.circle(x = unit(0.25, "npc"), y = unit(0.25, "npc"), r = unit(0.08, "npc"),
            gp = gpar(fill = venn_winter_color, col = venn_winter_color, alpha = 0.6))
grid.text("Winter", x = unit(0.55, "npc"), y = unit(0.25, "npc"),
          just = "left", gp = gpar(fontsize = 11, fontfamily = "sans"))
upViewport(2)

# Panel C
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
pushViewport(viewport(width = 0.95, height = 0.92))
grid.draw(ggplotGrob(p_ridge_summer))
upViewport()
grid.text("C", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "sans"))
upViewport()

# Panel D
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2:3))
pushViewport(viewport(width = 0.95, height = 0.92))
grid.draw(ggplotGrob(p_ridge_winter))
upViewport()
grid.text("D", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "sans"))
upViewport()

dev.off()

cat("\n✓ Saved: Fig1_combined_venn_ridgeline.pdf/png (12x10 inches)\n")

# Save individual ridgeline plots
ggsave("figures/Fig1C_ridgeline_host_summer.pdf", p_ridge_summer,
       width = 6, height = 4, dpi = 300)
ggsave("figures/Fig1D_ridgeline_winter_combined.pdf", p_ridge_winter,
       width = 6, height = 5, dpi = 300)

cat("✓ Saved individual ridgeline plots\n")

cat("\n=== Figure 1-2 Combined Summary ===\n")
cat("Layout: 2 rows × 3 columns (with legend)\n")
cat("  Row 1: A) Host Venn, B) Symbiont Venn, Season Legend\n")
cat("  Row 2: C) Host Summer L2FC, D) Winter Combined L2FC\n")
cat("\nChanges from previous version:\n")
cat("  - Removed 'Summer'/'Winter' text from Venn diagrams\n")
cat("  - Added Season legend panel\n")
cat("  - Matched fonts (sans family)\n")


# ==============================================================================
# UPDATED FIGURE 3: Heatmaps - Remove treatment legend from A and B
# ==============================================================================

cat("\n==============================================================================\n")
cat("Figure 3: Heatmaps (Updated - Legend only on Panel C)\n")
cat("==============================================================================\n\n")

# Color palettes for each panel
heatmap_colors_host_summer <- colorRampPalette(c("#FFF5EB", "#FDD49E", "#FDBB84",
                                                  "#FC8D59", "#EF6548", "#D7301F",
                                                  "#990000"))(100)

heatmap_colors_host_winter <- colorRampPalette(c("#F7FBFF", "#DEEBF7", "#C6DBEF",
                                                  "#9ECAE1", "#6BAED6", "#3182BD",
                                                  "#08519C"))(100)

heatmap_colors_sym_winter <- colorRampPalette(c("#F7FCF5", "#E5F5E0", "#C7E9C0",
                                                 "#A1D99B", "#74C476", "#31A354",
                                                 "#006D2C"))(100)

# Modified function with option to show/hide annotation legend
create_deg_heatmap_panel <- function(vsd_mat, deseq_results, color_palette, 
                                      max_genes = 250, show_legend = TRUE) {

    all_degs <- deseq_results %>%
        filter(!is.na(padj) & padj < 0.05) %>%
        arrange(desc(abs(log2FoldChange)))

    if (nrow(all_degs) < 5) {
        cat("  Warning: Fewer than 5 DEGs, skipping\n")
        return(NULL)
    }

    if (nrow(all_degs) > max_genes) {
        cat("  Note: Limiting to top", max_genes, "of", nrow(all_degs), "DEGs\n")
        all_degs <- head(all_degs, max_genes)
    } else {
        cat("  Using all", nrow(all_degs), "DEGs\n")
    }

    common_genes <- intersect(all_degs$gene_id, rownames(vsd_mat))
    if (length(common_genes) < 5) return(NULL)

    heatmap_mat <- vsd_mat[common_genes, , drop = FALSE]
    heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

    sample_names <- colnames(heatmap_mat)
    treatments <- ifelse(grepl("B", sample_names), "OA", "Ambient")
    anno_col <- data.frame(Treatment = treatments, row.names = sample_names)
    anno_colors <- list(Treatment = c("OA" = "#E69F00", "Ambient" = "#56B4E9"))

    # Create heatmap with optional legend
    ht <- pheatmap::pheatmap(heatmap_mat_scaled,
                              color = color_palette,
                              cluster_rows = TRUE, cluster_cols = TRUE,
                              show_rownames = FALSE, show_colnames = TRUE,
                              annotation_col = anno_col,
                              annotation_colors = anno_colors,
                              annotation_names_col = FALSE,
                              annotation_legend = show_legend,  # Control legend visibility
                              fontsize = 10,
                              border_color = NA,
                              silent = TRUE)

    return(list(heatmap = ht, n_genes = length(common_genes)))
}

# Generate heatmaps - NO legend for Host panels, YES legend for Symbiont
cat("Creating Host Summer heatmap (no legend)...\n")
ht_host_summer <- create_deg_heatmap_panel(host_summer_vsd, host_summer,
                                            heatmap_colors_host_summer, 
                                            max_genes = 250, show_legend = FALSE)

cat("Creating Host Winter heatmap (no legend)...\n")
ht_host_winter <- create_deg_heatmap_panel(host_winter_vsd, host_winter,
                                            heatmap_colors_host_winter, 
                                            max_genes = 250, show_legend = FALSE)

cat("Creating Symbiont Winter heatmap (with legend)...\n")
ht_sym_winter <- create_deg_heatmap_panel(sym_winter_vsd, sym_winter,
                                           heatmap_colors_sym_winter, 
                                           max_genes = 250, show_legend = TRUE)

cat("Note: Symbiont Summer has only", nrow(sym_summer_degs), "DEGs - excluding from combined figure\n")

# Save individual heatmaps
if (!is.null(ht_host_summer)) {
    pdf("figures/Fig3A_heatmap_host_summer.pdf", width = 8, height = 10)
    grid.draw(ht_host_summer$heatmap$gtable)
    dev.off()
    cat("  Saved: Fig3A_heatmap_host_summer.pdf with", ht_host_summer$n_genes, "DEGs (no legend)\n")
}

if (!is.null(ht_host_winter)) {
    pdf("figures/Fig3B_heatmap_host_winter.pdf", width = 8, height = 10)
    grid.draw(ht_host_winter$heatmap$gtable)
    dev.off()
    cat("  Saved: Fig3B_heatmap_host_winter.pdf with", ht_host_winter$n_genes, "DEGs (no legend)\n")
}

if (!is.null(ht_sym_winter)) {
    pdf("figures/Fig3C_heatmap_symbiont_winter.pdf", width = 8, height = 10)
    grid.draw(ht_sym_winter$heatmap$gtable)
    dev.off()
    cat("  Saved: Fig3C_heatmap_symbiont_winter.pdf with", ht_sym_winter$n_genes, "DEGs (with legend)\n")
}

# ==============================================================================
# Combined Heatmap Figure
# ==============================================================================

cat("\nCreating combined heatmap figure...\n")

gt_host_summer <- if (!is.null(ht_host_summer)) ht_host_summer$heatmap$gtable else NULL
gt_host_winter <- if (!is.null(ht_host_winter)) ht_host_winter$heatmap$gtable else NULL
gt_sym_winter <- if (!is.null(ht_sym_winter)) ht_sym_winter$heatmap$gtable else NULL

pdf("figures/Fig3_heatmap_combined.pdf", width = 18, height = 10)

grid.newpage()
vp_layout <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                                            heights = unit(c(0.05, 0.95), "npc"),
                                            widths = unit(c(1, 1, 1), "null")))
pushViewport(vp_layout)

# Row 1: Titles
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.text("A", x = 0.05, y = 0.5, just = "left", gp = gpar(fontsize = 18, fontface = "bold"))
grid.text(expression(italic("M. capitata")~"- Summer"), x = 0.5, y = 0.5,
          gp = gpar(fontsize = 12, fontface = "bold"))
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.text("B", x = 0.05, y = 0.5, just = "left", gp = gpar(fontsize = 18, fontface = "bold"))
grid.text(expression(italic("M. capitata")~"- Winter"), x = 0.5, y = 0.5,
          gp = gpar(fontsize = 12, fontface = "bold"))
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
grid.text("C", x = 0.05, y = 0.5, just = "left", gp = gpar(fontsize = 18, fontface = "bold"))
grid.text(expression(italic("D. trenchii")~"- Winter"), x = 0.5, y = 0.5,
          gp = gpar(fontsize = 12, fontface = "bold"))
upViewport()

# Row 2: Heatmaps
if (!is.null(gt_host_summer)) {
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
    pushViewport(viewport(width = 0.95, height = 0.95))
    grid.draw(gt_host_summer)
    upViewport(2)
}

if (!is.null(gt_host_winter)) {
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
    pushViewport(viewport(width = 0.95, height = 0.95))
    grid.draw(gt_host_winter)
    upViewport(2)
}

if (!is.null(gt_sym_winter)) {
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
    pushViewport(viewport(width = 0.95, height = 0.95))
    grid.draw(gt_sym_winter)
    upViewport(2)
}

dev.off()

# PNG version
png("figures/Fig3_heatmap_combined.png", width = 18, height = 10, units = "in", res = 300)

grid.newpage()
vp_layout <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                                            heights = unit(c(0.05, 0.95), "npc"),
                                            widths = unit(c(1, 1, 1), "null")))
pushViewport(vp_layout)

# Titles
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.text("A", x = 0.05, y = 0.5, just = "left", gp = gpar(fontsize = 18, fontface = "bold"))
grid.text(expression(italic("M. capitata")~"- Summer"), x = 0.5, y = 0.5,
          gp = gpar(fontsize = 12, fontface = "bold"))
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.text("B", x = 0.05, y = 0.5, just = "left", gp = gpar(fontsize = 18, fontface = "bold"))
grid.text(expression(italic("M. capitata")~"- Winter"), x = 0.5, y = 0.5,
          gp = gpar(fontsize = 12, fontface = "bold"))
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
grid.text("C", x = 0.05, y = 0.5, just = "left", gp = gpar(fontsize = 18, fontface = "bold"))
grid.text(expression(italic("D. trenchii")~"- Winter"), x = 0.5, y = 0.5,
          gp = gpar(fontsize = 12, fontface = "bold"))
upViewport()

# Heatmaps
if (!is.null(gt_host_summer)) {
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
    pushViewport(viewport(width = 0.95, height = 0.95))
    grid.draw(gt_host_summer)
    upViewport(2)
}

if (!is.null(gt_host_winter)) {
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
    pushViewport(viewport(width = 0.95, height = 0.95))
    grid.draw(gt_host_winter)
    upViewport(2)
}

if (!is.null(gt_sym_winter)) {
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
    pushViewport(viewport(width = 0.95, height = 0.95))
    grid.draw(gt_sym_winter)
    upViewport(2)
}

dev.off()

cat("✓ Saved: Fig3_heatmap_combined.pdf/png (18x10 inches)\n")

cat("\n=== Heatmap Summary ===\n")
cat("Panel A (Host Summer):", ifelse(!is.null(ht_host_summer),
    paste(ht_host_summer$n_genes, "DEGs, warm colors, NO legend"), "skipped"), "\n")
cat("Panel B (Host Winter):", ifelse(!is.null(ht_host_winter),
    paste(ht_host_winter$n_genes, "DEGs, cool colors, NO legend"), "skipped"), "\n")
cat("Panel C (Symbiont Winter):", ifelse(!is.null(ht_sym_winter),
    paste(ht_sym_winter$n_genes, "DEGs, green colors, WITH legend"), "skipped"), "\n")
cat("Note: Treatment legend only shown on Panel C (Symbiont Winter)\n")

# ==============================================================================
# FIGURE 4: Sankey Plot - Final with Priority Categories & Colored Flows
# Orange (host) and Green (symbiont) with flows matching source organism
# ==============================================================================

cat("\n==============================================================================\n")
cat("Figure 4: Sankey Plots - Complete Calcification Pathways\n")
cat("==============================================================================\n\n")

suppressPackageStartupMessages({
    library(networkD3)
    library(htmlwidgets)
})

# Configuration
CALC_KEYWORDS <- c(
    "calcium", "carbonate", "carbonic",
    "ion transport", "ion homeostasis", "metal ion",
    "proton", "ATPase", "pH",
    "solute", "biomineralization", "ossification",
    "cation", "anion", "sodium", "potassium",
    "phosphorylation", "phospholipid", "oxidative"
)

P_ADJ_CUTOFF <- 0.05

# How many top parent categories to show
SEASON_CONFIG <- list(
    summer = list(max_parent_categories = 10),
    winter = list(max_parent_categories = 8)  # Slightly more for new categories
)

# Load GO_MWU results
load_gomwu_results <- function(dir, organism, season, divisions = c("BP", "MF", "CC")) {
    results <- list()
    for (div in divisions) {
        file_prefix <- ifelse(organism == "host",
                             paste0(season, "_", div),
                             paste0("symbiont_", season, "_", div))
        file_path <- file.path(dir, paste0(file_prefix, "_MWU_results.csv"))

        if (file.exists(file_path)) {
            df <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE,
                           sep = "", quote = "\"", comment.char = "")
            names(df) <- gsub('^"', '', names(df))
            names(df) <- gsub('"$', '', names(df))
            names(df) <- trimws(names(df))

            if (!"p.adj" %in% names(df)) {
                if ("padj" %in% names(df)) {
                    names(df)[names(df) == "padj"] <- "p.adj"
                }
            }

            df$organism <- organism
            df$division <- div
            df$season <- season
            results[[paste(organism, div, sep="_")]] <- df
        }
    }
    return(bind_rows(results))
}

# Function to create Sankey for a given season
create_season_sankey <- function(season_name) {
    cat("\n--- Processing", toupper(season_name), "---\n")

    config <- SEASON_CONFIG[[season_name]]

    host_gomwu <- load_gomwu_results("../10_GO_MWU/output", "host", season_name)
    symbiont_gomwu <- load_gomwu_results("../11_symbiont_GO_MWU/output", "symbiont", season_name)
    all_gomwu <- bind_rows(host_gomwu, symbiont_gomwu)

    cat("  Loaded GO terms - Host:", nrow(host_gomwu), "| Symbiont:", nrow(symbiont_gomwu), "\n")

    # Filter for calcification terms
    calc_data <- all_gomwu %>%
        filter(p.adj < P_ADJ_CUTOFF) %>%
        filter(str_detect(name, regex(paste(CALC_KEYWORDS, collapse = "|"), ignore_case = TRUE)))

    cat("  Calcification terms at p.adj < 0.05:", nrow(calc_data), "\n")

    if (nrow(calc_data) == 0) {
        cat("  Warning: No calcification terms found\n")
        return(NULL)
    }

    # Assign parent categories - EXPANDED with priority categories
    assign_parent_category <- function(go_name) {
        go_name_lower <- tolower(go_name)
        case_when(
            # Priority 3: Oxidative Phosphorylation (check first - most specific)
            str_detect(go_name_lower, "oxidative phosphorylation") ~ "Oxidative Phosphorylation",

            # Priority 1: Specific Ion Transport (before general ion transport)
            str_detect(go_name_lower, "sodium.*transport|potassium.*transport") ~ "Na/K Transport",

            # Priority 4: Phospholipid/Membrane Transport
            str_detect(go_name_lower, "phospholipid.*transport|organophosphate.*transport") ~ "Phospholipid Transport",

            # Priority 2: Protein Phosphorylation (broad kinase activity)
            str_detect(go_name_lower, "protein phosphorylation|phosphotransferase.*alcohol|autophosphorylation") ~ "Protein Phosphorylation",

            # Original categories (more specific matches first)
            str_detect(go_name_lower, "calcium") ~ "Calcium Homeostasis",
            str_detect(go_name_lower, "carbonate|carbonic") ~ "Carbon/Carbonate",
            str_detect(go_name_lower, "proton.*transport|h\\+.*transport") ~ "Proton Transport",
            str_detect(go_name_lower, "metal ion") ~ "Metal Ion Homeostasis",
            str_detect(go_name_lower, "atpase") ~ "ATPase Activity",

            # Ion homeostasis & regulation (catch regulation terms)
            str_detect(go_name_lower, "ion homeostasis|regulation.*ion transport") ~ "Ion Homeostasis/Regulation",

            # General ion transport (last, catches remaining)
            str_detect(go_name_lower, "cation|anion|ion transport") ~ "Ion Transport",

            TRUE ~ NA_character_
        )
    }

    calc_data <- calc_data %>%
        mutate(
            parent_category = assign_parent_category(name),
            source_node = paste(organism, division, sep = "-")
        )

    # Diagnostic output
    cat("\n  === PARENT CATEGORY BREAKDOWN ===\n")

    parent_breakdown <- calc_data %>%
        filter(!is.na(parent_category)) %>%
        group_by(organism, parent_category) %>%
        summarise(n_terms = n(), .groups = "drop") %>%
        arrange(organism, desc(n_terms))

    cat("\n  HOST:\n")
    host_breakdown <- parent_breakdown %>% filter(organism == "host")
    if (nrow(host_breakdown) > 0) print(host_breakdown) else cat("    No host terms\n")

    cat("\n  SYMBIONT:\n")
    symbiont_breakdown <- parent_breakdown %>% filter(organism == "symbiont")
    if (nrow(symbiont_breakdown) > 0) print(symbiont_breakdown) else cat("    No symbiont terms\n")

    # Show excluded terms
    excluded_terms <- calc_data %>%
        filter(is.na(parent_category)) %>%
        select(organism, division, name, p.adj) %>%
        arrange(organism, p.adj)

    if (nrow(excluded_terms) > 0) {
        cat("\n  === EXCLUDED TERMS ===\n")
        cat("  Total excluded:", nrow(excluded_terms), "\n")
        excluded_file <- paste0("data/excluded_calcification_terms_", season_name, ".csv")
        write.csv(excluded_terms, excluded_file, row.names = FALSE)
        cat("  ✓ Saved to:", excluded_file, "\n")
    }

    n_before <- nrow(calc_data)
    calc_data <- calc_data %>% filter(!is.na(parent_category))
    n_after <- nrow(calc_data)
    cat("\n  Terms: ", n_before, " → ", n_after, " (excluded ", n_before - n_after, ")\n", sep = "")

    # Hierarchical filtering
    parent_ranking <- calc_data %>%
        group_by(parent_category) %>%
        summarise(total_significance = sum(-log10(p.adj + 1e-10)), n_terms = n(), .groups = "drop") %>%
        arrange(desc(total_significance))

    cat("\n  Parent category ranking:\n")
    print(parent_ranking)

    top_parents <- head(parent_ranking, config$max_parent_categories)$parent_category
    cat("\n  Selected top", config$max_parent_categories, "parents\n")

    calc_data <- calc_data %>%
        filter(parent_category %in% top_parents) %>%
        mutate(go_term_display = str_trunc(name, 60))

    cat("  Final GO terms:", nrow(calc_data), "\n")

    # Build nodes
    source_nodes <- calc_data %>%
        distinct(source_node, organism) %>%
        mutate(group = organism)

    parent_nodes <- calc_data %>%
        distinct(parent_category) %>%
        mutate(group = "parent", organism = "parent")
    names(parent_nodes)[1] <- "source_node"

    go_nodes <- calc_data %>%
        distinct(go_term_display, parent_category, organism) %>%
        mutate(group = organism)
    names(go_nodes)[1] <- "source_node"
    go_nodes <- go_nodes %>% select(source_node, group, organism)

    nodes <- bind_rows(
        source_nodes %>% select(source_node, group, organism),
        parent_nodes,
        go_nodes
    ) %>%
        distinct(source_node, .keep_all = TRUE) %>%
        mutate(node_id = row_number() - 1)

    # Build links WITH organism tracking for coloring
    link1 <- calc_data %>%
        group_by(source_node, parent_category, organism) %>%
        summarise(value = sum(-log10(p.adj + 1e-10)), .groups = "drop") %>%
        left_join(nodes %>% select(source_node, source_id = node_id), by = "source_node") %>%
        left_join(nodes %>% select(source_node, target_id = node_id), by = c("parent_category" = "source_node"))

    link2 <- calc_data %>%
        mutate(value = -log10(p.adj + 1e-10)) %>%
        left_join(nodes %>% select(source_node, source_id = node_id), by = c("parent_category" = "source_node")) %>%
        left_join(nodes %>% select(source_node, target_id = node_id), by = c("go_term_display" = "source_node"))

    links <- bind_rows(
        link1 %>% select(source = source_id, target = target_id, value, organism),
        link2 %>% select(source = source_id, target = target_id, value, organism)
    )

    cat("  Sankey: Nodes =", nrow(nodes), "| Links =", nrow(links), "\n")

    # Prepare for Sankey
    nodes_renamed <- nodes
    colnames(nodes_renamed)[colnames(nodes_renamed) == "source_node"] <- "name"
    nodes_for_sankey <- as.data.frame(nodes_renamed)
    links_for_sankey <- as.data.frame(links)

    # CRITICAL: Add link group for coloring flows by source organism
    links_for_sankey$group <- links_for_sankey$organism

    # Colors: Orange (host), Green (symbiont), Gray (parent)
    color_scale <- 'd3.scaleOrdinal()
        .domain(["host", "symbiont", "parent"])
        .range(["#E69F00", "#2ECC71", "#95A5A6"])'

    sankey <- sankeyNetwork(
        Links = links_for_sankey,
        Nodes = nodes_for_sankey,
        Source = "source",
        Target = "target",
        Value = "value",
        NodeID = "name",
        NodeGroup = "organism",
        LinkGroup = "group",  # Color links by source organism!
        fontSize = 11,
        nodeWidth = 25,
        nodePadding = 8,
        iterations = 100,
        sinksRight = FALSE,
        colourScale = color_scale,
        fontFamily = "Arial"
    )

    html_file <- paste0("figures/Fig4_calcification_sankey_", season_name, ".html")
    saveWidget(sankey, html_file, selfcontained = TRUE)
    cat("  ✓ Saved:", html_file, "\n")

    summary_calcs <- calc_data %>%
        group_by(organism, division, parent_category) %>%
        summarise(n_terms = n(), mean_padj = mean(p.adj), .groups = "drop") %>%
        arrange(organism, parent_category)

    write.csv(summary_calcs, paste0("data/calcification_summary_", season_name, ".csv"),
              row.names = FALSE)

    return(list(summary = summary_calcs, parent_ranking = parent_ranking))
}

# Create Sankies for both seasons
summer_result <- create_season_sankey("summer")
winter_result <- create_season_sankey("winter")

# ==============================================================================
# COMBINED SUMMER + WINTER SANKEY
# ==============================================================================

cat("\n==============================================================================\n")
cat("Creating Combined Summer + Winter Sankey\n")
cat("NOTE: Excludes Protein Phosphorylation, Na/K Transport, Phospholipid Transport\n")
cat("==============================================================================\n\n")

# Load both seasons together
host_summer_go <- load_gomwu_results("../10_GO_MWU/output", "host", "summer")
host_winter_go <- load_gomwu_results("../10_GO_MWU/output", "host", "winter")
symbiont_summer_go <- load_gomwu_results("../11_symbiont_GO_MWU/output", "symbiont", "summer")
symbiont_winter_go <- load_gomwu_results("../11_symbiont_GO_MWU/output", "symbiont", "winter")

all_combined <- bind_rows(host_summer_go, host_winter_go, symbiont_summer_go, symbiont_winter_go)

cat("Loaded combined data - Total GO terms:", nrow(all_combined), "\n")

# Filter for calcification
calc_combined <- all_combined %>%
    filter(p.adj < P_ADJ_CUTOFF) %>%
    filter(str_detect(name, regex(paste(CALC_KEYWORDS, collapse = "|"), ignore_case = TRUE)))

cat("Calcification terms at p.adj < 0.05:", nrow(calc_combined), "\n")

# Assign parent categories (same function)
assign_parent_category <- function(go_name) {
    go_name_lower <- tolower(go_name)
    case_when(
        str_detect(go_name_lower, "oxidative phosphorylation") ~ "Oxidative Phosphorylation",
        str_detect(go_name_lower, "sodium.*transport|potassium.*transport") ~ "Na/K Transport",
        str_detect(go_name_lower, "phospholipid.*transport|organophosphate.*transport") ~ "Phospholipid Transport",
        str_detect(go_name_lower, "protein phosphorylation|phosphotransferase.*alcohol|autophosphorylation") ~ "Protein Phosphorylation",
        str_detect(go_name_lower, "calcium") ~ "Calcium Homeostasis",
        str_detect(go_name_lower, "carbonate|carbonic") ~ "Carbon/Carbonate",
        str_detect(go_name_lower, "proton.*transport|h\\+.*transport") ~ "Proton Transport",
        str_detect(go_name_lower, "metal ion") ~ "Metal Ion Homeostasis",
        str_detect(go_name_lower, "atpase") ~ "ATPase Activity",
        str_detect(go_name_lower, "ion homeostasis|regulation.*ion transport") ~ "Ion Homeostasis/Regulation",
        str_detect(go_name_lower, "cation|anion|ion transport") ~ "Ion Transport",
        TRUE ~ NA_character_
    )
}

calc_combined <- calc_combined %>%
    mutate(
        parent_category = assign_parent_category(name),
        season_prefix = ifelse(season == "summer", "S", "W"),
        source_node = paste(season_prefix, organism, division, sep = "-"),
        organism_season = paste(organism, season, sep = "_")  # For 4-color scheme
    ) %>%
    filter(!is.na(parent_category)) %>%
    # EXCLUDE specific categories for combined plot
    filter(!parent_category %in% c("Protein Phosphorylation", "Na/K Transport", "Phospholipid Transport"))

cat("Terms after filtering (with category exclusions):", nrow(calc_combined), "\n")

# Hierarchical filtering - top 10 parent categories overall
parent_ranking_combined <- calc_combined %>%
    group_by(parent_category) %>%
    summarise(total_significance = sum(-log10(p.adj + 1e-10)), n_terms = n(), .groups = "drop") %>%
    arrange(desc(total_significance))

cat("\nCombined parent category ranking:\n")
print(parent_ranking_combined)

top_parents_combined <- head(parent_ranking_combined, 10)$parent_category

calc_combined <- calc_combined %>%
    filter(parent_category %in% top_parents_combined) %>%
    mutate(go_term_display = str_trunc(name, 60))

cat("\nFinal combined terms:", nrow(calc_combined), "\n")
cat("Breakdown:\n")
print(table(calc_combined$season, calc_combined$organism))

# Build nodes
source_nodes_combined <- calc_combined %>%
    distinct(source_node, organism_season) %>%
    mutate(group = organism_season)

parent_nodes_combined <- calc_combined %>%
    distinct(parent_category) %>%
    mutate(group = "parent", organism_season = "parent")
names(parent_nodes_combined)[1] <- "source_node"

go_nodes_combined <- calc_combined %>%
    distinct(go_term_display, parent_category, organism_season) %>%
    mutate(group = organism_season)
names(go_nodes_combined)[1] <- "source_node"
go_nodes_combined <- go_nodes_combined %>% select(source_node, group, organism_season)

nodes_combined <- bind_rows(
    source_nodes_combined %>% select(source_node, group, organism_season),
    parent_nodes_combined,
    go_nodes_combined
) %>%
    distinct(source_node, .keep_all = TRUE) %>%
    mutate(node_id = row_number() - 1)

# Build links
link1_combined <- calc_combined %>%
    group_by(source_node, parent_category, organism_season) %>%
    summarise(value = sum(-log10(p.adj + 1e-10)), .groups = "drop") %>%
    left_join(nodes_combined %>% select(source_node, source_id = node_id), by = "source_node") %>%
    left_join(nodes_combined %>% select(source_node, target_id = node_id), by = c("parent_category" = "source_node"))

link2_combined <- calc_combined %>%
    mutate(value = -log10(p.adj + 1e-10)) %>%
    left_join(nodes_combined %>% select(source_node, source_id = node_id), by = c("parent_category" = "source_node")) %>%
    left_join(nodes_combined %>% select(source_node, target_id = node_id), by = c("go_term_display" = "source_node"))

links_combined <- bind_rows(
    link1_combined %>% select(source = source_id, target = target_id, value, organism_season),
    link2_combined %>% select(source = source_id, target = target_id, value, organism_season)
)

cat("Combined Sankey: Nodes =", nrow(nodes_combined), "| Links =", nrow(links_combined), "\n")

# Prepare for Sankey
nodes_combined_renamed <- nodes_combined
colnames(nodes_combined_renamed)[colnames(nodes_combined_renamed) == "source_node"] <- "name"
nodes_combined_for_sankey <- as.data.frame(nodes_combined_renamed)
links_combined_for_sankey <- as.data.frame(links_combined)
links_combined_for_sankey$group <- links_combined_for_sankey$organism_season

# 4-color scheme:
# #ff671f = Host Summer, #ffb81c = Host Winter
# #006341 = Symbiont Summer, #8fe2b0 = Symbiont Winter
color_scale_combined <- 'd3.scaleOrdinal()
    .domain(["host_summer", "host_winter", "symbiont_summer", "symbiont_winter", "parent"])
    .range(["#ff671f", "#ffb81c", "#006341", "#8fe2b0", "#95A5A6"])'

sankey_combined <- sankeyNetwork(
    Links = links_combined_for_sankey,
    Nodes = nodes_combined_for_sankey,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    NodeGroup = "organism_season",
    LinkGroup = "group",
    fontSize = 10,
    nodeWidth = 25,
    nodePadding = 6,
    iterations = 100,
    sinksRight = FALSE,
    colourScale = color_scale_combined,
    fontFamily = "Arial"
)

html_combined <- "figures/Fig4_calcification_sankey_combined.html"
saveWidget(sankey_combined, html_combined, selfcontained = TRUE)
cat("\n✓ Saved combined plot:", html_combined, "\n")

# Summary
summary_combined <- calc_combined %>%
    group_by(season, organism, division, parent_category) %>%
    summarise(n_terms = n(), mean_padj = mean(p.adj), .groups = "drop") %>%
    arrange(season, organism, parent_category)

write.csv(summary_combined, "data/calcification_summary_combined.csv", row.names = FALSE)
cat("✓ Saved combined summary\n")

cat("\n==============================================================================\n")
cat("SUMMARY\n")
cat("==============================================================================\n")

if (!is.null(summer_result)) {
    cat("\n--- SUMMER ---\n")
    print(summer_result$summary)
}

if (!is.null(winter_result)) {
    cat("\n--- WINTER ---\n")
    print(winter_result$summary)
}

cat("\n--- COMBINED (Summer + Winter) ---\n")
print(summary_combined)

cat("\n✓ All Sankey plots complete!\n")
cat("  - Fig4_calcification_sankey_summer.html (Orange=Host, Green=Symbiont)\n")
cat("  - Fig4_calcification_sankey_winter.html (Orange=Host, Green=Symbiont)\n")
cat("  - Fig4_calcification_sankey_combined.html\n")
cat("      Host Summer=#ff671f, Host Winter=#ffb81c\n")
cat("      Symbiont Summer=#006341, Symbiont Winter=#8fe2b0\n")
cat("      Excludes: Protein Phosphorylation, Na/K Transport, Phospholipid Transport\n")

# ==============================================================================
# FIGURE 4D-E: GO_MWU Bubble Plots - VERSION 5 (Staggered Lanes)
# Y-axis: -log10(p.adj), X-axis: delta.rank (direction), Size: nseqs
# Faceted by Division (BP/MF/CC) × Season (Summer/Winter)
# Annotations: Top 10 significant terms + keyword matches (bolded)
# Labels positioned by delta rank sign: positive=RIGHT, negative=LEFT
# Staggered lanes: Bold=inner, Plain=outer
# ==============================================================================

cat("\n==============================================================================\n")
cat("Figure 4D-E: GO_MWU Bubble Plots (v5 - Staggered Lanes)\n")
cat("==============================================================================\n\n")

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

ANNOTATION_KEYWORDS <- c(
    "\\bion\\b",
    "ion transport",
    "ion channel",
    "ion homeostasis",
    "proton",
    "calcium",
    "carbonate",
    "carbonic",
    "bicarbonate",
    "calcium channel",
    "chemosensory",
    "plasmamembrane",
    "ossification",
    "biomineralization",
    "calcification",
    "endoplasmic reticulum"
)

BUBBLE_P_ADJ_CUTOFF <- 0.05
TOP_N_ANNOTATE <- 6

division_colors <- c(
    "BP" = "#4DAF4A",
    "CC" = "#E41A1C",
    "MF" = "#377EB8"
)

# -----------------------------------------------------------------------------
# GO Term Simplification Function
# -----------------------------------------------------------------------------

simplify_go_term <- function(term) {
    
    # Remove uninformative words
    remove_words <- c(
        "\\bobsolete\\s*",
        "\\bmonoatomic\\s*",
        "\\bprocess\\b",
        "\\bactivity\\b",
        "\\bpart\\b"
    )
    for (pattern in remove_words) {
        term <- str_replace_all(term, regex(pattern, ignore_case = TRUE), "")
    }
    
    # Ion symbols
    ion_replacements <- c(
        "\\bcalcium\\b" = "Ca2+",
        "\\bsodium\\b" = "Na+",
        "\\bpotassium\\b" = "K+",
        "\\bcadmium\\b" = "Cd2+",
        "\\blithium\\b" = "Li+",
        "\\biron\\b" = "Fe",
        "\\bmagnesium\\b" = "Mg2+",
        "\\bzinc\\b" = "Zn2+",
        "\\bcopper\\b" = "Cu2+",
        "\\bproton\\b" = "H+",
        "\\bhydrogen\\b" = "H+"
    )
    for (pattern in names(ion_replacements)) {
        term <- str_replace_all(term, regex(pattern, ignore_case = TRUE), 
                                ion_replacements[pattern])
    }
    
    # Multi-word phrase abbreviations
    phrase_abbrevs <- c(
        "endoplasmic reticulum" = "ER",
        "plasma membrane" = "PM",
        "cell membrane" = "membrane",
        "rough ER" = "rough ER",
        "ER-Golgi intermediate compartment" = "ER-Golgi",
        "positive regulation of" = "+reg.",
        "negative regulation of" = "-reg.",
        "regulation of" = "reg.",
        "response to " = "",
        "involved in " = "",
        "establishment of " = "",
        " involved in.*$" = "",
        "-containing complex" = "",
        "-transporting " = " ",
        "-driven active" = "-driven",
        "-dependent " = "-dep. ",
        "-mediated " = "-med. ",
        "-directed " = "-dir. ",
        "-transcribed " = "-txd ",
        "-type " = " ",
        "transmembrane transporter" = "TM transporter",
        "transmembrane transport" = "TM transport",
        "ion transmembrane" = "ion TM",
        "active transmembrane" = "active TM",
        "ATP synthase" = "ATP synth.",
        "two-sector ATPase" = "ATPase",
        "nervous system" = "NS",
        "central nervous system" = "CNS",
        "nucleic acid" = "nucl. acid",
        "amino acid" = "AA",
        "fatty acid" = "FA",
        "cell-cell" = "cell-cell",
        "protein-RNA" = "prot.-RNA"
    )
    for (pattern in names(phrase_abbrevs)) {
        term <- str_replace_all(term, regex(pattern, ignore_case = TRUE), 
                                phrase_abbrevs[pattern])
    }
    
    # Single word abbreviations
    word_abbrevs <- c(
        "\\bchemosensory\\b" = "chemosen.",
        "\\bprotein\\b" = "prot.",
        "\\bcomplex\\b" = "cplx",
        "\\bligand\\b" = "lig.",
        "\\borganic\\b" = "org.",
        "\\bmitochondrial\\b" = "mito.",
        "\\bmitochondrion\\b" = "mito.",
        "\\bcytoplasmic\\b" = "cyto.",
        "\\bcytosolic\\b" = "cyto.",
        "\\bcytoplasm\\b" = "cyto.",
        "\\bnuclear\\b" = "nucl.",
        "\\bribosomal\\b" = "ribo.",
        "\\bribosome\\b" = "ribo.",
        "\\bpolysomal\\b" = "polysomal",
        "\\bpolysome\\b" = "polysome",
        "\\bpreribosome\\b" = "preribo.",
        "\\bribonucleoprotein\\b" = "RNP",
        "\\bspliceosomal\\b" = "spliceo.",
        "\\btransmembrane\\b" = "TM",
        "\\btransporter\\b" = "transporter",
        "\\btransport\\b" = "transport",
        "\\bchannel\\b" = "ch.",
        "\\bvoltage-gated\\b" = "V-gated",
        "\\btransmitter-gated\\b" = "transmit.-gated",
        "\\borganelle\\b" = "organelle",
        "\\benvelope\\b" = "env.",
        "\\bmembrane\\b" = "memb.",
        "\\blumenal\\b" = "lumen.",
        "\\blumen\\b" = "lumen",
        "\\bcomponent\\b" = "comp.",
        "\\bcompartment\\b" = "compart.",
        "\\bintermediate\\b" = "interm.",
        "\\bintrinsic\\b" = "intrins.",
        "\\bbounded\\b" = "bound.",
        "\\bprojection\\b" = "proj.",
        "\\bjunction\\b" = "junct.",
        "\\banchoring\\b" = "anchor.",
        "\\bsynaptic\\b" = "synap.",
        "\\blocalization\\b" = "local.",
        "\\borganization\\b" = "org.",
        "\\bbiosynthetic\\b" = "biosynth.",
        "\\bcatabolic\\b" = "catabol.",
        "\\bmetabolic\\b" = "metab.",
        "\\bphosphorylation\\b" = "phosph.",
        "\\boxidative\\b" = "oxid.",
        "\\boxidoreductase\\b" = "oxidored.",
        "\\boxidoreduction\\b" = "redox",
        "\\btranslational\\b" = "transl.",
        "\\btranslation\\b" = "transl.",
        "\\btranscription\\b" = "txn",
        "\\binitiation\\b" = "init.",
        "\\bassembly\\b" = "assemb.",
        "\\bbiogenesis\\b" = "biogen.",
        "\\bdevelopment\\b" = "dev.",
        "\\bsignaling\\b" = "signal.",
        "\\bstructural\\b" = "struct.",
        "\\bconstituent\\b" = "const.",
        "\\bmolecule\\b" = "mol.",
        "\\bbinding\\b" = "bind.",
        "\\bregulator\\b" = "reg.",
        "\\bexcitatory\\b" = "excit.",
        "\\bextracellular\\b" = "EC",
        "\\bintracellular\\b" = "IC",
        "\\binorganic\\b" = "inorg.",
        "\\bmonovalent\\b" = "monoval.",
        "\\bmicrotubule\\b" = "MT",
        "\\bcytoskeletal\\b" = "cytoskel.",
        "\\blocomotory\\b" = "locomo.",
        "\\bmovement\\b" = "movmt",
        "\\bmotor\\b" = "motor",
        "\\bdynein\\b" = "dynein",
        "\\bpolypeptide\\b" = "polypep.",
        "\\bconformation\\b" = "conform.",
        "\\bdehydrogenase\\b" = "DH",
        "\\bendopeptidase\\b" = "endopep.",
        "\\bphospholipid\\b" = "phospholip.",
        "\\bphosphatidylinositol\\b" = "PI",
        "\\bbisphosphate\\b" = "bisP",
        "\\bcysteine\\b" = "Cys",
        "\\bhomeostasis\\b" = "homeo.",
        "\\bbehavior\\b" = "behav.",
        "\\bsensory\\b" = "sens.",
        "\\borgan\\b" = "organ"
    )
    for (pattern in names(word_abbrevs)) {
        term <- str_replace_all(term, regex(pattern, ignore_case = TRUE), 
                                word_abbrevs[pattern])
    }
    
    # Cleanup
    term <- str_replace_all(term, "\\s+", " ")
    term <- str_replace_all(term, "^\\s+|\\s+$", "")
    term <- str_replace_all(term, "\\s*-\\s*", "-")
    term <- str_replace_all(term, "\\s+,", ",")
    term <- str_replace_all(term, ",\\s*$", "")
    term <- str_replace_all(term, "\\s*\\.$", "")
    term <- str_replace_all(term, "^\\s*of\\s+", "")
    term <- str_replace_all(term, "\\s+of$", "")
    term <- str_replace_all(term, "\\s+to$", "")
    term <- str_replace_all(term, "\\s+or$", "")
    
    # Limit to 4 words
    words <- str_split(term, "\\s+")[[1]]
    words <- words[words != ""]
    if (length(words) > 4) {
        term <- paste(words[1:4], collapse = " ")
    } else {
        term <- paste(words, collapse = " ")
    }
    
    return(term)
}

# -----------------------------------------------------------------------------
# Load GO_MWU Results
# -----------------------------------------------------------------------------

load_gomwu_bubble_results <- function(base_dir, prefix = "", organism_label) {
    
    seasons <- c("summer", "winter")
    divisions <- c("BP", "MF", "CC")
    all_results <- list()
    
    for (season in seasons) {
        for (div in divisions) {
            if (prefix == "") {
                filename <- file.path(base_dir, "output",
                                     paste0(season, "_", div, "_MWU_results.csv"))
            } else {
                filename <- file.path(base_dir, "output",
                                     paste0(prefix, "_", season, "_", div, "_MWU_results.csv"))
            }
            
            if (file.exists(filename)) {
                df <- read.table(filename, header = TRUE, stringsAsFactors = FALSE,
                                sep = "", quote = "\"", comment.char = "")
                names(df) <- gsub('^"', '', names(df))
                names(df) <- gsub('"$', '', names(df))
                names(df) <- trimws(names(df))
                df$season <- season
                df$division <- div
                df$organism <- organism_label
                all_results[[paste(organism_label, season, div, sep = "_")]] <- df
            } else {
                cat("  Warning: File not found -", filename, "\n")
            }
        }
    }
    return(bind_rows(all_results))
}

# -----------------------------------------------------------------------------
# Process Data for Plotting
# -----------------------------------------------------------------------------

process_for_bubble <- function(gomwu_data, p_cutoff = 0.05, top_n = 6) {
    
    keyword_pattern <- paste(ANNOTATION_KEYWORDS, collapse = "|")
    
    df <- gomwu_data %>%
        filter(!is.na(p.adj)) %>%
        mutate(
            neg_log10_padj = pmin(-log10(p.adj), 30),
            x_value = delta.rank / 100,
            go_name = name,
            season_label = factor(
                ifelse(season == "summer", "Summer", "Winter"),
                levels = c("Summer", "Winter")
            ),
            division = factor(division, levels = c("BP", "CC", "MF")),
            matches_keyword = str_detect(tolower(go_name), regex(keyword_pattern, ignore_case = TRUE))
        )
    
    top_significant <- df %>%
        filter(p.adj < p_cutoff) %>%
        group_by(season_label, division) %>%
        slice_min(order_by = p.adj, n = top_n, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(is_top_significant = TRUE) %>%
        select(go_name, season_label, division, is_top_significant)
    
    df <- df %>%
        left_join(top_significant, by = c("go_name", "season_label", "division")) %>%
        mutate(
            is_top_significant = replace_na(is_top_significant, FALSE),
            should_annotate = (p.adj < p_cutoff & matches_keyword) | is_top_significant,
            annotate_label = ifelse(
                should_annotate, 
                sapply(go_name, simplify_go_term),
                NA_character_
            ),
            label_fontface = case_when(
                !should_annotate ~ NA_character_,
                matches_keyword ~ "bold",
                TRUE ~ "plain"
            ),
            label_side = case_when(
                !should_annotate ~ NA_character_,
                delta.rank >= 0 ~ "right",
                TRUE ~ "left"
            )
        )
    
    return(df)
}

# -----------------------------------------------------------------------------
# Create Bubble Plot Function - STAGGERED LANES
# -----------------------------------------------------------------------------

create_gomwu_bubble <- function(bubble_data, organism_name, org_color) {
    
    sig_line <- -log10(BUBBLE_P_ADJ_CUTOFF)
    
    # Separate data by side AND fontface
    label_data <- bubble_data %>%
        filter(!is.na(annotate_label))
    
    bold_left <- label_data %>% filter(label_fontface == "bold", label_side == "left")
    bold_right <- label_data %>% filter(label_fontface == "bold", label_side == "right")
    plain_left <- label_data %>% filter(label_fontface == "plain", label_side == "left")
    plain_right <- label_data %>% filter(label_fontface == "plain", label_side == "right")
    
    # STAGGERED nudge distances
    # Bold = 20% (inner lane, closer to points)
    # Plain = 38% (outer lane, further out)
    x_range <- range(bubble_data$x_value, na.rm = TRUE)
    nudge_bold <- diff(x_range) * 0.20   # Inner lane
    nudge_plain <- diff(x_range) * 0.38  # Outer lane
    
    p <- ggplot(bubble_data, aes(x = x_value, y = neg_log10_padj)) +
        
        # Points
        geom_point(aes(size = nseqs, fill = division),
                   shape = 21, alpha = 0.6, color = "gray30", stroke = 0.3) +
        
        # Significance threshold line
        geom_hline(yintercept = sig_line, color = "darkorange",
                   linetype = "solid", linewidth = 0.8) +
        
        # Vertical line at 0
        geom_vline(xintercept = 0, color = "gray50",
                   linetype = "dashed", linewidth = 0.5) +
        
        # BOLD LEFT labels (inner lane)
        geom_text_repel(
            data = bold_left,
            aes(label = annotate_label),
            size = 2.4,
            fontface = "bold",
            color = "black",
            xlim = c(-Inf, NA),
            ylim = c(0, NA),  # Keep labels within plot
            hjust = 1,
            direction = "y",
            nudge_x = -nudge_bold,
            segment.color = "gray40",
            segment.size = 0.4,
            segment.linetype = "solid",
            min.segment.length = 0,
            box.padding = 0.2,
            point.padding = 0.15,
            force = 4,
            force_pull = 0,
            max.overlaps = 50,
            seed = 42,
            na.rm = TRUE
        ) +
        
        # BOLD RIGHT labels (inner lane)
        geom_text_repel(
            data = bold_right,
            aes(label = annotate_label),
            size = 2.4,
            fontface = "bold",
            color = "black",
            xlim = c(NA, Inf),
            ylim = c(0, NA),
            hjust = 0,
            direction = "y",
            nudge_x = nudge_bold,
            segment.color = "gray40",
            segment.size = 0.4,
            segment.linetype = "solid",
            min.segment.length = 0,
            box.padding = 0.2,
            point.padding = 0.15,
            force = 4,
            force_pull = 0,
            max.overlaps = 50,
            seed = 42,
            na.rm = TRUE
        ) +
        
        # PLAIN LEFT labels (outer lane)
        geom_text_repel(
            data = plain_left,
            aes(label = annotate_label),
            size = 2.2,
            fontface = "plain",
            color = "gray30",
            xlim = c(-Inf, NA),
            ylim = c(0, NA),
            hjust = 1,
            direction = "y",
            nudge_x = -nudge_plain,
            segment.color = "gray60",
            segment.size = 0.3,
            segment.linetype = "dashed",
            min.segment.length = 0,
            box.padding = 0.2,
            point.padding = 0.15,
            force = 4,
            force_pull = 0,
            max.overlaps = 50,
            seed = 123,  # Different seed for plain labels
            na.rm = TRUE
        ) +
        
        # PLAIN RIGHT labels (outer lane)
        geom_text_repel(
            data = plain_right,
            aes(label = annotate_label),
            size = 2.2,
            fontface = "plain",
            color = "gray30",
            xlim = c(NA, Inf),
            ylim = c(0, NA),
            hjust = 0,
            direction = "y",
            nudge_x = nudge_plain,
            segment.color = "gray60",
            segment.size = 0.3,
            segment.linetype = "dashed",
            min.segment.length = 0,
            box.padding = 0.2,
            point.padding = 0.15,
            force = 4,
            force_pull = 0,
            max.overlaps = 50,
            seed = 123,
            na.rm = TRUE
        ) +
        
        # Scales
        scale_fill_manual(
            values = division_colors,
            name = "GO Division",
            labels = c("BP" = "Biological Process",
                      "CC" = "Cellular Component",
                      "MF" = "Molecular Function")
        ) +
        
        scale_size_continuous(
            name = "# Genes",
            range = c(1, 12),
            breaks = c(10, 50, 100, 200, 500)
        ) +
        
        # Axes - expanded with 8% extra at top
        scale_x_continuous(expand = expansion(mult = c(0.40, 0.40))) +
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
        
        # Facet
        facet_grid(season_label ~ division, scales = "free_x") +
        
        # Labels
        labs(
            x = "Delta Rank (Direction Score)",
            y = expression(-log[10]~italic(p)[adj]),
            title = paste0("GO Enrichment: ", organism_name),
            caption = "Bold = calcification/ion keywords | Plain = top 10 significant | Left = down | Right = up"
        ) +
        
        # Theme
        theme_bw(base_size = 11) +
        theme(
            plot.title = element_text(size = 14, face = "bold.italic",
                                      hjust = 0.5, color = org_color),
            plot.caption = element_text(size = 8, face = "italic", hjust = 0.5,
                                        color = "gray40"),
            axis.title = element_text(size = 11, face = "bold"),
            axis.text = element_text(size = 9),
            strip.text = element_text(size = 10, face = "bold"),
            strip.background = element_rect(fill = "gray95"),
            legend.position = "right",
            legend.box = "vertical",
            panel.grid.minor = element_blank(),
            panel.spacing = unit(1.2, "lines"),  # Increased spacing
            plot.margin = margin(10, 15, 10, 10)  # Extra right margin
        ) +
        
        guides(
            fill = guide_legend(override.aes = list(size = 5), order = 1),
            size = guide_legend(order = 2)
        )
    
    return(p)
}

# -----------------------------------------------------------------------------
# Generate Plots
# -----------------------------------------------------------------------------

cat("Loading Host GO_MWU results...\n")
host_gomwu <- load_gomwu_bubble_results("../10_GO_MWU", prefix = "", "Host")
cat("  Loaded", nrow(host_gomwu), "GO terms\n")

cat("Loading Symbiont GO_MWU results...\n")
symbiont_gomwu <- load_gomwu_bubble_results("../11_symbiont_GO_MWU", prefix = "symbiont", "Symbiont")
cat("  Loaded", nrow(symbiont_gomwu), "GO terms\n")

# Process data
host_bubble_data <- process_for_bubble(host_gomwu, BUBBLE_P_ADJ_CUTOFF, TOP_N_ANNOTATE)
symbiont_bubble_data <- process_for_bubble(symbiont_gomwu, BUBBLE_P_ADJ_CUTOFF, TOP_N_ANNOTATE)

# Show simplified terms
cat("\n=== Simplified Labels (Host Sample) ===\n")
host_sample <- host_bubble_data %>%
    filter(!is.na(annotate_label)) %>%
    select(go_name, annotate_label, label_fontface) %>%
    head(20)
for (i in 1:nrow(host_sample)) {
    cat(sprintf("  %-50s -> %s\n", 
                substr(host_sample$go_name[i], 1, 50),
                host_sample$annotate_label[i]))
}

cat("\n=== Simplified Labels (Symbiont Sample) ===\n")
sym_sample <- symbiont_bubble_data %>%
    filter(!is.na(annotate_label)) %>%
    select(go_name, annotate_label, label_fontface) %>%
    head(20)
for (i in 1:nrow(sym_sample)) {
    cat(sprintf("  %-50s -> %s\n", 
                substr(sym_sample$go_name[i], 1, 50),
                sym_sample$annotate_label[i]))
}

# Summary stats
cat("\nHost annotation breakdown:\n")
host_labels <- host_bubble_data %>% filter(!is.na(annotate_label))
cat("  Bold-Left:", sum(host_labels$label_fontface == "bold" & host_labels$label_side == "left"), "\n")
cat("  Bold-Right:", sum(host_labels$label_fontface == "bold" & host_labels$label_side == "right"), "\n")
cat("  Plain-Left:", sum(host_labels$label_fontface == "plain" & host_labels$label_side == "left"), "\n")
cat("  Plain-Right:", sum(host_labels$label_fontface == "plain" & host_labels$label_side == "right"), "\n")
cat("  Total:", nrow(host_labels), "\n")

cat("\nSymbiont annotation breakdown:\n")
sym_labels <- symbiont_bubble_data %>% filter(!is.na(annotate_label))
cat("  Bold-Left:", sum(sym_labels$label_fontface == "bold" & sym_labels$label_side == "left"), "\n")
cat("  Bold-Right:", sum(sym_labels$label_fontface == "bold" & sym_labels$label_side == "right"), "\n")
cat("  Plain-Left:", sum(sym_labels$label_fontface == "plain" & sym_labels$label_side == "left"), "\n")
cat("  Plain-Right:", sum(sym_labels$label_fontface == "plain" & sym_labels$label_side == "right"), "\n")
cat("  Total:", nrow(sym_labels), "\n")

cat("\nGenerating Host bubble plot...\n")
fig4d_host_bubble <- create_gomwu_bubble(
    host_bubble_data,
    "M. capitata (Host)",
    "#E69F00"
)

cat("Generating Symbiont bubble plot...\n")
fig4e_symbiont_bubble <- create_gomwu_bubble(
    symbiont_bubble_data,
    "D. trenchii (Symbiont)",
    "#56B4E9"
)

# -----------------------------------------------------------------------------
# Save Plots - 15 x 11 inches
# -----------------------------------------------------------------------------

ggsave("figures/Fig4D_bubble_host_GOMWU.pdf", fig4d_host_bubble,
       width = 15, height = 11, dpi = 300)
ggsave("figures/Fig4D_bubble_host_GOMWU.png", fig4d_host_bubble,
       width = 15, height = 11, dpi = 300)
cat("✓ Saved: Fig4D_bubble_host_GOMWU.pdf/png (15x11 inches)\n")

ggsave("figures/Fig4E_bubble_symbiont_GOMWU.pdf", fig4e_symbiont_bubble,
       width = 15, height = 11, dpi = 300)
ggsave("figures/Fig4E_bubble_symbiont_GOMWU.png", fig4e_symbiont_bubble,
       width = 15, height = 11, dpi = 300)
cat("✓ Saved: Fig4E_bubble_symbiont_GOMWU.pdf/png (15x11 inches)\n")

# -----------------------------------------------------------------------------
# Save Summary
# -----------------------------------------------------------------------------

annotated_summary <- bind_rows(
    host_bubble_data %>%
        filter(!is.na(annotate_label)) %>%
        mutate(organism = "Host",
               annotation_type = ifelse(matches_keyword, "Keyword (bold)", "Top significant")) %>%
        select(organism, season = season_label, division, go_name, annotate_label, p.adj, delta.rank, 
               nseqs, matches_keyword, is_top_significant, annotation_type, label_side),
    symbiont_bubble_data %>%
        filter(!is.na(annotate_label)) %>%
        mutate(organism = "Symbiont",
               annotation_type = ifelse(matches_keyword, "Keyword (bold)", "Top significant")) %>%
        select(organism, season = season_label, division, go_name, annotate_label, p.adj, delta.rank, 
               nseqs, matches_keyword, is_top_significant, annotation_type, label_side)
) %>%
    arrange(organism, season, division, label_side, p.adj)

write.csv(annotated_summary, "data/bubble_plot_annotated_terms.csv", row.names = FALSE)
cat("✓ Saved: bubble_plot_annotated_terms.csv\n")

cat("\n=== Summary by Side ===\n")
cat("\nHost:\n")
print(table(host_labels$label_side, host_labels$label_fontface))
cat("\nSymbiont:\n")
print(table(sym_labels$label_side, sym_labels$label_fontface))

cat("\n✓ Figure 4D-E bubble plots complete!\n")

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
            neg_log10_pval = pmin(-log10(pvalue), 50),
            log2FoldChange = pmax(pmin(log2FoldChange, 10), -10),
            significance = case_when(
                padj < padj_cutoff & log2FoldChange > lfc_cutoff ~ "Up",
                padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Down",
                padj < padj_cutoff ~ "Sig (|LFC|<1)",
                TRUE ~ "NS"
            ),
            significance = factor(significance, levels = c("Up", "Down", "Sig (|LFC|<1)", "NS"))
        )

    n_up <- sum(df$significance == "Up", na.rm = TRUE)
    n_down <- sum(df$significance == "Down", na.rm = TRUE)

    ggplot(df, aes(x = log2FoldChange, y = neg_log10_pval, color = significance)) +
        geom_point(alpha = 0.5, size = 1) +
        geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "gray50") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
        scale_color_manual(values = c("Up" = up_color, "Down" = down_color,
                                      "Sig (|LFC|<1)" = "gray50", "NS" = "gray80"),
                           name = "Expression") +
        labs(x = expression(Log[2]~Fold~Change), y = expression(-Log[10]~P-value),
             title = title_text, subtitle = paste0("Up: ", n_up, " | Down: ", n_down)) +
        theme_pub + theme(legend.position = "right")
}

vol_host_summer <- create_volcano(host_summer, "M. cap - Summer")
vol_host_winter <- create_volcano(host_winter, "M. cap - Winter")
vol_sym_summer <- create_volcano(sym_summer, "D. tre - Summer")
vol_sym_winter <- create_volcano(sym_winter, "D. tre - Winter")

ggsave("figures/Fig5A_volcano_host_summer.pdf", vol_host_summer, width = 8, height = 6)
ggsave("figures/Fig5B_volcano_host_winter.pdf", vol_host_winter, width = 8, height = 6)
ggsave("figures/Fig5C_volcano_symbiont_summer.pdf", vol_sym_summer, width = 8, height = 6)
ggsave("figures/Fig5D_volcano_symbiont_winter.pdf", vol_sym_winter, width = 8, height = 6)

vol_combined <- plot_grid(
    vol_host_summer + theme(legend.position = "none"),
    vol_host_winter + theme(legend.position = "none"),
    vol_sym_summer + theme(legend.position = "none"),
    vol_sym_winter + theme(legend.position = "none"),
    ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), label_size = 16
)
legend <- get_legend(vol_host_summer + theme(legend.position = "bottom"))
vol_combined_legend <- plot_grid(vol_combined, legend, ncol = 1, rel_heights = c(1, 0.1))

ggsave("figures/Fig5_volcano_combined.pdf", vol_combined_legend, width = 12, height = 10)
ggsave("figures/Fig5_volcano_combined.png", vol_combined_legend, width = 12, height = 10, dpi = 300)
cat("Saved: Fig5_volcano_*.pdf/png\n")

# ==============================================================================
# FIGURE 6: PCA Plots - Updated with heatmap colors and sample labels
# ==============================================================================

cat("\n==============================================================================\n")
cat("Figure 6: PCA Plots (Updated - Heatmap Colors + Sample Labels)\n")
cat("==============================================================================\n\n")

# Color palettes matching heatmaps
# Host Summer: Warm (orange-red)
# Host Winter: Cool (blue)
# Symbiont Winter: Green (matches heatmap)
# Symbiont Summer: Purple/violet (distinct - not in heatmaps)

pca_colors <- list(
    host_summer = c("OA" = "#D7301F", "Ambient" = "#FDBB84"),
    host_winter = c("OA" = "#08519C", "Ambient" = "#9ECAE1"),
    sym_summer = c("OA" = "#7B2D8E", "Ambient" = "#C9A0DC"),
    sym_winter = c("OA" = "#006D2C", "Ambient" = "#A1D99B")
)

create_pca_from_matrix <- function(vsd_mat, title_text, color_palette) {
    pca_result <- prcomp(t(vsd_mat), scale. = FALSE)
    percent_var <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)
    
    pca_data <- data.frame(
        PC1 = pca_result$x[, 1], 
        PC2 = pca_result$x[, 2],
        sample = rownames(pca_result$x),
        Treatment = ifelse(grepl("B", rownames(pca_result$x)), "OA", "Ambient")
    )
    
    ggplot(pca_data, aes(x = PC1, y = PC2, color = Treatment, shape = Treatment)) +
        geom_point(size = 5, alpha = 0.8) +
        geom_text_repel(aes(label = sample), 
                        size = 3, 
                        fontface = "bold",
                        box.padding = 0.5,
                        point.padding = 0.3,
                        min.segment.length = 0,
                        segment.color = "gray50",
                        show.legend = FALSE) +
        scale_color_manual(values = color_palette, name = "Treatment") +
        scale_shape_manual(values = c("OA" = 16, "Ambient" = 17), name = "Treatment") +
        labs(x = paste0("PC1: ", percent_var[1], "% variance"),
             y = paste0("PC2: ", percent_var[2], "% variance"),
             title = title_text) +
        theme_pub + 
        theme(legend.position = "right")
}

# Create individual PCA plots
pca_host_summer <- create_pca_from_matrix(host_summer_vsd, 
                                           expression(italic("M. capitata")~"- Summer"),
                                           pca_colors$host_summer)

pca_host_winter <- create_pca_from_matrix(host_winter_vsd, 
                                           expression(italic("M. capitata")~"- Winter"),
                                           pca_colors$host_winter)

pca_sym_summer <- create_pca_from_matrix(sym_summer_vsd, 
                                          expression(italic("D. trenchii")~"- Summer"),
                                          pca_colors$sym_summer)

pca_sym_winter <- create_pca_from_matrix(sym_winter_vsd, 
                                          expression(italic("D. trenchii")~"- Winter"),
                                          pca_colors$sym_winter)

# Save individual plots
ggsave("figures/Fig6A_pca_host_summer.pdf", pca_host_summer, width = 7, height = 6)
ggsave("figures/Fig6B_pca_host_winter.pdf", pca_host_winter, width = 7, height = 6)
ggsave("figures/Fig6C_pca_symbiont_summer.pdf", pca_sym_summer, width = 7, height = 6)
ggsave("figures/Fig6D_pca_symbiont_winter.pdf", pca_sym_winter, width = 7, height = 6)

cat("Saved individual PCA plots\n")

# ==============================================================================
# Combined PCA Figure - Manual Legend Approach
# ==============================================================================

# Create plots without legends for the grid
pca_A <- pca_host_summer + theme(legend.position = "none")
pca_B <- pca_host_winter + theme(legend.position = "none")
pca_C <- pca_sym_summer + theme(legend.position = "none")
pca_D <- pca_sym_winter + theme(legend.position = "none")

# Create a standalone legend plot
legend_data <- data.frame(
    x = c(1, 2),
    y = c(1, 2),
    Treatment = factor(c("OA", "Ambient"), levels = c("OA", "Ambient"))
)

legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = Treatment, shape = Treatment)) +
    geom_point(size = 5) +
    scale_color_manual(
        values = c("OA" = "gray30", "Ambient" = "gray70"),
        name = "Treatment",
        labels = c("OA" = "Ocean Acidification (OA)", "Ambient" = "Ambient")
    ) +
    scale_shape_manual(
        values = c("OA" = 16, "Ambient" = 17),
        name = "Treatment",
        labels = c("OA" = "Ocean Acidification (OA)", "Ambient" = "Ambient")
    ) +
    theme_void() +
    theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        legend.key.size = unit(1.2, "cm")
    ) +
    guides(
        color = guide_legend(override.aes = list(size = 5)),
        shape = guide_legend(override.aes = list(size = 5))
    )

# Extract legend using cowplot
pca_legend <- cowplot::get_legend(legend_plot)

# Combine panels
pca_grid <- plot_grid(
    pca_A, pca_B,
    pca_C, pca_D,
    ncol = 2, nrow = 2,
    labels = c("A", "B", "C", "D"),
    label_size = 16,
    align = "hv"
)

# Add legend below
pca_combined_final <- plot_grid(
    pca_grid,
    pca_legend,
    ncol = 1,
    rel_heights = c(1, 0.1)
)

# Save combined figure
ggsave("figures/Fig6_pca_combined.pdf", pca_combined_final, width = 12, height = 11)
ggsave("figures/Fig6_pca_combined.png", pca_combined_final, width = 12, height = 11, dpi = 300)

cat("Saved: Fig6_pca_combined.pdf/png (12x11 inches)\n")

# ==============================================================================
# Alternative: Use gridExtra for more control over legend
# ==============================================================================

# If cowplot legend still doesn't work, try this alternative approach:
library(gridExtra)

# Create the legend as a separate grob
legend_grob <- ggplotGrob(legend_plot)$grobs
legend_index <- which(sapply(legend_grob, function(x) x$name) == "guide-box")

if (length(legend_index) > 0) {
    legend_only <- legend_grob[[legend_index]]
    
    # Alternative combined figure using grid.arrange
    pdf("figures/Fig6_pca_combined_v2.pdf", width = 12, height = 11)
    grid.arrange(
        arrangeGrob(pca_A, pca_B, pca_C, pca_D, 
                    ncol = 2, 
                    top = NULL),
        legend_only,
        nrow = 2,
        heights = c(10, 1)
    )
    dev.off()
    
    png("figures/Fig6_pca_combined_v2.png", width = 12, height = 11, units = "in", res = 300)
    grid.arrange(
        arrangeGrob(pca_A, pca_B, pca_C, pca_D, 
                    ncol = 2,
                    top = NULL),
        legend_only,
        nrow = 2,
        heights = c(10, 1)
    )
    dev.off()
    
    cat("Saved alternative: Fig6_pca_combined_v2.pdf/png\n")
}

cat("\n=== PCA Summary ===\n")
cat("Color scheme (matches heatmaps where applicable):\n")
cat("  Panel A (Host Summer): OA = #D7301F (dark red), Ambient = #FDBB84 (light orange)\n")
cat("  Panel B (Host Winter): OA = #08519C (dark blue), Ambient = #9ECAE1 (light blue)\n")
cat("  Panel C (Sym Summer):  OA = #7B2D8E (dark purple), Ambient = #C9A0DC (light purple)\n")
cat("  Panel D (Sym Winter):  OA = #006D2C (dark green), Ambient = #A1D99B (light green)\n")
cat("\nFeatures:\n")
cat("  - Full sample codes (1BS, 2BW, etc.) as labels\n")
cat("  - No ellipses (too few points per group)\n")
cat("  - Shapes: OA = circle, Ambient = triangle\n")
cat("  - Combined legend at bottom (gray scale for universal reference)\n")

# ==============================================================================
# FIGURE 6E-F: Combined Host + Symbiont PCA per Season (UPDATED)
# Color = Treatment, Shape = Organism, Ellipses by Treatment
# ==============================================================================

cat("\n==============================================================================\n")
cat("Figure 6E-F: Combined Host + Symbiont PCA by Season (Updated)\n")
cat("==============================================================================\n\n")

# Function to run PCA and return data frame with metadata
run_pca_with_metadata <- function(vsd_mat, organism_label) {
    pca_result <- prcomp(t(vsd_mat), scale. = FALSE)
    percent_var <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)
    
    pca_data <- data.frame(
        PC1 = pca_result$x[, 1], 
        PC2 = pca_result$x[, 2],
        sample = rownames(pca_result$x),
        Treatment = ifelse(grepl("B", rownames(pca_result$x)), "OA", "Ambient"),
        Organism = organism_label
    )
    
    return(list(data = pca_data, percent_var = percent_var))
}

# Scale PC coordinates for comparability
scale_pca <- function(pca_data) {
    pca_data$PC1_scaled <- scale(pca_data$PC1)[,1]
    pca_data$PC2_scaled <- scale(pca_data$PC2)[,1]
    return(pca_data)
}

# ------------------------------------------------------------------------------
# Summer: Host + Symbiont Combined
# ------------------------------------------------------------------------------

cat("Creating Summer combined PCA...\n")

host_summer_pca <- run_pca_with_metadata(host_summer_vsd, "Host")
sym_summer_pca <- run_pca_with_metadata(sym_summer_vsd, "Symbiont")

host_summer_scaled <- scale_pca(host_summer_pca$data)
sym_summer_scaled <- scale_pca(sym_summer_pca$data)

summer_combined <- bind_rows(host_summer_scaled, sym_summer_scaled)
summer_combined$Organism <- factor(summer_combined$Organism, levels = c("Host", "Symbiont"))
summer_combined$Treatment <- factor(summer_combined$Treatment, levels = c("OA", "Ambient"))

# Treatment colors (consistent across both plots)
treatment_colors <- c("OA" = "#E69F00", "Ambient" = "#56B4E9")

pca_summer_combined <- ggplot(summer_combined, 
                               aes(x = PC1_scaled, y = PC2_scaled, 
                                   color = Treatment, shape = Organism)) +
    geom_point(size = 5, alpha = 0.8, stroke = 1.5) +
    stat_ellipse(aes(group = Treatment), level = 0.95, linetype = "dashed", linewidth = 0.8) +
    geom_text_repel(aes(label = sample),
                    size = 3,
                    fontface = "bold",
                    box.padding = 0.5,
                    point.padding = 0.3,
                    min.segment.length = 0,
                    segment.color = "gray50",
                    show.legend = FALSE) +
    scale_color_manual(values = treatment_colors, name = "Treatment") +
    scale_shape_manual(values = c("Host" = 16, "Symbiont" = 1), 
                       name = "Organism",
                       labels = c("Host" = expression(italic("M. capitata")),
                                  "Symbiont" = expression(italic("D. trenchii")))) +
    labs(x = "PC1 (scaled)",
         y = "PC2 (scaled)",
         title = "Summer") +
    theme_pub +
    theme(legend.position = "right",
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
    guides(color = guide_legend(override.aes = list(size = 4)),
           shape = guide_legend(override.aes = list(size = 4)))

# ------------------------------------------------------------------------------
# Winter: Host + Symbiont Combined
# ------------------------------------------------------------------------------

cat("Creating Winter combined PCA...\n")

host_winter_pca <- run_pca_with_metadata(host_winter_vsd, "Host")
sym_winter_pca <- run_pca_with_metadata(sym_winter_vsd, "Symbiont")

host_winter_scaled <- scale_pca(host_winter_pca$data)
sym_winter_scaled <- scale_pca(sym_winter_pca$data)

winter_combined <- bind_rows(host_winter_scaled, sym_winter_scaled)
winter_combined$Organism <- factor(winter_combined$Organism, levels = c("Host", "Symbiont"))
winter_combined$Treatment <- factor(winter_combined$Treatment, levels = c("OA", "Ambient"))

pca_winter_combined <- ggplot(winter_combined, 
                               aes(x = PC1_scaled, y = PC2_scaled, 
                                   color = Treatment, shape = Organism)) +
    geom_point(size = 5, alpha = 0.8, stroke = 1.5) +
    stat_ellipse(aes(group = Treatment), level = 0.95, linetype = "dashed", linewidth = 0.8) +
    geom_text_repel(aes(label = sample),
                    size = 3,
                    fontface = "bold",
                    box.padding = 0.5,
                    point.padding = 0.3,
                    min.segment.length = 0,
                    segment.color = "gray50",
                    show.legend = FALSE) +
    scale_color_manual(values = treatment_colors, name = "Treatment") +
    scale_shape_manual(values = c("Host" = 16, "Symbiont" = 1), 
                       name = "Organism",
                       labels = c("Host" = expression(italic("M. capitata")),
                                  "Symbiont" = expression(italic("D. trenchii")))) +
    labs(x = "PC1 (scaled)",
         y = "PC2 (scaled)",
         title = "Winter") +
    theme_pub +
    theme(legend.position = "right",
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
    guides(color = guide_legend(override.aes = list(size = 4)),
           shape = guide_legend(override.aes = list(size = 4)))

# ------------------------------------------------------------------------------
# Save individual combined plots
# ------------------------------------------------------------------------------

ggsave("figures/Fig6E_pca_summer_holobiont.pdf", pca_summer_combined, width = 9, height = 7)
ggsave("figures/Fig6E_pca_summer_holobiont.png", pca_summer_combined, width = 9, height = 7, dpi = 300)
cat("Saved: Fig6E_pca_summer_holobiont.pdf/png\n")

ggsave("figures/Fig6F_pca_winter_holobiont.pdf", pca_winter_combined, width = 9, height = 7)
ggsave("figures/Fig6F_pca_winter_holobiont.png", pca_winter_combined, width = 9, height = 7, dpi = 300)
cat("Saved: Fig6F_pca_winter_holobiont.pdf/png\n")

# ------------------------------------------------------------------------------
# Combined E+F figure
# ------------------------------------------------------------------------------

cat("Creating combined holobiont PCA figure...\n")

pca_E <- pca_summer_combined + theme(legend.position = "none")
pca_F <- pca_winter_combined + theme(legend.position = "none")

# Create unified legend
legend_data_holobiont <- expand.grid(
    Treatment = c("OA", "Ambient"),
    Organism = c("Host", "Symbiont")
)
legend_data_holobiont$x <- 1:4
legend_data_holobiont$y <- 1:4

legend_plot_holobiont <- ggplot(legend_data_holobiont, 
                                 aes(x = x, y = y, color = Treatment, shape = Organism)) +
    geom_point(size = 5, stroke = 1.5) +
    scale_color_manual(values = treatment_colors, name = "Treatment") +
    scale_shape_manual(values = c("Host" = 16, "Symbiont" = 1), 
                       name = "Organism",
                       labels = c("Host" = expression(italic("M. capitata")),
                                  "Symbiont" = expression(italic("D. trenchii")))) +
    theme_void() +
    theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.8, "cm"),
        legend.spacing.x = unit(1, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 4), order = 1),
           shape = guide_legend(override.aes = list(size = 4), order = 2))

# Extract legend
holobiont_legend <- ggplotGrob(legend_plot_holobiont)$grobs
legend_idx <- which(sapply(holobiont_legend, function(x) x$name) == "guide-box")

if (length(legend_idx) > 0) {
    holobiont_legend_grob <- holobiont_legend[[legend_idx]]
    
    # Add panel labels A and B
    pca_E_labeled <- arrangeGrob(pca_E, top = textGrob("A", x = unit(0.02, "npc"), 
                                                        just = "left",
                                                        gp = gpar(fontsize = 16, fontface = "bold")))
    pca_F_labeled <- arrangeGrob(pca_F, top = textGrob("B", x = unit(0.02, "npc"), 
                                                        just = "left",
                                                        gp = gpar(fontsize = 16, fontface = "bold")))
    
    pdf("figures/Fig6EF_pca_holobiont_combined.pdf", width = 14, height = 7)
    grid.arrange(
        arrangeGrob(pca_E_labeled, pca_F_labeled, ncol = 2),
        holobiont_legend_grob,
        nrow = 2,
        heights = c(10, 1)
    )
    dev.off()
    
    png("figures/Fig6EF_pca_holobiont_combined.png", width = 14, height = 7, units = "in", res = 300)
    grid.arrange(
        arrangeGrob(pca_E_labeled, pca_F_labeled, ncol = 2),
        holobiont_legend_grob,
        nrow = 2,
        heights = c(10, 1)
    )
    dev.off()
    
    cat("Saved: Fig6EF_pca_holobiont_combined.pdf/png (14x7 inches)\n")
}

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------

cat("\n=== Holobiont PCA Summary ===\n")
cat("Method: Separate PCAs run on host and symbiont (different gene sets),\n")
cat("        then PC coordinates scaled (z-score) for visual comparison.\n")
cat("\nVisual encoding:\n")
cat("  Color = Treatment (OA = orange, Ambient = blue)\n
  Shape = Organism (Host = filled circle, Symbiont = open circle)\n")
cat("  Ellipses = 95% confidence by treatment (n=6 per treatment)\n")
cat("\nOutput files:\n")
cat("  - Fig6E_pca_summer_holobiont.pdf/png\n")
cat("  - Fig6F_pca_winter_holobiont.pdf/png\n")
cat("  - Fig6EF_pca_holobiont_combined.pdf/png\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n==============================================================================\n")
cat("Summary Statistics\n")
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

cat("\n==============================================================================\n")
cat("Analysis Complete!\n")
cat("==============================================================================\n\n")
sessionInfo()
