#!/usr/bin/env Rscript
# ==============================================================================
# Publication-Ready Figures: M. capitata and D. trenchii OA Response
# Unified script combining figure generation, dendrograms, and chord diagrams
#
# Figure scheme:
#   Fig 1  - Calcification barplot (generated separately in 13_cal_chem)
#   Fig 2  - Venn diagrams + Ridgeline plots (DEG overview)
#   Fig 3  - Combined heatmaps (Host Summer, Host Winter, Symbiont Winter)
#   Fig 4  - Holobiont PCA (Host + Symbiont combined per season)
#   Fig S1 - Individual PCA 4-panel grid (supplementary)
#   Fig 5  - Summer representative GO dendrogram (Host + Symbiont)
#   Fig 6  - Winter representative GO dendrogram (Host + Symbiont)
#   Fig 7  - Chord diagram (shared GO terms between Host and Symbiont)
#
# Author: David Armstrong
# ==============================================================================

setwd("/home/darmstrong4/mc_rework/12_publication_figures")

# ==============================================================================
# LOAD LIBRARIES
# ==============================================================================

cat("Loading required packages...\n")

suppressPackageStartupMessages({
    library(tidyverse)
    library(conflicted)
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
    library(ape)
    library(gt)
})

# Resolve namespace conflicts
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::count)
conflicts_prefer(base::intersect)
conflicts_prefer(base::setdiff)
conflicts_prefer(base::union)

# ==============================================================================
# GLOBAL THEME AND COLOR PALETTES
# ==============================================================================

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

host_color <- "#E69F00"
symbiont_color <- "#56B4E9"
summer_color <- "#D55E00"
winter_color <- "#0072B2"
up_color <- "#D73027"
down_color <- "#4575B4"

# ==============================================================================
# CREATE OUTPUT DIRECTORIES
# ==============================================================================

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)
dir.create("data", showWarnings = FALSE, recursive = TRUE)

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


# ##############################################################################
# FIGURE 2: Venn Diagrams + Ridgeline Plots
# ##############################################################################

cat("\n==============================================================================\n")
cat("Figure 2: Venn Diagrams + Ridgeline Plots\n")
cat("==============================================================================\n\n")

# --- Color Palettes ---
host_summer_colors <- c("1BS" = "#FDBB84", "2BS" = "#EF6548", "3BS" = "#D7301F")
host_winter_colors <- c("1BW" = "#9ECAE1", "2BW" = "#4292C6", "3BW" = "#08519C")
sym_winter_colors <- c("1BW" = "#A1D99B", "2BW" = "#41AB5D", "3BW" = "#006D2C")
venn_summer_color <- "#EF6548"
venn_winter_color <- "#4292C6"

# --- DEG IDs ---
host_summer_ids <- host_summer_degs$gene_id
host_winter_ids <- host_winter_degs$gene_id
host_overlap <- length(intersect(host_summer_ids, host_winter_ids))
sym_summer_ids <- sym_summer_degs$gene_id
sym_winter_ids <- sym_winter_degs$gene_id
sym_overlap <- length(intersect(sym_summer_ids, sym_winter_ids))

# --- Ridgeline data preparation ---
prepare_ridgeline_l2fc <- function(vsd_mat, deg_ids, season = "S") {
    common_genes <- intersect(deg_ids, rownames(vsd_mat))
    cat("  Using", length(common_genes), "DEGs\n")
    if (season == "S") {
        oa_pattern <- "BS$"; ambient_pattern <- "DS$"
    } else {
        oa_pattern <- "BW$"; ambient_pattern <- "DW$"
    }
    oa_samples <- grep(oa_pattern, colnames(vsd_mat), value = TRUE)
    ambient_samples <- grep(ambient_pattern, colnames(vsd_mat), value = TRUE)
    ambient_mean <- rowMeans(vsd_mat[common_genes, ambient_samples, drop = FALSE])
    expr_long <- lapply(oa_samples, function(samp) {
        l2fc <- vsd_mat[common_genes, samp] - ambient_mean
        data.frame(gene_id = common_genes, sample = samp, log2FC = l2fc,
                   sample_num = gsub("[^0-9]", "", samp),
                   sample_label = paste0(gsub("[^0-9]", "", samp), "B", season))
    }) %>% bind_rows()
    return(expr_long)
}

cat("Preparing ridgeline data...\n")
ridge_host_summer <- prepare_ridgeline_l2fc(host_summer_vsd, host_summer_ids, "S")
ridge_host_summer$organism <- "M. capitata"; ridge_host_summer$group <- "Host Summer"

ridge_host_winter <- prepare_ridgeline_l2fc(host_winter_vsd, host_winter_ids, "W")
ridge_host_winter$organism <- "M. capitata"; ridge_host_winter$group <- "Host Winter"

ridge_sym_winter <- prepare_ridgeline_l2fc(sym_winter_vsd, sym_winter_ids, "W")
ridge_sym_winter$organism <- "D. trenchii"; ridge_sym_winter$group <- "Symbiont Winter"

# --- Panel C: Host Summer ridgeline ---
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
    theme(legend.position = "none",
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, color = "#D7301F"),
          axis.text.y = element_text(size = 11, face = "bold"),
          axis.title.x = element_text(size = 11), panel.grid.minor = element_blank())

# --- Panel D: Winter combined ridgeline ---
ridge_winter_combined <- bind_rows(
    ridge_host_winter %>% mutate(group_label = paste0(sample_label, "\n(Host)")),
    ridge_sym_winter %>% mutate(group_label = paste0(sample_label, "\n(Sym)"))
)

ridge_winter_combined <- ridge_winter_combined %>%
    mutate(
        display_label = case_when(
            organism == "D. trenchii" ~ paste0(sample_label, " (Sym)"),
            organism == "M. capitata" ~ paste0(sample_label, " (Host)")
        ),
        display_label = factor(display_label,
                               levels = c("3BW (Sym)", "2BW (Sym)", "1BW (Sym)",
                                          "3BW (Host)", "2BW (Host)", "1BW (Host)"))
    )

winter_combined_colors <- c(
    "1BW (Host)" = "#9ECAE1", "2BW (Host)" = "#4292C6", "3BW (Host)" = "#08519C",
    "1BW (Sym)" = "#A1D99B", "2BW (Sym)" = "#41AB5D", "3BW (Sym)" = "#006D2C"
)

p_ridge_winter <- ggplot(ridge_winter_combined,
                          aes(x = log2FC, y = display_label, fill = display_label)) +
    geom_density_ridges(alpha = 0.8, scale = 1.2, rel_min_height = 0.01) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray60", linewidth = 0.4) +
    scale_fill_manual(values = winter_combined_colors) +
    scale_x_continuous(limits = c(-4, 4), breaks = seq(-4, 4, 1)) +
    labs(x = expression(Log[2]~Fold~Change), y = "", title = "Winter") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title.x = element_text(size = 11), panel.grid.minor = element_blank())

# --- Combined Figure 2 ---
cat("Creating combined Figure 2...\n")

for (ext in c("pdf", "png")) {
    if (ext == "pdf") {
        pdf("figures/Fig2_combined_venn_ridgeline.pdf", width = 12, height = 10)
    } else {
        png("figures/Fig2_combined_venn_ridgeline.png", width = 12, height = 10, units = "in", res = 300)
    }

    grid.newpage()
    pushViewport(viewport(layout = grid.layout(
        nrow = 2, ncol = 3,
        heights = unit(c(0.45, 0.55), "npc"),
        widths = unit(c(0.42, 0.42, 0.16), "npc")
    )))

    # Panel A: Host Venn
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
    pushViewport(viewport(width = 0.9, height = 0.9))
    grid.rect(gp = gpar(col = NA, fill = "white"))
    venn_host <- draw.pairwise.venn(
        area1 = length(host_summer_ids), area2 = length(host_winter_ids),
        cross.area = host_overlap, category = c("", ""),
        fill = c(venn_summer_color, venn_winter_color), alpha = 0.6,
        col = c(venn_summer_color, venn_winter_color), lwd = 2, cex = 1.8,
        fontface = "bold", fontfamily = "sans", cat.cex = 0, margin = 0.05, ind = FALSE)
    grid.draw(venn_host)
    upViewport()
    grid.text("A", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
              just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "sans"))
    grid.text(expression(italic("M. capitata")), x = unit(0.5, "npc"), y = unit(0.97, "npc"),
              just = "center", gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "sans"))
    upViewport()

    # Panel B: Symbiont Venn
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
    pushViewport(viewport(width = 0.9, height = 0.9))
    grid.rect(gp = gpar(col = NA, fill = "white"))
    venn_sym <- draw.pairwise.venn(
        area1 = length(sym_summer_ids), area2 = length(sym_winter_ids),
        cross.area = sym_overlap, category = c("", ""),
        fill = c(venn_summer_color, venn_winter_color), alpha = 0.6,
        col = c(venn_summer_color, venn_winter_color), lwd = 2, cex = 1.8,
        fontface = "bold", fontfamily = "sans", cat.cex = 0, margin = 0.05, ind = FALSE)
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

    # Panel C: Host Summer Ridgeline
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
    pushViewport(viewport(width = 0.95, height = 0.92))
    grid.draw(ggplotGrob(p_ridge_summer))
    upViewport()
    grid.text("C", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
              just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "sans"))
    upViewport()

    # Panel D: Winter Combined Ridgeline
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2:3))
    pushViewport(viewport(width = 0.95, height = 0.92))
    grid.draw(ggplotGrob(p_ridge_winter))
    upViewport()
    grid.text("D", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
              just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "sans"))
    upViewport()

    dev.off()
}

ggsave("figures/Fig2C_ridgeline_host_summer.pdf", p_ridge_summer, width = 6, height = 4)
ggsave("figures/Fig2D_ridgeline_winter_combined.pdf", p_ridge_winter, width = 6, height = 5)

cat("Saved: Fig2_combined_venn_ridgeline.pdf/png\n")


# ##############################################################################
# FIGURE 3: Combined Heatmaps
# ##############################################################################

cat("\n==============================================================================\n")
cat("Figure 3: Combined Heatmaps\n")
cat("==============================================================================\n\n")

heatmap_colors_host_summer <- colorRampPalette(c("#FFF5EB", "#FDD49E", "#FDBB84",
                                                  "#FC8D59", "#EF6548", "#D7301F", "#990000"))(100)
heatmap_colors_host_winter <- colorRampPalette(c("#F7FBFF", "#DEEBF7", "#C6DBEF",
                                                  "#9ECAE1", "#6BAED6", "#3182BD", "#08519C"))(100)
heatmap_colors_sym_winter <- colorRampPalette(c("#F7FCF5", "#E5F5E0", "#C7E9C0",
                                                 "#A1D99B", "#74C476", "#31A354", "#006D2C"))(100)

create_deg_heatmap_panel <- function(vsd_mat, deseq_results, color_palette,
                                      max_genes = 250, show_legend = TRUE) {
    all_degs <- deseq_results %>%
        filter(!is.na(padj) & padj < 0.05) %>%
        arrange(desc(abs(log2FoldChange)))
    if (nrow(all_degs) < 5) { cat("  Warning: Fewer than 5 DEGs\n"); return(NULL) }
    if (nrow(all_degs) > max_genes) {
        cat("  Note: Limiting to top", max_genes, "of", nrow(all_degs), "DEGs\n")
        all_degs <- head(all_degs, max_genes)
    } else { cat("  Using all", nrow(all_degs), "DEGs\n") }
    common_genes <- intersect(all_degs$gene_id, rownames(vsd_mat))
    if (length(common_genes) < 5) return(NULL)
    heatmap_mat <- vsd_mat[common_genes, , drop = FALSE]
    heatmap_mat_scaled <- t(scale(t(heatmap_mat)))
    sample_names <- colnames(heatmap_mat)
    treatments <- ifelse(grepl("B", sample_names), "OA", "Ambient")
    anno_col <- data.frame(Treatment = treatments, row.names = sample_names)
    anno_colors <- list(Treatment = c("OA" = "#E69F00", "Ambient" = "#56B4E9"))
    ht <- pheatmap::pheatmap(heatmap_mat_scaled, color = color_palette,
                              cluster_rows = TRUE, cluster_cols = TRUE,
                              show_rownames = FALSE, show_colnames = TRUE,
                              annotation_col = anno_col, annotation_colors = anno_colors,
                              annotation_names_col = FALSE, annotation_legend = show_legend,
                              fontsize = 10, border_color = NA, silent = TRUE)
    return(list(heatmap = ht, n_genes = length(common_genes)))
}

cat("Creating heatmaps...\n")
ht_host_summer <- create_deg_heatmap_panel(host_summer_vsd, host_summer,
                                            heatmap_colors_host_summer, 250, show_legend = FALSE)
ht_host_winter <- create_deg_heatmap_panel(host_winter_vsd, host_winter,
                                            heatmap_colors_host_winter, 250, show_legend = FALSE)
ht_sym_winter <- create_deg_heatmap_panel(sym_winter_vsd, sym_winter,
                                           heatmap_colors_sym_winter, 250, show_legend = TRUE)

# Individual panels
if (!is.null(ht_host_summer)) {
    pdf("figures/Fig3A_heatmap_host_summer.pdf", width = 8, height = 10)
    grid.draw(ht_host_summer$heatmap$gtable); dev.off()
}
if (!is.null(ht_host_winter)) {
    pdf("figures/Fig3B_heatmap_host_winter.pdf", width = 8, height = 10)
    grid.draw(ht_host_winter$heatmap$gtable); dev.off()
}
if (!is.null(ht_sym_winter)) {
    pdf("figures/Fig3C_heatmap_symbiont_winter.pdf", width = 8, height = 10)
    grid.draw(ht_sym_winter$heatmap$gtable); dev.off()
}

# Combined
gt_hs <- if (!is.null(ht_host_summer)) ht_host_summer$heatmap$gtable else NULL
gt_hw <- if (!is.null(ht_host_winter)) ht_host_winter$heatmap$gtable else NULL
gt_sw <- if (!is.null(ht_sym_winter)) ht_sym_winter$heatmap$gtable else NULL

for (ext in c("pdf", "png")) {
    if (ext == "pdf") pdf("figures/Fig3_heatmap_combined.pdf", width = 18, height = 10)
    else png("figures/Fig3_heatmap_combined.png", width = 18, height = 10, units = "in", res = 300)

    grid.newpage()
    vp_layout <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                                                heights = unit(c(0.05, 0.95), "npc"),
                                                widths = unit(c(1, 1, 1), "null")))
    pushViewport(vp_layout)

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid.text("A", x = 0.05, y = 0.5, just = "left", gp = gpar(fontsize = 18, fontface = "bold"))
    grid.text(expression(italic("M. capitata")~"- Summer"), x = 0.5, y = 0.5,
              gp = gpar(fontsize = 12, fontface = "bold")); upViewport()

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
    grid.text("B", x = 0.05, y = 0.5, just = "left", gp = gpar(fontsize = 18, fontface = "bold"))
    grid.text(expression(italic("M. capitata")~"- Winter"), x = 0.5, y = 0.5,
              gp = gpar(fontsize = 12, fontface = "bold")); upViewport()

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
    grid.text("C", x = 0.05, y = 0.5, just = "left", gp = gpar(fontsize = 18, fontface = "bold"))
    grid.text(expression(italic("D. trenchii")~"- Winter"), x = 0.5, y = 0.5,
              gp = gpar(fontsize = 12, fontface = "bold")); upViewport()

    if (!is.null(gt_hs)) {
        pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
        pushViewport(viewport(width = 0.95, height = 0.95))
        grid.draw(gt_hs); upViewport(2)
    }
    if (!is.null(gt_hw)) {
        pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
        pushViewport(viewport(width = 0.95, height = 0.95))
        grid.draw(gt_hw); upViewport(2)
    }
    if (!is.null(gt_sw)) {
        pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
        pushViewport(viewport(width = 0.95, height = 0.95))
        grid.draw(gt_sw); upViewport(2)
    }
    dev.off()
}

cat("Saved: Fig3_heatmap_combined.pdf/png\n")


# ##############################################################################
# FIGURE 4: Holobiont PCA (Host + Symbiont combined per season)
# ##############################################################################

cat("\n==============================================================================\n")
cat("Figure 4: Holobiont PCA\n")
cat("==============================================================================\n\n")

run_pca_with_metadata <- function(vsd_mat, organism_label) {
    pca_result <- prcomp(t(vsd_mat), scale. = FALSE)
    percent_var <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)
    pca_data <- data.frame(
        PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2],
        sample = rownames(pca_result$x),
        Treatment = ifelse(grepl("B", rownames(pca_result$x)), "OA", "Ambient"),
        Organism = organism_label)
    return(list(data = pca_data, percent_var = percent_var))
}

scale_pca <- function(pca_data) {
    pca_data$PC1_scaled <- scale(pca_data$PC1)[,1]
    pca_data$PC2_scaled <- scale(pca_data$PC2)[,1]
    return(pca_data)
}

treatment_colors <- c("OA" = "#E69F00", "Ambient" = "#56B4E9")

# Summer combined
host_summer_pca <- run_pca_with_metadata(host_summer_vsd, "Host")
sym_summer_pca <- run_pca_with_metadata(sym_summer_vsd, "Symbiont")
summer_combined <- bind_rows(scale_pca(host_summer_pca$data), scale_pca(sym_summer_pca$data))
summer_combined$Organism <- factor(summer_combined$Organism, levels = c("Host", "Symbiont"))
summer_combined$Treatment <- factor(summer_combined$Treatment, levels = c("OA", "Ambient"))

pca_summer_combined <- ggplot(summer_combined,
                               aes(x = PC1_scaled, y = PC2_scaled, color = Treatment, shape = Organism)) +
    geom_point(size = 5, alpha = 0.8, stroke = 1.5) +
    stat_ellipse(aes(group = Treatment), level = 0.9, linetype = "dashed", linewidth = 0.8) +
    geom_text_repel(aes(label = sample), size = 3, fontface = "bold", box.padding = 0.5,
                    point.padding = 0.3, min.segment.length = 0, segment.color = "gray50", show.legend = FALSE) +
    scale_color_manual(values = treatment_colors, name = "Treatment") +
    scale_shape_manual(values = c("Host" = 16, "Symbiont" = 1), name = "Organism",
                       labels = c("Host" = expression(italic("M. capitata")),
                                  "Symbiont" = expression(italic("D. trenchii")))) +
    labs(x = "PC1 (scaled)", y = "PC2 (scaled)", title = "Summer") +
    theme_pub + theme(legend.position = "right", plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

# Winter combined
host_winter_pca <- run_pca_with_metadata(host_winter_vsd, "Host")
sym_winter_pca <- run_pca_with_metadata(sym_winter_vsd, "Symbiont")
winter_combined <- bind_rows(scale_pca(host_winter_pca$data), scale_pca(sym_winter_pca$data))
winter_combined$Organism <- factor(winter_combined$Organism, levels = c("Host", "Symbiont"))
winter_combined$Treatment <- factor(winter_combined$Treatment, levels = c("OA", "Ambient"))

pca_winter_combined <- ggplot(winter_combined,
                               aes(x = PC1_scaled, y = PC2_scaled, color = Treatment, shape = Organism)) +
    geom_point(size = 5, alpha = 0.8, stroke = 1.5) +
    stat_ellipse(aes(group = Treatment), level = 0.9, linetype = "dashed", linewidth = 0.8) +
    geom_text_repel(aes(label = sample), size = 3, fontface = "bold", box.padding = 0.5,
                    point.padding = 0.3, min.segment.length = 0, segment.color = "gray50", show.legend = FALSE) +
    scale_color_manual(values = treatment_colors, name = "Treatment") +
    scale_shape_manual(values = c("Host" = 16, "Symbiont" = 1), name = "Organism",
                       labels = c("Host" = expression(italic("M. capitata")),
                                  "Symbiont" = expression(italic("D. trenchii")))) +
    labs(x = "PC1 (scaled)", y = "PC2 (scaled)", title = "Winter") +
    theme_pub + theme(legend.position = "right", plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

# Save individual
ggsave("figures/Fig4A_pca_summer_holobiont.pdf", pca_summer_combined, width = 9, height = 7)
ggsave("figures/Fig4A_pca_summer_holobiont.png", pca_summer_combined, width = 9, height = 7, dpi = 300)
ggsave("figures/Fig4B_pca_winter_holobiont.pdf", pca_winter_combined, width = 9, height = 7)
ggsave("figures/Fig4B_pca_winter_holobiont.png", pca_winter_combined, width = 9, height = 7, dpi = 300)

# Combined E+F
pca_E <- pca_summer_combined + theme(legend.position = "none")
pca_F <- pca_winter_combined + theme(legend.position = "none")

legend_data_holobiont <- expand.grid(Treatment = c("OA", "Ambient"), Organism = c("Host", "Symbiont"))
legend_data_holobiont$x <- 1:4; legend_data_holobiont$y <- 1:4

legend_plot_holobiont <- ggplot(legend_data_holobiont,
                                 aes(x = x, y = y, color = Treatment, shape = Organism)) +
    geom_point(size = 5, stroke = 1.5) +
    scale_color_manual(values = treatment_colors, name = "Treatment") +
    scale_shape_manual(values = c("Host" = 16, "Symbiont" = 1), name = "Organism",
                       labels = c("Host" = expression(italic("M. capitata")),
                                  "Symbiont" = expression(italic("D. trenchii")))) +
    theme_void() +
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal",
          legend.title = element_text(size = 11, face = "bold"), legend.text = element_text(size = 10),
          legend.key.size = unit(0.8, "cm"), legend.spacing.x = unit(1, "cm"))

holobiont_legend <- ggplotGrob(legend_plot_holobiont)$grobs
legend_idx <- which(sapply(holobiont_legend, function(x) x$name) == "guide-box")

if (length(legend_idx) > 0) {
    holobiont_legend_grob <- holobiont_legend[[legend_idx]]
    pca_E_labeled <- arrangeGrob(pca_E, top = textGrob("A", x = unit(0.02, "npc"),
                                                        just = "left", gp = gpar(fontsize = 16, fontface = "bold")))
    pca_F_labeled <- arrangeGrob(pca_F, top = textGrob("B", x = unit(0.02, "npc"),
                                                        just = "left", gp = gpar(fontsize = 16, fontface = "bold")))

    for (ext in c("pdf", "png")) {
        if (ext == "pdf") pdf("figures/Fig4_pca_holobiont_combined.pdf", width = 14, height = 7)
        else png("figures/Fig4_pca_holobiont_combined.png", width = 14, height = 7, units = "in", res = 300)
        grid.arrange(arrangeGrob(pca_E_labeled, pca_F_labeled, ncol = 2),
                     holobiont_legend_grob, nrow = 2, heights = c(10, 1))
        dev.off()
    }
}

cat("Saved: Fig4_pca_holobiont_combined.pdf/png\n")


# ##############################################################################
# FIGURE S1: Individual PCA 4-panel Grid (Supplementary)
# ##############################################################################

cat("\n==============================================================================\n")
cat("Figure S1: Individual PCA 4-panel Grid\n")
cat("==============================================================================\n\n")

pca_colors <- list(
    host_summer = c("OA" = "#D7301F", "Ambient" = "#FDBB84"),
    host_winter = c("OA" = "#08519C", "Ambient" = "#9ECAE1"),
    sym_summer = c("OA" = "#7B2D8E", "Ambient" = "#C9A0DC"),
    sym_winter = c("OA" = "#006D2C", "Ambient" = "#A1D99B")
)

create_pca_from_matrix <- function(vsd_mat, title_text, color_palette) {
    pca_result <- prcomp(t(vsd_mat), scale. = FALSE)
    percent_var <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)
    pca_data <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2],
                           sample = rownames(pca_result$x),
                           Treatment = ifelse(grepl("B", rownames(pca_result$x)), "OA", "Ambient"))
    ggplot(pca_data, aes(x = PC1, y = PC2, color = Treatment, shape = Treatment)) +
        geom_point(size = 5, alpha = 0.8) +
        geom_text_repel(aes(label = sample), size = 3, fontface = "bold", box.padding = 0.5,
                        point.padding = 0.3, min.segment.length = 0, segment.color = "gray50", show.legend = FALSE) +
        scale_color_manual(values = color_palette, name = "Treatment") +
        scale_shape_manual(values = c("OA" = 16, "Ambient" = 17), name = "Treatment") +
        labs(x = paste0("PC1: ", percent_var[1], "% variance"),
             y = paste0("PC2: ", percent_var[2], "% variance"), title = title_text) +
        theme_pub + theme(legend.position = "right")
}

pca_A <- create_pca_from_matrix(host_summer_vsd, expression(italic("M. capitata")~"- Summer"), pca_colors$host_summer) + theme(legend.position = "none")
pca_B <- create_pca_from_matrix(host_winter_vsd, expression(italic("M. capitata")~"- Winter"), pca_colors$host_winter) + theme(legend.position = "none")
pca_C <- create_pca_from_matrix(sym_summer_vsd, expression(italic("D. trenchii")~"- Summer"), pca_colors$sym_summer) + theme(legend.position = "none")
pca_D <- create_pca_from_matrix(sym_winter_vsd, expression(italic("D. trenchii")~"- Winter"), pca_colors$sym_winter) + theme(legend.position = "none")

legend_data <- data.frame(x = c(1, 2), y = c(1, 2),
                          Treatment = factor(c("OA", "Ambient"), levels = c("OA", "Ambient")))
legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = Treatment, shape = Treatment)) +
    geom_point(size = 5) +
    scale_color_manual(values = c("OA" = "gray30", "Ambient" = "gray70"), name = "Treatment",
                       labels = c("OA" = "Ocean Acidification (OA)", "Ambient" = "Ambient")) +
    scale_shape_manual(values = c("OA" = 16, "Ambient" = 17), name = "Treatment",
                       labels = c("OA" = "Ocean Acidification (OA)", "Ambient" = "Ambient")) +
    theme_void() + theme(legend.position = "bottom", legend.direction = "horizontal",
                         legend.title = element_text(size = 12, face = "bold"),
                         legend.text = element_text(size = 11), legend.key.size = unit(1.2, "cm"))

pca_legend <- cowplot::get_legend(legend_plot)
pca_grid <- plot_grid(pca_A, pca_B, pca_C, pca_D, ncol = 2, nrow = 2,
                      labels = c("A", "B", "C", "D"), label_size = 16, align = "hv")
pca_combined_final <- plot_grid(pca_grid, pca_legend, ncol = 1, rel_heights = c(1, 0.1))

ggsave("figures/FigS1_pca_individual_combined.pdf", pca_combined_final, width = 12, height = 11)
ggsave("figures/FigS1_pca_individual_combined.png", pca_combined_final, width = 12, height = 11, dpi = 300)

cat("Saved: FigS1_pca_individual_combined.pdf/png\n")


# ##############################################################################
# FIGURES 5 & 6: Representative GO Dendrograms (from script 14)
# ##############################################################################

cat("\n==============================================================================\n")
cat("Figures 5 & 6: Representative GO Dendrograms\n")
cat("==============================================================================\n\n")

# --- Configuration ---
host_go_dir <- "/home/darmstrong4/mc_rework/10_host_GO_MWU/output"
symbiont_go_dir <- "/home/darmstrong4/mc_rework/11_symbiont_GO_MWU/output"

level1 <- 0.1; level2 <- 0.05; level3 <- 0.01
dendro_colors <- c("dodgerblue2", "firebrick1", "skyblue2", "lightcoral")
division_labels <- c("BP" = "Biological Process", "CC" = "Cellular Component", "MF" = "Molecular Function")

# --- Helper functions ---

load_division_data <- function(rep_file, mwu_file, dissim_file, division) {
    if (!file.exists(rep_file) || !file.exists(mwu_file) || !file.exists(dissim_file)) return(NULL)
    rep_gos <- read.csv(rep_file, stringsAsFactors = FALSE)
    if (nrow(rep_gos) < 1) return(NULL)
    mwu <- read.table(mwu_file, header = TRUE, stringsAsFactors = FALSE, quote = "\"", fill = TRUE)

    rep_gos$go_id <- NA; rep_gos$delta_rank <- NA; rep_gos$nseqs_total <- NA
    for (i in 1:nrow(rep_gos)) {
        match_idx <- which(mwu$name == rep_gos$term[i])
        if (length(match_idx) == 0) match_idx <- which(tolower(mwu$name) == tolower(rep_gos$term[i]))
        if (length(match_idx) > 0) {
            rep_gos$go_id[i] <- mwu$term[match_idx[1]]
            rep_gos$delta_rank[i] <- mwu$delta.rank[match_idx[1]]
            rep_gos$nseqs_total[i] <- mwu$nseqs[match_idx[1]]
        }
    }
    rep_gos <- rep_gos[!is.na(rep_gos$go_id), ]
    if (nrow(rep_gos) < 2) return(NULL)

    diss <- read.table(dissim_file, sep = "\t", header = TRUE, check.names = FALSE)
    if (!any(grepl("^GO:", rownames(diss)))) rownames(diss) <- colnames(diss)

    available_ids <- c(); id_mapping <- list()
    for (go_id in rep_gos$go_id) {
        if (is.na(go_id)) next
        if (go_id %in% colnames(diss)) { available_ids <- c(available_ids, go_id); id_mapping[[go_id]] <- go_id
        } else {
            compound_match <- grep(paste0("(^|;)", go_id, "(;|$)"), colnames(diss), value = TRUE)
            if (length(compound_match) > 0) { available_ids <- c(available_ids, go_id); id_mapping[[go_id]] <- compound_match[1] }
        }
    }
    if (length(available_ids) < 2) return(NULL)
    rep_gos <- rep_gos[rep_gos$go_id %in% available_ids, ]
    matrix_ids <- sapply(rep_gos$go_id, function(x) id_mapping[[x]])
    diss_sub <- diss[matrix_ids, matrix_ids]
    rownames(diss_sub) <- rep_gos$go_id; colnames(diss_sub) <- rep_gos$go_id
    diss_sub <- as.matrix(diss_sub); mode(diss_sub) <- "numeric"
    if (any(is.na(diss_sub))) diss_sub[is.na(diss_sub)] <- 1.0
    hc <- hclust(as.dist(diss_sub), method = "average")
    rep_gos$label <- paste0(rep_gos$nseqs, "/", rep_gos$nseqs_total, " ", rep_gos$term)
    label_map <- setNames(rep_gos$label, rep_gos$go_id)
    hc$labels <- label_map[hc$labels]
    ordered_labels <- hc$labels[hc$order]
    plot_data <- data.frame(label = ordered_labels, stringsAsFactors = FALSE)
    for (i in 1:nrow(plot_data)) {
        match_idx <- which(rep_gos$label == plot_data$label[i])
        if (length(match_idx) > 0) {
            plot_data$pval[i] <- rep_gos$pval[match_idx[1]]
            plot_data$direction[i] <- rep_gos$direction[match_idx[1]]
            plot_data$delta_rank[i] <- rep_gos$delta_rank[match_idx[1]]
        }
    }
    if ("direction" %in% names(plot_data) && !all(is.na(plot_data$direction)))
        plot_data$is_up <- plot_data$direction == "up"
    else plot_data$is_up <- plot_data$delta_rank > 0
    plot_data$color <- NA
    for (i in 1:nrow(plot_data)) {
        p <- plot_data$pval[i]; up <- plot_data$is_up[i]
        if (p < level2) plot_data$color[i] <- ifelse(up, dendro_colors[2], dendro_colors[1])
        else plot_data$color[i] <- ifelse(up, dendro_colors[4], dendro_colors[3])
    }
    plot_data$font <- 3; plot_data$font[plot_data$pval < level2] <- 1; plot_data$font[plot_data$pval < level3] <- 2
    plot_data$cex <- 0.72; plot_data$cex[plot_data$pval < level3] <- 0.85
    return(list(hc = hc, plot_data = plot_data, n_terms = nrow(plot_data), division = division))
}

load_organism_data <- function(data_dir, prefix, season) {
    divisions <- c("BP", "CC", "MF"); div_data <- list()
    for (div in divisions) {
        if (prefix == "") {
            rep_file <- file.path(data_dir, paste0(season, "_", div, "_representative_GOs.csv"))
            mwu_file <- file.path(data_dir, paste0(season, "_", div, "_MWU_results.csv"))
            dissim_file <- file.path(data_dir, paste0(season, "_", div, "_dissimilarity.csv"))
        } else {
            rep_file <- file.path(data_dir, paste0(prefix, "_", season, "_", div, "_representative_GOs.csv"))
            mwu_file <- file.path(data_dir, paste0(prefix, "_", season, "_", div, "_MWU_results.csv"))
            dissim_file <- file.path(data_dir, paste0(prefix, "_", season, "_", div, "_dissimilarity.csv"))
        }
        div_data[[div]] <- load_division_data(rep_file, mwu_file, dissim_file, div)
    }
    return(div_data)
}

draw_division_dendrogram <- function(hc, plot_data, y_offset, tree_x_start, tree_x_end, label_x) {
    n_terms <- nrow(plot_data)
    y_positions <- y_offset + seq(0.5, n_terms - 0.5, by = 1)
    obs_y <- setNames(y_positions, hc$labels[hc$order])
    max_height <- max(hc$height); if (max_height == 0) max_height <- 1
    tree_width <- tree_x_end - tree_x_start
    node_y <- numeric(nrow(hc$merge)); node_x <- numeric(nrow(hc$merge))
    for (j in 1:nrow(hc$merge)) {
        left <- hc$merge[j, 1]; right <- hc$merge[j, 2]
        if (left < 0) { left_y <- obs_y[hc$labels[-left]]; left_x <- tree_x_end
        } else { left_y <- node_y[left]; left_x <- node_x[left] }
        if (right < 0) { right_y <- obs_y[hc$labels[-right]]; right_x <- tree_x_end
        } else { right_y <- node_y[right]; right_x <- node_x[right] }
        node_y[j] <- (left_y + right_y) / 2
        node_x[j] <- tree_x_end - (hc$height[j] / max_height) * tree_width
        segments(node_x[j], left_y, left_x, left_y, col = "gray30", lwd = 0.6)
        segments(node_x[j], right_y, right_x, right_y, col = "gray30", lwd = 0.6)
        segments(node_x[j], left_y, node_x[j], right_y, col = "gray30", lwd = 0.6)
    }
    for (i in 1:n_terms) {
        label <- sub(" activity$", "", plot_data$label[i])
        text(label_x, y_positions[i], label, font = plot_data$font[i],
             cex = plot_data$cex[i], col = plot_data$color[i], adj = c(0, 0.5))
    }
    return(c(y_offset, y_offset + n_terms))
}

draw_organism_column <- function(div_data, x_offset, x_width, total_height, gap_size) {
    divisions <- c("BP", "CC", "MF")
    valid_divs <- divisions[!sapply(div_data[divisions], is.null)]
    if (length(valid_divs) == 0) return(NULL)
    tree_x_start <- x_offset + 0.01 * x_width
    tree_x_end <- x_offset + 0.12 * x_width
    label_x <- x_offset + 0.13 * x_width
    y_cursor <- total_height - 0.5
    for (i in seq_along(valid_divs)) {
        div <- valid_divs[i]; d <- div_data[[div]]; n_terms <- d$n_terms
        y_offset <- y_cursor - n_terms
        text(tree_x_start, y_cursor + 0.15, division_labels[div],
             font = 2, cex = 0.62, adj = c(0, 0), col = "gray25")
        draw_division_dendrogram(d$hc, d$plot_data, y_offset, tree_x_start, tree_x_end, label_x)
        y_cursor <- y_offset - gap_size
    }
    return(TRUE)
}

plot_combined_season <- function(season, output_file, fig_label) {
    cat("\n=== Creating", fig_label, "for", season, "===\n")
    host_data <- load_organism_data(host_go_dir, "", season)
    symbiont_data <- load_organism_data(symbiont_go_dir, "symbiont", season)
    divisions <- c("BP", "CC", "MF")
    host_valid <- divisions[!sapply(host_data[divisions], is.null)]
    symbiont_valid <- divisions[!sapply(symbiont_data[divisions], is.null)]
    host_terms <- sum(sapply(host_data[host_valid], function(x) x$n_terms))
    symbiont_terms <- sum(sapply(symbiont_data[symbiont_valid], function(x) x$n_terms))
    cat("  Host terms:", host_terms, "| Symbiont terms:", symbiont_terms, "\n")
    max_terms <- max(host_terms, symbiont_terms)
    n_gaps <- max(length(host_valid), length(symbiont_valid)) - 1
    gap_size <- 1.2; total_height <- max_terms + n_gaps * gap_size + 2
    line_height <- 0.18; fig_height <- max(4, total_height * line_height); fig_width <- 13

    pdf(output_file, width = fig_width, height = fig_height)
    par(mar = c(0.5, 0.5, 0.5, 0.5))
    plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, total_height), axes = FALSE, xlab = "", ylab = "")
    col_width <- 0.46; col1_start <- 0.01; col2_start <- 0.47
    title_y <- total_height - 0.3
    text(col1_start, title_y, "A", font = 2, cex = 1.3, adj = c(0, 0.5))
    text(col1_start + 0.03, title_y, expression(italic("M. capitata")~"(host)"),
         font = 1, cex = 1.0, adj = c(0, 0.5))
    text(col2_start, title_y, "B", font = 2, cex = 1.3, adj = c(0, 0.5))
    text(col2_start + 0.03, title_y, expression(italic("D. trenchii")~"(symbiont)"),
         font = 1, cex = 1.0, adj = c(0, 0.5))
    draw_height <- total_height - 1.5
    draw_organism_column(host_data, col1_start, col_width, draw_height, gap_size)
    draw_organism_column(symbiont_data, col2_start, col_width, draw_height, gap_size)
    legend_x <- 0.88; legend_y <- total_height - 0.5; line_spacing <- 0.4
    text(legend_x, legend_y, expression(bold("p < 0.01")), cex = 0.6, adj = c(0, 0.5))
    text(legend_x, legend_y - line_spacing, "p < 0.05", cex = 0.6, adj = c(0, 0.5), font = 1)
    text(legend_x, legend_y - 2*line_spacing, expression(italic("p < 0.1")), cex = 0.6, adj = c(0, 0.5), col = "grey50")
    dev.off()
    cat("  Saved:", basename(output_file), "\n")
}

plot_combined_season("summer", "figures/Fig5_Summer_Representative_GO_Dendrogram.pdf", "Figure 5")
plot_combined_season("winter", "figures/Fig6_Winter_Representative_GO_Dendrogram.pdf", "Figure 6")


# ##############################################################################
# FIGURE 7: Chord Diagram (from script 15)
# ##############################################################################

cat("\n==============================================================================\n")
cat("Figure 7: Chord Diagram - Shared GO Terms\n")
cat("==============================================================================\n\n")

sig_threshold <- 0.05

interaction_colors <- c(
    "Synergistic Up" = "#907AD6", "Synergistic Down" = "#7DCFB6",
    "Antagonistic (Host Up)" = "#F79256", "Antagonistic (Host Down)" = "#6BBF59"
)

division_names <- c("BP" = "Biological Process", "CC" = "Cellular Component", "MF" = "Molecular Function")

load_significant_terms <- function(data_dir, prefix, season, sig_threshold) {
    divisions <- c("BP", "CC", "MF"); all_terms <- list()
    for (div in divisions) {
        if (prefix == "") mwu_file <- file.path(data_dir, paste0(season, "_", div, "_MWU_results.csv"))
        else mwu_file <- file.path(data_dir, paste0(prefix, "_", season, "_", div, "_MWU_results.csv"))
        if (!file.exists(mwu_file)) next
        mwu <- read.table(mwu_file, header = TRUE, stringsAsFactors = FALSE, quote = "\"", fill = TRUE)
        sig_terms <- mwu %>%
            filter(p.adj < sig_threshold) %>%
            mutate(division = div, direction = ifelse(delta.rank > 0, "up", "down"),
                   go_id = term,
                   go_id_clean = sapply(strsplit(as.character(term), ";"), function(x) {
                       id <- x[1]; if (!grepl("^GO:", id)) id <- paste0("GO:", id); return(id)
                   })) %>%
            select(go_id, go_id_clean, name, division, direction, delta.rank, pval, p.adj, nseqs)
        all_terms[[div]] <- sig_terms
    }
    result <- bind_rows(all_terms)
    if (nrow(result) == 0) {
        return(data.frame(go_id = character(), go_id_clean = character(), name = character(),
                          division = character(), direction = character(), delta.rank = numeric(),
                          pval = numeric(), p.adj = numeric(), nseqs = integer()))
    }
    return(result)
}

classify_shared_terms <- function(host_terms, symbiont_terms) {
    if (nrow(host_terms) == 0 || nrow(symbiont_terms) == 0) return(data.frame())
    shared_ids <- intersect(host_terms$go_id, symbiont_terms$go_id)
    if (length(shared_ids) == 0) return(data.frame())
    shared_data <- data.frame(go_id = shared_ids, stringsAsFactors = FALSE) %>%
        left_join(host_terms %>% select(go_id, go_id_clean, name, division, direction, p.adj, nseqs) %>%
                      rename(host_direction = direction, host_division = division,
                             term_name = name, host_padj = p.adj, host_nseqs = nseqs), by = "go_id") %>%
        left_join(symbiont_terms %>% select(go_id, direction, division, p.adj, nseqs) %>%
                      rename(symbiont_direction = direction, symbiont_division = division,
                             symbiont_padj = p.adj, symbiont_nseqs = nseqs), by = "go_id") %>%
        mutate(
            interaction_type = case_when(
                host_direction == "up" & symbiont_direction == "up" ~ "Synergistic Up",
                host_direction == "down" & symbiont_direction == "down" ~ "Synergistic Down",
                host_direction == "up" & symbiont_direction == "down" ~ "Antagonistic (Host Up)",
                host_direction == "down" & symbiont_direction == "up" ~ "Antagonistic (Host Down)"),
            interaction_order = case_when(
                interaction_type == "Synergistic Up" ~ 1, interaction_type == "Synergistic Down" ~ 2,
                interaction_type == "Antagonistic (Host Up)" ~ 3, interaction_type == "Antagonistic (Host Down)" ~ 4),
            division = host_division, total_genes = host_nseqs + symbiont_nseqs
        ) %>% arrange(interaction_order, division, go_id)
    return(shared_data)
}

create_supplementary_table <- function(shared_terms, season, output_dir) {
    table_data <- shared_terms %>%
        select(`GO ID` = go_id_clean, `Term Name` = term_name, `Interaction Type` = interaction_type,
               `Host Division` = host_division, `Host Direction` = host_direction,
               `Host DEGs` = host_nseqs, `Host p.adj` = host_padj,
               `Symbiont Division` = symbiont_division, `Symbiont Direction` = symbiont_direction,
               `Symbiont DEGs` = symbiont_nseqs, `Symbiont p.adj` = symbiont_padj)

    gt_table <- table_data %>% gt() %>%
        tab_header(title = paste0(tools::toTitleCase(season), " - Shared GO Terms Between Host and Symbiont"),
                   subtitle = paste0("n = ", nrow(table_data), " shared significant GO terms (p.adj < ", sig_threshold, ")")) %>%
        tab_spanner(label = "Host (M. capitata)", columns = c(`Host Division`, `Host Direction`, `Host DEGs`, `Host p.adj`)) %>%
        tab_spanner(label = "Symbiont (D. trenchii)", columns = c(`Symbiont Division`, `Symbiont Direction`, `Symbiont DEGs`, `Symbiont p.adj`)) %>%
        fmt_scientific(columns = c(`Host p.adj`, `Symbiont p.adj`), decimals = 2) %>%
        tab_options(table.font.size = px(11), column_labels.font.weight = "bold")

    gtsave(gt_table, file.path(output_dir, paste0("Table_S_", season, "_Shared_GO_Terms.html")))
    gtsave(gt_table, file.path(output_dir, paste0("Table_S_", season, "_Shared_GO_Terms.docx")))
    write.csv(table_data, file.path(output_dir, paste0("Table_S_", season, "_Shared_GO_Terms.csv")), row.names = FALSE)
    cat("  Saved supplementary tables for", season, "\n")
}

prepare_season_data <- function(season) {
    cat("  Loading", season, "data...\n")
    host_terms <- load_significant_terms(host_go_dir, "", season, sig_threshold)
    symbiont_terms <- load_significant_terms(symbiont_go_dir, "symbiont", season, sig_threshold)
    shared_terms <- classify_shared_terms(host_terms, symbiont_terms)
    cat("    Host:", nrow(host_terms), "| Symbiont:", nrow(symbiont_terms), "| Shared:", nrow(shared_terms), "\n")
    if (nrow(shared_terms) > 0) create_supplementary_table(shared_terms, season, "tables")
    return(list(host_terms = host_terms, symbiont_terms = symbiont_terms, shared_terms = shared_terms))
}

draw_chord_diagram <- function(shared_terms, season_label, show_legend = FALSE) {
    if (nrow(shared_terms) == 0) { plot.new(); text(0.5, 0.5, paste0(season_label, "\nNo shared terms"), cex = 1.2); return() }
    n_terms <- nrow(shared_terms)
    div_stats <- shared_terms %>% group_by(division) %>%
        summarise(n_terms = n(), total_width = sum(sqrt(total_genes + 1)), .groups = "drop")
    division_sectors <- c("BP", "CC", "MF")
    division_sectors <- division_sectors[division_sectors %in% div_stats$division]
    go_sectors <- shared_terms$go_id
    all_sectors <- c(go_sectors, division_sectors)
    go_sizes <- setNames(sqrt(shared_terms$total_genes + 1), go_sectors)
    div_sizes <- setNames(div_stats$total_width[match(division_sectors, div_stats$division)], division_sectors)
    sector_sizes <- c(go_sizes, div_sizes)

    circos.clear()
    n_go <- length(go_sectors); n_div <- length(division_sectors)
    go_gap <- ifelse(n_terms > 50, 0.9, 1.9)
    gaps <- c(rep(go_gap, n_go - 1), 20, rep(6, n_div - 1), 20)
    circos.par(start.degree = 270, gap.degree = gaps, track.margin = c(0.002, 0.002),
               cell.padding = c(0, 0, 0, 0), canvas.xlim = c(-0.85, 1.2), canvas.ylim = c(-0.9, 0.9))
    xlim_matrix <- cbind(rep(0, length(all_sectors)), sector_sizes); rownames(xlim_matrix) <- all_sectors
    circos.initialize(factors = factor(all_sectors, levels = all_sectors), xlim = xlim_matrix)

    circos.track(ylim = c(0, 1), track.height = 0.25, bg.border = NA, panel.fun = function(x, y) {
        sector.name <- CELL_META$sector.index
        if (sector.name %in% go_sectors) {
            clean_id <- shared_terms$go_id_clean[shared_terms$go_id == sector.name]
            label_cex <- ifelse(n_terms > 60, 0.58, ifelse(n_terms > 50, 0.45, 0.6))
            circos.text(CELL_META$xcenter, 0.1, as.character(clean_id), facing = "clockwise",
                        niceFacing = TRUE, cex = label_cex, adj = c(0, 0.5), font = 1)
        } else {
            short_names <- c("BP" = "Biological Process", "CC" = "Cellular Component", "MF" = "Molecular Function")
            circos.text(CELL_META$xcenter, 0.5, as.character(short_names[sector.name]),
                        facing = "bending.inside", niceFacing = TRUE, cex = 0.75, font = 2)
        }
    })

    div_positions <- setNames(rep(0, length(division_sectors)), division_sectors)
    for (i in 1:n_terms) {
        term <- shared_terms[i, ]; go_id <- term$go_id; div <- term$division
        chord_col <- adjustcolor(interaction_colors[term$interaction_type], alpha.f = 0.55)
        chord_border <- adjustcolor(interaction_colors[term$interaction_type], alpha.f = 0.75)
        chord_width <- sqrt(term$total_genes + 1); div_pos <- div_positions[[div]]
        circos.link(go_id, c(0, sector_sizes[go_id]), div, c(div_pos, div_pos + chord_width),
                    col = chord_col, border = chord_border, lwd = 0.2)
        div_positions[[div]] <- div_pos + chord_width
    }
    text(0, 0.02, season_label, cex = 1.1, font = 2)
    text(0, -0.12, paste0("n = ", n_terms), cex = 0.65)

    if (show_legend) {
        legend_x <- 0.87; legend_y <- 0.80; legend_spacing <- 0.05
        text(legend_x, legend_y + 0.08, "Interaction Type", cex = 0.6, font = 2, adj = c(0, 0.5))
        short_labels <- c("Synergistic Up" = "Synergy All Up", "Synergistic Down" = "Synergy All Down",
                          "Antagonistic (Host Up)" = "Antag. (Host Up / Sym Down)",
                          "Antagonistic (Host Down)" = "Antag. (Sym. Up / Host Down)")
        for (i in seq_along(interaction_colors)) {
            y_pos <- legend_y - (i - 1) * legend_spacing
            points(legend_x, y_pos, pch = 15, col = interaction_colors[i], cex = 1.3)
            text(legend_x + 0.08, y_pos, short_labels[names(interaction_colors)[i]], cex = 0.5, adj = c(0, 0.5))
        }
    }
    circos.clear()
}

summer_chord_data <- prepare_season_data("summer")
winter_chord_data <- prepare_season_data("winter")

pdf("figures/Fig7_Shared_GO_Chord_Diagram.pdf", width = 8, height = 14)
layout(matrix(c(1, 2), nrow = 2, ncol = 1), heights = c(1, 1))
par(mar = c(0.5, 0.5, 0.5, 0.5))
draw_chord_diagram(summer_chord_data$shared_terms, "Summer", show_legend = TRUE)
par(mar = c(0.5, 0.5, 0.5, 0.5))
draw_chord_diagram(winter_chord_data$shared_terms, "Winter", show_legend = FALSE)
dev.off()

cat("Saved: Fig7_Shared_GO_Chord_Diagram.pdf\n")


# ##############################################################################
# SUMMARY
# ##############################################################################

cat("\n==============================================================================\n")
cat("Summary Statistics\n")
cat("==============================================================================\n\n")

summary_stats <- data.frame(
    Organism = c("Host", "Host", "Symbiont", "Symbiont"),
    Season = c("Summer", "Winter", "Summer", "Winter"),
    Total_Genes = c(nrow(host_summer), nrow(host_winter), nrow(sym_summer), nrow(sym_winter)),
    DEGs_total = c(nrow(host_summer_degs), nrow(host_winter_degs), nrow(sym_summer_degs), nrow(sym_winter_degs)),
    DEGs_up = c(sum(host_summer_degs$log2FoldChange > 0), sum(host_winter_degs$log2FoldChange > 0),
                sum(sym_summer_degs$log2FoldChange > 0), sum(sym_winter_degs$log2FoldChange > 0)),
    DEGs_down = c(sum(host_summer_degs$log2FoldChange < 0), sum(host_winter_degs$log2FoldChange < 0),
                  sum(sym_summer_degs$log2FoldChange < 0), sum(sym_winter_degs$log2FoldChange < 0))
)

write.csv(summary_stats, "data/DEG_summary_statistics.csv", row.names = FALSE)
print(summary_stats)

cat("\n==============================================================================\n")
cat("All Figures Complete!\n")
cat("==============================================================================\n\n")

cat("Main figures:\n")
cat("  Fig 2: figures/Fig2_combined_venn_ridgeline.pdf/png\n")
cat("  Fig 3: figures/Fig3_heatmap_combined.pdf/png\n")
cat("  Fig 4: figures/Fig4_pca_holobiont_combined.pdf/png\n")
cat("  Fig 5: figures/Fig5_Summer_Representative_GO_Dendrogram.pdf\n")
cat("  Fig 6: figures/Fig6_Winter_Representative_GO_Dendrogram.pdf\n")
cat("  Fig 7: figures/Fig7_Shared_GO_Chord_Diagram.pdf\n")
cat("\nSupplementary:\n")
cat("  Fig S1: figures/FigS1_pca_individual_combined.pdf/png\n")
cat("  Tables: tables/Table_S_*_Shared_GO_Terms.html/docx/csv\n")

cat("\n")
sessionInfo()
