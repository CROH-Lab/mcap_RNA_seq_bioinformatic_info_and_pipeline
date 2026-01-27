#!/usr/bin/env Rscript
# ==============================================================================
# Publication-Ready Figures: M. capitata and D. trenchii OA Response
# Version 3 - Fixed chord diagram
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

# VSD matrices
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
# FIGURE 1: Density Ridgeline Plot (±2 log2FC, clean labels)
# ==============================================================================

cat("\n==============================================================================\n")
cat("Figure 1: Density Ridgeline Plot\n")
cat("==============================================================================\n\n")

ridgeline_data <- bind_rows(
    host_summer %>% 
        filter(!is.na(log2FoldChange)) %>%
        mutate(Organism = "M. cap", Season = "S", Group_Label = "M. cap - S") %>%
        select(log2FoldChange, Organism, Season, Group_Label),
    host_winter %>% 
        filter(!is.na(log2FoldChange)) %>%
        mutate(Organism = "M. cap", Season = "W", Group_Label = "M. cap - W") %>%
        select(log2FoldChange, Organism, Season, Group_Label),
    sym_summer %>% 
        filter(!is.na(log2FoldChange)) %>%
        mutate(Organism = "D. tre", Season = "S", Group_Label = "D. tre - S") %>%
        select(log2FoldChange, Organism, Season, Group_Label),
    sym_winter %>% 
        filter(!is.na(log2FoldChange)) %>%
        mutate(Organism = "D. tre", Season = "W", Group_Label = "D. tre - W") %>%
        select(log2FoldChange, Organism, Season, Group_Label)
)

ridgeline_data <- ridgeline_data %>%
    mutate(Group_Label = factor(Group_Label, levels = c(
        "D. tre - W",
        "D. tre - S",
        "M. cap - W",
        "M. cap - S"
    )))

fig1_ridgeline <- ggplot(ridgeline_data, aes(x = log2FoldChange, y = Group_Label, 
                                              fill = interaction(Organism, Season))) +
    geom_density_ridges(alpha = 0.8, scale = 1.5, rel_min_height = 0.01) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray60", linewidth = 0.5) +
    scale_fill_manual(values = c(
        "M. cap.S" = "#E69F00",
        "M. cap.W" = "#F0E442",
        "D. tre.S" = "#56B4E9",
        "D. tre.W" = "#009E73"
    )) +
    scale_x_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 0.5)) +
    labs(
        x = expression(Log[2]~Fold~Change~(OA~vs~Ambient)),
        y = ""
    ) +
    theme_pub +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 12, face = "bold.italic"))

ggsave("figures/Fig1_ridgeline_log2FC.pdf", fig1_ridgeline, width = 10, height = 6, dpi = 300)
ggsave("figures/Fig1_ridgeline_log2FC.png", fig1_ridgeline, width = 10, height = 6, dpi = 300)
cat("Saved: Fig1_ridgeline_log2FC.pdf/png\n")

# ==============================================================================
# FIGURE 2: Venn Diagrams (no titles, color legend only)
# ==============================================================================

cat("\n==============================================================================\n")
cat("Figure 2: Venn Diagrams\n")
cat("==============================================================================\n\n")

host_summer_ids <- host_summer_degs$gene_id
host_winter_ids <- host_winter_degs$gene_id
host_overlap <- length(intersect(host_summer_ids, host_winter_ids))

pdf("figures/Fig2A_host_venn.pdf", width = 6, height = 6)
grid.newpage()
venn_host <- draw.pairwise.venn(
    area1 = length(host_summer_ids),
    area2 = length(host_winter_ids),
    cross.area = host_overlap,
    category = c("", ""),
    fill = c(summer_color, winter_color),
    alpha = 0.6,
    col = c(summer_color, winter_color),
    lwd = 2,
    cex = 2,
    fontface = "bold",
    cat.cex = 0,
    margin = 0.1
)
dev.off()

sym_summer_ids <- sym_summer_degs$gene_id
sym_winter_ids <- sym_winter_degs$gene_id
sym_overlap <- length(intersect(sym_summer_ids, sym_winter_ids))

pdf("figures/Fig2B_symbiont_venn.pdf", width = 6, height = 6)
grid.newpage()
venn_sym <- draw.pairwise.venn(
    area1 = length(sym_summer_ids),
    area2 = length(sym_winter_ids),
    cross.area = sym_overlap,
    category = c("", ""),
    fill = c(summer_color, winter_color),
    alpha = 0.6,
    col = c(summer_color, winter_color),
    lwd = 2,
    cex = 2,
    fontface = "bold",
    cat.cex = 0,
    margin = 0.1
)
dev.off()

legend_df <- data.frame(
    Season = c("Summer", "Winter"),
    color = c(summer_color, winter_color)
)

legend_plot <- ggplot(legend_df, aes(x = 1, y = Season, fill = Season)) +
    geom_tile(width = 0.3, height = 0.8) +
    scale_fill_manual(values = c("Summer" = summer_color, "Winter" = winter_color)) +
    theme_void() +
    theme(legend.position = "right",
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12))

ggsave("figures/Fig2_legend.pdf", legend_plot, width = 3, height = 2)

cat("Saved: Fig2A_host_venn.pdf, Fig2B_symbiont_venn.pdf, Fig2_legend.pdf\n")
cat("  Host overlap:", host_overlap, "DEGs\n")
cat("  Symbiont overlap:", sym_overlap, "DEGs\n")

# ==============================================================================
# FIGURE 3: Heatmaps (ALL DEGs, no title)
# ==============================================================================

cat("\n==============================================================================\n")
cat("Figure 3: Heatmaps (All DEGs)\n")
cat("==============================================================================\n\n")

create_deg_heatmap_all <- function(vsd_mat, deseq_results, filename, max_genes = 500) {
    
    all_degs <- deseq_results %>%
        filter(!is.na(padj) & padj < 0.05) %>%
        arrange(desc(abs(log2FoldChange)))
    
    if (nrow(all_degs) < 5) {
        cat("  Warning: Fewer than 5 DEGs, skipping\n")
        return(NULL)
    }
    
    if (nrow(all_degs) > max_genes) {
        cat("  Note: Limiting to top", max_genes, "DEGs by |log2FC| for visualization\n")
        all_degs <- head(all_degs, max_genes)
    }
    
    common_genes <- intersect(all_degs$gene_id, rownames(vsd_mat))
    if (length(common_genes) < 5) {
        cat("  Warning: Fewer than 5 matching genes, skipping\n")
        return(NULL)
    }
    
    heatmap_mat <- vsd_mat[common_genes, , drop = FALSE]
    heatmap_mat_scaled <- t(scale(t(heatmap_mat)))
    
    sample_names <- colnames(heatmap_mat)
    treatments <- ifelse(grepl("B", sample_names), "OA", "Ambient")
    
    anno_col <- data.frame(
        Treatment = treatments,
        row.names = sample_names
    )
    
    anno_colors <- list(Treatment = c("OA" = "#E69F00", "Ambient" = "#56B4E9"))
    
    plot_height <- max(8, min(20, nrow(heatmap_mat_scaled) / 20))
    
    pdf(filename, width = 8, height = plot_height)
    pheatmap::pheatmap(heatmap_mat_scaled,
             color = colorRampPalette(c(down_color, "white", up_color))(100),
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             show_rownames = FALSE,
             show_colnames = TRUE,
             annotation_col = anno_col,
             annotation_colors = anno_colors,
             fontsize = 10,
             border_color = NA,
             main = "")
    dev.off()
    
    cat("  Saved:", filename, "with", length(common_genes), "DEGs\n")
    return(length(common_genes))
}

cat("Creating heatmaps with ALL DEGs...\n")

create_deg_heatmap_all(host_summer_vsd, host_summer, "figures/Fig3A_heatmap_host_summer.pdf")
create_deg_heatmap_all(host_winter_vsd, host_winter, "figures/Fig3B_heatmap_host_winter.pdf")
create_deg_heatmap_all(sym_summer_vsd, sym_summer, "figures/Fig3C_heatmap_symbiont_summer.pdf")
create_deg_heatmap_all(sym_winter_vsd, sym_winter, "figures/Fig3D_heatmap_symbiont_winter.pdf")

# ==============================================================================
# FIGURE 4: Chord Diagram - Calcification/Ion Transport
# ==============================================================================

cat("\n==============================================================================\n")
cat("Figure 4: Chord Diagram - Calcification/Ion Transport\n")
cat("==============================================================================\n\n")

calc_keywords <- c("calcif", "calcium", "carbonate", "bicarbon", "ion.transport", 
                   "homeostasis", "channel", "pump", "exchanger", "ATPase", 
                   "voltage.gated", "sodium", "potassium", "chloride", "proton",
                   "pH", "acid", "carbonic.anhydrase", "SLC", "transporter")

keyword_pattern <- paste(calc_keywords, collapse = "|")

# Function to extract category - returns character, not list
extract_category <- function(desc) {
    if (is.na(desc) || desc == "") return("Other Transport")
    desc_lower <- tolower(desc)
    if (grepl("calcium.channel|voltage.gated.*calcium|voltage.dependent.*calcium", desc_lower)) return("Ca2+ Channel")
    if (grepl("calcium.pump|calcium.transport.*atpase|plasma membrane calcium", desc_lower)) return("Ca2+ ATPase")
    if (grepl("calcium.bind|calmodulin|ef.hand|calcium.dependent.protein", desc_lower)) return("Ca2+ Signaling")
    if (grepl("sodium.*potassium|na\\+.*k\\+|nka|sodium.potassium", desc_lower)) return("Na+/K+ ATPase")
    if (grepl("sodium|na\\+", desc_lower)) return("Na+ Transport")
    if (grepl("potassium|k\\+", desc_lower)) return("K+ Channel")
    if (grepl("chloride|anion|clc", desc_lower)) return("Cl- Transport")
    if (grepl("carbonic.anhydrase", desc_lower)) return("Carbonic Anhydrase")
    if (grepl("exchanger|antiporter|nhe|ncx", desc_lower)) return("Ion Exchanger")
    if (grepl("proton|h\\+.*atpase|v.type|vacuolar", desc_lower)) return("H+ ATPase")
    if (grepl("bicarbonate|hco3|slc4|slc26", desc_lower)) return("HCO3- Transport")
    if (grepl("aquaporin|water.channel", desc_lower)) return("Aquaporin")
    return("Other Transport")
}

# Process each group separately and ensure category is character
process_calc_degs <- function(deseq_degs, annot_df, organism_name, season_name) {
    
    if (nrow(deseq_degs) == 0) {
        return(data.frame())
    }
    
    result <- deseq_degs %>%
        left_join(annot_df, by = "gene_id") %>%
        filter(grepl(keyword_pattern, description, ignore.case = TRUE))
    
    if (nrow(result) == 0) {
        return(data.frame())
    }
    
    result <- result %>%
        mutate(
            organism = organism_name,
            season = season_name,
            # Use vapply to ensure character output
            category = vapply(description, extract_category, character(1)),
            short_name = ifelse(!is.na(uniprot_entry), uniprot_entry, 
                               sub(".*g(\\d+)\\.?t?.*", "g\\1", gene_id))
        )
    
    return(result)
}

host_calc_summer <- process_calc_degs(host_summer_degs, host_annot, "Host", "Summer")
host_calc_winter <- process_calc_degs(host_winter_degs, host_annot, "Host", "Winter")
sym_calc_summer <- process_calc_degs(sym_summer_degs, sym_annot, "Symbiont", "Summer")
sym_calc_winter <- process_calc_degs(sym_winter_degs, sym_annot, "Symbiont", "Winter")

cat("Calcification-related DEGs:\n")
cat("  Host Summer:", nrow(host_calc_summer), "\n")
cat("  Host Winter:", nrow(host_calc_winter), "\n")
cat("  Symbiont Summer:", nrow(sym_calc_summer), "\n")
cat("  Symbiont Winter:", nrow(sym_calc_winter), "\n")

# Combine only non-empty dataframes
calc_dfs <- list(host_calc_summer, host_calc_winter, sym_calc_summer, sym_calc_winter)
calc_dfs <- calc_dfs[sapply(calc_dfs, function(x) nrow(x) > 0)]

if (length(calc_dfs) > 0) {
    all_calc_degs <- bind_rows(calc_dfs)
    cat("  Total:", nrow(all_calc_degs), "\n")
    
    # Save data
    write.csv(all_calc_degs, "data/all_calcification_degs.csv", row.names = FALSE)
    
    if (nrow(all_calc_degs) >= 5) {
        
        # Order: Host first (by abs log2FC descending), then Symbiont
        all_calc_degs <- all_calc_degs %>%
            arrange(desc(organism == "Host"), desc(abs(log2FoldChange)))
        
        # Make unique short names by appending organism prefix
        all_calc_degs$gene_label <- paste0(
            ifelse(all_calc_degs$organism == "Host", "H:", "S:"),
            all_calc_degs$short_name
        )
        
        # Handle duplicate labels
        all_calc_degs <- all_calc_degs %>%
            group_by(gene_label) %>%
            mutate(gene_label = ifelse(n() > 1, 
                                        paste0(gene_label, ".", row_number()),
                                        gene_label)) %>%
            ungroup()
        
        # Get unique categories present in data
        categories <- unique(all_calc_degs$category)
        cat("  Categories:", paste(categories, collapse = ", "), "\n")
        
        # Create data frame for chord connections
        chord_df <- all_calc_degs %>%
            select(gene_label, category, log2FoldChange, season, organism) %>%
            as.data.frame()
        
        # Create the chord diagram
        pdf("figures/Fig4_chord_calcification.pdf", width = 16, height = 14)
        
        circos.clear()
        circos.par(start.degree = 90, gap.after = c(rep(2, nrow(chord_df)-1), 15,
                                                     rep(2, length(categories)-1), 15))
        
        # Create color mapping for log2FC
        lfc_range <- range(chord_df$log2FoldChange, na.rm = TRUE)
        lfc_col_fun <- colorRamp2(c(min(-2, lfc_range[1]), 0, max(2, lfc_range[2])), 
                                   c(down_color, "white", up_color))
        
        # Season colors for links
        season_cols <- c("Summer" = summer_color, "Winter" = winter_color)
        
        # Category colors
        cat_cols <- setNames(
            colorRampPalette(brewer.pal(min(length(categories), 8), "Set2"))(length(categories)),
            categories
        )
        
        # Build the adjacency list for chordDiagram
        # Format: from, to, value
        adj_list <- chord_df %>%
            mutate(value = 1) %>%
            select(from = gene_label, to = category, value)
        
        # Set sector colors: genes colored by organism, categories by Set2
        gene_colors <- setNames(
            ifelse(chord_df$organism == "Host", host_color, symbiont_color),
            chord_df$gene_label
        )
        
        all_sector_colors <- c(gene_colors, cat_cols)
        
        # Create link colors based on season
        link_colors <- sapply(1:nrow(chord_df), function(i) {
            adjustcolor(season_cols[chord_df$season[i]], alpha.f = 0.6)
        })
        
        # Draw chord diagram
        chordDiagram(
            adj_list,
            order = c(chord_df$gene_label, categories),
            grid.col = all_sector_colors,
            col = link_colors,
            transparency = 0.3,
            annotationTrack = c("grid"),
            preAllocateTracks = list(
                track.height = 0.15,
                track.margin = c(0.01, 0)
            ),
            annotationTrackHeight = c(0.03)
        )
        
        # Add gene labels (rotated text)
        circos.track(track.index = 1, panel.fun = function(x, y) {
            sector.name <- get.cell.meta.data("sector.index")
            xlim <- get.cell.meta.data("xlim")
            ylim <- get.cell.meta.data("ylim")
            
            # Only label genes, not categories
            if (sector.name %in% chord_df$gene_label) {
                circos.text(mean(xlim), ylim[1] + 0.5, sector.name,
                           facing = "clockwise", niceFacing = TRUE,
                           adj = c(0, 0.5), cex = 0.4)
            } else {
                # Category labels - larger and bending
                circos.text(mean(xlim), ylim[1] + 0.5, sector.name,
                           facing = "bending.inside", niceFacing = TRUE,
                           adj = c(0.5, 0), cex = 0.8, font = 2)
            }
        }, bg.border = NA)
        
        # Add log2FC color track for genes only
        circos.track(track.index = 2, panel.fun = function(x, y) {
            sector.name <- get.cell.meta.data("sector.index")
            if (sector.name %in% chord_df$gene_label) {
                lfc <- chord_df$log2FoldChange[chord_df$gene_label == sector.name][1]
                xlim <- get.cell.meta.data("xlim")
                circos.rect(xlim[1], 0, xlim[2], 1, col = lfc_col_fun(lfc), border = NA)
            }
        }, ylim = c(0, 1), track.height = 0.05, bg.border = NA)
        
        # Add legends
        lgd_lfc <- Legend(
            col_fun = lfc_col_fun,
            title = "Log2FC",
            title_position = "topleft",
            legend_height = unit(3, "cm")
        )
        
        lgd_season <- Legend(
            labels = c("Summer", "Winter"),
            legend_gp = gpar(fill = c(summer_color, winter_color)),
            title = "Season",
            title_position = "topleft"
        )
        
        lgd_org <- Legend(
            labels = c("Host (H:)", "Symbiont (S:)"),
            legend_gp = gpar(fill = c(host_color, symbiont_color)),
            title = "Organism",
            title_position = "topleft"
        )
        
        draw(lgd_lfc, x = unit(0.92, "npc"), y = unit(0.8, "npc"))
        draw(lgd_season, x = unit(0.92, "npc"), y = unit(0.55, "npc"))
        draw(lgd_org, x = unit(0.92, "npc"), y = unit(0.35, "npc"))
        
        circos.clear()
        dev.off()
        
        cat("Saved: Fig4_chord_calcification.pdf\n")
        
        # Save category summary
        cat_summary <- all_calc_degs %>%
            group_by(category, organism, season) %>%
            summarise(count = n(), avg_lfc = mean(log2FoldChange, na.rm = TRUE), .groups = "drop")
        
        write.csv(cat_summary, "data/calcification_category_summary.csv", row.names = FALSE)
        
    } else {
        cat("Fewer than 5 calcification DEGs - skipping chord diagram\n")
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
    
    pca_result <- prcomp(t(vsd_mat), scale. = FALSE)
    percent_var <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)
    
    pca_data <- data.frame(
        PC1 = pca_result$x[, 1],
        PC2 = pca_result$x[, 2],
        sample = rownames(pca_result$x)
    )
    
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

pca_host_summer <- create_pca_from_matrix(host_summer_vsd, "M. cap - Summer")
pca_host_winter <- create_pca_from_matrix(host_winter_vsd, "M. cap - Winter")
pca_sym_summer <- create_pca_from_matrix(sym_summer_vsd, "D. tre - Summer")
pca_sym_winter <- create_pca_from_matrix(sym_winter_vsd, "D. tre - Winter")

ggsave("figures/Fig6A_pca_host_summer.pdf", pca_host_summer, width = 7, height = 6)
ggsave("figures/Fig6B_pca_host_winter.pdf", pca_host_winter, width = 7, height = 6)
ggsave("figures/Fig6C_pca_symbiont_summer.pdf", pca_sym_summer, width = 7, height = 6)
ggsave("figures/Fig6D_pca_symbiont_winter.pdf", pca_sym_winter, width = 7, height = 6)

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

cat("\n==============================================================================\n")
cat("Analysis Complete!\n")
cat("==============================================================================\n\n")

sessionInfo()
