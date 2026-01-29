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
# FIGURE 1: Density Ridgeline Plot
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
        "D. tre - W", "D. tre - S", "M. cap - W", "M. cap - S"
    )))

fig1_ridgeline <- ggplot(ridgeline_data, aes(x = log2FoldChange, y = Group_Label,
                                              fill = interaction(Organism, Season))) +
    geom_density_ridges(alpha = 0.8, scale = 1.5, rel_min_height = 0.01) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray60", linewidth = 0.5) +
    scale_fill_manual(values = c(
        "M. cap.S" = "#E69F00", "M. cap.W" = "#F0E442",
        "D. tre.S" = "#56B4E9", "D. tre.W" = "#009E73"
    )) +
    scale_x_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 0.5)) +
    labs(x = expression(Log[2]~Fold~Change~(OA~vs~Ambient)), y = "") +
    theme_pub +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 12, face = "bold.italic"))

ggsave("figures/Fig1_ridgeline_log2FC.pdf", fig1_ridgeline, width = 10, height = 6, dpi = 300)
ggsave("figures/Fig1_ridgeline_log2FC.png", fig1_ridgeline, width = 10, height = 6, dpi = 300)
cat("Saved: Fig1_ridgeline_log2FC.pdf/png\n")

# ==============================================================================
# FIGURE 2: Venn Diagrams
# ==============================================================================

cat("\n==============================================================================\n")
cat("Figure 2: Venn Diagrams\n")
cat("==============================================================================\n\n")

host_summer_ids <- host_summer_degs$gene_id
host_winter_ids <- host_winter_degs$gene_id
host_overlap <- length(intersect(host_summer_ids, host_winter_ids))

pdf("figures/Fig2A_host_venn.pdf", width = 6, height = 6)
grid.newpage()
draw.pairwise.venn(
    area1 = length(host_summer_ids), area2 = length(host_winter_ids),
    cross.area = host_overlap, category = c("", ""),
    fill = c(summer_color, winter_color), alpha = 0.6,
    col = c(summer_color, winter_color), lwd = 2, cex = 2,
    fontface = "bold", cat.cex = 0, margin = 0.1
)
dev.off()

sym_summer_ids <- sym_summer_degs$gene_id
sym_winter_ids <- sym_winter_degs$gene_id
sym_overlap <- length(intersect(sym_summer_ids, sym_winter_ids))

pdf("figures/Fig2B_symbiont_venn.pdf", width = 6, height = 6)
grid.newpage()
draw.pairwise.venn(
    area1 = length(sym_summer_ids), area2 = length(sym_winter_ids),
    cross.area = sym_overlap, category = c("", ""),
    fill = c(summer_color, winter_color), alpha = 0.6,
    col = c(summer_color, winter_color), lwd = 2, cex = 2,
    fontface = "bold", cat.cex = 0, margin = 0.1
)
dev.off()

legend_plot <- ggplot(data.frame(Season = c("Summer", "Winter")),
                      aes(x = 1, y = Season, fill = Season)) +
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
# FIGURE 3: Heatmaps
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
        cat("  Note: Limiting to top", max_genes, "DEGs\n")
        all_degs <- head(all_degs, max_genes)
    }

    common_genes <- intersect(all_degs$gene_id, rownames(vsd_mat))
    if (length(common_genes) < 5) return(NULL)

    heatmap_mat <- vsd_mat[common_genes, , drop = FALSE]
    heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

    sample_names <- colnames(heatmap_mat)
    treatments <- ifelse(grepl("B", sample_names), "OA", "Ambient")
    anno_col <- data.frame(Treatment = treatments, row.names = sample_names)
    anno_colors <- list(Treatment = c("OA" = "#E69F00", "Ambient" = "#56B4E9"))

    plot_height <- max(8, min(20, nrow(heatmap_mat_scaled) / 20))

    pdf(filename, width = 8, height = plot_height)
    pheatmap::pheatmap(heatmap_mat_scaled,
             color = colorRampPalette(c(down_color, "white", up_color))(100),
             cluster_rows = TRUE, cluster_cols = TRUE,
             show_rownames = FALSE, show_colnames = TRUE,
             annotation_col = anno_col, annotation_colors = anno_colors,
             fontsize = 10, border_color = NA, main = "")
    dev.off()

    cat("  Saved:", filename, "with", length(common_genes), "DEGs\n")
    return(length(common_genes))
}

cat("Creating heatmaps...\n")
create_deg_heatmap_all(host_summer_vsd, host_summer, "figures/Fig3A_heatmap_host_summer.pdf")
create_deg_heatmap_all(host_winter_vsd, host_winter, "figures/Fig3B_heatmap_host_winter.pdf")
create_deg_heatmap_all(sym_summer_vsd, sym_summer, "figures/Fig3C_heatmap_symbiont_summer.pdf")
create_deg_heatmap_all(sym_winter_vsd, sym_winter, "figures/Fig3D_heatmap_symbiont_winter.pdf")

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
# FIGURE 4D-E: GO_MWU Bubble Plots - Host and Symbiont
# Y-axis: -log10(p.adj), X-axis: delta.rank (direction), Size: nseqs
# Faceted by Division (BP/MF/CC) × Season (Summer/Winter)
# Only annotate terms matching calcification/transport keywords
# ==============================================================================

cat("\n==============================================================================\n")
cat("Figure 4D-E: GO_MWU Bubble Plots\n")
cat("==============================================================================\n\n")

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Keywords for annotation (case-insensitive)
ANNOTATION_KEYWORDS <- c(
    "ion", "transport", "atpase", "proton", "calcium", "carbon", 
    "channel", "acid", "chemosensory", "respiration", "membrane", 
    "potential", "ossification", "adhesion", "homeostasis", "reticulum"
)

# Significance cutoff for plotting (orange line)
BUBBLE_P_ADJ_CUTOFF <- 0.05

# Colors matching reference image (BP=green, CC=red, MF=blue)
division_colors <- c(
    "BP" = "#4DAF4A",  # Green
    "CC" = "#E41A1C",  # Red  
    "MF" = "#377EB8"   # Blue
)

# -----------------------------------------------------------------------------
# Load GO_MWU Results
# -----------------------------------------------------------------------------

load_gomwu_bubble_results <- function(base_dir, prefix = "", organism_label) {
    
    seasons <- c("summer", "winter")
    divisions <- c("BP", "MF", "CC")
    
    all_results <- list()
    
    for (season in seasons) {
        for (div in divisions) {
            
            # Construct filename
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
                
                # Clean column names
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

cat("Loading Host GO_MWU results...\n")
host_gomwu <- load_gomwu_bubble_results("../10_GO_MWU", prefix = "", "Host")
cat("  Loaded", nrow(host_gomwu), "GO terms\n")

cat("Loading Symbiont GO_MWU results...\n")
symbiont_gomwu <- load_gomwu_bubble_results("../11_symbiont_GO_MWU", prefix = "symbiont", "Symbiont")
cat("  Loaded", nrow(symbiont_gomwu), "GO terms\n")

# -----------------------------------------------------------------------------
# Process Data for Plotting
# -----------------------------------------------------------------------------

process_for_bubble <- function(gomwu_data, p_cutoff = 0.05) {
    
    df <- gomwu_data %>%
        filter(!is.na(p.adj)) %>%
        mutate(
            # Y-axis: -log10(p.adj), cap at reasonable value
            neg_log10_padj = pmin(-log10(p.adj), 30),
            
            # X-axis: delta.rank (already signed for direction)
            # Normalize to z-score-like scale for better visualization
            x_value = delta.rank / 100,  # Scale down for readability
            
            # Clean up GO term name
            go_name = name,
            
            # Season labels
            season_label = factor(
                ifelse(season == "summer", "Summer", "Winter"),
                levels = c("Summer", "Winter")
            ),
            
            # Division factor
            division = factor(division, levels = c("BP", "CC", "MF"))
        )
    
    # Check if term matches any annotation keyword
    keyword_pattern <- paste(ANNOTATION_KEYWORDS, collapse = "|")
    df <- df %>%
        mutate(
            matches_keyword = str_detect(tolower(go_name), regex(keyword_pattern, ignore_case = TRUE)),
            # Only annotate if significant AND matches keyword
            annotate_label = ifelse(p.adj < p_cutoff & matches_keyword, 
                                    str_trunc(go_name, 35), NA)
        )
    
    return(df)
}

host_bubble_data <- process_for_bubble(host_gomwu, BUBBLE_P_ADJ_CUTOFF)
symbiont_bubble_data <- process_for_bubble(symbiont_gomwu, BUBBLE_P_ADJ_CUTOFF)

cat("\nHost data for plotting:", nrow(host_bubble_data), "terms\n")
cat("  Significant (p.adj <", BUBBLE_P_ADJ_CUTOFF, "):", 
    sum(host_bubble_data$p.adj < BUBBLE_P_ADJ_CUTOFF, na.rm = TRUE), "\n")
cat("  To be annotated:", sum(!is.na(host_bubble_data$annotate_label)), "\n")

cat("\nSymbiont data for plotting:", nrow(symbiont_bubble_data), "terms\n")
cat("  Significant (p.adj <", BUBBLE_P_ADJ_CUTOFF, "):", 
    sum(symbiont_bubble_data$p.adj < BUBBLE_P_ADJ_CUTOFF, na.rm = TRUE), "\n")
cat("  To be annotated:", sum(!is.na(symbiont_bubble_data$annotate_label)), "\n")

# -----------------------------------------------------------------------------
# Create Bubble Plot Function
# -----------------------------------------------------------------------------

create_gomwu_bubble <- function(bubble_data, organism_name, org_color) {
    
    # Significance threshold line
    sig_line <- -log10(BUBBLE_P_ADJ_CUTOFF)
    
    p <- ggplot(bubble_data, aes(x = x_value, y = neg_log10_padj)) +
        
        # Points - all terms
        geom_point(aes(size = nseqs, fill = division), 
                   shape = 21, alpha = 0.6, color = "gray30", stroke = 0.3) +
        
        # Significance threshold line
        geom_hline(yintercept = sig_line, color = "darkorange", 
                   linetype = "solid", linewidth = 0.8) +
        
        # Vertical line at 0 (no change)
        geom_vline(xintercept = 0, color = "gray50", 
                   linetype = "dashed", linewidth = 0.5) +
        
        # Labels for keyword-matching significant terms
        geom_text_repel(
            aes(label = annotate_label),
            size = 2.5,
            max.overlaps = 20,
            box.padding = 0.3,
            point.padding = 0.2,
            segment.color = "gray50",
            segment.size = 0.3,
            min.segment.length = 0.2,
            na.rm = TRUE
        ) +
        
        # Color scale
        scale_fill_manual(
            values = division_colors,
            name = "GO Division",
            labels = c("BP" = "Biological Process", 
                      "CC" = "Cellular Component", 
                      "MF" = "Molecular Function")
        ) +
        
        # Size scale
        scale_size_continuous(
            name = "# Genes",
            range = c(1, 12),
            breaks = c(10, 50, 100, 200, 500)
        ) +
        
        # Facet: columns = division, rows = season
        facet_grid(season_label ~ division, scales = "free_x") +
        
        # Labels
        labs(
            x = "Delta Rank (Direction Score)",
            y = expression(-log[10]~italic(p)[adj]),
            title = paste0("GO Enrichment: ", organism_name)
        ) +
        
        # Theme
        theme_bw(base_size = 11) +
        theme(
            plot.title = element_text(size = 14, face = "bold.italic", 
                                      hjust = 0.5, color = org_color),
            axis.title = element_text(size = 11, face = "bold"),
            axis.text = element_text(size = 9),
            strip.text = element_text(size = 10, face = "bold"),
            strip.background = element_rect(fill = "gray95"),
            legend.position = "right",
            legend.box = "vertical",
            panel.grid.minor = element_blank(),
            panel.spacing = unit(0.8, "lines")
        ) +
        
        # Guides
        guides(
            fill = guide_legend(override.aes = list(size = 5), order = 1),
            size = guide_legend(order = 2)
        )
    
    return(p)
}

# -----------------------------------------------------------------------------
# Generate Plots
# -----------------------------------------------------------------------------

cat("\nGenerating Host bubble plot...\n")
fig4d_host_bubble <- create_gomwu_bubble(
    host_bubble_data, 
    "M. capitata (Host)", 
    host_color
)

cat("Generating Symbiont bubble plot...\n")
fig4e_symbiont_bubble <- create_gomwu_bubble(
    symbiont_bubble_data, 
    "D. trenchii (Symbiont)", 
    symbiont_color
)

# -----------------------------------------------------------------------------
# Save Plots
# -----------------------------------------------------------------------------

ggsave("figures/Fig4D_bubble_host_GOMWU.pdf", fig4d_host_bubble, 
       width = 14, height = 8, dpi = 300)
ggsave("figures/Fig4D_bubble_host_GOMWU.png", fig4d_host_bubble, 
       width = 14, height = 8, dpi = 300)
cat("✓ Saved: Fig4D_bubble_host_GOMWU.pdf/png\n")

ggsave("figures/Fig4E_bubble_symbiont_GOMWU.pdf", fig4e_symbiont_bubble, 
       width = 14, height = 8, dpi = 300)
ggsave("figures/Fig4E_bubble_symbiont_GOMWU.png", fig4e_symbiont_bubble, 
       width = 14, height = 8, dpi = 300)
cat("✓ Saved: Fig4E_bubble_symbiont_GOMWU.pdf/png\n")

# -----------------------------------------------------------------------------
# Save Summary Statistics
# -----------------------------------------------------------------------------

# Summary of annotated terms
annotated_summary <- bind_rows(
    host_bubble_data %>% 
        filter(!is.na(annotate_label)) %>%
        select(organism, season, division, go_name, p.adj, delta.rank, nseqs),
    symbiont_bubble_data %>%
        filter(!is.na(annotate_label)) %>%
        select(organism, season, division, go_name, p.adj, delta.rank, nseqs)
) %>%
    arrange(organism, season, division, p.adj)

write.csv(annotated_summary, "data/bubble_plot_annotated_terms.csv", row.names = FALSE)
cat("✓ Saved: bubble_plot_annotated_terms.csv\n")

# Print summary
cat("\n=== Annotated Terms Summary ===\n")
print(table(annotated_summary$organism, annotated_summary$season))

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
# FIGURE 6: PCA Plots
# ==============================================================================

cat("\n==============================================================================\n")
cat("Figure 6: PCA Plots\n")
cat("==============================================================================\n\n")

create_pca_from_matrix <- function(vsd_mat, title_text) {
    pca_result <- prcomp(t(vsd_mat), scale. = FALSE)
    percent_var <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)

    pca_data <- data.frame(
        PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2],
        sample = rownames(pca_result$x),
        Treatment = ifelse(grepl("B", rownames(pca_result$x)), "OA", "Ambient")
    )

    ggplot(pca_data, aes(x = PC1, y = PC2, color = Treatment, shape = Treatment)) +
        geom_point(size = 4, alpha = 0.8) +
        stat_ellipse(level = 0.95, linetype = "dashed") +
        scale_color_manual(values = c("OA" = "#E69F00", "Ambient" = "#56B4E9"), name = "Treatment") +
        scale_shape_manual(values = c("OA" = 16, "Ambient" = 17), name = "Treatment") +
        labs(x = paste0("PC1: ", percent_var[1], "% variance"),
             y = paste0("PC2: ", percent_var[2], "% variance"),
             title = title_text) +
        theme_pub + theme(legend.position = "right")
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
    ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), label_size = 16
)
pca_legend <- get_legend(pca_host_summer + theme(legend.position = "bottom"))
pca_combined_legend <- plot_grid(pca_combined, pca_legend, ncol = 1, rel_heights = c(1, 0.1))

ggsave("figures/Fig6_pca_combined.pdf", pca_combined_legend, width = 12, height = 10)
ggsave("figures/Fig6_pca_combined.png", pca_combined_legend, width = 12, height = 10, dpi = 300)
cat("Saved: Fig6_pca_*.pdf/png\n")

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
