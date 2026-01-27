#!/usr/bin/env Rscript
# ==============================================================================
# Refined Chord Diagram - Calcification/Ion Transport DEGs
# ==============================================================================

setwd("/home/darmstrong4/mc_rework/12_publication_figures")

suppressPackageStartupMessages({
    library(tidyverse)
    library(circlize)
    library(RColorBrewer)
    library(ComplexHeatmap)
})

# Color palettes
host_color <- "#E69F00"
symbiont_color <- "#56B4E9"
summer_color <- "#D55E00"
winter_color <- "#0072B2"
up_color <- "#D73027"
down_color <- "#4575B4"

cat("Loading data...\n")

# Load the full calcification DEGs
calc_degs <- read.csv("data/all_calcification_degs.csv")

cat("Original DEGs:", nrow(calc_degs), "\n")

# ==============================================================================
# REFINED FILTERING - Only truly relevant genes
# ==============================================================================

# Define truly relevant descriptions for calcification
calcification_relevant <- c(
    "carbonic anhydrase",
    "calcium.*channel", "calcium.*atpase", "calcium.*pump", "calcium.*transport",
    "calcium.*dependent.*kinase", "calcium.*binding",
    "skeletal", "bone morphogenetic", "collagen",
    "alkaline phosphatase",
    "transient receptor potential", "TRP",
    "sodium.*channel", "potassium.*channel", "chloride.*channel",
    "monocarboxylate transporter",
    "phosphate.*translocator", "sugar phosphate",
    "bicarbonate", "HCO3",
    "exchanger", "antiporter",
    "aquaporin"
)

# False positives to EXCLUDE
false_positives <- c(
    "ubiquitin", "ligase",
    "glyceraldehyde", "GAPDH",
    "tubulin", "actin",
    "phosphatase.*regulatory", "protein phosphatase [0-9]",
    "steroid", "hydroxylase",
    "methyltransferase",
    "lyase",
    "dehydrogenase",
    "phosphodiesterase",
    "phospholipase",
    "GABA.*receptor",
    "adrenergic receptor",
    "lymphotoxin", "laminin",
    "mannose.*receptor",
    "mannosyl",
    "josephin",
    "endonuclease", "exonuclease",
    "peridinin", "chlorophyll", "fucoxanthin", "photosystem",
    "elongation factor",
    "amylase",
    "aminotransferase"
)

relevant_pattern <- paste(calcification_relevant, collapse = "|")
exclude_pattern <- paste(false_positives, collapse = "|")

# Filter for truly relevant genes
truly_relevant <- calc_degs %>%
    filter(
        grepl(relevant_pattern, description, ignore.case = TRUE) |
        category %in% c("Ca2+ Channel", "Ca2+ ATPase", "Ca2+ Signaling", 
                        "Carbonic Anhydrase", "H+ ATPase", "HCO3- Transport",
                        "Na+/K+ ATPase", "Ion Exchanger", "Aquaporin")
    ) %>%
    filter(!grepl(exclude_pattern, description, ignore.case = TRUE))

cat("Refined DEGs:", nrow(truly_relevant), "\n")

# Reclassify with refined categories
truly_relevant <- truly_relevant %>%
    mutate(refined_category = case_when(
        grepl("skeletal|SARP", description, ignore.case = TRUE) ~ "Biomineralization Matrix",
        grepl("alkaline phosphatase", description, ignore.case = TRUE) ~ "Biomineralization Enzyme",
        grepl("bone morphogenetic|BMP", description, ignore.case = TRUE) ~ "BMP Signaling",
        grepl("collagen", description, ignore.case = TRUE) ~ "Collagen/ECM",
        grepl("transient receptor potential|TRP", description, ignore.case = TRUE) ~ "TRP Channel",
        grepl("monocarboxylate", description, ignore.case = TRUE) ~ "MCT (pH Regulation)",
        grepl("zinc transporter|ZIP", description, ignore.case = TRUE) ~ "Zn2+ Transport",
        grepl("sugar phosphate|phosphate.*translocator", description, ignore.case = TRUE) ~ "Phosphate Transport",
        grepl("sodium.*channel", description, ignore.case = TRUE) ~ "Na+ Channel",
        grepl("potassium.*channel|potassium.*voltage", description, ignore.case = TRUE) ~ "K+ Channel",
        grepl("chloride.*channel", description, ignore.case = TRUE) ~ "Cl- Channel",
        grepl("ABC transporter", description, ignore.case = TRUE) ~ "ABC Transporter",
        grepl("carbonic anhydrase", description, ignore.case = TRUE) ~ "Carbonic Anhydrase",
        grepl("calcium.*channel|voltage.*calcium", description, ignore.case = TRUE) ~ "Ca2+ Channel",
        grepl("calcium.*atpase|calcium.*pump|plasma membrane calcium", description, ignore.case = TRUE) ~ "Ca2+ ATPase",
        grepl("calcium.*kinase", description, ignore.case = TRUE) ~ "Ca2+ Signaling (CDPK)",
        TRUE ~ category
    ))

# Save refined list
write.csv(truly_relevant, "data/calcification_degs_refined.csv", row.names = FALSE)

# Summary
cat("\n--- REFINED CATEGORY COUNTS ---\n")
truly_relevant %>%
    count(refined_category, organism) %>%
    pivot_wider(names_from = organism, values_from = n, values_fill = 0) %>%
    arrange(desc(Host + Symbiont)) %>%
    as.data.frame() %>%
    print()

# ==============================================================================
# CREATE CHORD DIAGRAM
# ==============================================================================

cat("\nCreating chord diagram...\n")

# Order: Host first (by abs log2FC descending), then Symbiont
truly_relevant <- truly_relevant %>%
    arrange(desc(organism == "Host"), desc(abs(log2FoldChange)))

# Create gene labels with organism prefix and short name
truly_relevant$gene_label <- paste0(
    ifelse(truly_relevant$organism == "Host", "H:", "S:"),
    truly_relevant$short_name
)

# Make unique labels
truly_relevant <- truly_relevant %>%
    group_by(gene_label) %>%
    mutate(gene_label = ifelse(n() > 1, 
                                paste0(gene_label, ".", row_number()),
                                gene_label)) %>%
    ungroup()

# Get categories present in data
categories <- unique(truly_relevant$refined_category)
n_cats <- length(categories)

cat("Categories:", n_cats, "\n")
cat("Genes:", nrow(truly_relevant), "\n")

# Prepare chord data
chord_df <- truly_relevant %>%
    select(gene_label, refined_category, log2FoldChange, season, organism) %>%
    as.data.frame()

# Define sectors: genes first, then categories
all_sectors <- c(chord_df$gene_label, categories)

# Sector colors - genes by organism
gene_colors <- setNames(
    ifelse(chord_df$organism == "Host", host_color, symbiont_color),
    chord_df$gene_label
)

# Category colors - use a nice palette
category_colors <- setNames(
    colorRampPalette(brewer.pal(min(n_cats, 12), "Set3"))(n_cats),
    categories
)

all_sector_colors <- c(gene_colors, category_colors)

# Link colors by season
season_cols <- c("Summer" = summer_color, "Winter" = winter_color)
link_colors <- adjustcolor(season_cols[chord_df$season], alpha.f = 0.5)

# Log2FC color function
lfc_range <- range(chord_df$log2FoldChange, na.rm = TRUE)
lfc_col_fun <- colorRamp2(c(min(-3, lfc_range[1]), 0, max(3, lfc_range[2])), 
                           c(down_color, "white", up_color))

# Create adjacency list
adj_list <- chord_df %>%
    mutate(value = 1) %>%
    select(from = gene_label, to = refined_category, value)

# ==============================================================================
# DRAW CHORD DIAGRAM
# ==============================================================================

pdf("figures/Fig4_chord_calcification.pdf", width = 14, height = 12)

circos.clear()
circos.par(
    start.degree = 90, 
    gap.degree = 2,
    track.margin = c(0.01, 0.01)
)

chordDiagram(
    adj_list,
    order = all_sectors,
    grid.col = all_sector_colors,
    col = link_colors,
    transparency = 0.3,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.18)
)

# Add labels track
circos.track(track.index = 1, panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    
    if (sector.name %in% chord_df$gene_label) {
        # Gene labels - rotated, smaller
        circos.text(mean(xlim), ylim[1] + 0.4, sector.name,
                   facing = "clockwise", niceFacing = TRUE,
                   adj = c(0, 0.5), cex = 0.45)
    } else {
        # Category labels - bending, larger, bold
        circos.text(mean(xlim), ylim[1] + 0.4, sector.name,
                   facing = "bending.inside", niceFacing = TRUE,
                   adj = c(0.5, 0), cex = 0.75, font = 2)
    }
}, bg.border = NA)

# Add log2FC color highlight for genes
for (i in 1:nrow(chord_df)) {
    gene <- chord_df$gene_label[i]
    lfc <- chord_df$log2FoldChange[i]
    highlight.sector(gene, track.index = 1, 
                   col = lfc_col_fun(lfc), 
                   border = NA)
}

# Draw legends
lgd_lfc <- Legend(
    col_fun = lfc_col_fun, 
    title = "Log2FC",
    title_position = "topleft", 
    legend_height = unit(3, "cm"),
    title_gp = gpar(fontsize = 11, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
)

lgd_season <- Legend(
    labels = c("Summer", "Winter"),
    legend_gp = gpar(fill = c(summer_color, winter_color)),
    title = "Season", 
    title_position = "topleft",
    title_gp = gpar(fontsize = 11, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
)

lgd_org <- Legend(
    labels = c("Host (H:)", "Symbiont (S:)"),
    legend_gp = gpar(fill = c(host_color, symbiont_color)),
    title = "Organism", 
    title_position = "topleft",
    title_gp = gpar(fontsize = 11, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
)

draw(lgd_lfc, x = unit(0.92, "npc"), y = unit(0.82, "npc"))
draw(lgd_season, x = unit(0.92, "npc"), y = unit(0.55, "npc"))
draw(lgd_org, x = unit(0.92, "npc"), y = unit(0.35, "npc"))

circos.clear()
dev.off()

cat("Saved: figures/Fig4_chord_calcification.pdf\n")

# Also save as PNG
png("figures/Fig4_chord_calcification.png", width = 14, height = 12, units = "in", res = 300)

circos.clear()
circos.par(start.degree = 90, gap.degree = 2, track.margin = c(0.01, 0.01))

chordDiagram(
    adj_list, order = all_sectors, grid.col = all_sector_colors,
    col = link_colors, transparency = 0.3,
    annotationTrack = "grid", preAllocateTracks = list(track.height = 0.18)
)

circos.track(track.index = 1, panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    
    if (sector.name %in% chord_df$gene_label) {
        circos.text(mean(xlim), ylim[1] + 0.4, sector.name,
                   facing = "clockwise", niceFacing = TRUE,
                   adj = c(0, 0.5), cex = 0.45)
    } else {
        circos.text(mean(xlim), ylim[1] + 0.4, sector.name,
                   facing = "bending.inside", niceFacing = TRUE,
                   adj = c(0.5, 0), cex = 0.75, font = 2)
    }
}, bg.border = NA)

for (i in 1:nrow(chord_df)) {
    highlight.sector(chord_df$gene_label[i], track.index = 1, 
                   col = lfc_col_fun(chord_df$log2FoldChange[i]), border = NA)
}

draw(lgd_lfc, x = unit(0.92, "npc"), y = unit(0.82, "npc"))
draw(lgd_season, x = unit(0.92, "npc"), y = unit(0.55, "npc"))
draw(lgd_org, x = unit(0.92, "npc"), y = unit(0.35, "npc"))

circos.clear()
dev.off()

cat("Saved: figures/Fig4_chord_calcification.png\n")

# Save category summary
cat_summary <- truly_relevant %>%
    group_by(refined_category, organism, season) %>%
    summarise(count = n(), avg_lfc = mean(log2FoldChange, na.rm = TRUE), .groups = "drop")
write.csv(cat_summary, "data/calcification_category_summary_refined.csv", row.names = FALSE)

cat("\nDone!\n")
