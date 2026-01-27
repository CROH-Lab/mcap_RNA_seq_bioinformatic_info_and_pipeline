#!/usr/bin/env Rscript
# ==============================================================================
# Refined Chord Diagrams - Calcification/Ion Transport DEGs
# Version 6 - Combined figure with A/B panels and shared legend
# ==============================================================================

setwd("/home/darmstrong4/mc_rework/12_publication_figures")
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("gridBase")

suppressPackageStartupMessages({
    library(tidyverse)
    library(circlize)
    library(RColorBrewer)
    library(ComplexHeatmap)
    library(grid)
    library(gridBase)
})

# Color palettes
host_color <- "#E69F00"
symbiont_color <- "#56B4E9"
summer_color <- "#D55E00"
winter_color <- "#0072B2"
up_color <- "#D73027"
down_color <- "#4575B4"

cat("Loading data...\n")

calc_degs <- read.csv("data/all_calcification_degs.csv")

# ==============================================================================
# REFINED FILTERING
# ==============================================================================

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

truly_relevant <- calc_degs %>%
    filter(
        grepl(relevant_pattern, description, ignore.case = TRUE) |
        category %in% c("Ca2+ Channel", "Ca2+ ATPase", "Ca2+ Signaling", 
                        "Carbonic Anhydrase", "H+ ATPase", "HCO3- Transport",
                        "Na+/K+ ATPase", "Ion Exchanger", "Aquaporin")
    ) %>%
    filter(!grepl(exclude_pattern, description, ignore.case = TRUE))

cat("Refined DEGs:", nrow(truly_relevant), "\n")

# ==============================================================================
# BROADER CATEGORY GROUPINGS
# ==============================================================================

truly_relevant <- truly_relevant %>%
    mutate(broad_category = case_when(
        grepl("calcium.*channel|voltage.*calcium|TRP|transient receptor", description, ignore.case = TRUE) ~ "Ca2+ Transport",
        grepl("calcium.*atpase|calcium.*pump|plasma membrane calcium", description, ignore.case = TRUE) ~ "Ca2+ Transport",
        grepl("calcium.*kinase|CDPK", description, ignore.case = TRUE) ~ "Ca2+ Signaling",
        category %in% c("Ca2+ Channel", "Ca2+ ATPase", "Ca2+ Signaling") ~ "Ca2+ Transport",
        
        grepl("potassium.*channel|potassium.*voltage|K\\+", description, ignore.case = TRUE) ~ "K+ Channels",
        grepl("sodium.*channel|Na\\+.*channel", description, ignore.case = TRUE) ~ "Na+ Channels",
        grepl("chloride.*channel|Cl-", description, ignore.case = TRUE) ~ "Cl- Channels",
        category == "K+ Channel" ~ "K+ Channels",
        category == "Na+ Transport" ~ "Na+ Channels",
        category == "Cl- Transport" ~ "Cl- Channels",
        
        grepl("Na\\+/K\\+|sodium.*potassium.*exchanger|NCKX", description, ignore.case = TRUE) ~ "Ion Pumps",
        grepl("H\\+.*ATPase|proton.*transport|V-type", description, ignore.case = TRUE) ~ "Ion Pumps",
        category %in% c("Na+/K+ ATPase", "H+ ATPase", "Ion Exchanger") ~ "Ion Pumps",
        
        grepl("carbonic anhydrase", description, ignore.case = TRUE) ~ "pH/CO2 Regulation",
        grepl("monocarboxylate", description, ignore.case = TRUE) ~ "pH/CO2 Regulation",
        category == "Carbonic Anhydrase" ~ "pH/CO2 Regulation",
        
        grepl("skeletal|SARP|alkaline phosphatase|bone morphogenetic|BMP|collagen", description, ignore.case = TRUE) ~ "Biomineralization",
        
        grepl("sugar phosphate|phosphate.*translocator", description, ignore.case = TRUE) ~ "Metabolite Transport",
        
        TRUE ~ "Other"
    ))

truly_relevant <- truly_relevant %>% filter(broad_category != "Other")

write.csv(truly_relevant, "data/calcification_degs_refined.csv", row.names = FALSE)

# ==============================================================================
# PREPARE DATA FOR BOTH SEASONS
# ==============================================================================

prepare_season_data <- function(data, season_name) {
    
    season_data <- data %>% filter(season == season_name)
    
    if (nrow(season_data) < 3) return(NULL)
    
    season_data <- season_data %>%
        arrange(desc(organism == "Host"), desc(abs(log2FoldChange)))
    
    # Clean gene labels
    season_data$gene_label <- paste0(
        ifelse(season_data$organism == "Host", "H:", "S:"),
        sub("_.*", "", season_data$short_name)
    )
    
    season_data <- season_data %>%
        group_by(gene_label) %>%
        mutate(gene_label = ifelse(n() > 1, 
                                    paste0(gene_label, ".", row_number()),
                                    gene_label)) %>%
        ungroup()
    
    return(season_data)
}

summer_data <- prepare_season_data(truly_relevant, "Summer")
winter_data <- prepare_season_data(truly_relevant, "Winter")

cat("Summer DEGs:", nrow(summer_data), "\n")
cat("Winter DEGs:", nrow(winter_data), "\n")

# Global LFC range for consistent color scale
all_lfc <- c(summer_data$log2FoldChange, winter_data$log2FoldChange)
lfc_min <- min(-3, min(all_lfc, na.rm = TRUE))
lfc_max <- max(3, max(all_lfc, na.rm = TRUE))
lfc_col_fun <- colorRamp2(c(lfc_min, 0, lfc_max), c(down_color, "white", up_color))

# ==============================================================================
# FUNCTION TO DRAW ONE CHORD DIAGRAM
# ==============================================================================

draw_chord <- function(season_data, label_text) {
    
    chord_df <- season_data %>%
        select(gene_label, broad_category, log2FoldChange, organism) %>%
        as.data.frame()
    
    categories <- unique(chord_df$broad_category)
    all_sectors <- c(chord_df$gene_label, categories)
    
    # Gene colors by LFC
    gene_colors <- setNames(
        sapply(chord_df$log2FoldChange, function(x) lfc_col_fun(x)),
        chord_df$gene_label
    )
    
    # Category colors - gray
    category_colors <- setNames(rep("gray85", length(categories)), categories)
    all_sector_colors <- c(gene_colors, category_colors)
    
    # Link colors by organism
    link_colors <- adjustcolor(
        ifelse(chord_df$organism == "Host", host_color, symbiont_color), 
        alpha.f = 0.6
    )
    
    # Adjacency list
    adj_list <- chord_df %>%
        mutate(value = 1) %>%
        select(from = gene_label, to = broad_category, value)
    
    circos.clear()
    circos.par(
        start.degree = 90, 
        gap.degree = 3,
        track.margin = c(0.01, 0.01)
    )
    
    chordDiagram(
        adj_list,
        order = all_sectors,
        grid.col = all_sector_colors,
        col = link_colors,
        transparency = 0.2,
        annotationTrack = "grid",
        preAllocateTracks = list(list(track.height = 0.18))
    )
    
    circos.track(track.index = 1, panel.fun = function(x, y) {
        sector.name <- get.cell.meta.data("sector.index")
        xlim <- get.cell.meta.data("xlim")
        ylim <- get.cell.meta.data("ylim")
        circos.text(mean(xlim), ylim[1] + 0.3, sector.name,
                   facing = "clockwise", niceFacing = TRUE,
                   adj = c(0, 0.5), cex = 0.7, font = 2, col = "black")
    }, bg.border = NA)
    
    circos.clear()
}

# ==============================================================================
# CREATE COMBINED PDF WITH STACKED PANELS
# ==============================================================================

cat("\nCreating combined figure...\n")

pdf("figures/Fig4_chord_combined.pdf", width = 10, height = 18)

# Layout: 2 rows for plots, space at bottom for legend
layout(matrix(c(1, 2, 3, 3), nrow = 2, ncol = 2, byrow = TRUE),
       widths = c(1, 1),
       heights = c(1, 1, 0.25, 0.25))

# Adjusted layout for stacked plots with legend on right
layout(matrix(c(1, 3,
                2, 3), nrow = 2, ncol = 2, byrow = TRUE),
       widths = c(4, 1),
       heights = c(1, 1))

# Panel A - Summer
par(mar = c(1, 1, 2, 1))
draw_chord(summer_data, "A")
mtext("A", side = 3, line = 0, adj = 0, cex = 2, font = 2)

# Panel B - Winter
par(mar = c(1, 1, 2, 1))
draw_chord(winter_data, "B")
mtext("B", side = 3, line = 0, adj = 0, cex = 2, font = 2)

# Panel 3 - Legends
par(mar = c(2, 0, 2, 0))
plot.new()

# Draw legends using grid
pushViewport(viewport())

# Log2FC legend
lgd_lfc <- Legend(
    col_fun = lfc_col_fun, 
    title = "Log2FC",
    title_position = "topleft", 
    legend_height = unit(4, "cm"),
    legend_width = unit(1, "cm"),
    title_gp = gpar(fontsize = 14, fontface = "bold"),
    labels_gp = gpar(fontsize = 12),
    grid_width = unit(0.8, "cm")
)

# Organism legend
lgd_org <- Legend(
    labels = c("Host (H:)", "Symbiont (S:)"),
    legend_gp = gpar(fill = c(host_color, symbiont_color)),
    title = "Organism\n(chord color)", 
    title_position = "topleft",
    title_gp = gpar(fontsize = 14, fontface = "bold"),
    labels_gp = gpar(fontsize = 12),
    grid_height = unit(0.8, "cm"),
    grid_width = unit(0.8, "cm")
)

draw(lgd_lfc, x = unit(0.5, "npc"), y = unit(0.7, "npc"))
draw(lgd_org, x = unit(0.5, "npc"), y = unit(0.25, "npc"))

popViewport()

dev.off()

cat("Saved: figures/Fig4_chord_combined.pdf\n")

# ==============================================================================
# CREATE COMBINED PNG - Better approach with viewport
# ==============================================================================

png("figures/Fig4_chord_combined.png", width = 10, height = 18, units = "in", res = 300)

layout(matrix(c(1, 3,
                2, 3), nrow = 2, ncol = 2, byrow = TRUE),
       widths = c(4, 1),
       heights = c(1, 1))

# Panel A - Summer
par(mar = c(1, 1, 2, 1))
draw_chord(summer_data, "A")
mtext("A", side = 3, line = 0, adj = 0, cex = 2, font = 2)

# Panel B - Winter
par(mar = c(1, 1, 2, 1))
draw_chord(winter_data, "B")
mtext("B", side = 3, line = 0, adj = 0, cex = 2, font = 2)

# Panel 3 - Legends
par(mar = c(2, 0, 2, 0))
plot.new()

pushViewport(viewport())
draw(lgd_lfc, x = unit(0.5, "npc"), y = unit(0.7, "npc"))
draw(lgd_org, x = unit(0.5, "npc"), y = unit(0.25, "npc"))
popViewport()

dev.off()

cat("Saved: figures/Fig4_chord_combined.png\n")

# ==============================================================================
# ALSO SAVE INDIVIDUAL PLOTS (without legends, for flexibility)
# ==============================================================================

# Summer only
pdf("figures/Fig4A_chord_summer.pdf", width = 10, height = 10)
par(mar = c(1, 1, 1, 1))
draw_chord(summer_data, "A")
dev.off()

png("figures/Fig4A_chord_summer.png", width = 10, height = 10, units = "in", res = 300)
par(mar = c(1, 1, 1, 1))
draw_chord(summer_data, "A")
dev.off()

# Winter only
pdf("figures/Fig4B_chord_winter.pdf", width = 10, height = 10)
par(mar = c(1, 1, 1, 1))
draw_chord(winter_data, "B")
dev.off()

png("figures/Fig4B_chord_winter.png", width = 10, height = 10, units = "in", res = 300)
par(mar = c(1, 1, 1, 1))
draw_chord(winter_data, "B")
dev.off()

cat("\nSaved individual plots as well.\n")

# Save summary
cat_summary <- truly_relevant %>%
    group_by(broad_category, organism, season) %>%
    summarise(count = n(), avg_lfc = mean(log2FoldChange, na.rm = TRUE), .groups = "drop")
write.csv(cat_summary, "data/calcification_category_summary_refined.csv", row.names = FALSE)

cat("\n=== COMPLETE ===\n")
