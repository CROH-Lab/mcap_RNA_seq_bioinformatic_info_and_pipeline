#!/usr/bin/env Rscript
# ==============================================================================
# Refined Chord Diagrams - Calcification/Ion Transport DEGs
# Version 5 - Clean gene names, no titles, matching font sizes
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

cat("\n--- COUNTS BY SEASON ---\n")
truly_relevant %>%
    count(season, organism) %>%
    pivot_wider(names_from = organism, values_from = n, values_fill = 0) %>%
    as.data.frame() %>%
    print()

# ==============================================================================
# FUNCTION TO CREATE CHORD DIAGRAM FOR ONE SEASON
# ==============================================================================

create_seasonal_chord <- function(data, season_name, output_prefix) {
    
    season_data <- data %>% filter(season == season_name)
    
    if (nrow(season_data) < 3) {
        cat("  Skipping", season_name, "- too few genes (", nrow(season_data), ")\n")
        return(NULL)
    }
    
    cat("\nCreating", season_name, "chord diagram with", nrow(season_data), "genes...\n")
    
    # Order: Host first, then Symbiont, by abs(log2FC)
    season_data <- season_data %>%
        arrange(desc(organism == "Host"), desc(abs(log2FoldChange)))
    
    # Create gene labels - CLEAN: remove species notation (everything after _)
    season_data$gene_label <- paste0(
        ifelse(season_data$organism == "Host", "H:", "S:"),
        sub("_.*", "", season_data$short_name)  # Remove _ and everything after
    )
    
    # Make unique labels
    season_data <- season_data %>%
        group_by(gene_label) %>%
        mutate(gene_label = ifelse(n() > 1, 
                                    paste0(gene_label, ".", row_number()),
                                    gene_label)) %>%
        ungroup()
    
    categories <- unique(season_data$broad_category)
    n_cats <- length(categories)
    
    cat("  Categories:", paste(categories, collapse = ", "), "\n")
    
    chord_df <- season_data %>%
        select(gene_label, broad_category, log2FoldChange, organism) %>%
        as.data.frame()
    
    all_sectors <- c(chord_df$gene_label, categories)
    
    # Log2FC color function
    lfc_range <- range(chord_df$log2FoldChange, na.rm = TRUE)
    lfc_min <- min(-3, lfc_range[1])
    lfc_max <- max(3, lfc_range[2])
    lfc_col_fun <- colorRamp2(c(lfc_min, 0, lfc_max), c(down_color, "white", up_color))
    
    # Gene sector colors - BY LOG2FC
    gene_colors <- setNames(
        sapply(chord_df$log2FoldChange, function(x) lfc_col_fun(x)),
        chord_df$gene_label
    )
    
    # Category colors - gray/neutral
    category_colors <- setNames(rep("gray85", n_cats), categories)
    
    all_sector_colors <- c(gene_colors, category_colors)
    
    # Link/chord colors by organism
    link_colors <- adjustcolor(
        ifelse(chord_df$organism == "Host", host_color, symbiont_color), 
        alpha.f = 0.6
    )
    
    # Adjacency list
    adj_list <- chord_df %>%
        mutate(value = 1) %>%
        select(from = gene_label, to = broad_category, value)
    
    # ==============================================================================
    # DRAW PDF
    # ==============================================================================
    
    pdf(paste0("figures/", output_prefix, ".pdf"), width = 12, height = 11)
    
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
        preAllocateTracks = list(
            list(track.height = 0.18)
        )
    )
    
    # Single track for labels - all perpendicular, bold, same size
    circos.track(track.index = 1, panel.fun = function(x, y) {
        sector.name <- get.cell.meta.data("sector.index")
        xlim <- get.cell.meta.data("xlim")
        ylim <- get.cell.meta.data("ylim")
        
        # Both genes and categories: perpendicular, bold, same size
        circos.text(mean(xlim), ylim[1] + 0.3, sector.name,
                   facing = "clockwise", niceFacing = TRUE,
                   adj = c(0, 0.5), cex = 0.8, font = 2, col = "black")
    }, bg.border = NA)
    
    # NO TITLE
    
    # Legends
    lgd_lfc <- Legend(
        col_fun = lfc_col_fun, 
        title = "Log2FC",
        title_position = "topleft", 
        legend_height = unit(2.5, "cm"),
        title_gp = gpar(fontsize = 11, fontface = "bold"),
        labels_gp = gpar(fontsize = 9)
    )
    
    lgd_org <- Legend(
        labels = c("Host (H:)", "Symbiont (S:)"),
        legend_gp = gpar(fill = c(host_color, symbiont_color)),
        title = "Organism\n(chord color)", 
        title_position = "topleft",
        title_gp = gpar(fontsize = 11, fontface = "bold"),
        labels_gp = gpar(fontsize = 9)
    )
    
    draw(lgd_lfc, x = unit(0.92, "npc"), y = unit(0.75, "npc"))
    draw(lgd_org, x = unit(0.92, "npc"), y = unit(0.45, "npc"))
    
    circos.clear()
    dev.off()
    
    cat("  Saved:", paste0("figures/", output_prefix, ".pdf"), "\n")
    
    # ==============================================================================
    # DRAW PNG
    # ==============================================================================
    
    png(paste0("figures/", output_prefix, ".png"), width = 12, height = 11, units = "in", res = 300)
    
    circos.clear()
    circos.par(start.degree = 90, gap.degree = 3, track.margin = c(0.01, 0.01))
    
    chordDiagram(
        adj_list, order = all_sectors, grid.col = all_sector_colors,
        col = link_colors, transparency = 0.2,
        annotationTrack = "grid",
        preAllocateTracks = list(list(track.height = 0.18))
    )
    
    circos.track(track.index = 1, panel.fun = function(x, y) {
        sector.name <- get.cell.meta.data("sector.index")
        xlim <- get.cell.meta.data("xlim")
        ylim <- get.cell.meta.data("ylim")
        circos.text(mean(xlim), ylim[1] + 0.3, sector.name,
                   facing = "clockwise", niceFacing = TRUE,
                   adj = c(0, 0.5), cex = 0.8, font = 2, col = "black")
    }, bg.border = NA)
    
    draw(lgd_lfc, x = unit(0.92, "npc"), y = unit(0.75, "npc"))
    draw(lgd_org, x = unit(0.92, "npc"), y = unit(0.45, "npc"))
    
    circos.clear()
    dev.off()
    
    cat("  Saved:", paste0("figures/", output_prefix, ".png"), "\n")
    
    return(season_data)
}

# ==============================================================================
# CREATE BOTH SEASONAL CHORD DIAGRAMS
# ==============================================================================

summer_data <- create_seasonal_chord(truly_relevant, "Summer", "Fig4A_chord_summer")
winter_data <- create_seasonal_chord(truly_relevant, "Winter", "Fig4B_chord_winter")

# Save summary
cat_summary <- truly_relevant %>%
    group_by(broad_category, organism, season) %>%
    summarise(count = n(), avg_lfc = mean(log2FoldChange, na.rm = TRUE), .groups = "drop")
write.csv(cat_summary, "data/calcification_category_summary_refined.csv", row.names = FALSE)

cat("\n=== SUMMARY ===\n")
cat("Summer DEGs:", ifelse(is.null(summer_data), 0, nrow(summer_data)), "\n")
cat("Winter DEGs:", ifelse(is.null(winter_data), 0, nrow(winter_data)), "\n")

cat("\nDone!\n")
