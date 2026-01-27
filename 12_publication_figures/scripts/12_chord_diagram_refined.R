#!/usr/bin/env Rscript
# ==============================================================================
# Refined Chord Diagrams - Calcification/Ion Transport DEGs
# Version 3 - Separate diagrams for Summer and Winter
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
        # Ca2+ Transport & Signaling
        grepl("calcium.*channel|voltage.*calcium|TRP|transient receptor", description, ignore.case = TRUE) ~ "Ca2+ Transport",
        grepl("calcium.*atpase|calcium.*pump|plasma membrane calcium", description, ignore.case = TRUE) ~ "Ca2+ Transport",
        grepl("calcium.*kinase|CDPK", description, ignore.case = TRUE) ~ "Ca2+ Signaling",
        category %in% c("Ca2+ Channel", "Ca2+ ATPase", "Ca2+ Signaling") ~ "Ca2+ Transport",
        
        # Other Ion Channels
        grepl("potassium.*channel|potassium.*voltage|K\\+", description, ignore.case = TRUE) ~ "K+ Channels",
        grepl("sodium.*channel|Na\\+.*channel", description, ignore.case = TRUE) ~ "Na+ Channels",
        grepl("chloride.*channel|Cl-", description, ignore.case = TRUE) ~ "Cl- Channels",
        category == "K+ Channel" ~ "K+ Channels",
        category == "Na+ Transport" ~ "Na+ Channels",
        category == "Cl- Transport" ~ "Cl- Channels",
        
        # Ion Pumps & Exchangers
        grepl("Na\\+/K\\+|sodium.*potassium.*exchanger|NCKX", description, ignore.case = TRUE) ~ "Ion Pumps",
        grepl("H\\+.*ATPase|proton.*transport|V-type", description, ignore.case = TRUE) ~ "Ion Pumps",
        category %in% c("Na+/K+ ATPase", "H+ ATPase", "Ion Exchanger") ~ "Ion Pumps",
        
        # pH/CO2 Regulation
        grepl("carbonic anhydrase", description, ignore.case = TRUE) ~ "pH/CO2 Regulation",
        grepl("monocarboxylate", description, ignore.case = TRUE) ~ "pH/CO2 Regulation",
        category == "Carbonic Anhydrase" ~ "pH/CO2 Regulation",
        
        # Biomineralization & ECM
        grepl("skeletal|SARP|alkaline phosphatase|bone morphogenetic|BMP|collagen", description, ignore.case = TRUE) ~ "Biomineralization",
        
        # Metabolite Transport
        grepl("sugar phosphate|phosphate.*translocator", description, ignore.case = TRUE) ~ "Metabolite Transport",
        
        TRUE ~ "Other"
    ))

# Remove "Other" if any slipped through
truly_relevant <- truly_relevant %>% filter(broad_category != "Other")

# Save refined list
write.csv(truly_relevant, "data/calcification_degs_refined.csv", row.names = FALSE)

# Summary by season
cat("\n--- COUNTS BY SEASON ---\n")
truly_relevant %>%
    count(season, organism) %>%
    pivot_wider(names_from = organism, values_from = n, values_fill = 0) %>%
    as.data.frame() %>%
    print()

cat("\n--- BROAD CATEGORY COUNTS ---\n")
truly_relevant %>%
    count(broad_category, season, organism) %>%
    pivot_wider(names_from = c(organism, season), values_from = n, values_fill = 0) %>%
    as.data.frame() %>%
    print()

# ==============================================================================
# FUNCTION TO CREATE CHORD DIAGRAM FOR ONE SEASON
# ==============================================================================

create_seasonal_chord <- function(data, season_name, output_prefix) {
    
    # Filter for this season
    season_data <- data %>% filter(season == season_name)
    
    if (nrow(season_data) < 3) {
        cat("  Skipping", season_name, "- too few genes (", nrow(season_data), ")\n")
        return(NULL)
    }
    
    cat("\nCreating", season_name, "chord diagram with", nrow(season_data), "genes...\n")
    
    # Order: Host first, then Symbiont, by abs(log2FC)
    season_data <- season_data %>%
        arrange(desc(organism == "Host"), desc(abs(log2FoldChange)))
    
    # Create gene labels
    season_data$gene_label <- paste0(
        ifelse(season_data$organism == "Host", "H:", "S:"),
        season_data$short_name
    )
    
    # Make unique labels
    season_data <- season_data %>%
        group_by(gene_label) %>%
        mutate(gene_label = ifelse(n() > 1, 
                                    paste0(gene_label, ".", row_number()),
                                    gene_label)) %>%
        ungroup()
    
    # Get categories present in this season's data
    categories <- unique(season_data$broad_category)
    n_cats <- length(categories)
    
    cat("  Categories:", paste(categories, collapse = ", "), "\n")
    
    # Prepare chord data
    chord_df <- season_data %>%
        select(gene_label, broad_category, log2FoldChange, organism) %>%
        as.data.frame()
    
    # Define sectors
    all_sectors <- c(chord_df$gene_label, categories)
    
    # Sector colors - genes by organism
    gene_colors <- setNames(
        ifelse(chord_df$organism == "Host", host_color, symbiont_color),
        chord_df$gene_label
    )
    
    # Category colors
    cat_color_palette <- c(
        "Ca2+ Transport" = "#E41A1C",
        "Ca2+ Signaling" = "#FC8D62", 
        "K+ Channels" = "#377EB8",
        "Na+ Channels" = "#4DAF4A",
        "Cl- Channels" = "#984EA3",
        "Ion Pumps" = "#FF7F00",
        "pH/CO2 Regulation" = "#A65628",
        "Biomineralization" = "#F781BF",
        "Metabolite Transport" = "#999999"
    )
    
    category_colors <- cat_color_palette[categories]
    all_sector_colors <- c(gene_colors, category_colors)
    
    # Link colors by organism for this season
    link_colors <- adjustcolor(
        ifelse(chord_df$organism == "Host", host_color, symbiont_color), 
        alpha.f = 0.5
    )
    
    # Log2FC color function
    lfc_range <- range(chord_df$log2FoldChange, na.rm = TRUE)
    lfc_min <- min(-3, lfc_range[1])
    lfc_max <- max(3, lfc_range[2])
    lfc_col_fun <- colorRamp2(c(lfc_min, 0, lfc_max), c(down_color, "white", up_color))
    
    # Adjacency list
    adj_list <- chord_df %>%
        mutate(value = 1) %>%
        select(from = gene_label, to = broad_category, value)
    
    # Season-specific colors for title
    season_title_color <- ifelse(season_name == "Summer", summer_color, winter_color)
    
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
        transparency = 0.3,
        annotationTrack = "grid",
        preAllocateTracks = list(
            list(track.height = 0.05),
            list(track.height = 0.18)
        )
    )
    
    # Track 1: LFC color bar
    circos.track(track.index = 1, panel.fun = function(x, y) {
        sector.name <- get.cell.meta.data("sector.index")
        xlim <- get.cell.meta.data("xlim")
        if (sector.name %in% chord_df$gene_label) {
            lfc <- chord_df$log2FoldChange[chord_df$gene_label == sector.name][1]
            circos.rect(xlim[1], 0, xlim[2], 1, col = lfc_col_fun(lfc), border = NA)
        }
    }, bg.border = NA)
    
    # Track 2: Labels
    circos.track(track.index = 2, panel.fun = function(x, y) {
        sector.name <- get.cell.meta.data("sector.index")
        xlim <- get.cell.meta.data("xlim")
        ylim <- get.cell.meta.data("ylim")
        if (sector.name %in% chord_df$gene_label) {
            circos.text(mean(xlim), ylim[1] + 0.3, sector.name,
                       facing = "clockwise", niceFacing = TRUE,
                       adj = c(0, 0.5), cex = 0.55, col = "black")
        } else {
            circos.text(mean(xlim), ylim[1] + 0.3, sector.name,
                       facing = "bending.inside", niceFacing = TRUE,
                       adj = c(0.5, 0), cex = 1.0, font = 2, col = "black")
        }
    }, bg.border = NA)
    
    # Title
    title(main = paste0(season_name, " - Calcification/Ion Transport DEGs"), 
          cex.main = 1.5, col.main = season_title_color, font.main = 2)
    
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
        title = "Organism", 
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
        col = link_colors, transparency = 0.3,
        annotationTrack = "grid",
        preAllocateTracks = list(
            list(track.height = 0.05),
            list(track.height = 0.18)
        )
    )
    
    circos.track(track.index = 1, panel.fun = function(x, y) {
        sector.name <- get.cell.meta.data("sector.index")
        xlim <- get.cell.meta.data("xlim")
        if (sector.name %in% chord_df$gene_label) {
            lfc <- chord_df$log2FoldChange[chord_df$gene_label == sector.name][1]
            circos.rect(xlim[1], 0, xlim[2], 1, col = lfc_col_fun(lfc), border = NA)
        }
    }, bg.border = NA)
    
    circos.track(track.index = 2, panel.fun = function(x, y) {
        sector.name <- get.cell.meta.data("sector.index")
        xlim <- get.cell.meta.data("xlim")
        ylim <- get.cell.meta.data("ylim")
        if (sector.name %in% chord_df$gene_label) {
            circos.text(mean(xlim), ylim[1] + 0.3, sector.name,
                       facing = "clockwise", niceFacing = TRUE,
                       adj = c(0, 0.5), cex = 0.55, col = "black")
        } else {
            circos.text(mean(xlim), ylim[1] + 0.3, sector.name,
                       facing = "bending.inside", niceFacing = TRUE,
                       adj = c(0.5, 0), cex = 1.0, font = 2, col = "black")
        }
    }, bg.border = NA)
    
    title(main = paste0(season_name, " - Calcification/Ion Transport DEGs"), 
          cex.main = 1.5, col.main = season_title_color, font.main = 2)
    
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

# ==============================================================================
# SAVE SUMMARIES
# ==============================================================================

# Category summary by season
cat_summary <- truly_relevant %>%
    group_by(broad_category, organism, season) %>%
    summarise(count = n(), avg_lfc = mean(log2FoldChange, na.rm = TRUE), .groups = "drop")
write.csv(cat_summary, "data/calcification_category_summary_refined.csv", row.names = FALSE)

cat("\n=== SUMMARY ===\n")
cat("Summer DEGs:", ifelse(is.null(summer_data), 0, nrow(summer_data)), "\n")
cat("Winter DEGs:", ifelse(is.null(winter_data), 0, nrow(winter_data)), "\n")

cat("\nDone!\n")
