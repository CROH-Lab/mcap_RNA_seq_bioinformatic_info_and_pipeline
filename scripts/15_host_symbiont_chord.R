#!/usr/bin/env Rscript
# =============================================================================
# Chord Diagram: Host-Symbiont Shared GO Term Analysis
# =============================================================================
# Visualizes shared and unique significant GO terms between host and symbiont
# - Left half: M. capitata (host)
# - Right half: D. trenchii (symbiont)
# - Chords connect matching GO terms by ID
# - Colors indicate expression coordination:
#   - Red: Synergistic up (both up-regulated)
#   - Blue: Synergistic down (both down-regulated)
#   - Light green: Antagonistic (host up, symbiont down)
#   - Dark green: Antagonistic (host down, symbiont up)
#
# Author: Claude (for David Armstrong)
# Date: 2025-02-06
# =============================================================================

# Load libraries
library(tidyverse)
library(circlize)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Directories
host_dir <- "/home/darmstrong4/mc_rework/10_GO_MWU/output"
symbiont_dir <- "/home/darmstrong4/mc_rework/11_symbiont_GO_MWU/output"
output_dir <- "/home/darmstrong4/mc_rework/12_publication_figures/figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Significance threshold for including terms
sig_threshold <- 0.1  # p.adj < 0.1

# Colors
color_synergistic_up <- "#E41A1C"      # Red - both up
color_synergistic_down <- "#377EB8"    # Blue - both down
color_antag_host_up <- "#90EE90"       # Light green - host up, symbiont down
color_antag_host_down <- "#228B22"     # Dark/forest green - host down, symbiont up

# Sector colors by division
division_colors <- list(
    host = c("BP" = "#FDB462", "CC" = "#80B1D3", "MF" = "#FB8072"),
    symbiont = c("BP" = "#B3DE69", "CC" = "#BEBADA", "MF" = "#FCCDE5")
)

# =============================================================================
# HELPER: Load significant GO terms from MWU results
# =============================================================================

load_significant_terms <- function(data_dir, prefix, season, sig_threshold) {
    
    divisions <- c("BP", "CC", "MF")
    all_terms <- list()
    
    for (div in divisions) {
        if (prefix == "") {
            mwu_file <- file.path(data_dir, paste0(season, "_", div, "_MWU_results.csv"))
        } else {
            mwu_file <- file.path(data_dir, paste0(prefix, "_", season, "_", div, "_MWU_results.csv"))
        }
        
        if (!file.exists(mwu_file)) {
            cat("  Warning: Missing", mwu_file, "\n")
            next
        }
        
        # Read MWU results (space/tab-separated)
        mwu <- read.table(mwu_file, header = TRUE, stringsAsFactors = FALSE,
                          quote = "\"", fill = TRUE)
        
        # Filter significant terms
        sig_terms <- mwu %>%
            filter(p.adj < sig_threshold) %>%
            mutate(
                division = div,
                direction = ifelse(delta.rank > 0, "up", "down"),
                go_id = term
            ) %>%
            select(go_id, name, division, direction, delta.rank, pval, p.adj)
        
        all_terms[[div]] <- sig_terms
    }
    
    # Combine all divisions
    combined <- bind_rows(all_terms)
    return(combined)
}

# =============================================================================
# HELPER: Classify shared terms by expression pattern
# =============================================================================

classify_shared_terms <- function(host_terms, symbiont_terms) {
    
    # Find shared GO IDs
    shared_ids <- intersect(host_terms$go_id, symbiont_terms$go_id)
    
    if (length(shared_ids) == 0) {
        return(data.frame())
    }
    
    # Get direction for each shared term
    shared_data <- data.frame(go_id = shared_ids, stringsAsFactors = FALSE)
    
    shared_data <- shared_data %>%
        left_join(host_terms %>% select(go_id, name, division, direction) %>% 
                      rename(host_direction = direction, host_division = division, term_name = name),
                  by = "go_id") %>%
        left_join(symbiont_terms %>% select(go_id, direction, division) %>%
                      rename(symbiont_direction = direction, symbiont_division = division),
                  by = "go_id")
    
    # Classify expression pattern
    shared_data <- shared_data %>%
        mutate(
            pattern = case_when(
                host_direction == "up" & symbiont_direction == "up" ~ "synergistic_up",
                host_direction == "down" & symbiont_direction == "down" ~ "synergistic_down",
                host_direction == "up" & symbiont_direction == "down" ~ "antagonistic_host_up",
                host_direction == "down" & symbiont_direction == "up" ~ "antagonistic_host_down"
            ),
            chord_color = case_when(
                pattern == "synergistic_up" ~ color_synergistic_up,
                pattern == "synergistic_down" ~ color_synergistic_down,
                pattern == "antagonistic_host_up" ~ color_antag_host_up,
                pattern == "antagonistic_host_down" ~ color_antag_host_down
            )
        )
    
    return(shared_data)
}

# =============================================================================
# MAIN: Create chord diagram for one season
# =============================================================================

create_chord_diagram <- function(season, output_file) {
    
    cat("\n=== Creating chord diagram for", season, "===\n")
    
    # Load significant terms
    cat("  Loading host terms...\n")
    host_terms <- load_significant_terms(host_dir, "", season, sig_threshold)
    cat("    Significant terms:", nrow(host_terms), "\n")
    
    cat("  Loading symbiont terms...\n")
    symbiont_terms <- load_significant_terms(symbiont_dir, "symbiont", season, sig_threshold)
    cat("    Significant terms:", nrow(symbiont_terms), "\n")
    
    # Classify shared terms
    cat("  Identifying shared terms...\n")
    shared_terms <- classify_shared_terms(host_terms, symbiont_terms)
    cat("    Shared terms:", nrow(shared_terms), "\n")
    
    if (nrow(shared_terms) > 0) {
        cat("    - Synergistic up:", sum(shared_terms$pattern == "synergistic_up"), "\n")
        cat("    - Synergistic down:", sum(shared_terms$pattern == "synergistic_down"), "\n")
        cat("    - Antagonistic (host up):", sum(shared_terms$pattern == "antagonistic_host_up"), "\n")
        cat("    - Antagonistic (host down):", sum(shared_terms$pattern == "antagonistic_host_down"), "\n")
    }
    
    # Create sector data
    # Each sector represents a division within an organism
    # Sector size proportional to number of significant terms
    
    host_sectors <- host_terms %>%
        group_by(division) %>%
        summarise(n_terms = n(), .groups = "drop") %>%
        mutate(
            sector = paste0("Host_", division),
            organism = "host"
        )
    
    symbiont_sectors <- symbiont_terms %>%
        group_by(division) %>%
        summarise(n_terms = n(), .groups = "drop") %>%
        mutate(
            sector = paste0("Symbiont_", division),
            organism = "symbiont"
        )
    
    all_sectors <- bind_rows(host_sectors, symbiont_sectors)
    
    # Create sector colors
    sector_colors <- c(
        setNames(division_colors$host, paste0("Host_", names(division_colors$host))),
        setNames(division_colors$symbiont, paste0("Symbiont_", names(division_colors$symbiont)))
    )
    
    # Define sector order: Host on left (BP, CC, MF), Symbiont on right (MF, CC, BP)
    # This creates a mirror layout
    sector_order <- c("Host_BP", "Host_CC", "Host_MF", "Symbiont_MF", "Symbiont_CC", "Symbiont_BP")
    sector_order <- sector_order[sector_order %in% all_sectors$sector]
    
    # Create the chord diagram
    cat("  Generating chord diagram...\n")
    
    pdf(output_file, width = 10, height = 10)
    
    # Set up circos parameters
    circos.clear()
    circos.par(
        start.degree = 90,
        gap.degree = c(rep(2, 2), 8, rep(2, 2), 8),  # Larger gap between organisms
        track.margin = c(0.01, 0.01)
    )
    
    # Initialize sectors with sizes
    sector_sizes <- setNames(all_sectors$n_terms, all_sectors$sector)
    xlim_matrix <- matrix(0, nrow = length(sector_order), ncol = 2)
    rownames(xlim_matrix) <- sector_order
    for (i in seq_along(sector_order)) {
        xlim_matrix[i, 2] <- sector_sizes[sector_order[i]]
    }
    
    circos.initialize(factors = factor(sector_order, levels = sector_order),
                      xlim = xlim_matrix)
    
    # Add track for sector labels
    circos.track(
        ylim = c(0, 1),
        track.height = 0.08,
        bg.border = NA,
        panel.fun = function(x, y) {
            sector.name <- CELL_META$sector.index
            division <- sub(".*_", "", sector.name)
            circos.text(
                CELL_META$xcenter, 0.5,
                as.character(division),
                facing = "bending.inside",
                niceFacing = TRUE,
                cex = 0.9,
                font = 2
            )
        }
    )
    
    # Add main track (colored sectors)
    circos.track(
        ylim = c(0, 1),
        track.height = 0.15,
        bg.col = sector_colors[sector_order],
        bg.border = "grey30",
        panel.fun = function(x, y) {
            # Optional: add term count
            sector.name <- CELL_META$sector.index
            n <- sector_sizes[sector.name]
            circos.text(
                CELL_META$xcenter, 0.5,
                as.character(n),
                facing = "bending.inside",
                niceFacing = TRUE,
                cex = 0.7,
                col = "grey30"
            )
        }
    )
    
    # Add chords for shared terms
    if (nrow(shared_terms) > 0) {
        
        # Track position within each sector for placing chord endpoints
        host_positions <- list(BP = 0, CC = 0, MF = 0)
        symbiont_positions <- list(BP = 0, CC = 0, MF = 0)
        
        for (i in 1:nrow(shared_terms)) {
            term <- shared_terms[i, ]
            
            # Get sectors
            host_sector <- paste0("Host_", term$host_division)
            symbiont_sector <- paste0("Symbiont_", term$symbiont_division)
            
            # Get positions (increment for each term)
            host_pos <- host_positions[[term$host_division]]
            symbiont_pos <- symbiont_positions[[term$symbiont_division]]
            
            # Draw chord
            circos.link(
                host_sector, c(host_pos, host_pos + 1),
                symbiont_sector, c(symbiont_pos, symbiont_pos + 1),
                col = adjustcolor(term$chord_color, alpha.f = 0.6),
                border = adjustcolor(term$chord_color, alpha.f = 0.8),
                lwd = 0.5
            )
            
            # Increment positions
            host_positions[[term$host_division]] <- host_pos + 1
            symbiont_positions[[term$symbiont_division]] <- symbiont_pos + 1
        }
    }
    
    # Add organism labels
    # Host label (left side)
    text(-0.85, 0, expression(italic("M. capitata")), cex = 1.1, font = 2)
    text(-0.85, -0.08, "(host)", cex = 0.9)
    
    # Symbiont label (right side)
    text(0.85, 0, expression(italic("D. trenchii")), cex = 1.1, font = 2)
    text(0.85, -0.08, "(symbiont)", cex = 0.9)
    
    # Add legend
    legend_x <- -0.15
    legend_y <- -0.75
    legend_spacing <- 0.08
    
    # Legend title
    text(legend_x, legend_y + 0.12, "Expression Pattern", cex = 0.85, font = 2, adj = c(0, 0.5))
    
    # Legend items
    points(legend_x, legend_y, pch = 15, col = color_synergistic_up, cex = 1.5)
    text(legend_x + 0.03, legend_y, "Synergistic up (both up)", cex = 0.7, adj = c(0, 0.5))
    
    points(legend_x, legend_y - legend_spacing, pch = 15, col = color_synergistic_down, cex = 1.5)
    text(legend_x + 0.03, legend_y - legend_spacing, "Synergistic down (both down)", cex = 0.7, adj = c(0, 0.5))
    
    points(legend_x, legend_y - 2*legend_spacing, pch = 15, col = color_antag_host_up, cex = 1.5)
    text(legend_x + 0.03, legend_y - 2*legend_spacing, "Antagonistic (host up, symbiont down)", cex = 0.7, adj = c(0, 0.5))
    
    points(legend_x, legend_y - 3*legend_spacing, pch = 15, col = color_antag_host_down, cex = 1.5)
    text(legend_x + 0.03, legend_y - 3*legend_spacing, "Antagonistic (host down, symbiont up)", cex = 0.7, adj = c(0, 0.5))
    
    # Season title
    title(main = paste0(tools::toTitleCase(season), " - Shared GO Term Expression Patterns"),
          line = -1, cex.main = 1.3)
    
    circos.clear()
    dev.off()
    
    cat("  SUCCESS:", basename(output_file), "\n")
    
    # Return summary statistics
    return(list(
        season = season,
        host_total = nrow(host_terms),
        symbiont_total = nrow(symbiont_terms),
        shared_total = nrow(shared_terms),
        shared_details = if(nrow(shared_terms) > 0) shared_terms else NULL
    ))
}

# =============================================================================
# GENERATE CHORD DIAGRAMS
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("GENERATING HOST-SYMBIONT CHORD DIAGRAMS\n")
cat(strrep("=", 60), "\n")

# Summer
summer_results <- create_chord_diagram(
    season = "summer",
    output_file = file.path(output_dir, "Fig_Summer_Host_Symbiont_GO_Chord.pdf")
)

# Winter
winter_results <- create_chord_diagram(
    season = "winter",
    output_file = file.path(output_dir, "Fig_Winter_Host_Symbiont_GO_Chord.pdf")
)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("SUMMARY\n")
cat(strrep("=", 60), "\n")

cat("\nSummer:\n")
cat("  Host significant terms:", summer_results$host_total, "\n")
cat("  Symbiont significant terms:", summer_results$symbiont_total, "\n")
cat("  Shared terms:", summer_results$shared_total, "\n")

cat("\nWinter:\n")
cat("  Host significant terms:", winter_results$host_total, "\n")
cat("  Symbiont significant terms:", winter_results$symbiont_total, "\n")
cat("  Shared terms:", winter_results$shared_total, "\n")

# Save shared terms to CSV for reference
if (!is.null(summer_results$shared_details)) {
    write.csv(summer_results$shared_details, 
              file.path(output_dir, "Summer_Shared_GO_Terms.csv"),
              row.names = FALSE)
    cat("\nSaved: Summer_Shared_GO_Terms.csv\n")
}

if (!is.null(winter_results$shared_details)) {
    write.csv(winter_results$shared_details,
              file.path(output_dir, "Winter_Shared_GO_Terms.csv"),
              row.names = FALSE)
    cat("Saved: Winter_Shared_GO_Terms.csv\n")
}

cat("\nOutput directory:", output_dir, "\n")
