#!/usr/bin/env Rscript
# =============================================================================
# Chord Diagram: Host-Symbiont Shared GO Terms (Redesigned)
# =============================================================================
# Layout:
# - LEFT: Individual GO terms (labeled by GO ID), colored by division
# - RIGHT: 6 nodes (Host_BP, Host_CC, Host_MF, Symbiont_BP, Symbiont_CC, Symbiont_MF)
# - OUTER RING: Interaction type grouping (Synergistic Up/Down, Antagonistic)
# - Each GO term has 2 chords: one to host division, one to symbiont division
#
# Author: Claude (for David Armstrong)
# Date: 2025-02-06
# =============================================================================

library(tidyverse)
library(circlize)
library(gt)

# =============================================================================
# CONFIGURATION
# =============================================================================

host_dir <- "/home/darmstrong4/mc_rework/10_GO_MWU/output"
symbiont_dir <- "/home/darmstrong4/mc_rework/11_symbiont_GO_MWU/output"
output_dir <- "/home/darmstrong4/mc_rework/12_publication_figures/figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

sig_threshold <- 0.1

# Interaction type colors (for outer ring)
interaction_colors <- c(
    "Synergistic Up" = "#E41A1C",
    "Synergistic Down" = "#377EB8",
    "Antagonistic (Host Up)" = "#90EE90",
    "Antagonistic (Host Down)" = "#228B22"
)

# Division colors
division_colors <- c(
    "BP" = "#FDB462",
    "CC" = "#80B1D3",
    "MF" = "#FB8072"
)

# Organism-division colors (for right side nodes)
org_division_colors <- c(
    "Host_BP" = "#FDB462",
    "Host_CC" = "#80B1D3", 
    "Host_MF" = "#FB8072",
    "Symbiont_BP" = "#B3DE69",
    "Symbiont_CC" = "#BEBADA",
    "Symbiont_MF" = "#FCCDE5"
)

# =============================================================================
# HELPER: Load significant GO terms
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
        
        if (!file.exists(mwu_file)) next
        
        mwu <- read.table(mwu_file, header = TRUE, stringsAsFactors = FALSE,
                          quote = "\"", fill = TRUE)
        
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
    
    bind_rows(all_terms)
}

# =============================================================================
# HELPER: Classify shared terms
# =============================================================================

classify_shared_terms <- function(host_terms, symbiont_terms) {
    
    shared_ids <- intersect(host_terms$go_id, symbiont_terms$go_id)
    if (length(shared_ids) == 0) return(data.frame())
    
    shared_data <- data.frame(go_id = shared_ids, stringsAsFactors = FALSE)
    
    shared_data <- shared_data %>%
        left_join(host_terms %>% select(go_id, name, division, direction, p.adj) %>% 
                      rename(host_direction = direction, host_division = division, 
                             term_name = name, host_padj = p.adj),
                  by = "go_id") %>%
        left_join(symbiont_terms %>% select(go_id, direction, division, p.adj) %>%
                      rename(symbiont_direction = direction, symbiont_division = division,
                             symbiont_padj = p.adj),
                  by = "go_id")
    
    shared_data <- shared_data %>%
        mutate(
            interaction_type = case_when(
                host_direction == "up" & symbiont_direction == "up" ~ "Synergistic Up",
                host_direction == "down" & symbiont_direction == "down" ~ "Synergistic Down",
                host_direction == "up" & symbiont_direction == "down" ~ "Antagonistic (Host Up)",
                host_direction == "down" & symbiont_direction == "up" ~ "Antagonistic (Host Down)"
            ),
            # Order factor for sorting
            interaction_order = case_when(
                interaction_type == "Synergistic Up" ~ 1,
                interaction_type == "Synergistic Down" ~ 2,
                interaction_type == "Antagonistic (Host Up)" ~ 3,
                interaction_type == "Antagonistic (Host Down)" ~ 4
            )
        ) %>%
        arrange(interaction_order, host_division, go_id)
    
    return(shared_data)
}

# =============================================================================
# HELPER: Create GT table for supplementary
# =============================================================================

create_supplementary_table <- function(shared_terms, season, output_dir) {
    
    table_data <- shared_terms %>%
        select(
            `GO ID` = go_id,
            `Term Name` = term_name,
            `Interaction Type` = interaction_type,
            `Host Division` = host_division,
            `Host Direction` = host_direction,
            `Host p.adj` = host_padj,
            `Symbiont Division` = symbiont_division,
            `Symbiont Direction` = symbiont_direction,
            `Symbiont p.adj` = symbiont_padj
        )
    
    gt_table <- table_data %>%
        gt() %>%
        tab_header(
            title = paste0(tools::toTitleCase(season), " - Shared GO Terms Between Host and Symbiont"),
            subtitle = paste0("n = ", nrow(table_data), " shared significant GO terms (p.adj < ", sig_threshold, ")")
        ) %>%
        tab_spanner(
            label = "Host (M. capitata)",
            columns = c(`Host Division`, `Host Direction`, `Host p.adj`)
        ) %>%
        tab_spanner(
            label = "Symbiont (D. trenchii)",
            columns = c(`Symbiont Division`, `Symbiont Direction`, `Symbiont p.adj`)
        ) %>%
        fmt_scientific(
            columns = c(`Host p.adj`, `Symbiont p.adj`),
            decimals = 2
        ) %>%
        data_color(
            columns = `Interaction Type`,
            colors = scales::col_factor(
                palette = c("#E41A1C", "#377EB8", "#90EE90", "#228B22"),
                domain = c("Synergistic Up", "Synergistic Down", 
                          "Antagonistic (Host Up)", "Antagonistic (Host Down)")
            )
        ) %>%
        data_color(
            columns = `Host Division`,
            colors = scales::col_factor(
                palette = c("#FDB462", "#80B1D3", "#FB8072"),
                domain = c("BP", "CC", "MF")
            )
        ) %>%
        data_color(
            columns = `Symbiont Division`,
            colors = scales::col_factor(
                palette = c("#FDB462", "#80B1D3", "#FB8072"),
                domain = c("BP", "CC", "MF")
            )
        ) %>%
        tab_style(
            style = cell_text(weight = "bold"),
            locations = cells_column_labels()
        ) %>%
        tab_footnote(
            footnote = "BP = Biological Process, CC = Cellular Component, MF = Molecular Function",
            locations = cells_column_labels(columns = `Host Division`)
        ) %>%
        tab_options(
            table.font.size = 11,
            heading.title.font.size = 14,
            heading.subtitle.font.size = 12
        )
    
    # Save as HTML
    html_file <- file.path(output_dir, paste0("Table_S_", season, "_Shared_GO_Terms.html"))
    gtsave(gt_table, html_file)
    cat("  Saved:", basename(html_file), "\n")
    
    # Save as CSV
    csv_file <- file.path(output_dir, paste0("Table_S_", season, "_Shared_GO_Terms.csv"))
    write.csv(table_data, csv_file, row.names = FALSE)
    cat("  Saved:", basename(csv_file), "\n")
    
    return(gt_table)
}

# =============================================================================
# MAIN: Create chord diagram
# =============================================================================

create_chord_diagram <- function(season, output_file) {
    
    cat("\n=== Creating chord diagram for", season, "===\n")
    
    # Load data
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
    
    if (nrow(shared_terms) == 0) {
        cat("  No shared terms - skipping\n")
        return(NULL)
    }
    
    # Print summary
    cat("    - Synergistic Up:", sum(shared_terms$interaction_type == "Synergistic Up"), "\n")
    cat("    - Synergistic Down:", sum(shared_terms$interaction_type == "Synergistic Down"), "\n")
    cat("    - Antagonistic (Host Up):", sum(shared_terms$interaction_type == "Antagonistic (Host Up)"), "\n")
    cat("    - Antagonistic (Host Down):", sum(shared_terms$interaction_type == "Antagonistic (Host Down)"), "\n")
    
    # Create supplementary table
    cat("  Creating supplementary table...\n")
    create_supplementary_table(shared_terms, season, output_dir)
    
    # ==========================================================================
    # Build chord diagram data
    # ==========================================================================
    
    n_terms <- nrow(shared_terms)
    
    # Sectors: GO terms on left, organism-divisions on right
    # GO term sectors are small and equal-sized
    # Division sectors are sized by how many connections they receive
    
    # Count connections to each division
    host_div_counts <- shared_terms %>% count(host_division) %>% 
        mutate(sector = paste0("Host_", host_division))
    symbiont_div_counts <- shared_terms %>% count(symbiont_division) %>% 
        mutate(sector = paste0("Symbiont_", symbiont_division))
    
    # Create sector sizes
    go_term_size <- 1  # Each GO term gets size 1
    
    # Sector order: GO terms (ordered by interaction type), then divisions
    go_sectors <- shared_terms$go_id
    
    # Right side: Host divisions then Symbiont divisions
    right_sectors <- c("Host_BP", "Host_CC", "Host_MF", "Symbiont_MF", "Symbiont_CC", "Symbiont_BP")
    right_sectors <- right_sectors[right_sectors %in% c(host_div_counts$sector, symbiont_div_counts$sector)]
    
    all_sectors <- c(go_sectors, right_sectors)
    
    # Sector sizes
    sector_sizes <- c(
        setNames(rep(go_term_size, n_terms), go_sectors),
        setNames(host_div_counts$n, host_div_counts$sector),
        setNames(symbiont_div_counts$n, symbiont_div_counts$sector)
    )
    sector_sizes <- sector_sizes[all_sectors]
    
    # Sector colors
    go_colors <- setNames(division_colors[shared_terms$host_division], go_sectors)
    sector_colors <- c(go_colors, org_division_colors[right_sectors])
    
    # ==========================================================================
    # Create the diagram
    # ==========================================================================
    
    cat("  Generating chord diagram...\n")
    
    pdf(output_file, width = 12, height = 10)
    
    circos.clear()
    
    # Calculate gaps
    # Small gaps between GO terms, larger gap between GO terms and divisions
    n_go <- length(go_sectors)
    n_right <- length(right_sectors)
    
    gaps <- c(
        rep(0.5, n_go - 1),  # Small gaps between GO terms
        15,                   # Large gap between GO terms and divisions
        rep(3, n_right - 1), # Medium gaps between divisions
        15                    # Large gap back to GO terms
    )
    
    circos.par(
        start.degree = 90,
        gap.degree = gaps,
        track.margin = c(0.01, 0.01),
        cell.padding = c(0, 0, 0, 0)
    )
    
    # Initialize
    xlim_matrix <- cbind(rep(0, length(all_sectors)), sector_sizes)
    rownames(xlim_matrix) <- all_sectors
    
    circos.initialize(
        factors = factor(all_sectors, levels = all_sectors),
        xlim = xlim_matrix
    )
    
    # --------------------------------------------------------------------------
    # Track 1: Outer ring - Interaction type (only for GO terms)
    # --------------------------------------------------------------------------
    
    circos.track(
        ylim = c(0, 1),
        track.height = 0.06,
        bg.border = NA,
        panel.fun = function(x, y) {
            sector.name <- CELL_META$sector.index
            
            # Only draw for GO term sectors
            if (sector.name %in% go_sectors) {
                interaction <- shared_terms$interaction_type[shared_terms$go_id == sector.name]
                bg_col <- interaction_colors[interaction]
                
                circos.rect(
                    CELL_META$xlim[1], 0,
                    CELL_META$xlim[2], 1,
                    col = bg_col,
                    border = bg_col
                )
            }
        }
    )
    
    # --------------------------------------------------------------------------
    # Track 2: Main sectors with labels
    # --------------------------------------------------------------------------
    
    circos.track(
        ylim = c(0, 1),
        track.height = 0.08,
        bg.col = sector_colors[all_sectors],
        bg.border = "grey50",
        bg.lwd = 0.5,
        panel.fun = function(x, y) {
            sector.name <- CELL_META$sector.index
            
            # Label GO terms with GO ID
            if (sector.name %in% go_sectors) {
                # Short label - just the number part
                short_id <- sub("GO:", "", sector.name)
                circos.text(
                    CELL_META$xcenter, 0.5,
                    as.character(short_id),
                    facing = "clockwise",
                    niceFacing = TRUE,
                    cex = 0.35,
                    adj = c(0, 0.5)
                )
            } else {
                # Division labels
                div_name <- sub(".*_", "", sector.name)
                org_name <- sub("_.*", "", sector.name)
                circos.text(
                    CELL_META$xcenter, 0.5,
                    as.character(div_name),
                    facing = "bending.inside",
                    niceFacing = TRUE,
                    cex = 0.8,
                    font = 2
                )
            }
        }
    )
    
    # --------------------------------------------------------------------------
    # Track 3: Organism labels for right side
    # --------------------------------------------------------------------------
    
    circos.track(
        ylim = c(0, 1),
        track.height = 0.05,
        bg.border = NA,
        panel.fun = function(x, y) {
            sector.name <- CELL_META$sector.index
            
            if (!(sector.name %in% go_sectors)) {
                org_name <- sub("_.*", "", sector.name)
                if (org_name == "Host") {
                    label <- "Host"
                } else {
                    label <- "Symbiont"
                }
                circos.text(
                    CELL_META$xcenter, 0.5,
                    as.character(label),
                    facing = "bending.inside",
                    niceFacing = TRUE,
                    cex = 0.6,
                    col = "grey30"
                )
            }
        }
    )
    
    # --------------------------------------------------------------------------
    # Draw chords
    # --------------------------------------------------------------------------
    
    # Track positions for right side sectors
    right_positions <- list()
    for (s in right_sectors) {
        right_positions[[s]] <- 0
    }
    
    for (i in 1:n_terms) {
        term <- shared_terms[i, ]
        go_id <- term$go_id
        
        host_sector <- paste0("Host_", term$host_division)
        symbiont_sector <- paste0("Symbiont_", term$symbiont_division)
        
        # Chord color based on interaction type (with transparency)
        chord_col <- adjustcolor(interaction_colors[term$interaction_type], alpha.f = 0.5)
        chord_border <- adjustcolor(interaction_colors[term$interaction_type], alpha.f = 0.7)
        
        # Draw chord to host division
        if (host_sector %in% right_sectors) {
            host_pos <- right_positions[[host_sector]]
            circos.link(
                go_id, c(0, go_term_size),
                host_sector, c(host_pos, host_pos + 1),
                col = chord_col,
                border = chord_border,
                lwd = 0.5
            )
            right_positions[[host_sector]] <- host_pos + 1
        }
        
        # Draw chord to symbiont division
        if (symbiont_sector %in% right_sectors) {
            symbiont_pos <- right_positions[[symbiont_sector]]
            circos.link(
                go_id, c(0, go_term_size),
                symbiont_sector, c(symbiont_pos, symbiont_pos + 1),
                col = chord_col,
                border = chord_border,
                lwd = 0.5
            )
            right_positions[[symbiont_sector]] <- symbiont_pos + 1
        }
    }
    
    # --------------------------------------------------------------------------
    # Legends and annotations
    # --------------------------------------------------------------------------
    
    # Title
    title(main = paste0(tools::toTitleCase(season), " - Shared GO Terms"),
          line = -2, cex.main = 1.4, font.main = 2)
    
    # Interaction type legend (top left)
    legend_x <- -1.15
    legend_y <- 0.95
    legend_spacing <- 0.08
    
    text(legend_x, legend_y + 0.08, "Interaction Type", cex = 0.75, font = 2, adj = c(0, 0.5))
    
    for (i in seq_along(interaction_colors)) {
        y_pos <- legend_y - (i - 1) * legend_spacing
        points(legend_x, y_pos, pch = 15, col = interaction_colors[i], cex = 1.5)
        text(legend_x + 0.05, y_pos, names(interaction_colors)[i], cex = 0.65, adj = c(0, 0.5))
    }
    
    # Division legend (bottom left)
    legend_y2 <- legend_y - 5 * legend_spacing
    text(legend_x, legend_y2 + 0.08, "GO Division", cex = 0.75, font = 2, adj = c(0, 0.5))
    
    for (i in seq_along(division_colors)) {
        y_pos <- legend_y2 - (i - 1) * legend_spacing
        points(legend_x, y_pos, pch = 15, col = division_colors[i], cex = 1.5)
        text(legend_x + 0.05, y_pos, names(division_colors)[i], cex = 0.65, adj = c(0, 0.5))
    }
    
    # Organism labels
    text(0.85, 0.6, expression(italic("M. capitata")), cex = 0.9, font = 2)
    text(0.85, 0.55, "(Host)", cex = 0.75)
    
    text(0.85, -0.55, expression(italic("D. trenchii")), cex = 0.9, font = 2)
    text(0.85, -0.6, "(Symbiont)", cex = 0.75)
    
    # N terms annotation
    text(0, -1.15, paste0("n = ", n_terms, " shared GO terms"), cex = 0.8)
    
    circos.clear()
    dev.off()
    
    cat("  SUCCESS:", basename(output_file), "\n")
    
    return(list(
        season = season,
        host_total = nrow(host_terms),
        symbiont_total = nrow(symbiont_terms),
        shared_total = n_terms,
        shared_details = shared_terms
    ))
}

# =============================================================================
# GENERATE DIAGRAMS
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("GENERATING CHORD DIAGRAMS AND SUPPLEMENTARY TABLES\n")
cat(strrep("=", 60), "\n")

summer_results <- create_chord_diagram(
    season = "summer",
    output_file = file.path(output_dir, "Fig_Summer_Host_Symbiont_GO_Chord.pdf")
)

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

cat("\nOutput files in:", output_dir, "\n")
