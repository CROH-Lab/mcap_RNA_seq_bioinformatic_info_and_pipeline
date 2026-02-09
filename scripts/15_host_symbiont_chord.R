#!/usr/bin/env Rscript
# =============================================================================
# Chord Diagram: Host-Symbiont Shared GO Terms (Simplified)
# =============================================================================
# Layout:
# - LEFT: Individual GO terms (labeled by GO ID), colored by interaction type
# - RIGHT: 3 nodes (BP, MF, CC)
# - Chords connect GO terms to their division, colored by interaction type
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

# Interaction type colors
interaction_colors <- c(
    "Synergistic Up" = "#E41A1C",
    "Synergistic Down" = "#377EB8",
    "Antagonistic (Host Up)" = "#90EE90",
    "Antagonistic (Host Down)" = "#228B22"
)

# Division colors (neutral grey tones for the 3 nodes)
division_colors <- c(
    "BP" = "#D9D9D9",
    "CC" = "#BDBDBD",
    "MF" = "#969696"
)

# Division full names
division_names <- c(
    "BP" = "Biological\nProcess",
    "CC" = "Cellular\nComponent",
    "MF" = "Molecular\nFunction"
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
            interaction_order = case_when(
                interaction_type == "Synergistic Up" ~ 1,
                interaction_type == "Synergistic Down" ~ 2,
                interaction_type == "Antagonistic (Host Up)" ~ 3,
                interaction_type == "Antagonistic (Host Down)" ~ 4
            ),
            # Use host division as the primary division for connecting
            division = host_division
        ) %>%
        arrange(interaction_order, division, go_id)
    
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
    
    html_file <- file.path(output_dir, paste0("Table_S_", season, "_Shared_GO_Terms.html"))
    gtsave(gt_table, html_file)
    cat("  Saved:", basename(html_file), "\n")
    
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
    
    cat("  Identifying shared terms...\n")
    shared_terms <- classify_shared_terms(host_terms, symbiont_terms)
    cat("    Shared terms:", nrow(shared_terms), "\n")
    
    if (nrow(shared_terms) == 0) {
        cat("  No shared terms - skipping\n")
        return(NULL)
    }
    
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
    
    # Count connections to each division
    div_counts <- shared_terms %>% count(division)
    
    # Sector sizes
    go_term_size <- 1
    
    # Sector order: Divisions on LEFT (BP, CC, MF), GO terms on RIGHT
    # But we want GO terms on LEFT - so we reverse the start degree
    division_sectors <- c("BP", "CC", "MF")
    division_sectors <- division_sectors[division_sectors %in% div_counts$division]
    go_sectors <- shared_terms$go_id
    
    # Order: GO terms first (will be on left with start.degree = 270), then divisions
    all_sectors <- c(go_sectors, division_sectors)
    
    # Sector sizes
    sector_sizes <- c(
        setNames(rep(go_term_size, n_terms), go_sectors),
        setNames(div_counts$n[match(division_sectors, div_counts$division)], division_sectors)
    )
    
    # Sector colors - GO terms colored by interaction type
    go_colors <- setNames(
        interaction_colors[shared_terms$interaction_type],
        go_sectors
    )
    sector_colors <- c(go_colors, division_colors[division_sectors])
    
    # ==========================================================================
    # Create the diagram
    # ==========================================================================
    
    cat("  Generating chord diagram...\n")
    
    pdf(output_file, width = 11, height = 9)
    
    circos.clear()
    
    # Gaps: small between GO terms, large between GO terms and divisions
    n_go <- length(go_sectors)
    n_div <- length(division_sectors)
    
    gaps <- c(
        rep(0.3, n_go - 1),   # Tiny gaps between GO terms
        20,                    # Large gap between GO terms and divisions
        rep(4, n_div - 1),    # Medium gaps between divisions
        20                     # Large gap back to GO terms
    )
    
    circos.par(
        start.degree = 270,    # Start from bottom, GO terms will be on left
        gap.degree = gaps,
        track.margin = c(0.005, 0.005),
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
    # Track 1: Main colored sectors
    # --------------------------------------------------------------------------
    
    circos.track(
        ylim = c(0, 1),
        track.height = 0.08,
        bg.col = sector_colors[all_sectors],
        bg.border = "grey40",
        bg.lwd = 0.5,
        panel.fun = function(x, y) {
            sector.name <- CELL_META$sector.index
            
            if (sector.name %in% go_sectors) {
                # GO ID label (just the number)
                short_id <- sub("GO:", "", sector.name)
                circos.text(
                    CELL_META$xcenter, 0.5,
                    as.character(short_id),
                    facing = "clockwise",
                    niceFacing = TRUE,
                    cex = 0.3,
                    adj = c(0, 0.5)
                )
            }
        }
    )
    
    # --------------------------------------------------------------------------
    # Track 2: Division labels (only for division sectors)
    # --------------------------------------------------------------------------
    
    circos.track(
        ylim = c(0, 1),
        track.height = 0.12,
        bg.border = NA,
        panel.fun = function(x, y) {
            sector.name <- CELL_META$sector.index
            
            if (sector.name %in% division_sectors) {
                circos.text(
                    CELL_META$xcenter, 0.5,
                    as.character(division_names[sector.name]),
                    facing = "bending.inside",
                    niceFacing = TRUE,
                    cex = 0.85,
                    font = 2
                )
            }
        }
    )
    
    # --------------------------------------------------------------------------
    # Draw chords
    # --------------------------------------------------------------------------
    
    # Track positions for division sectors
    div_positions <- setNames(rep(0, length(division_sectors)), division_sectors)
    
    for (i in 1:n_terms) {
        term <- shared_terms[i, ]
        go_id <- term$go_id
        div <- term$division
        
        # Chord color based on interaction type
        chord_col <- adjustcolor(interaction_colors[term$interaction_type], alpha.f = 0.6)
        chord_border <- adjustcolor(interaction_colors[term$interaction_type], alpha.f = 0.8)
        
        # Draw chord from GO term to division
        div_pos <- div_positions[[div]]
        circos.link(
            go_id, c(0, go_term_size),
            div, c(div_pos, div_pos + 1),
            col = chord_col,
            border = chord_border,
            lwd = 0.3
        )
        div_positions[[div]] <- div_pos + 1
    }
    
    # --------------------------------------------------------------------------
    # Legend
    # --------------------------------------------------------------------------
    
    # Interaction type legend
    legend_x <- -1.2
    legend_y <- 0.5
    legend_spacing <- 0.1
    
    text(legend_x, legend_y + 0.15, "Interaction Type", cex = 0.85, font = 2, adj = c(0, 0.5))
    
    for (i in seq_along(interaction_colors)) {
        y_pos <- legend_y - (i - 1) * legend_spacing
        points(legend_x, y_pos, pch = 15, col = interaction_colors[i], cex = 2)
        text(legend_x + 0.06, y_pos, names(interaction_colors)[i], cex = 0.7, adj = c(0, 0.5))
    }
    
    # Title
    title(main = paste0(tools::toTitleCase(season), " - Shared GO Terms Between Host and Symbiont"),
          line = -1.5, cex.main = 1.3, font.main = 2)
    
    # Subtitle with count
    mtext(paste0("n = ", n_terms, " shared GO terms"), side = 1, line = -1.5, cex = 0.9)
    
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
