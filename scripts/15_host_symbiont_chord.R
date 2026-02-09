#!/usr/bin/env Rscript
# =============================================================================
# Chord Diagram: Host-Symbiont Shared GO Terms (Final Version)
# =============================================================================
# Layout:
# - Summer and Winter stacked vertically on one plot
# - GO terms on left, divisions (labels only) on right
# - Chord width proportional to significant gene count
# - GO IDs outside nodes, larger text
# - Single shared legend
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
    "Synergistic Up" = "#FBD1A2",
    "Synergistic Down" = "#7DCFB6",
    "Antagonistic (Host Up)" = "#F79256",
    "Antagonistic (Host Down)" = "#6BBF59"
)

# Division full names
division_names <- c(
    "BP" = "Biological Process",
    "CC" = "Cellular Component",
    "MF" = "Molecular Function"
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
                go_id = term,
                # Clean GO ID - take first ID if compound, ensure GO: prefix
                go_id_clean = sapply(strsplit(as.character(term), ";"), function(x) {
                    id <- x[1]
                    if (!grepl("^GO:", id)) id <- paste0("GO:", id)
                    return(id)
                })
            ) %>%
            select(go_id, go_id_clean, name, division, direction, delta.rank, pval, p.adj, nseqs)
        
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
        left_join(host_terms %>% select(go_id, go_id_clean, name, division, direction, p.adj, nseqs) %>% 
                      rename(host_direction = direction, host_division = division, 
                             term_name = name, host_padj = p.adj, host_nseqs = nseqs),
                  by = "go_id") %>%
        left_join(symbiont_terms %>% select(go_id, direction, division, p.adj, nseqs) %>%
                      rename(symbiont_direction = direction, symbiont_division = division,
                             symbiont_padj = p.adj, symbiont_nseqs = nseqs),
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
            division = host_division,
            # Total gene count for chord width
            total_genes = host_nseqs + symbiont_nseqs
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
            `GO ID` = go_id_clean,
            `Term Name` = term_name,
            `Interaction Type` = interaction_type,
            `Host Division` = host_division,
            `Host Direction` = host_direction,
            `Host DEGs` = host_nseqs,
            `Host p.adj` = host_padj,
            `Symbiont Division` = symbiont_division,
            `Symbiont Direction` = symbiont_direction,
            `Symbiont DEGs` = symbiont_nseqs,
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
            columns = c(`Host Division`, `Host Direction`, `Host DEGs`, `Host p.adj`)
        ) %>%
        tab_spanner(
            label = "Symbiont (D. trenchii)",
            columns = c(`Symbiont Division`, `Symbiont Direction`, `Symbiont DEGs`, `Symbiont p.adj`)
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
            footnote = "BP = Biological Process, CC = Cellular Component, MF = Molecular Function; DEGs = Differentially Expressed Genes",
            locations = cells_column_labels(columns = `Host Division`)
        ) %>%
        tab_options(
            table.font.size = 11,
            heading.title.font.size = 14,
            heading.subtitle.font.size = 12
        )
    
    html_file <- file.path(output_dir, paste0("Table_S_", season, "_Shared_GO_Terms.html"))
    gtsave(gt_table, html_file)
    cat("    Saved:", basename(html_file), "\n")
    
    csv_file <- file.path(output_dir, paste0("Table_S_", season, "_Shared_GO_Terms.csv"))
    write.csv(table_data, csv_file, row.names = FALSE)
    cat("    Saved:", basename(csv_file), "\n")
    
    return(gt_table)
}

# =============================================================================
# HELPER: Prepare data for one season
# =============================================================================

prepare_season_data <- function(season) {
    
    cat("  Loading", season, "data...\n")
    
    host_terms <- load_significant_terms(host_dir, "", season, sig_threshold)
    symbiont_terms <- load_significant_terms(symbiont_dir, "symbiont", season, sig_threshold)
    shared_terms <- classify_shared_terms(host_terms, symbiont_terms)
    
    cat("    Host terms:", nrow(host_terms), "\n")
    cat("    Symbiont terms:", nrow(symbiont_terms), "\n")
    cat("    Shared terms:", nrow(shared_terms), "\n")
    
    if (nrow(shared_terms) > 0) {
        cat("    - Synergistic Up:", sum(shared_terms$interaction_type == "Synergistic Up"), "\n")
        cat("    - Synergistic Down:", sum(shared_terms$interaction_type == "Synergistic Down"), "\n")
        cat("    - Antagonistic (Host Up):", sum(shared_terms$interaction_type == "Antagonistic (Host Up)"), "\n")
        cat("    - Antagonistic (Host Down):", sum(shared_terms$interaction_type == "Antagonistic (Host Down)"), "\n")
        
        # Create supplementary table
        create_supplementary_table(shared_terms, season, output_dir)
    }
    
    return(list(
        host_terms = host_terms,
        symbiont_terms = symbiont_terms,
        shared_terms = shared_terms
    ))
}

# =============================================================================
# HELPER: Draw single chord diagram
# =============================================================================

draw_chord_diagram <- function(shared_terms, season_label, show_legend = FALSE) {
    
    if (nrow(shared_terms) == 0) {
        plot.new()
        text(0.5, 0.5, paste0(season_label, "\nNo shared terms"), cex = 1.2)
        return()
    }
    
    n_terms <- nrow(shared_terms)
    
    # Count connections and total genes to each division
    div_stats <- shared_terms %>%
        group_by(division) %>%
        summarise(
            n_terms = n(),
            total_width = sum(sqrt(total_genes + 1)),
            .groups = "drop"
        )
    
    # Sector order: GO terms first, then divisions
    division_sectors <- c("BP", "CC", "MF")
    division_sectors <- division_sectors[division_sectors %in% div_stats$division]
    go_sectors <- shared_terms$go_id
    
    all_sectors <- c(go_sectors, division_sectors)
    
    # Sector sizes
    go_sizes <- setNames(sqrt(shared_terms$total_genes + 1), go_sectors)
    div_sizes <- setNames(
        div_stats$total_width[match(division_sectors, div_stats$division)],
        division_sectors
    )
    sector_sizes <- c(go_sizes, div_sizes)
    
    # GO term colors
    go_colors <- setNames(
        interaction_colors[shared_terms$interaction_type],
        go_sectors
    )
    
    # Initialize circos
    circos.clear()
    
    # Gaps
    n_go <- length(go_sectors)
    n_div <- length(division_sectors)
    
    go_gap <- ifelse(n_terms > 50, 0.9, 1.9)
    
    gaps <- c(
        rep(go_gap, n_go - 1),
        20,
        rep(6, n_div - 1),
        20
    )
    
    circos.par(
        start.degree = 270,
        gap.degree = gaps,
        track.margin = c(0.002, 0.002),
        cell.padding = c(0, 0, 0, 0),
        canvas.xlim = c(-0.85, 1.2),
        canvas.ylim = c(-0.9, 0.9)
    )
    
    # Initialize sectors
    xlim_matrix <- cbind(rep(0, length(all_sectors)), sector_sizes)
    rownames(xlim_matrix) <- all_sectors
    
    circos.initialize(
        factors = factor(all_sectors, levels = all_sectors),
        xlim = xlim_matrix
    )
    
    # --------------------------------------------------------------------------
    # Single track: Labels only (no nodes)
    # --------------------------------------------------------------------------
    
    circos.track(
        ylim = c(0, 1),
        track.height = 0.25,
        bg.border = NA,
        panel.fun = function(x, y) {
            sector.name <- CELL_META$sector.index
            
            if (sector.name %in% go_sectors) {
                # GO ID label
                clean_id <- shared_terms$go_id_clean[shared_terms$go_id == sector.name]
                label_cex <- ifelse(n_terms > 80, 0.48, ifelse(n_terms > 50, 0.45, 0.6))
                circos.text(
                    CELL_META$xcenter, 0.1,
                    as.character(clean_id),
                    facing = "clockwise",
                    niceFacing = TRUE,
                    cex = label_cex,
                    adj = c(0, 0.5),
                    font = 1
                )
            } else {
                # Division labels
                short_names <- c(
                    "BP" = "Biological Process",
                    "CC" = "Cellular Component", 
                    "MF" = "Molecular Function"
                )
                circos.text(
                    CELL_META$xcenter, 0.5,
                    as.character(short_names[sector.name]),
                    facing = "bending.inside",
                    niceFacing = TRUE,
                    cex = 0.75,
                    font = 2
                )
            }
        }
    )
    
    # --------------------------------------------------------------------------
    # Draw chords
    # --------------------------------------------------------------------------
    
    div_positions <- setNames(rep(0, length(division_sectors)), division_sectors)
    
    for (i in 1:n_terms) {
        term <- shared_terms[i, ]
        go_id <- term$go_id
        div <- term$division
        
        chord_col <- adjustcolor(interaction_colors[term$interaction_type], alpha.f = 0.55)
        chord_border <- adjustcolor(interaction_colors[term$interaction_type], alpha.f = 0.75)
        
        chord_width <- sqrt(term$total_genes + 1)
        div_pos <- div_positions[[div]]
        
        circos.link(
            go_id, c(0, sector_sizes[go_id]),
            div, c(div_pos, div_pos + chord_width),
            col = chord_col,
            border = chord_border,
            lwd = 0.2
        )
        
        div_positions[[div]] <- div_pos + chord_width
    }
    
    # Season label (center)
    text(0, 0.02, season_label, cex = 1.1, font = 2)
    text(0, -0.12, paste0("n = ", n_terms), cex = 0.65)
    
    # Legend (only if requested)
    if (show_legend) {
        legend_x <- 0.87
        legend_y <- 0.80
        legend_spacing <- 0.05
        
        text(legend_x, legend_y + 0.08, "Interaction Type", cex = 0.6, font = 2, adj = c(0, 0.5))
        
        short_labels <- c(
            "Synergistic Up" = "Synergy All Up",
            "Synergistic Down" = "Synergy All Down", 
            "Antagonistic (Host Up)" = "Antag. (Host Up / Sym Down)",
            "Antagonistic (Host Down)" = "Antag. (Host Down / Sym. Up)"
        )
        
        for (i in seq_along(interaction_colors)) {
            y_pos <- legend_y - (i - 1) * legend_spacing
            points(legend_x, y_pos, pch = 15, col = interaction_colors[i], cex = 1.3)
            text(legend_x + 0.08, y_pos, short_labels[names(interaction_colors)[i]], cex = 0.5, adj = c(0, 0.5))
        }
    }
    
    circos.clear()
}

# =============================================================================
# MAIN: Create combined stacked diagram
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("GENERATING COMBINED CHORD DIAGRAM\n")
cat(strrep("=", 60), "\n")

# Prepare data for both seasons
summer_data <- prepare_season_data("summer")
winter_data <- prepare_season_data("winter")

# Create combined figure
output_file <- file.path(output_dir, "Fig_Shared_GO_Chord_Diagram.pdf")

cat("\n  Creating combined figure...\n")

pdf(output_file, width = 8, height = 14)

# Layout: 2 rows stacked vertically
layout(matrix(c(1, 2), nrow = 2, ncol = 1), heights = c(1, 1))

# Summer diagram (with legend)
par(mar = c(0.5, 0.5, 0.5, 0.5))
draw_chord_diagram(summer_data$shared_terms, "Summer", show_legend = TRUE)

# Winter diagram (no legend)
par(mar = c(0.5, 0.5, 0.5, 0.5))
draw_chord_diagram(winter_data$shared_terms, "Winter", show_legend = FALSE)

dev.off()

cat("  SUCCESS:", basename(output_file), "\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("SUMMARY\n")
cat(strrep("=", 60), "\n")

cat("\nSummer:\n")
cat("  Host significant terms:", nrow(summer_data$host_terms), "\n")
cat("  Symbiont significant terms:", nrow(summer_data$symbiont_terms), "\n")
cat("  Shared terms:", nrow(summer_data$shared_terms), "\n")

cat("\nWinter:\n")
cat("  Host significant terms:", nrow(winter_data$host_terms), "\n")
cat("  Symbiont significant terms:", nrow(winter_data$symbiont_terms), "\n")
cat("  Shared terms:", nrow(winter_data$shared_terms), "\n")

cat("\nOutput files:\n")
cat("  ", basename(output_file), "\n")
cat("  Table_S_summer_Shared_GO_Terms.html/csv\n")
cat("  Table_S_winter_Shared_GO_Terms.html/csv\n")
cat("\nOutput directory:", output_dir, "\n")
