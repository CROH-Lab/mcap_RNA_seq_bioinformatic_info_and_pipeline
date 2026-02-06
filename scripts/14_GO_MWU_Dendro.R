#!/usr/bin/env Rscript
# =============================================================================
# Plot Combined Representative GO Term Dendrograms
# =============================================================================
# Creates clean, readable dendrograms with BP, CC, MF stacked vertically
# One figure per season/organism combination (4 total)
# 
# Author: Claude (for David Armstrong)
# Date: 2025-02-06
# =============================================================================

library(ape)
library(tidyverse)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Directories
host_dir <- "/home/darmstrong4/mc_rework/10_GO_MWU/output"
symbiont_dir <- "/home/darmstrong4/mc_rework/11_symbiont_GO_MWU/output"

# Output directory
output_dir <- "/home/darmstrong4/mc_rework/12_publication_figures/figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Significance thresholds (match GO_MWU defaults)
level1 <- 0.1   # Italic text
level2 <- 0.05  # Normal text  
level3 <- 0.01  # Bold text

# Colors (GO_MWU defaults)
colors <- c("dodgerblue2", "firebrick1", "skyblue2", "lightcoral")

# Division labels
division_labels <- c("BP" = "Biological Process", 
                     "CC" = "Cellular Component", 
                     "MF" = "Molecular Function")

# =============================================================================
# HELPER FUNCTION: Load and prepare data for one division
# =============================================================================

load_division_data <- function(rep_file, mwu_file, dissim_file, division) {
    
    # Check files exist
    if (!file.exists(rep_file) || !file.exists(mwu_file) || !file.exists(dissim_file)) {
        cat("    Missing files for", division, "\n")
        return(NULL)
    }
    
    # Load representative GOs
    rep_gos <- read.csv(rep_file, stringsAsFactors = FALSE)
    if (nrow(rep_gos) < 1) {
        cat("    No representative GOs for", division, "\n")
        return(NULL)
    }
    
    # Load MWU results
    mwu <- read.table(mwu_file, header = TRUE, stringsAsFactors = FALSE,
                      quote = "\"", fill = TRUE)
    
    # Match terms to get GO IDs
    rep_gos$go_id <- NA
    rep_gos$delta_rank <- NA
    rep_gos$nseqs_total <- NA
    
    for (i in 1:nrow(rep_gos)) {
        match_idx <- which(mwu$name == rep_gos$term[i])
        if (length(match_idx) == 0) {
            match_idx <- which(tolower(mwu$name) == tolower(rep_gos$term[i]))
        }
        if (length(match_idx) > 0) {
            rep_gos$go_id[i] <- mwu$term[match_idx[1]]
            rep_gos$delta_rank[i] <- mwu$delta.rank[match_idx[1]]
            rep_gos$nseqs_total[i] <- mwu$nseqs[match_idx[1]]
        }
    }
    
    rep_gos <- rep_gos[!is.na(rep_gos$go_id), ]
    if (nrow(rep_gos) < 2) {
        cat("    Less than 2 matched terms for", division, "\n")
        return(NULL)
    }
    
    # Load dissimilarity matrix
    diss <- read.table(dissim_file, sep = "\t", header = TRUE, check.names = FALSE)
    if (!any(grepl("^GO:", rownames(diss)))) {
        rownames(diss) <- colnames(diss)
    }
    
    # Find matching GO IDs (handle compound IDs)
    available_ids <- c()
    id_mapping <- list()
    
    for (go_id in rep_gos$go_id) {
        if (is.na(go_id)) next
        if (go_id %in% colnames(diss)) {
            available_ids <- c(available_ids, go_id)
            id_mapping[[go_id]] <- go_id
        } else {
            compound_match <- grep(paste0("(^|;)", go_id, "(;|$)"), colnames(diss), value = TRUE)
            if (length(compound_match) > 0) {
                available_ids <- c(available_ids, go_id)
                id_mapping[[go_id]] <- compound_match[1]
            }
        }
    }
    
    if (length(available_ids) < 2) {
        cat("    Less than 2 terms in dissimilarity matrix for", division, "\n")
        return(NULL)
    }
    
    rep_gos <- rep_gos[rep_gos$go_id %in% available_ids, ]
    matrix_ids <- sapply(rep_gos$go_id, function(x) id_mapping[[x]])
    
    diss_sub <- diss[matrix_ids, matrix_ids]
    rownames(diss_sub) <- rep_gos$go_id
    colnames(diss_sub) <- rep_gos$go_id
    
    # Ensure numeric
    diss_sub <- as.matrix(diss_sub)
    mode(diss_sub) <- "numeric"
    if (any(is.na(diss_sub))) diss_sub[is.na(diss_sub)] <- 1.0
    
    # Cluster
    hc <- hclust(as.dist(diss_sub), method = "average")
    
    # Prepare labels
    rep_gos$label <- paste0(rep_gos$nseqs, "/", rep_gos$nseqs_total, " ", rep_gos$term)
    label_map <- setNames(rep_gos$label, rep_gos$go_id)
    hc$labels <- label_map[hc$labels]
    
    ordered_labels <- hc$labels[hc$order]
    
    # Create plot data
    plot_data <- data.frame(label = ordered_labels, stringsAsFactors = FALSE)
    
    for (i in 1:nrow(plot_data)) {
        match_idx <- which(rep_gos$label == plot_data$label[i])
        if (length(match_idx) > 0) {
            plot_data$pval[i] <- rep_gos$pval[match_idx[1]]
            plot_data$direction[i] <- rep_gos$direction[match_idx[1]]
            plot_data$delta_rank[i] <- rep_gos$delta_rank[match_idx[1]]
        }
    }
    
    # Direction and colors
    if ("direction" %in% names(plot_data) && !all(is.na(plot_data$direction))) {
        plot_data$is_up <- plot_data$direction == "up"
    } else {
        plot_data$is_up <- plot_data$delta_rank > 0
    }
    
    plot_data$color <- NA
    for (i in 1:nrow(plot_data)) {
        p <- plot_data$pval[i]
        up <- plot_data$is_up[i]
        if (p < level2) {
            plot_data$color[i] <- ifelse(up, colors[2], colors[1])
        } else {
            plot_data$color[i] <- ifelse(up, colors[4], colors[3])
        }
    }
    
    # Font styling
    plot_data$font <- 3
    plot_data$font[plot_data$pval < level2] <- 1
    plot_data$font[plot_data$pval < level3] <- 2
    
    plot_data$cex <- 0.7
    plot_data$cex[plot_data$pval < level3] <- 0.85
    
    return(list(
        hc = hc,
        plot_data = plot_data,
        n_terms = nrow(plot_data),
        division = division
    ))
}

# =============================================================================
# MAIN FUNCTION: Plot combined dendrogram (BP, CC, MF stacked)
# =============================================================================

plot_combined_dendrogram <- function(data_dir, prefix, season, organism, output_file) {
    
    cat("\n=== Creating combined plot:", organism, season, "===\n")
    
    # Load data for each division
    divisions <- c("BP", "CC", "MF")
    div_data <- list()
    
    for (div in divisions) {
        if (prefix == "") {
            rep_file <- file.path(data_dir, paste0(season, "_", div, "_representative_GOs.csv"))
            mwu_file <- file.path(data_dir, paste0(season, "_", div, "_MWU_results.csv"))
            dissim_file <- file.path(data_dir, paste0(season, "_", div, "_dissimilarity.csv"))
        } else {
            rep_file <- file.path(data_dir, paste0(prefix, "_", season, "_", div, "_representative_GOs.csv"))
            mwu_file <- file.path(data_dir, paste0(prefix, "_", season, "_", div, "_MWU_results.csv"))
            dissim_file <- file.path(data_dir, paste0(prefix, "_", season, "_", div, "_dissimilarity.csv"))
        }
        
        cat("  Loading", div, "...\n")
        div_data[[div]] <- load_division_data(rep_file, mwu_file, dissim_file, div)
        
        if (!is.null(div_data[[div]])) {
            cat("    Terms:", div_data[[div]]$n_terms, "\n")
        }
    }
    
    # Count valid divisions
    valid_divs <- divisions[!sapply(div_data, is.null)]
    if (length(valid_divs) == 0) {
        cat("  No valid divisions - skipping\n")
        return(NULL)
    }
    
    # Calculate total terms and figure dimensions
    total_terms <- sum(sapply(div_data[valid_divs], function(x) x$n_terms))
    cat("  Total terms:", total_terms, "\n")
    
    # Figure dimensions - condensed spacing
    line_height <- 0.22  # inches per term (condensed)
    fig_height <- max(4, total_terms * line_height + 1.5)  # Add space for legend
    fig_width <- 9
    
    cat("  Figure size:", fig_width, "x", fig_height, "\n")
    
    # Create PDF
    pdf(output_file, width = fig_width, height = fig_height)
    
    # Set up layout: tree | labels
    # Tree panel is narrower now that trees are more compact
    layout(matrix(c(1, 2), 1, 2, byrow = TRUE), widths = c(0.28, 0.72))
    
    # Collect all data for plotting
    all_labels <- c()
    all_colors <- c()
    all_fonts <- c()
    all_cex <- c()
    div_boundaries <- c()  # Track where divisions start
    div_names <- c()
    
    current_pos <- 0
    for (div in valid_divs) {
        d <- div_data[[div]]
        div_boundaries <- c(div_boundaries, current_pos)
        div_names <- c(div_names, div)
        
        all_labels <- c(all_labels, d$plot_data$label)
        all_colors <- c(all_colors, d$plot_data$color)
        all_fonts <- c(all_fonts, d$plot_data$font)
        all_cex <- c(all_cex, d$plot_data$cex)
        
        current_pos <- current_pos + d$n_terms
    }
    div_boundaries <- c(div_boundaries, current_pos)  # End boundary
    
    n_total <- length(all_labels)
    
    # --- Panel 1: Dendrograms (stacked) ---
    par(mar = c(0.5, 0.5, 0.5, 0))
    
    # Create blank plot for dendrograms
    plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, n_total),
         axes = FALSE, xlab = "", ylab = "")
    
    # Draw each dendrogram
    for (i in seq_along(valid_divs)) {
        div <- valid_divs[i]
        d <- div_data[[div]]
        
        y_start <- div_boundaries[i]
        y_end <- div_boundaries[i + 1]
        n_terms <- d$n_terms
        
        # Convert hclust to dendrogram and plot
        dend <- as.dendrogram(d$hc)
        
        # Calculate y positions for this division
        y_positions <- seq(y_start + 0.5, y_end - 0.5, length.out = n_terms)
        
        # Draw the dendrogram using segments
        # Get dendrogram coordinates
        dend_data <- dendrapply(dend, function(n) {
            if (is.leaf(n)) {
                attr(n, "members")
            }
            n
        })
        
        # Simple approach: draw phylo tree horizontally
        hc <- d$hc
        
        # Tree drawing parameters - keep tree compact on left side
        tree_left <- 0.15    # Left margin for tree
        tree_right <- 0.65   # Right edge of tree (leaves connect here)
        tree_width <- tree_right - tree_left
        
        # Normalize merge heights
        max_height <- max(hc$height)
        if (max_height == 0) max_height <- 1
        
        # Calculate y position for each original observation
        obs_y <- setNames(y_positions[hc$order], hc$labels[hc$order])
        
        # Track node y positions
        node_y <- numeric(nrow(hc$merge))
        
        for (j in 1:nrow(hc$merge)) {
            left <- hc$merge[j, 1]
            right <- hc$merge[j, 2]
            
            # Get y positions of children
            if (left < 0) {
                left_y <- obs_y[hc$labels[-left]]
            } else {
                left_y <- node_y[left]
            }
            
            if (right < 0) {
                right_y <- obs_y[hc$labels[-right]]
            } else {
                right_y <- node_y[right]
            }
            
            # Node y is midpoint of children
            node_y[j] <- (left_y + right_y) / 2
            
            # X position based on height (root at left, leaves at right)
            x_pos <- tree_right - (hc$height[j] / max_height) * tree_width
            
            # Get x positions of children
            if (left < 0) {
                left_x <- tree_right
            } else {
                left_x <- tree_right - (hc$height[left] / max_height) * tree_width
            }
            
            if (right < 0) {
                right_x <- tree_right
            } else {
                right_x <- tree_right - (hc$height[right] / max_height) * tree_width
            }
            
            # Draw horizontal lines to children
            segments(x_pos, left_y, left_x, left_y, col = "gray40", lwd = 0.8)
            segments(x_pos, right_y, right_x, right_y, col = "gray40", lwd = 0.8)
            
            # Draw vertical connector
            segments(x_pos, left_y, x_pos, right_y, col = "gray40", lwd = 0.8)
        }
        
        # Add division label at top left, outside the tree area
        text(0.02, y_end - 0.5, division_labels[div], 
             font = 2, cex = 0.7, adj = c(0, 1), col = "gray30")
        
        # Add subtle separator line between divisions (except after last)
        if (i < length(valid_divs)) {
            segments(0, y_end, 1, y_end, col = "gray80", lty = 2, lwd = 0.5)
        }
    }
    
    # --- Panel 2: Labels ---
    par(mar = c(0.5, 0, 0.5, 0.5))
    
    plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(0, n_total),
         axes = FALSE, xlab = "", ylab = "")
    
    # Draw labels
    for (i in 1:n_total) {
        y_pos <- i - 0.5
        label <- sub(" activity$", "", all_labels[i])
        
        text(0.1, y_pos, label,
             font = all_fonts[i],
             cex = all_cex[i],
             col = all_colors[i],
             adj = c(0, 0.5))
    }
    
    # Add legend at bottom
    legend_y <- -0.8
    text(0.1, legend_y, expression(bold("p < 0.01")), cex = 0.65, adj = c(0, 0.5))
    text(2.5, legend_y, "p < 0.05", cex = 0.65, adj = c(0, 0.5), font = 1)
    text(4.5, legend_y, expression(italic("p < 0.1")), cex = 0.65, adj = c(0, 0.5), col = "grey50")
    
    dev.off()
    
    cat("  SUCCESS: Saved to", basename(output_file), "\n")
    return(TRUE)
}

# =============================================================================
# GENERATE ALL FOUR COMBINED PLOTS
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("GENERATING COMBINED REPRESENTATIVE GO DENDROGRAMS\n")
cat(strrep("=", 60), "\n")

# Host Summer
plot_combined_dendrogram(
    data_dir = host_dir,
    prefix = "",
    season = "summer",
    organism = "Host",
    output_file = file.path(output_dir, "Fig4_Host_Summer_Representative_GO_Dendrogram.pdf")
)

# Host Winter
plot_combined_dendrogram(
    data_dir = host_dir,
    prefix = "",
    season = "winter",
    organism = "Host",
    output_file = file.path(output_dir, "Fig4_Host_Winter_Representative_GO_Dendrogram.pdf")
)

# Symbiont Summer
plot_combined_dendrogram(
    data_dir = symbiont_dir,
    prefix = "symbiont",
    season = "summer",
    organism = "Symbiont",
    output_file = file.path(output_dir, "Fig4_Symbiont_Summer_Representative_GO_Dendrogram.pdf")
)

# Symbiont Winter
plot_combined_dendrogram(
    data_dir = symbiont_dir,
    prefix = "symbiont",
    season = "winter",
    organism = "Symbiont",
    output_file = file.path(output_dir, "Fig4_Symbiont_Winter_Representative_GO_Dendrogram.pdf")
)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("COMPLETE\n")
cat(strrep("=", 60), "\n")

output_files <- list.files(output_dir, pattern = "Representative_GO_Dendrogram\\.pdf$", full.names = TRUE)
cat("\nGenerated figures:\n")
for (f in output_files) {
    cat("  ", basename(f), "\n")
}

cat("\nOutput directory:", output_dir, "\n")
