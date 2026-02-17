#!/usr/bin/env Rscript
# =============================================================================
# Plot Combined Host + Symbiont Representative GO Dendrograms
# =============================================================================
# Creates side-by-side dendrograms for host and symbiont with:
# - Column headers: "M. capitata (host)" and "D. trenchii (symbiont)"
# - Panel letters (A, B)
# - BP, CC, MF stacked vertically in each column
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
host_dir <- "/home/darmstrong4/mc_rework/10_host_GO_MWU/output"
symbiont_dir <- "/home/darmstrong4/mc_rework/11_symbiont_GO_MWU/output"

# Output directory
output_dir <- "/home/darmstrong4/mc_rework/12_publication_figures/figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Significance thresholds
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
# HELPER: Load division data
# =============================================================================

load_division_data <- function(rep_file, mwu_file, dissim_file, division) {
    
    if (!file.exists(rep_file) || !file.exists(mwu_file) || !file.exists(dissim_file)) {
        return(NULL)
    }
    
    rep_gos <- read.csv(rep_file, stringsAsFactors = FALSE)
    if (nrow(rep_gos) < 1) return(NULL)
    
    mwu <- read.table(mwu_file, header = TRUE, stringsAsFactors = FALSE,
                      quote = "\"", fill = TRUE)
    
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
    if (nrow(rep_gos) < 2) return(NULL)
    
    diss <- read.table(dissim_file, sep = "\t", header = TRUE, check.names = FALSE)
    if (!any(grepl("^GO:", rownames(diss)))) {
        rownames(diss) <- colnames(diss)
    }
    
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
    
    if (length(available_ids) < 2) return(NULL)
    
    rep_gos <- rep_gos[rep_gos$go_id %in% available_ids, ]
    matrix_ids <- sapply(rep_gos$go_id, function(x) id_mapping[[x]])
    
    diss_sub <- diss[matrix_ids, matrix_ids]
    rownames(diss_sub) <- rep_gos$go_id
    colnames(diss_sub) <- rep_gos$go_id
    
    diss_sub <- as.matrix(diss_sub)
    mode(diss_sub) <- "numeric"
    if (any(is.na(diss_sub))) diss_sub[is.na(diss_sub)] <- 1.0
    
    hc <- hclust(as.dist(diss_sub), method = "average")
    
    rep_gos$label <- paste0(rep_gos$nseqs, "/", rep_gos$nseqs_total, " ", rep_gos$term)
    label_map <- setNames(rep_gos$label, rep_gos$go_id)
    hc$labels <- label_map[hc$labels]
    
    ordered_labels <- hc$labels[hc$order]
    
    plot_data <- data.frame(label = ordered_labels, stringsAsFactors = FALSE)
    
    for (i in 1:nrow(plot_data)) {
        match_idx <- which(rep_gos$label == plot_data$label[i])
        if (length(match_idx) > 0) {
            plot_data$pval[i] <- rep_gos$pval[match_idx[1]]
            plot_data$direction[i] <- rep_gos$direction[match_idx[1]]
            plot_data$delta_rank[i] <- rep_gos$delta_rank[match_idx[1]]
        }
    }
    
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
    
    plot_data$font <- 3
    plot_data$font[plot_data$pval < level2] <- 1
    plot_data$font[plot_data$pval < level3] <- 2
    
    plot_data$cex <- 0.72
    plot_data$cex[plot_data$pval < level3] <- 0.85
    
    return(list(
        hc = hc,
        plot_data = plot_data,
        n_terms = nrow(plot_data),
        division = division
    ))
}

# =============================================================================
# HELPER: Load all divisions for one organism/season
# =============================================================================

load_organism_data <- function(data_dir, prefix, season) {
    
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
        
        div_data[[div]] <- load_division_data(rep_file, mwu_file, dissim_file, div)
    }
    
    return(div_data)
}

# =============================================================================
# HELPER: Draw dendrogram for one division
# =============================================================================

draw_division_dendrogram <- function(hc, plot_data, y_offset, tree_x_start, tree_x_end, label_x) {
    
    n_terms <- nrow(plot_data)
    
    y_positions <- y_offset + seq(0.5, n_terms - 0.5, by = 1)
    obs_y <- setNames(y_positions, hc$labels[hc$order])
    
    max_height <- max(hc$height)
    if (max_height == 0) max_height <- 1
    
    tree_width <- tree_x_end - tree_x_start
    
    node_y <- numeric(nrow(hc$merge))
    node_x <- numeric(nrow(hc$merge))
    
    for (j in 1:nrow(hc$merge)) {
        left <- hc$merge[j, 1]
        right <- hc$merge[j, 2]
        
        if (left < 0) {
            left_y <- obs_y[hc$labels[-left]]
            left_x <- tree_x_end
        } else {
            left_y <- node_y[left]
            left_x <- node_x[left]
        }
        
        if (right < 0) {
            right_y <- obs_y[hc$labels[-right]]
            right_x <- tree_x_end
        } else {
            right_y <- node_y[right]
            right_x <- node_x[right]
        }
        
        node_y[j] <- (left_y + right_y) / 2
        node_x[j] <- tree_x_end - (hc$height[j] / max_height) * tree_width
        
        segments(node_x[j], left_y, left_x, left_y, col = "gray30", lwd = 0.6)
        segments(node_x[j], right_y, right_x, right_y, col = "gray30", lwd = 0.6)
        segments(node_x[j], left_y, node_x[j], right_y, col = "gray30", lwd = 0.6)
    }
    
    for (i in 1:n_terms) {
        label <- sub(" activity$", "", plot_data$label[i])
        
        text(label_x, y_positions[i], label,
             font = plot_data$font[i],
             cex = plot_data$cex[i],
             col = plot_data$color[i],
             adj = c(0, 0.5))
    }
    
    return(c(y_offset, y_offset + n_terms))
}

# =============================================================================
# HELPER: Draw one column (organism)
# =============================================================================

draw_organism_column <- function(div_data, x_offset, x_width, total_height, gap_size) {
    
    divisions <- c("BP", "CC", "MF")
    valid_divs <- divisions[!sapply(div_data[divisions], is.null)]
    
    if (length(valid_divs) == 0) return(NULL)
    
    # Tree coordinates within this column
    tree_x_start <- x_offset + 0.01 * x_width
    tree_x_end <- x_offset + 0.12 * x_width
    label_x <- x_offset + 0.13 * x_width
    
    # Draw divisions from top to bottom
    y_cursor <- total_height - 0.5
    
    for (i in seq_along(valid_divs)) {
        div <- valid_divs[i]
        d <- div_data[[div]]
        n_terms <- d$n_terms
        
        y_offset <- y_cursor - n_terms
        
        # Division label
        text(tree_x_start, y_cursor + 0.15,
             division_labels[div],
             font = 2, cex = 0.62, adj = c(0, 0), col = "gray25")
        
        # Draw dendrogram
        draw_division_dendrogram(
            hc = d$hc,
            plot_data = d$plot_data,
            y_offset = y_offset,
            tree_x_start = tree_x_start,
            tree_x_end = tree_x_end,
            label_x = label_x
        )
        
        y_cursor <- y_offset - gap_size
    }
    
    return(TRUE)
}

# =============================================================================
# MAIN: Plot combined host + symbiont dendrogram
# =============================================================================

plot_combined_season <- function(season, output_file) {
    
    cat("\n=== Creating combined plot for", season, "===\n")
    
    # Load data for both organisms
    cat("  Loading Host data...\n")
    host_data <- load_organism_data(host_dir, "", season)
    
    cat("  Loading Symbiont data...\n")
    symbiont_data <- load_organism_data(symbiont_dir, "symbiont", season)
    
    # Count terms for each organism
    divisions <- c("BP", "CC", "MF")
    
    host_valid <- divisions[!sapply(host_data[divisions], is.null)]
    symbiont_valid <- divisions[!sapply(symbiont_data[divisions], is.null)]
    
    host_terms <- sum(sapply(host_data[host_valid], function(x) x$n_terms))
    symbiont_terms <- sum(sapply(symbiont_data[symbiont_valid], function(x) x$n_terms))
    
    cat("  Host terms:", host_terms, "\n")
    cat("  Symbiont terms:", symbiont_terms, "\n")
    
    # Use max terms to determine height
    max_terms <- max(host_terms, symbiont_terms)
    n_gaps <- max(length(host_valid), length(symbiont_valid)) - 1
    gap_size <- 1.2
    total_height <- max_terms + n_gaps * gap_size + 2  # +2 for title space
    
    # Figure dimensions
    line_height <- 0.18
    fig_height <- max(4, total_height * line_height)
    fig_width <- 13  # Slightly narrower
    
    cat("  Figure size:", fig_width, "x", fig_height, "\n")
    
    # Create PDF
    pdf(output_file, width = fig_width, height = fig_height)
    
    par(mar = c(0.5, 0.5, 0.5, 0.5))
    
    # Set up plot
    plot(0, 0, type = "n",
         xlim = c(0, 1),
         ylim = c(0, total_height),
         axes = FALSE, xlab = "", ylab = "")
    
    # Column widths - closer together
    col_width <- 0.46
    col1_start <- 0.01
    col2_start <- 0.47  # Moved closer to host column
    
    # --- Column titles ---
    title_y <- total_height - 0.3
    
    # Panel A - Host
    text(col1_start, title_y, "A", font = 2, cex = 1.3, adj = c(0, 0.5))
    text(col1_start + 0.03, title_y, expression(italic("M. capitata")~"(host)"), 
         font = 1, cex = 1.0, adj = c(0, 0.5))
    
    # Panel B - Symbiont
    text(col2_start, title_y, "B", font = 2, cex = 1.3, adj = c(0, 0.5))
    text(col2_start + 0.03, title_y, expression(italic("D. trenchii")~"(symbiont)"), 
         font = 1, cex = 1.0, adj = c(0, 0.5))
    
    # Adjust total_height for drawing (below titles)
    draw_height <- total_height - 1.5
    
    # --- Draw Host column ---
    cat("  Drawing Host column...\n")
    draw_organism_column(
        div_data = host_data,
        x_offset = col1_start,
        x_width = col_width,
        total_height = draw_height,
        gap_size = gap_size
    )
    
    # --- Draw Symbiont column ---
    cat("  Drawing Symbiont column...\n")
    draw_organism_column(
        div_data = symbiont_data,
        x_offset = col2_start,
        x_width = col_width,
        total_height = draw_height,
        gap_size = gap_size
    )
    
    # --- Legend at top right (tighter to content) ---
    legend_x <- 0.88
    legend_y <- total_height - 0.5
    line_spacing <- 0.4
    
    text(legend_x, legend_y, expression(bold("p < 0.01")), cex = 0.6, adj = c(0, 0.5))
    text(legend_x, legend_y - line_spacing, "p < 0.05", cex = 0.6, adj = c(0, 0.5), font = 1)
    text(legend_x, legend_y - 2*line_spacing, expression(italic("p < 0.1")), cex = 0.6, adj = c(0, 0.5), col = "grey50")
    
    dev.off()
    
    cat("  SUCCESS:", basename(output_file), "\n")
    return(TRUE)
}

# =============================================================================
# GENERATE COMBINED PLOTS
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("GENERATING COMBINED HOST + SYMBIONT DENDROGRAMS\n")
cat(strrep("=", 60), "\n")

# Summer
plot_combined_season(
    season = "summer",
    output_file = file.path(output_dir, "Fig4_Summer_Representative_GO_Dendrogram.pdf")
)

# Winter
plot_combined_season(
    season = "winter",
    output_file = file.path(output_dir, "Fig4_Winter_Representative_GO_Dendrogram.pdf")
)

cat("\n", strrep("=", 60), "\n")
cat("COMPLETE\n")
cat(strrep("=", 60), "\n")

output_files <- list.files(output_dir, pattern = "Fig4_.*_Representative_GO_Dendrogram\\.pdf$", full.names = TRUE)
cat("\nGenerated figures:\n")
for (f in output_files) cat("  ", basename(f), "\n")
