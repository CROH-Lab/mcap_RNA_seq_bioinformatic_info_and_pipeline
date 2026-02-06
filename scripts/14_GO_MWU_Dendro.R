#!/usr/bin/env Rscript
# =============================================================================
# Plot Representative GO Term Dendrograms
# =============================================================================
# Creates clean, readable dendrograms showing ONLY representative GO terms
# Uses identical styling to GO_MWU's gomwuPlot() function
# 
# Author: David Armstrong
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

# Output directory for new figures
output_dir <- "/home/darmstrong4/mc_rework/12_publication_figures/figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Significance thresholds (match GO_MWU defaults)
level1 <- 0.05   # Italic text
level2 <- 0.01  # Normal text  
level3 <- 0.001  # Bold text

# Colors (GO_MWU defaults)
colors <- c("dodgerblue2", "firebrick1", "skyblue2", "lightcoral")
# colors[1] = significant DOWN, colors[2] = significant UP
# colors[3] = less significant DOWN, colors[4] = less significant UP

# =============================================================================
# MAIN FUNCTION: Plot Representative GO Dendrogram
# =============================================================================

plot_representative_dendrogram <- function(
    rep_go_file,      # Path to representative GOs CSV
    mwu_file,         # Path to MWU results CSV
    dissim_file,      # Path to dissimilarity matrix CSV
    output_file,      # Output PDF path
    title = NULL,     # Optional title
    txtsize = 1.0,    # Text size multiplier
    treeHeight = 0.5  # Relative width of dendrogram
) {
    
    cat("\n=== Processing:", basename(rep_go_file), "===\n")
    
    # -------------------------------------------------------------------------
    # 1. Load representative GOs
    # -------------------------------------------------------------------------
    rep_gos <- read.csv(rep_go_file, stringsAsFactors = FALSE)
    cat("  Representative GOs:", nrow(rep_gos), "\n")
    
    if (nrow(rep_gos) < 2) {
        cat("  SKIPPING: Need at least 2 terms for dendrogram\n")
        return(NULL)
    }
    
    # -------------------------------------------------------------------------
    # 2. Load MWU results to get GO IDs and full statistics
    # -------------------------------------------------------------------------
    # MWU results are space/tab-separated, not CSV!
    mwu <- read.table(mwu_file, header = TRUE, stringsAsFactors = FALSE, 
                      quote = "\"", fill = TRUE)
    cat("  MWU results loaded:", nrow(mwu), "terms\n")
    cat("  MWU columns:", paste(names(mwu), collapse = ", "), "\n")
    
    # Match representative terms to MWU results by name
    # Need to handle potential slight differences in naming
    rep_gos$matched_term <- NA
    rep_gos$go_id <- NA
    rep_gos$delta_rank <- NA
    rep_gos$nseqs_total <- NA
    
    for (i in 1:nrow(rep_gos)) {
        term_name <- rep_gos$term[i]
        
        # Try exact match first
        match_idx <- which(mwu$name == term_name)
        
        if (length(match_idx) == 0) {
            # Try case-insensitive match
            match_idx <- which(tolower(mwu$name) == tolower(term_name))
        }
        
        if (length(match_idx) > 0) {
            # Take first match (should be unique)
            match_idx <- match_idx[1]
            rep_gos$matched_term[i] <- mwu$name[match_idx]
            rep_gos$go_id[i] <- mwu$term[match_idx]
            rep_gos$delta_rank[i] <- mwu$delta.rank[match_idx]
            rep_gos$nseqs_total[i] <- mwu$nseqs[match_idx]
        } else {
            cat("  WARNING: Could not match term:", term_name, "\n")
        }
    }
    
    # Remove unmatched terms
    rep_gos <- rep_gos[!is.na(rep_gos$go_id), ]
    cat("  Matched terms:", nrow(rep_gos), "\n")
    
    if (nrow(rep_gos) < 2) {
        cat("  SKIPPING: Less than 2 matched terms\n")
        return(NULL)
    }
    
    # -------------------------------------------------------------------------
    # 3. Load and subset dissimilarity matrix
    # -------------------------------------------------------------------------
    cat("  Loading dissimilarity matrix...\n")
    # Dissimilarity matrix is tab-separated
    # Column names are GO IDs (some are compound: GO:0001;GO:0002)
    # Row structure mirrors columns (symmetric matrix)
    diss <- read.table(dissim_file, sep = "\t", header = TRUE, 
                       check.names = FALSE)
    
    # Set row names from first column or from column names
    # The matrix should be symmetric with GO IDs as both row and column names
    if (!any(grepl("^GO:", rownames(diss)))) {
        # Row names are just numbers - use column names
        rownames(diss) <- colnames(diss)
    }
    
    cat("  Full matrix size:", nrow(diss), "x", ncol(diss), "\n")
    
    # Find matching GO IDs - need to handle compound IDs
    # A compound ID like "GO:0001;GO:0002" contains multiple terms
    available_ids <- c()
    id_mapping <- list()  # Maps our GO ID to the matrix column name
    
    for (go_id in rep_gos$go_id) {
        if (is.na(go_id)) next
        
        # First try direct match
        if (go_id %in% colnames(diss)) {
            available_ids <- c(available_ids, go_id)
            id_mapping[[go_id]] <- go_id
        } else {
            # Try to find as part of compound ID
            compound_match <- grep(paste0("(^|;)", go_id, "(;|$)"), colnames(diss), value = TRUE)
            if (length(compound_match) > 0) {
                available_ids <- c(available_ids, go_id)
                id_mapping[[go_id]] <- compound_match[1]  # Use first match
            }
        }
    }
    
    cat("  GO IDs found in dissimilarity matrix:", length(available_ids), "\n")
    
    if (length(available_ids) < 2) {
        cat("  SKIPPING: Less than 2 terms found in dissimilarity matrix\n")
        return(NULL)
    }
    
    # Subset to representative terms only
    rep_gos <- rep_gos[rep_gos$go_id %in% available_ids, ]
    
    # Get the actual column names to use for subsetting
    matrix_ids <- sapply(rep_gos$go_id, function(x) id_mapping[[x]])
    
    # Subset dissimilarity matrix
    diss_sub <- diss[matrix_ids, matrix_ids]
    
    # Use simple GO IDs as row/column names for clarity
    rownames(diss_sub) <- rep_gos$go_id
    colnames(diss_sub) <- rep_gos$go_id
    
    # -------------------------------------------------------------------------
    # 4. Hierarchical clustering
    # -------------------------------------------------------------------------
    # Ensure matrix is numeric
    diss_sub <- as.matrix(diss_sub)
    mode(diss_sub) <- "numeric"
    
    # Check for any NA values
    if (any(is.na(diss_sub))) {
        cat("  WARNING: NA values in distance matrix, replacing with 1.0\n")
        diss_sub[is.na(diss_sub)] <- 1.0
    }
    
    hc <- hclust(as.dist(diss_sub), method = "average")
    
    # -------------------------------------------------------------------------
    # 5. Prepare data for plotting (GO_MWU style)
    # -------------------------------------------------------------------------
    
    # Create labels in GO_MWU format: "nGood/nTotal termName"
    # nGood = genes with |L2FC| > threshold (we use nseqs from rep file)
    # nTotal = total genes annotated to term
    rep_gos$label <- paste0(
        rep_gos$nseqs, "/", rep_gos$nseqs_total, " ", rep_gos$term
    )
    
    # Map labels to clustering order
    label_map <- setNames(rep_gos$label, rep_gos$go_id)
    hc$labels <- label_map[hc$labels]
    
    # Get ordered labels and data
    ordered_labels <- hc$labels[hc$order]
    
    # Create plotting data frame in dendrogram order
    plot_data <- data.frame(
        label = ordered_labels,
        stringsAsFactors = FALSE
    )
    
    # Match back to rep_gos data
    for (i in 1:nrow(plot_data)) {
        match_idx <- which(rep_gos$label == plot_data$label[i])
        if (length(match_idx) > 0) {
            plot_data$pval[i] <- rep_gos$pval[match_idx[1]]
            plot_data$direction[i] <- rep_gos$direction[match_idx[1]]
            plot_data$delta_rank[i] <- rep_gos$delta_rank[match_idx[1]]
        }
    }
    
    # Determine direction from delta_rank or direction column
    if ("direction" %in% names(plot_data) && !all(is.na(plot_data$direction))) {
        plot_data$is_up <- plot_data$direction == "up"
    } else {
        plot_data$is_up <- plot_data$delta_rank > 0
    }
    
    # Assign colors based on significance and direction
    plot_data$color <- NA
    for (i in 1:nrow(plot_data)) {
        p <- plot_data$pval[i]
        up <- plot_data$is_up[i]
        
        if (p < level2) {
            # Significant (use strong colors)
            plot_data$color[i] <- ifelse(up, colors[2], colors[1])
        } else {
            # Less significant (use pale colors)
            plot_data$color[i] <- ifelse(up, colors[4], colors[3])
        }
    }
    
    # Assign font based on significance
    plot_data$font <- 3  # italic (p < 0.1)
    plot_data$font[plot_data$pval < level2] <- 1  # normal (p < 0.05)
    plot_data$font[plot_data$pval < level3] <- 2  # bold (p < 0.01)
    
    # Text size adjustment
    plot_data$cex <- 0.8 * txtsize
    plot_data$cex[plot_data$pval < level3] <- 1.0 * txtsize
    
    # -------------------------------------------------------------------------
    # 6. Create plot (GO_MWU style)
    # -------------------------------------------------------------------------
    
    # Calculate figure height based on number of terms
    n_terms <- nrow(plot_data)
    fig_height <- max(4, n_terms * 0.3 + 2)
    fig_width <- 10
    
    cat("  Creating figure:", output_file, "\n")
    pdf(output_file, width = fig_width, height = fig_height)
    
    # Layout: dendrogram | labels | legend
    old.par <- par(no.readonly = TRUE)
    layout(matrix(c(1,2,3), 1, 3, byrow = TRUE), 
           widths = c(treeHeight, 3, 0.8), heights = 1)
    
    # --- Panel 1: Dendrogram ---
    par(mar = c(2, 2, 1, 0))
    plot(as.phylo(hc), show.tip.label = FALSE, cex = 0.0001)
    
    # --- Panel 2: Labels ---
    par(mar = c(0, 0, 0.5, 0))
    step <- 100
    top <- step * (2 + n_terms)
    
    plot(c(1:top) ~ c(1:top), type = "n", axes = FALSE, xlab = "", ylab = "")
    
    # Add title if provided
    if (!is.null(title)) {
        text(1, top - step * 0.5, title, font = 2, cex = 1.2 * txtsize, adj = c(0, 0))
    }
    
    # Plot labels from bottom to top (to match dendrogram)
    for (i in n_terms:1) {
        ypos <- top - step * (n_terms - i + 2)
        
        # Clean up label (remove " activity" suffix like GO_MWU does)
        label <- sub(" activity$", "", plot_data$label[i])
        
        text(1, ypos, label,
             font = plot_data$font[i],
             cex = plot_data$cex[i],
             col = plot_data$color[i],
             adj = c(0, 0))
    }
    
    # --- Panel 3: Legend ---
    par(mar = c(3, 1, 1, 0))
    plot(c(1:top) ~ c(1:top), type = "n", axes = FALSE, xlab = "", ylab = "")
    
    text(1, top - step * 2, paste("p <", level3), font = 2, cex = 1.0 * txtsize, adj = c(0, 0))
    text(1, top - step * 3, paste("p <", level2), font = 1, cex = 0.8 * txtsize, adj = c(0, 0))
    text(1, top - step * 4, paste("p <", level1), font = 3, col = "grey50", cex = 0.8 * txtsize, adj = c(0, 0))
    
    # Color legend
    text(1, top - step * 6, "Up-regulated", col = colors[2], cex = 0.8 * txtsize, adj = c(0, 0))
    text(1, top - step * 7, "Down-regulated", col = colors[1], cex = 0.8 * txtsize, adj = c(0, 0))
    
    par(old.par)
    dev.off()
    
    cat("  SUCCESS: Figure saved\n")
    cat("  Terms displayed:", n_terms, "\n")
    
    return(plot_data)
}

# =============================================================================
# PROCESS ALL HOST GO ANALYSES
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("PROCESSING HOST REPRESENTATIVE GO DENDROGRAMS\n")
cat(strrep("=", 60), "\n")

host_analyses <- list(
    list(season = "summer", division = "BP"),
    list(season = "summer", division = "MF"),
    list(season = "summer", division = "CC"),
    list(season = "winter", division = "BP"),
    list(season = "winter", division = "MF"),
    list(season = "winter", division = "CC")
)

for (analysis in host_analyses) {
    season <- analysis$season
    division <- analysis$division
    
    rep_file <- file.path(host_dir, paste0(season, "_", division, "_representative_GOs.csv"))
    mwu_file <- file.path(host_dir, paste0(season, "_", division, "_MWU_results.csv"))
    dissim_file <- file.path(host_dir, paste0(season, "_", division, "_dissimilarity.csv"))
    output_file <- file.path(output_dir, paste0("Fig4_Host_", season, "_", division, "_Representative_Dendrogram.pdf"))
    
    # Check if all input files exist
    if (!file.exists(rep_file)) {
        cat("\nSKIPPING: Missing", rep_file, "\n")
        next
    }
    if (!file.exists(mwu_file)) {
        cat("\nSKIPPING: Missing", mwu_file, "\n")
        next
    }
    if (!file.exists(dissim_file)) {
        cat("\nSKIPPING: Missing", dissim_file, "\n")
        next
    }
    
    title <- paste0("Host ", tools::toTitleCase(season), " - ", division)
    
    tryCatch({
        plot_representative_dendrogram(
            rep_go_file = rep_file,
            mwu_file = mwu_file,
            dissim_file = dissim_file,
            output_file = output_file,
            title = title
        )
    }, error = function(e) {
        cat("\nERROR processing", season, division, ":", conditionMessage(e), "\n")
    })
}

# =============================================================================
# PROCESS ALL SYMBIONT GO ANALYSES
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("PROCESSING SYMBIONT REPRESENTATIVE GO DENDROGRAMS\n")
cat(strrep("=", 60), "\n")

symbiont_analyses <- list(
    list(season = "summer", division = "BP"),
    list(season = "summer", division = "MF"),
    list(season = "summer", division = "CC"),
    list(season = "winter", division = "BP"),
    list(season = "winter", division = "MF"),
    list(season = "winter", division = "CC")
)

for (analysis in symbiont_analyses) {
    season <- analysis$season
    division <- analysis$division
    
    rep_file <- file.path(symbiont_dir, paste0("symbiont_", season, "_", division, "_representative_GOs.csv"))
    mwu_file <- file.path(symbiont_dir, paste0("symbiont_", season, "_", division, "_MWU_results.csv"))
    dissim_file <- file.path(symbiont_dir, paste0("symbiont_", season, "_", division, "_dissimilarity.csv"))
    output_file <- file.path(output_dir, paste0("Fig4_Symbiont_", season, "_", division, "_Representative_Dendrogram.pdf"))
    
    # Check if all input files exist
    if (!file.exists(rep_file)) {
        cat("\nSKIPPING: Missing", rep_file, "\n")
        next
    }
    if (!file.exists(mwu_file)) {
        cat("\nSKIPPING: Missing", mwu_file, "\n")
        next
    }
    if (!file.exists(dissim_file)) {
        cat("\nSKIPPING: Missing", dissim_file, "\n")
        next
    }
    
    title <- paste0("Symbiont ", tools::toTitleCase(season), " - ", division)
    
    tryCatch({
        plot_representative_dendrogram(
            rep_go_file = rep_file,
            mwu_file = mwu_file,
            dissim_file = dissim_file,
            output_file = output_file,
            title = title
        )
    }, error = function(e) {
        cat("\nERROR processing symbiont", season, division, ":", conditionMessage(e), "\n")
    })
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("COMPLETE\n")
cat(strrep("=", 60), "\n")

output_files <- list.files(output_dir, pattern = "Representative_Dendrogram\\.pdf$", full.names = TRUE)
cat("\nGenerated figures:\n")
for (f in output_files) {
    cat("  ", basename(f), "\n")
}

cat("\nOutput directory:", output_dir, "\n")
