#!/usr/bin/env Rscript
# =============================================================================
# Simplified Host-Symbiont Cellular Interaction Schematic
# =============================================================================
# Two-panel figure (Summer | Winter) showing shared GO term regulation
# mapped to cellular compartments.
#
# Required packages: ggplot2, ggforce, patchwork
# Install if needed:
#   conda install -c conda-forge r-ggplot2 r-ggforce r-patchwork
#   # or within R:
#   install.packages(c("ggforce", "patchwork"))
#
# Author: Claude (for David Armstrong)
# Date: 2025-02-10
# =============================================================================

library(ggplot2)
library(ggforce)
library(patchwork)

# =============================================================================
# CONFIGURATION
# =============================================================================

output_dir <- "/home/darmstrong4/mc_rework/12_publication_figures/figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Interaction colors (matching chord diagram)
col_syn_up     <- "#907AD6"
col_syn_down   <- "#7DCFB6"
col_antag_hup  <- "#F79256"
col_antag_hdown <- "#6BBF59"

# Organelle colors
col_host_fill   <- "#E8F0F8"
col_host_border <- "#4A7BA7"
col_sym_fill    <- "#E8F5E2"
col_sym_border  <- "#4A8A3A"
col_symbiosome  <- "#F5EDD0"
col_symbio_bdr  <- "#C4A84A"
col_nuc_fill    <- "#F0EAF5"
col_nuc_border  <- "#8A6AAA"
col_mito_fill   <- "#FCEAEA"
col_mito_border <- "#C46A6A"
col_er_fill     <- "#EAF0FC"
col_er_border   <- "#6A8AC4"
col_ribo        <- "#6A8AC4"
col_membrane    <- "#4A7BA7"

# =============================================================================
# HELPER: Build one panel
# =============================================================================

build_panel <- function(season = "Summer") {

    is_summer <- (season == "Summer")

    # --- Base canvas ---
    p <- ggplot() +
        coord_fixed(xlim = c(-0.05, 1.05), ylim = c(-0.05, 1.05), clip = "off") +
        theme_void() +
        theme(
            plot.margin = margin(5, 5, 5, 5),
            plot.title = element_text(hjust = 0.5, size = 13, face = "bold",
                                      margin = margin(b = 2)),
            plot.subtitle = element_text(hjust = 0.5, size = 8, color = "grey50",
                                         margin = margin(b = 6))
        ) +
        labs(
            title = season,
            subtitle = if (is_summer) {
                "18 shared terms \u2022 17 synergistic down"
            } else {
                "63 shared terms \u2022 60 antagonistic"
            }
        )

    # =========================================================================
    # ORGANELLES
    # =========================================================================

    # --- Host cell (rounded rectangle via geom_shape) ---
    host_x <- c(0.03, 0.97, 0.97, 0.03)
    host_y <- c(0.03, 0.03, 0.93, 0.93)
    host_df <- data.frame(x = host_x, y = host_y)

    p <- p + geom_shape(data = host_df, aes(x = x, y = y),
                        fill = col_host_fill, color = col_host_border,
                        linewidth = 1.2, radius = unit(0.3, "cm"))

    # Host label
    p <- p + annotate("text", x = 0.13, y = 0.89, label = "italic('M. capitata')",
                      parse = TRUE, size = 3.8, fontface = "bold", color = col_host_border) +
             annotate("text", x = 0.11, y = 0.855, label = "(Host cell)",
                      size = 2.8, color = col_host_border)

    # --- Plasma membrane label ---
    p <- p + annotate("text", x = 0.88, y = 0.895, label = "Plasma membrane",
                      size = 2.3, color = col_membrane, fontface = "italic")

    # --- Nucleus ---
    p <- p + geom_ellipse(aes(x0 = 0.25, y0 = 0.62, a = 0.14, b = 0.11, angle = 0),
                          fill = col_nuc_fill, color = col_nuc_border, linewidth = 0.8) +
             annotate("text", x = 0.25, y = 0.62, label = "Nucleus",
                      size = 2.8, fontface = "bold", color = col_nuc_border)

    # --- Rough ER (simplified wavy lines near nucleus) ---
    er_y_vals <- seq(0.46, 0.54, length.out = 50)
    er_x_vals1 <- 0.44 + 0.02 * sin(er_y_vals * 60)
    er_x_vals2 <- 0.48 + 0.02 * sin(er_y_vals * 60 + 1)
    er_df1 <- data.frame(x = er_x_vals1, y = er_y_vals)
    er_df2 <- data.frame(x = er_x_vals2, y = er_y_vals)

    p <- p + geom_path(data = er_df1, aes(x = x, y = y),
                       color = col_er_border, linewidth = 0.7) +
             geom_path(data = er_df2, aes(x = x, y = y),
                       color = col_er_border, linewidth = 0.7) +
             annotate("text", x = 0.46, y = 0.44, label = "ER",
                      size = 2.3, color = col_er_border, fontface = "bold")

    # --- Ribosomes (small dots) ---
    ribo_df <- data.frame(
        x = c(0.33, 0.36, 0.30, 0.38, 0.34, 0.31),
        y = c(0.42, 0.40, 0.38, 0.43, 0.36, 0.44)
    )
    p <- p + geom_point(data = ribo_df, aes(x = x, y = y),
                        size = 1.2, color = col_ribo, alpha = 0.6) +
             annotate("text", x = 0.34, y = 0.33, label = "Ribosomes",
                      size = 2.3, color = col_ribo, fontface = "bold")

    # --- Mitochondria ---
    mito_positions <- if (is_summer) {
        data.frame(x0 = c(0.18, 0.28), y0 = c(0.20, 0.15),
                   a = c(0.05, 0.04), b = c(0.025, 0.02))
    } else {
        data.frame(x0 = c(0.15, 0.25, 0.35), y0 = c(0.20, 0.14, 0.20),
                   a = c(0.05, 0.04, 0.045), b = c(0.025, 0.02, 0.022))
    }

    for (i in 1:nrow(mito_positions)) {
        p <- p + geom_ellipse(
            aes(x0 = .x0, y0 = .y0, a = .a, b = .b, angle = 0),
            data = data.frame(.x0 = mito_positions$x0[i],
                              .y0 = mito_positions$y0[i],
                              .a = mito_positions$a[i],
                              .b = mito_positions$b[i]),
            fill = col_mito_fill, color = col_mito_border, linewidth = 0.6
        )
    }
    p <- p + annotate("text", x = 0.25, y = 0.07, label = "Mitochondria",
                      size = 2.3, color = col_mito_border, fontface = "bold")

    # --- Symbiosome membrane ---
    p <- p + geom_ellipse(aes(x0 = 0.70, y0 = 0.45, a = 0.20, b = 0.22, angle = 0),
                          fill = col_symbiosome, color = col_symbio_bdr,
                          linewidth = 0.8, linetype = "dashed") +
             annotate("text", x = 0.70, y = 0.695, label = "Symbiosome",
                      size = 2.2, color = col_symbio_bdr, fontface = "italic")

    # --- Symbiont cell ---
    p <- p + geom_ellipse(aes(x0 = 0.70, y0 = 0.45, a = 0.15, b = 0.17, angle = 0),
                          fill = col_sym_fill, color = col_sym_border, linewidth = 1.0) +
             annotate("text", x = 0.70, y = 0.47, label = "italic('D. trenchii')",
                      parse = TRUE, size = 3.0, fontface = "bold", color = col_sym_border) +
             annotate("text", x = 0.70, y = 0.43, label = "(Symbiont)",
                      size = 2.5, color = col_sym_border)

    # =========================================================================
    # REGULATION ANNOTATIONS
    # =========================================================================

    # Helper: draw a regulation box with arrows
    # Returns a list of ggplot layers
    add_reg_box <- function(p, x, y, label, n_terms,
                            host_dir, sym_dir, box_color,
                            label_size = 2.6, box_w = 0.22, box_h = 0.055) {

        # Background box
        p <- p + annotate("rect",
                          xmin = x - box_w/2, xmax = x + box_w/2,
                          ymin = y - box_h/2, ymax = y + box_h/2,
                          fill = box_color, alpha = 0.15,
                          color = box_color, linewidth = 0.5)

        # Label text
        p <- p + annotate("text", x = x, y = y + 0.005,
                          label = label,
                          size = label_size, fontface = "bold",
                          color = "grey20")

        # N terms
        p <- p + annotate("text", x = x, y = y - box_h/2 - 0.018,
                          label = paste0("(", n_terms, " terms)"),
                          size = 2.0, color = "grey50")

        # Host arrow
        h_arrow <- if (host_dir == "up") "\u2191" else "\u2193"
        p <- p + annotate("text", x = x - box_w/2 - 0.025, y = y + 0.012,
                          label = "H", size = 2.0, color = "grey40") +
                 annotate("text", x = x - box_w/2 - 0.025, y = y - 0.012,
                          label = h_arrow, size = 3.0, fontface = "bold",
                          color = box_color)

        # Symbiont arrow
        s_arrow <- if (sym_dir == "up") "\u2191" else "\u2193"
        p <- p + annotate("text", x = x + box_w/2 + 0.025, y = y + 0.012,
                          label = "S", size = 2.0, color = "grey40") +
                 annotate("text", x = x + box_w/2 + 0.025, y = y - 0.012,
                          label = s_arrow, size = 3.0, fontface = "bold",
                          color = box_color)

        return(p)
    }

    # =========================================================================
    # SEASON-SPECIFIC REGULATION BOXES
    # =========================================================================

    if (is_summer) {

        # Translation & Ribosome suppression (main signal)
        p <- add_reg_box(p, 0.34, 0.77, "Translation & Ribosomes", "15",
                         "down", "down", col_syn_down, box_w = 0.26)

        # Protein-to-membrane targeting
        p <- add_reg_box(p, 0.75, 0.80, "Protein \u2192 Membrane", "2",
                         "down", "down", col_syn_down, box_w = 0.22)

        # Microtubule motors (only syn up)
        p <- add_reg_box(p, 0.50, 0.10, "Microtubule Motors", "1",
                         "up", "up", col_syn_up, box_w = 0.22)

    } else {

        # Translation - HOST UP, SYM DOWN
        p <- add_reg_box(p, 0.34, 0.80, "Translation & Ribosomes", "11",
                         "up", "down", col_antag_hup, box_w = 0.26)

        # Protein folding
        p <- add_reg_box(p, 0.75, 0.83, "Protein Folding", "3",
                         "up", "down", col_antag_hup, box_w = 0.20)

        # Mitochondrial energy
        p <- add_reg_box(p, 0.25, 0.27, "Mito. Energy", "8",
                         "up", "down", col_antag_hup, box_w = 0.18)

        # Nuclear transport / splicing
        p <- add_reg_box(p, 0.11, y = 0.50, "Nuclear\nTransport", "4",
                         "up", "down", col_antag_hup,
                         box_w = 0.13, box_h = 0.065)

        # Ion transport - HOST DOWN, SYM UP
        p <- add_reg_box(p, 0.50, 0.06, "Ion Channels", "4",
                         "down", "up", col_antag_hdown, box_w = 0.18)

        # Membrane & epithelial remodeling
        p <- add_reg_box(p, 0.85, 0.55, "Membrane\nRemodeling", "16",
                         "down", "up", col_antag_hdown,
                         box_w = 0.15, box_h = 0.065)

        # Signaling / neural-like
        p <- add_reg_box(p, 0.85, 0.35, "Signaling\nPathways", "6",
                         "down", "up", col_antag_hdown,
                         box_w = 0.15, box_h = 0.065)
    }

    return(p)
}

# =============================================================================
# BUILD LEGEND
# =============================================================================

build_legend <- function() {

    legend_df <- data.frame(
        x = c(1, 2, 3, 4),
        y = rep(0, 4),
        label = c("Synergistic Up", "Synergistic Down",
                   "Antagonistic\n(Host Up / Sym Down)",
                   "Antagonistic\n(Host Down / Sym Up)"),
        color = c(col_syn_up, col_syn_down, col_antag_hup, col_antag_hdown),
        stringsAsFactors = FALSE
    )

    p <- ggplot() +
        coord_fixed(xlim = c(0.3, 4.7), ylim = c(-0.3, 0.3), clip = "off") +
        theme_void()

    for (i in 1:4) {
        p <- p + annotate("rect",
                          xmin = legend_df$x[i] - 0.08,
                          xmax = legend_df$x[i] + 0.08,
                          ymin = -0.08, ymax = 0.08,
                          fill = legend_df$color[i], color = legend_df$color[i],
                          linewidth = 0.3) +
                 annotate("text", x = legend_df$x[i], y = -0.20,
                          label = legend_df$label[i],
                          size = 2.3, color = "grey30", lineheight = 0.85)
    }

    return(p)
}

# =============================================================================
# COMBINE AND SAVE
# =============================================================================

cat("\nBuilding schematic...\n")

p_summer <- build_panel("Summer")
p_winter <- build_panel("Winter")
p_legend <- build_legend()

# Combine: two panels side by side, legend below
combined <- (p_summer | p_winter) / p_legend +
    plot_layout(heights = c(10, 1)) +
    plot_annotation(
        title = "Host\u2013Symbiont Cellular Interactions Under Ocean Acidification",
        subtitle = "Shared significant GO terms between M. capitata and D. trenchii (p.adj < 0.1)",
        theme = theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey50",
                                         margin = margin(b = 8))
        )
    )

# Save
output_file <- file.path(output_dir, "Fig_Host_Symbiont_Schematic.pdf")
ggsave(output_file, combined, width = 12, height = 7, device = "pdf")
cat("Saved:", output_file, "\n")

# Also save PNG for quick preview
output_png <- file.path(output_dir, "Fig_Host_Symbiont_Schematic.png")
ggsave(output_png, combined, width = 12, height = 7, dpi = 300)
cat("Saved:", output_png, "\n")

cat("\nDone! Check output in:", output_dir, "\n")
