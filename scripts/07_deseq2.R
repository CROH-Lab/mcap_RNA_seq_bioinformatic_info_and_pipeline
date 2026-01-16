# ============================================================================
# DESeq2 Analysis - Within-Season Treatment Effect (B vs D)
# ============================================================================
#
# Design: ~ treatment (B vs D) within each season
#   - Summer: B vs D (n=6, 3 vs 3)
#   - Winter: B vs D (n=6, 3 vs 3)
#
# Treatment B = High CO2 (acidified)
# Treatment D = Control
#
# Revision notes:
#   - Added saveRDS() for fitted dds objects
#   - Added VST transformation (blind=FALSE) for downstream visualization
#   - Explicit factor releveling after dds creation
#   - Results extracted at both alpha=0.01 and alpha=0.05
#   - Added PCA visualization
#   - Added session info for reproducibility
#
# ============================================================================

library(DESeq2)
library(tidyverse)
library(qvalue)

# Install cowplot if not available (needed for combined PCA panel)
if (!requireNamespace("cowplot", quietly = TRUE)) {
    install.packages("cowplot", repos = "https://cloud.r-project.org")
}
library(cowplot)

cat("============================================\n")
cat("DESeq2 Analysis - Within-Season B vs D\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n\n")

# --- Configuration ---
BASE_DIR <- "/home/darmstrong4/mc_rework"
COUNTS_DIR <- file.path(BASE_DIR, "06_featurecounts")
OUT_DIR <- file.path(BASE_DIR, "07_deseq2")

dir.create(file.path(OUT_DIR, "results"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "objects"), showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# SECTION 1: Load and Prepare Data
# ============================================================================

cat("=== SECTION 1: Loading Data ===\n\n")

sample_info <- read.delim(file.path(BASE_DIR, "sample_info_with_chemistry.txt"))
rownames(sample_info) <- sample_info$sample_id

# Subset to treatments B and D only (Ocean acidification analysis)
sample_info_BD <- sample_info %>% filter(treatment %in% c("B", "D"))

# Set treatment as factor with D as reference (control)
sample_info_BD$treatment <- factor(sample_info_BD$treatment, levels = c("D", "B"))
sample_info_BD$season <- factor(sample_info_BD$season, levels = c("Summer", "Winter"))

cat("Samples:", nrow(sample_info_BD), "\n")
cat("  Treatment D (control):", sum(sample_info_BD$treatment == "D"), "\n")
cat("  Treatment B (high CO2):", sum(sample_info_BD$treatment == "B"), "\n\n")

# Show design
cat("=== Experimental Design ===\n\n")
print(table(sample_info_BD$treatment, sample_info_BD$season))

cat("\n\n=== Chemistry by Group ===\n\n")
sample_info_BD %>%
    group_by(treatment, season) %>%
    summarise(
        n = n(),
        H_nmol = round(mean(H_nmol), 2),
        pH = round(-log10(mean(H_nmol) * 1e-9), 2),
        DIC = round(mean(DIC), 0),
        .groups = "drop"
    ) %>%
    print()

# Load counts
counts_mcap <- read.delim(file.path(COUNTS_DIR, "Mcap_counts.txt"),
                          row.names = 1, check.names = FALSE)
counts_dtre <- read.delim(file.path(COUNTS_DIR, "Dtre_counts.txt"),
                          row.names = 1, check.names = FALSE)

# Subset to B+D samples
BD_samples <- rownames(sample_info_BD)
counts_mcap_BD <- counts_mcap[, BD_samples]
counts_dtre_BD <- counts_dtre[, BD_samples]

cat("\n\nCount matrices (B+D only):\n")
cat("  M. capitata:", nrow(counts_mcap_BD), "genes ×", ncol(counts_mcap_BD), "samples\n")
cat("  D. trenchii:", nrow(counts_dtre_BD), "genes ×", ncol(counts_dtre_BD), "samples\n\n")

# ============================================================================
# SECTION 2: Within-Season Analysis (B vs D)
# ============================================================================

cat("============================================\n")
cat("WITHIN-SEASON ANALYSIS: ~ treatment (B vs D)\n")
cat("============================================\n\n")

run_season_analysis <- function(counts, sample_info, season_name, organism, out_dir) {

    cat("--- ", organism, " - ", season_name, " ---\n\n", sep = "")

    # Create clean organism name for filenames
    org_short <- gsub("\\. ", "", organism)  # "M. capitata" -> "Mcapitata"
    org_short <- gsub(" ", "_", org_short)

    # Subset to season
    season_samples <- rownames(sample_info)[sample_info$season == season_name]
    counts_season <- counts[, season_samples]
    info_season <- sample_info[season_samples, ]

    cat("Samples:", nrow(info_season), "\n")
    cat("  Treatment D (control):", sum(info_season$treatment == "D"), "\n")
    cat("  Treatment B (high CO2):", sum(info_season$treatment == "B"), "\n\n")

    # Create DESeq2 dataset
    dds <- DESeqDataSetFromMatrix(
        countData = counts_season,
        colData = info_season,
        design = ~ treatment
    )

    # Explicit factor releveling (best practice)
    dds$treatment <- relevel(dds$treatment, ref = "D")

    # Pre-filtering
    keep <- rowSums(counts(dds) >= 10) >= 3
    dds <- dds[keep, ]
    cat("Genes after filtering:", nrow(dds), "\n")

    # Save pre-fit dds object
    saveRDS(dds, file.path(out_dir, "objects",
                           paste0(org_short, "_", season_name, "_dds_prefit.rds")))

    # Run DESeq2
    cat("Running DESeq2...\n")
    dds <- DESeq(dds)
    cat("Complete.\n\n")

    # Save fitted dds object
    saveRDS(dds, file.path(out_dir, "objects",
                           paste0(org_short, "_", season_name, "_dds_fitted.rds")))

    # --- VST Transformation (blind=FALSE for design-aware transformation) ---
    cat("Computing VST transformation (blind=FALSE)...\n")
    vsd <- vst(dds, blind = FALSE)
    vsd_mat <- assay(vsd)

    # Save VST matrix
    saveRDS(vsd_mat, file.path(out_dir, "objects",
                               paste0(org_short, "_", season_name, "_vsd.rds")))
    cat("VST matrix saved.\n\n")

    # --- Extract Results at Two Alpha Thresholds ---

    # Alpha = 0.05 (standard)
    res_05 <- results(dds, name = "treatment_B_vs_D", alpha = 0.05)
    res_05 <- res_05[order(res_05$padj), ]

    # Alpha = 0.01 (conservative)
    res_01 <- results(dds, name = "treatment_B_vs_D", alpha = 0.01)
    res_01 <- res_01[order(res_01$padj), ]

    cat("Treatment effect (B vs D) - Alpha = 0.05:\n")
    cat("  p < 0.01:", sum(res_05$pvalue < 0.01, na.rm = TRUE), "\n")
    cat("  p < 0.05:", sum(res_05$pvalue < 0.05, na.rm = TRUE), "\n")
    cat("  padj < 0.05:", sum(res_05$padj < 0.05, na.rm = TRUE), "\n")
    cat("  padj < 0.10:", sum(res_05$padj < 0.10, na.rm = TRUE), "\n")

    if (sum(res_05$padj < 0.05, na.rm = TRUE) > 0) {
        cat("    Up in B:", sum(res_05$padj < 0.05 & res_05$log2FoldChange > 0, na.rm = TRUE), "\n")
        cat("    Down in B:", sum(res_05$padj < 0.05 & res_05$log2FoldChange < 0, na.rm = TRUE), "\n")
    }

    cat("\nTreatment effect (B vs D) - Alpha = 0.01:\n")
    cat("  padj < 0.01:", sum(res_01$padj < 0.01, na.rm = TRUE), "\n")
    cat("  padj < 0.05:", sum(res_01$padj < 0.05, na.rm = TRUE), "\n")

    n_genes <- sum(!is.na(res_05$pvalue))
    cat("\nGenes tested:", n_genes, "\n")
    cat("Expected by chance: p<0.01 =", round(n_genes * 0.01),
        ", p<0.05 =", round(n_genes * 0.05), "\n\n")

    return(list(
        dds = dds,
        res_05 = res_05,
        res_01 = res_01,
        vsd = vsd_mat,
        info = info_season
    ))
}

# M. capitata
mcap_summer <- run_season_analysis(counts_mcap_BD, sample_info_BD, "Summer", "M. capitata", OUT_DIR)
mcap_winter <- run_season_analysis(counts_mcap_BD, sample_info_BD, "Winter", "M. capitata", OUT_DIR)

# D. trenchii
dtre_summer <- run_season_analysis(counts_dtre_BD, sample_info_BD, "Summer", "D. trenchii", OUT_DIR)
dtre_winter <- run_season_analysis(counts_dtre_BD, sample_info_BD, "Winter", "D. trenchii", OUT_DIR)

# ============================================================================
# SECTION 3: Save Results
# ============================================================================

cat("============================================\n")
cat("Saving Results\n")
cat("============================================\n\n")

# M. capitata - Alpha 0.05
write.csv(as.data.frame(mcap_summer$res_05),
          file.path(OUT_DIR, "results", "Mcap_Summer_BvsD_alpha05.csv"))
write.csv(as.data.frame(mcap_winter$res_05),
          file.path(OUT_DIR, "results", "Mcap_Winter_BvsD_alpha05.csv"))

# M. capitata - Alpha 0.01
write.csv(as.data.frame(mcap_summer$res_01),
          file.path(OUT_DIR, "results", "Mcap_Summer_BvsD_alpha01.csv"))
write.csv(as.data.frame(mcap_winter$res_01),
          file.path(OUT_DIR, "results", "Mcap_Winter_BvsD_alpha01.csv"))

# D. trenchii - Alpha 0.05
write.csv(as.data.frame(dtre_summer$res_05),
          file.path(OUT_DIR, "results", "Dtre_Summer_BvsD_alpha05.csv"))
write.csv(as.data.frame(dtre_winter$res_05),
          file.path(OUT_DIR, "results", "Dtre_Winter_BvsD_alpha05.csv"))

# D. trenchii - Alpha 0.01
write.csv(as.data.frame(dtre_summer$res_01),
          file.path(OUT_DIR, "results", "Dtre_Summer_BvsD_alpha01.csv"))
write.csv(as.data.frame(dtre_winter$res_01),
          file.path(OUT_DIR, "results", "Dtre_Winter_BvsD_alpha01.csv"))

cat("Saved CSV results to: 07_deseq2/results/\n")
cat("  - *_alpha05.csv: Results with alpha=0.05 (standard)\n")
cat("  - *_alpha01.csv: Results with alpha=0.01 (conservative)\n\n")

# ============================================================================
# SECTION 4: Pi0 Estimation
# ============================================================================

cat("============================================\n")
cat("Pi0 Estimation\n")
cat("============================================\n\n")

estimate_pi0 <- function(res, label) {
    res_df <- as.data.frame(res) %>% filter(!is.na(pvalue))
    if (nrow(res_df) > 100) {
        qobj <- tryCatch(qvalue(res_df$pvalue), error = function(e) NULL)
        if (!is.null(qobj)) {
            cat(label, ":\n", sep = "")
            cat("  Pi0:", round(qobj$pi0, 3), "\n")
            cat("  Estimated true DEGs:", round((1 - qobj$pi0) * nrow(res_df)), "\n\n")
            return(qobj$pi0)
        }
    }
    return(NA)
}

pi0_mcap_summer <- estimate_pi0(mcap_summer$res_05, "M. capitata Summer")
pi0_mcap_winter <- estimate_pi0(mcap_winter$res_05, "M. capitata Winter")
pi0_dtre_summer <- estimate_pi0(dtre_summer$res_05, "D. trenchii Summer")
pi0_dtre_winter <- estimate_pi0(dtre_winter$res_05, "D. trenchii Winter")

# ============================================================================
# SECTION 5: Summary Table
# ============================================================================

cat("============================================\n")
cat("SUMMARY\n")
cat("============================================\n\n")

summary_df <- data.frame(
    Organism = c("M. capitata", "M. capitata", "D. trenchii", "D. trenchii"),
    Season = c("Summer", "Winter", "Summer", "Winter"),
    n = c(6, 6, 6, 6),
    Genes_Tested = c(
        sum(!is.na(mcap_summer$res_05$pvalue)),
        sum(!is.na(mcap_winter$res_05$pvalue)),
        sum(!is.na(dtre_summer$res_05$pvalue)),
        sum(!is.na(dtre_winter$res_05$pvalue))
    ),
    p_01 = c(
        sum(mcap_summer$res_05$pvalue < 0.01, na.rm = TRUE),
        sum(mcap_winter$res_05$pvalue < 0.01, na.rm = TRUE),
        sum(dtre_summer$res_05$pvalue < 0.01, na.rm = TRUE),
        sum(dtre_winter$res_05$pvalue < 0.01, na.rm = TRUE)
    ),
    p_05 = c(
        sum(mcap_summer$res_05$pvalue < 0.05, na.rm = TRUE),
        sum(mcap_winter$res_05$pvalue < 0.05, na.rm = TRUE),
        sum(dtre_summer$res_05$pvalue < 0.05, na.rm = TRUE),
        sum(dtre_winter$res_05$pvalue < 0.05, na.rm = TRUE)
    ),
    padj_05 = c(
        sum(mcap_summer$res_05$padj < 0.05, na.rm = TRUE),
        sum(mcap_winter$res_05$padj < 0.05, na.rm = TRUE),
        sum(dtre_summer$res_05$padj < 0.05, na.rm = TRUE),
        sum(dtre_winter$res_05$padj < 0.05, na.rm = TRUE)
    ),
    padj_10 = c(
        sum(mcap_summer$res_05$padj < 0.10, na.rm = TRUE),
        sum(mcap_winter$res_05$padj < 0.10, na.rm = TRUE),
        sum(dtre_summer$res_05$padj < 0.10, na.rm = TRUE),
        sum(dtre_winter$res_05$padj < 0.10, na.rm = TRUE)
    ),
    Pi0 = c(pi0_mcap_summer, pi0_mcap_winter, pi0_dtre_summer, pi0_dtre_winter)
)

print(summary_df, row.names = FALSE)

write.csv(summary_df,
          file.path(OUT_DIR, "results", "DEG_summary.csv"),
          row.names = FALSE)

# ============================================================================
# SECTION 6: PCA Visualization
# ============================================================================

cat("\n============================================\n")
cat("Creating PCA Plots\n")
cat("============================================\n\n")

make_pca <- function(vsd_mat, sample_info, title, filename, out_dir) {

    # Run PCA
    pca <- prcomp(t(vsd_mat))
    var_explained <- round(100 * summary(pca)$importance[2, 1:2], 1)

    # Create plot data
    pca_df <- as.data.frame(pca$x[, 1:2])
    pca_df$sample <- rownames(pca_df)
    pca_df$treatment <- sample_info[pca_df$sample, "treatment"]

    # Plot
    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = treatment)) +
        geom_point(size = 4, alpha = 0.8) +
        geom_text(aes(label = sample), vjust = -1, size = 2.5, show.legend = FALSE) +
        scale_color_manual(
            values = c("D" = "#4575b4", "B" = "#d73027"),
            labels = c("D" = "Control (D)", "B" = "High CO2 (B)")
        ) +
        labs(
            title = title,
            x = paste0("PC1 (", var_explained[1], "%)"),
            y = paste0("PC2 (", var_explained[2], "%)"),
            color = "Treatment"
        ) +
        theme_bw(base_size = 12) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "bottom",
            panel.grid.minor = element_blank()
        )

    ggsave(file.path(out_dir, "figures", filename), p, width = 7, height = 6)
    cat("Saved:", filename, "\n")

    # Return PCA results for potential downstream use
    return(list(pca = pca, var_explained = var_explained))
}

# M. capitata PCAs
pca_mcap_summer <- make_pca(mcap_summer$vsd, mcap_summer$info,
                            "M. capitata Summer - PCA", "pca_Mcap_Summer.pdf", OUT_DIR)
pca_mcap_winter <- make_pca(mcap_winter$vsd, mcap_winter$info,
                            "M. capitata Winter - PCA", "pca_Mcap_Winter.pdf", OUT_DIR)

# D. trenchii PCAs
pca_dtre_summer <- make_pca(dtre_summer$vsd, dtre_summer$info,
                            "D. trenchii Summer - PCA", "pca_Dtre_Summer.pdf", OUT_DIR)
pca_dtre_winter <- make_pca(dtre_winter$vsd, dtre_winter$info,
                            "D. trenchii Winter - PCA", "pca_Dtre_Winter.pdf", OUT_DIR)

# --- Combined PCA panel ---
cat("\nCreating combined PCA panel...\n")

make_pca_panel <- function(results_list, out_dir) {

    plot_list <- list()
    names_list <- c("M. capitata Summer", "M. capitata Winter",
                    "D. trenchii Summer", "D. trenchii Winter")

    for (i in seq_along(results_list)) {
        res <- results_list[[i]]
        vsd_mat <- res$vsd
        sample_info <- res$info

        pca <- prcomp(t(vsd_mat))
        var_explained <- round(100 * summary(pca)$importance[2, 1:2], 1)

        pca_df <- as.data.frame(pca$x[, 1:2])
        pca_df$treatment <- sample_info[rownames(pca_df), "treatment"]

        p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = treatment)) +
            geom_point(size = 3) +
            scale_color_manual(
                values = c("D" = "#4575b4", "B" = "#d73027"),
                labels = c("D" = "Control", "B" = "High CO2")
            ) +
            labs(
                title = names_list[i],
                x = paste0("PC1 (", var_explained[1], "%)"),
                y = paste0("PC2 (", var_explained[2], "%)")
            ) +
            theme_bw(base_size = 10) +
            theme(
                plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
                legend.position = "none"
            )

        plot_list[[i]] <- p
    }

    # Combine plots
    combined <- plot_grid(
        plot_list[[1]], plot_list[[2]],
        plot_list[[3]], plot_list[[4]],
        nrow = 2, ncol = 2,
        labels = c("A", "B", "C", "D"),
        label_size = 12
    )

    # Add shared legend
    legend_plot <- ggplot(data.frame(x = 1:2, y = 1:2, treatment = c("D", "B")),
                          aes(x, y, color = treatment)) +
        geom_point(size = 3) +
        scale_color_manual(
            values = c("D" = "#4575b4", "B" = "#d73027"),
            labels = c("D" = "Control (D)", "B" = "High CO2 (B)")
        ) +
        theme_void() +
        theme(legend.position = "bottom") +
        guides(color = guide_legend(title = "Treatment"))

    legend <- get_legend(legend_plot)

    final_plot <- plot_grid(combined, legend, ncol = 1, rel_heights = c(1, 0.08))

    ggsave(file.path(out_dir, "figures", "pca_combined_panel.pdf"),
           final_plot, width = 10, height = 9)

    cat("Saved: pca_combined_panel.pdf\n")
}

make_pca_panel(list(mcap_summer, mcap_winter, dtre_summer, dtre_winter), OUT_DIR)

# ============================================================================
# SECTION 7: P-value Histograms
# ============================================================================

cat("\n=== Creating P-value Histograms ===\n\n")

pdf(file.path(OUT_DIR, "figures", "pvalue_histograms.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

hist(mcap_summer$res_05$pvalue, breaks = 50, main = "M. capitata Summer - B vs D",
     xlab = "p-value", col = "steelblue", border = "white")
abline(h = sum(!is.na(mcap_summer$res_05$pvalue))/50, col = "red", lty = 2)

hist(mcap_winter$res_05$pvalue, breaks = 50, main = "M. capitata Winter - B vs D",
     xlab = "p-value", col = "steelblue", border = "white")
abline(h = sum(!is.na(mcap_winter$res_05$pvalue))/50, col = "red", lty = 2)

hist(dtre_summer$res_05$pvalue, breaks = 50, main = "D. trenchii Summer - B vs D",
     xlab = "p-value", col = "darkgreen", border = "white")
abline(h = sum(!is.na(dtre_summer$res_05$pvalue))/50, col = "red", lty = 2)

hist(dtre_winter$res_05$pvalue, breaks = 50, main = "D. trenchii Winter - B vs D",
     xlab = "p-value", col = "darkgreen", border = "white")
abline(h = sum(!is.na(dtre_winter$res_05$pvalue))/50, col = "red", lty = 2)

dev.off()
cat("Saved: figures/pvalue_histograms.pdf\n")

# ============================================================================
# SECTION 8: Volcano Plots
# ============================================================================

cat("\n=== Creating Volcano Plots ===\n\n")

make_volcano <- function(res, title, filename, out_dir) {
    df <- as.data.frame(res) %>%
        rownames_to_column("gene") %>%
        filter(!is.na(padj)) %>%
        mutate(
            direction = case_when(
                padj < 0.05 & log2FoldChange > 0 ~ "Up",
                padj < 0.05 & log2FoldChange < 0 ~ "Down",
                TRUE ~ "NS"
            )
        )

    n_up <- sum(df$direction == "Up")
    n_down <- sum(df$direction == "Down")

    p <- ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue))) +
        geom_point(aes(color = direction), alpha = 0.5, size = 1) +
        scale_color_manual(values = c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "grey70")) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
        labs(
            title = title,
            subtitle = paste0("DEGs (padj < 0.05): ", n_up, " up, ", n_down, " down"),
            x = "log2 Fold Change (B vs D)",
            y = "-log10(p-value)"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5)
        )

    ggsave(file.path(out_dir, "figures", filename), p, width = 8, height = 6)
}

make_volcano(mcap_summer$res_05, "M. capitata Summer - B vs D", "volcano_Mcap_Summer.pdf", OUT_DIR)
make_volcano(mcap_winter$res_05, "M. capitata Winter - B vs D", "volcano_Mcap_Winter.pdf", OUT_DIR)
make_volcano(dtre_summer$res_05, "D. trenchii Summer - B vs D", "volcano_Dtre_Summer.pdf", OUT_DIR)
make_volcano(dtre_winter$res_05, "D. trenchii Winter - B vs D", "volcano_Dtre_Winter.pdf", OUT_DIR)

cat("Saved volcano plots to: figures/\n")

# ============================================================================
# SECTION 9: Sample Distance Heatmap
# ============================================================================

cat("\n=== Creating Sample Distance Heatmaps ===\n\n")

make_sample_heatmap <- function(vsd_mat, sample_info, title, filename, out_dir) {

    # Calculate sample distances
    sampleDists <- dist(t(vsd_mat))
    sampleDistMatrix <- as.matrix(sampleDists)

    # Set row/column names
    rownames(sampleDistMatrix) <- paste(sample_info$treatment,
                                         rownames(sample_info), sep = "_")
    colnames(sampleDistMatrix) <- NULL

    # Color palette
    colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)

    # Annotation
    annotation_row <- data.frame(
        Treatment = sample_info$treatment,
        row.names = rownames(sampleDistMatrix)
    )

    # Plot
    pdf(file.path(out_dir, "figures", filename), width = 8, height = 7)
    pheatmap::pheatmap(
        sampleDistMatrix,
        clustering_distance_rows = sampleDists,
        clustering_distance_cols = sampleDists,
        col = colors,
        annotation_row = annotation_row,
        annotation_colors = list(Treatment = c("D" = "#4575b4", "B" = "#d73027")),
        main = title,
        fontsize = 10
    )
    dev.off()

    cat("Saved:", filename, "\n")
}

make_sample_heatmap(mcap_summer$vsd, mcap_summer$info,
                    "M. capitata Summer - Sample Distances",
                    "heatmap_Mcap_Summer.pdf", OUT_DIR)
make_sample_heatmap(mcap_winter$vsd, mcap_winter$info,
                    "M. capitata Winter - Sample Distances",
                    "heatmap_Mcap_Winter.pdf", OUT_DIR)
make_sample_heatmap(dtre_summer$vsd, dtre_summer$info,
                    "D. trenchii Summer - Sample Distances",
                    "heatmap_Dtre_Summer.pdf", OUT_DIR)
make_sample_heatmap(dtre_winter$vsd, dtre_winter$info,
                    "D. trenchii Winter - Sample Distances",
                    "heatmap_Dtre_Winter.pdf", OUT_DIR)

# ============================================================================
# SECTION 10: Top DEGs
# ============================================================================

cat("\n=== Top DEGs ===\n\n")

print_top_degs <- function(res, label, n = 10) {
    n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
    if (n_sig > 0) {
        cat(label, " (", n_sig, " DEGs at padj < 0.05):\n\n", sep = "")
        top <- as.data.frame(res) %>%
            rownames_to_column("gene") %>%
            filter(padj < 0.05) %>%
            arrange(padj) %>%
            head(n)
        print(top[, c("gene", "log2FoldChange", "pvalue", "padj")], row.names = FALSE)
        cat("\n")
    } else {
        cat(label, ": No DEGs at padj < 0.05\n\n", sep = "")
    }
}

print_top_degs(mcap_summer$res_05, "M. capitata Summer")
print_top_degs(mcap_winter$res_05, "M. capitata Winter")
print_top_degs(dtre_summer$res_05, "D. trenchii Summer")
print_top_degs(dtre_winter$res_05, "D. trenchii Winter")

cat("Positive LFC = higher expression in B (high CO2/acidified)\n")
cat("Negative LFC = lower expression in B (high CO2/acidified)\n")

# ============================================================================
# SECTION 11: Save R Objects
# ============================================================================

cat("\n=== Saving R Objects ===\n\n")

# Save complete analysis objects
save(mcap_summer, mcap_winter,
     file = file.path(OUT_DIR, "Mcap_DESeq2_complete.RData"))

save(dtre_summer, dtre_winter,
     file = file.path(OUT_DIR, "Dtre_DESeq2_complete.RData"))

cat("Saved complete analysis objects:\n")
cat("  - Mcap_DESeq2_complete.RData\n")
cat("  - Dtre_DESeq2_complete.RData\n\n")

cat("Individual objects saved to objects/:\n")
cat("  - *_dds_prefit.rds  (DESeq2 object before fitting)\n")
cat("  - *_dds_fitted.rds  (DESeq2 object after fitting)\n")
cat("  - *_vsd.rds         (VST-transformed counts)\n")

# ============================================================================
# SECTION 12: Session Info
# ============================================================================

cat("\n============================================\n")
cat("Session Info\n")
cat("============================================\n\n")

session_info <- capture.output(sessionInfo())
writeLines(session_info, file.path(OUT_DIR, "session_info.txt"))
cat("Saved: session_info.txt\n")

# Print to console as well
cat("\n")
sessionInfo()

# ============================================================================
# SECTION 13: Output Summary
# ============================================================================

cat("\n============================================\n")
cat("Analysis Complete!\n")
cat("============================================\n\n")

cat("Output directory:", OUT_DIR, "\n\n")

cat("Results files:\n")
cat("  results/*_alpha05.csv  - DESeq2 results (alpha=0.05)\n")
cat("  results/*_alpha01.csv  - DESeq2 results (alpha=0.01)\n")
cat("  results/DEG_summary.csv\n\n")

cat("Figure files:\n")
cat("  figures/pca_*.pdf           - PCA plots\n")
cat("  figures/pca_combined_panel.pdf\n")
cat("  figures/volcano_*.pdf       - Volcano plots\n")
cat("  figures/heatmap_*.pdf       - Sample distance heatmaps\n")
cat("  figures/pvalue_histograms.pdf\n\n")

cat("R objects:\n")
cat("  objects/*_dds_prefit.rds    - Pre-fit DESeq2 objects\n")
cat("  objects/*_dds_fitted.rds    - Fitted DESeq2 objects\n")
cat("  objects/*_vsd.rds           - VST matrices\n")
cat("  *_DESeq2_complete.RData     - Complete analysis lists\n\n")

cat("Reproducibility:\n")
cat("  session_info.txt\n")
