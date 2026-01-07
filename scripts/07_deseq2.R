# ============================================================================
# M. capitata Holobiont RNA-seq Pipeline - Step 7: DESeq2 Analysis
#
# Differential expression analysis using DESeq2.
# Analyzes host and symbiont gene expression separately.
#
# Input:  Count matrices from featureCounts (06_counts/)
# Output: DE results, figures, normalized counts (07_deseq2/)
#
# Usage:
#   Interactive: Open in RStudio and run sections
#   Command line: Rscript scripts/07_deseq2.R
# ============================================================================

# --- Load Libraries ---
suppressPackageStartupMessages({
    library(DESeq2)
    library(tidyverse)
    library(pheatmap)
    library(RColorBrewer)
    library(ggrepel)
})

cat("============================================\n")
cat("DESeq2 Differential Expression Analysis\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n\n")

# --- Configuration ---
BASE_DIR <- "/home/darmstrong4/mc_rework"
COUNTS_DIR <- file.path(BASE_DIR, "06_counts")
OUT_DIR <- file.path(BASE_DIR, "07_deseq2")

# Create output directories
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "results"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figures"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "normalized_counts"), showWarnings = FALSE)

# ============================================================================
# SECTION 1: Load Data
# ============================================================================
cat("=== Loading Data ===\n")

# Load sample metadata
sample_info <- read.delim(file.path(BASE_DIR, "sample_info.txt"), 
                          stringsAsFactors = FALSE)
rownames(sample_info) <- sample_info$sample

cat("Sample info loaded:", nrow(sample_info), "samples\n")
print(table(sample_info$treatment, sample_info$season))
cat("\n")

# Load count matrices
counts_all <- read.delim(file.path(COUNTS_DIR, "gene_counts.txt"), 
                         row.names = 1)
counts_mcap <- read.delim(file.path(COUNTS_DIR, "gene_counts_Mcap.txt"), 
                          row.names = 1)
counts_cgor <- read.delim(file.path(COUNTS_DIR, "gene_counts_Cgor.txt"), 
                          row.names = 1)
counts_dtre <- read.delim(file.path(COUNTS_DIR, "gene_counts_Dtre.txt"), 
                          row.names = 1)

cat("Count matrices loaded:\n")
cat("  All genes:", nrow(counts_all), "genes x", ncol(counts_all), "samples\n")
cat("  M. capitata:", nrow(counts_mcap), "genes\n")
cat("  C. goreaui:", nrow(counts_cgor), "genes\n")
cat("  D. trenchii:", nrow(counts_dtre), "genes\n\n")

# Ensure sample order matches
sample_info <- sample_info[colnames(counts_all), ]

# Convert factors
sample_info$treatment <- factor(sample_info$treatment, 
                                 levels = c("D", "C", "B", "A"))  # D = control
sample_info$season <- factor(sample_info$season, 
                              levels = c("Summer", "Winter"))
sample_info$genotype <- factor(sample_info$genotype)

# ============================================================================
# SECTION 2: DESeq2 Analysis - Host (M. capitata)
# ============================================================================
cat("=== DESeq2 Analysis: Host (M. capitata) ===\n")

# Create DESeq2 dataset
dds_mcap <- DESeqDataSetFromMatrix(
    countData = counts_mcap,
    colData = sample_info,
    design = ~ genotype + season + treatment
)

# Filter low-count genes (at least 10 reads in at least 3 samples)
keep <- rowSums(counts(dds_mcap) >= 10) >= 3
dds_mcap <- dds_mcap[keep, ]
cat("Genes after filtering:", nrow(dds_mcap), "\n")

# Run DESeq2
dds_mcap <- DESeq(dds_mcap)
cat("DESeq2 analysis complete\n\n")

# --- Results: Treatment effect (A vs D, B vs D, C vs D) ---
cat("Extracting treatment contrasts (vs control D)...\n")

# Treatment A vs D (control)
res_mcap_AvD <- results(dds_mcap, contrast = c("treatment", "A", "D"), 
                        alpha = 0.05)
res_mcap_AvD <- res_mcap_AvD[order(res_mcap_AvD$padj), ]
cat("  A vs D: ", sum(res_mcap_AvD$padj < 0.05, na.rm = TRUE), "DEGs (padj < 0.05)\n")

# Treatment B vs D (control)
res_mcap_BvD <- results(dds_mcap, contrast = c("treatment", "B", "D"), 
                        alpha = 0.05)
res_mcap_BvD <- res_mcap_BvD[order(res_mcap_BvD$padj), ]
cat("  B vs D: ", sum(res_mcap_BvD$padj < 0.05, na.rm = TRUE), "DEGs (padj < 0.05)\n")

# Treatment C vs D (control)
res_mcap_CvD <- results(dds_mcap, contrast = c("treatment", "C", "D"), 
                        alpha = 0.05)
res_mcap_CvD <- res_mcap_CvD[order(res_mcap_CvD$padj), ]
cat("  C vs D: ", sum(res_mcap_CvD$padj < 0.05, na.rm = TRUE), "DEGs (padj < 0.05)\n")

# Season effect
res_mcap_season <- results(dds_mcap, contrast = c("season", "Winter", "Summer"), 
                           alpha = 0.05)
res_mcap_season <- res_mcap_season[order(res_mcap_season$padj), ]
cat("  Winter vs Summer: ", sum(res_mcap_season$padj < 0.05, na.rm = TRUE), "DEGs\n\n")

# Save results
write.csv(as.data.frame(res_mcap_AvD), 
          file.path(OUT_DIR, "results", "Mcap_AvD_results.csv"))
write.csv(as.data.frame(res_mcap_BvD), 
          file.path(OUT_DIR, "results", "Mcap_BvD_results.csv"))
write.csv(as.data.frame(res_mcap_CvD), 
          file.path(OUT_DIR, "results", "Mcap_CvD_results.csv"))
write.csv(as.data.frame(res_mcap_season), 
          file.path(OUT_DIR, "results", "Mcap_WinterVsSummer_results.csv"))

# Get normalized counts
norm_counts_mcap <- counts(dds_mcap, normalized = TRUE)
write.csv(norm_counts_mcap, 
          file.path(OUT_DIR, "normalized_counts", "Mcap_normalized_counts.csv"))

# ============================================================================
# SECTION 3: DESeq2 Analysis - Symbionts (if sufficient reads)
# ============================================================================
cat("=== DESeq2 Analysis: Symbionts ===\n")

# --- C. goreaui (Clade C) ---
cat("\nC. goreaui (Clade C):\n")

# Check which samples have enough Clade C reads
cgor_totals <- colSums(counts_cgor)
cgor_samples <- names(cgor_totals[cgor_totals > 100000])  # Min 100k reads
cat("  Samples with >100k Clade C reads:", length(cgor_samples), "\n")

if (length(cgor_samples) >= 6) {
    dds_cgor <- DESeqDataSetFromMatrix(
        countData = counts_cgor[, cgor_samples],
        colData = sample_info[cgor_samples, ],
        design = ~ season + treatment
    )
    keep <- rowSums(counts(dds_cgor) >= 10) >= 3
    dds_cgor <- dds_cgor[keep, ]
    cat("  Genes after filtering:", nrow(dds_cgor), "\n")
    
    dds_cgor <- DESeq(dds_cgor)
    norm_counts_cgor <- counts(dds_cgor, normalized = TRUE)
    write.csv(norm_counts_cgor, 
              file.path(OUT_DIR, "normalized_counts", "Cgor_normalized_counts.csv"))
    cat("  Analysis complete\n")
} else {
    cat("  Insufficient samples for Clade C analysis\n")
}

# --- D. trenchii (Clade D) ---
cat("\nD. trenchii (Clade D):\n")

# Check which samples have enough Clade D reads
dtre_totals <- colSums(counts_dtre)
dtre_samples <- names(dtre_totals[dtre_totals > 100000])
cat("  Samples with >100k Clade D reads:", length(dtre_samples), "\n")

if (length(dtre_samples) >= 6) {
    dds_dtre <- DESeqDataSetFromMatrix(
        countData = counts_dtre[, dtre_samples],
        colData = sample_info[dtre_samples, ],
        design = ~ season + treatment
    )
    keep <- rowSums(counts(dds_dtre) >= 10) >= 3
    dds_dtre <- dds_dtre[keep, ]
    cat("  Genes after filtering:", nrow(dds_dtre), "\n")
    
    dds_dtre <- DESeq(dds_dtre)
    norm_counts_dtre <- counts(dds_dtre, normalized = TRUE)
    write.csv(norm_counts_dtre, 
              file.path(OUT_DIR, "normalized_counts", "Dtre_normalized_counts.csv"))
    cat("  Analysis complete\n")
} else {
    cat("  Insufficient samples for Clade D analysis\n")
}

cat("\n")

# ============================================================================
# SECTION 4: Visualization - Host
# ============================================================================
cat("=== Creating Visualizations ===\n")

# --- PCA Plot ---
cat("Creating PCA plot...\n")
vsd_mcap <- vst(dds_mcap, blind = FALSE)

pca_data <- plotPCA(vsd_mcap, intgroup = c("treatment", "season"), 
                    returnData = TRUE)
pca_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, 
                                  color = treatment, shape = season)) +
    geom_point(size = 4, alpha = 0.8) +
    scale_color_brewer(palette = "Set1") +
    labs(
        title = "PCA: M. capitata Host Gene Expression",
        x = paste0("PC1 (", pca_var[1], "% variance)"),
        y = paste0("PC2 (", pca_var[2], "% variance)")
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
    )

ggsave(file.path(OUT_DIR, "figures", "PCA_Mcap_host.pdf"), 
       pca_plot, width = 8, height = 6)
ggsave(file.path(OUT_DIR, "figures", "PCA_Mcap_host.png"), 
       pca_plot, width = 8, height = 6, dpi = 300)

# --- Sample Distance Heatmap ---
cat("Creating sample distance heatmap...\n")
sample_dists <- dist(t(assay(vsd_mcap)))
sample_dist_matrix <- as.matrix(sample_dists)

rownames(sample_dist_matrix) <- paste(vsd_mcap$treatment, vsd_mcap$season, 
                                       vsd_mcap$genotype, sep = "_")
colnames(sample_dist_matrix) <- NULL

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pdf(file.path(OUT_DIR, "figures", "sample_distance_heatmap.pdf"), 
    width = 10, height = 8)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colors,
         main = "Sample Distance Matrix - M. capitata Host")
dev.off()

# --- Volcano Plot (Treatment A vs D) ---
cat("Creating volcano plot (A vs D)...\n")

volcano_data <- as.data.frame(res_mcap_AvD) %>%
    rownames_to_column("gene") %>%
    mutate(
        significance = case_when(
            padj < 0.05 & log2FoldChange > 1 ~ "Up",
            padj < 0.05 & log2FoldChange < -1 ~ "Down",
            TRUE ~ "NS"
        )
    )

volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), 
                                          color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray60")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
    labs(
        title = "Volcano Plot: Treatment A vs Control D",
        subtitle = "M. capitata Host",
        x = "log2 Fold Change",
        y = "-log10(adjusted p-value)"
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
    )

ggsave(file.path(OUT_DIR, "figures", "volcano_Mcap_AvD.pdf"), 
       volcano_plot, width = 8, height = 6)
ggsave(file.path(OUT_DIR, "figures", "volcano_Mcap_AvD.png"), 
       volcano_plot, width = 8, height = 6, dpi = 300)

# --- Top DEGs Heatmap ---
cat("Creating top DEGs heatmap...\n")

# Get top 50 DEGs by adjusted p-value
top_genes <- head(rownames(res_mcap_AvD[order(res_mcap_AvD$padj), ]), 50)
top_genes <- top_genes[!is.na(top_genes)]

if (length(top_genes) > 5) {
    mat <- assay(vsd_mcap)[top_genes, ]
    mat <- t(scale(t(mat)))  # Z-score normalize
    
    annotation_col <- data.frame(
        Treatment = sample_info$treatment,
        Season = sample_info$season,
        row.names = colnames(mat)
    )
    
    pdf(file.path(OUT_DIR, "figures", "heatmap_top50_DEGs_AvD.pdf"), 
        width = 12, height = 10)
    pheatmap(mat,
             annotation_col = annotation_col,
             show_rownames = FALSE,
             cluster_cols = TRUE,
             cluster_rows = TRUE,
             main = "Top 50 DEGs: Treatment A vs Control D")
    dev.off()
}

# ============================================================================
# SECTION 5: Summary Statistics
# ============================================================================
cat("\n=== Summary Statistics ===\n")

# Create summary table
summary_df <- data.frame(
    Comparison = c("A vs D", "B vs D", "C vs D", "Winter vs Summer"),
    Total_Tested = c(
        sum(!is.na(res_mcap_AvD$padj)),
        sum(!is.na(res_mcap_BvD$padj)),
        sum(!is.na(res_mcap_CvD$padj)),
        sum(!is.na(res_mcap_season$padj))
    ),
    DEGs_padj05 = c(
        sum(res_mcap_AvD$padj < 0.05, na.rm = TRUE),
        sum(res_mcap_BvD$padj < 0.05, na.rm = TRUE),
        sum(res_mcap_CvD$padj < 0.05, na.rm = TRUE),
        sum(res_mcap_season$padj < 0.05, na.rm = TRUE)
    ),
    Upregulated = c(
        sum(res_mcap_AvD$padj < 0.05 & res_mcap_AvD$log2FoldChange > 0, na.rm = TRUE),
        sum(res_mcap_BvD$padj < 0.05 & res_mcap_BvD$log2FoldChange > 0, na.rm = TRUE),
        sum(res_mcap_CvD$padj < 0.05 & res_mcap_CvD$log2FoldChange > 0, na.rm = TRUE),
        sum(res_mcap_season$padj < 0.05 & res_mcap_season$log2FoldChange > 0, na.rm = TRUE)
    ),
    Downregulated = c(
        sum(res_mcap_AvD$padj < 0.05 & res_mcap_AvD$log2FoldChange < 0, na.rm = TRUE),
        sum(res_mcap_BvD$padj < 0.05 & res_mcap_BvD$log2FoldChange < 0, na.rm = TRUE),
        sum(res_mcap_CvD$padj < 0.05 & res_mcap_CvD$log2FoldChange < 0, na.rm = TRUE),
        sum(res_mcap_season$padj < 0.05 & res_mcap_season$log2FoldChange < 0, na.rm = TRUE)
    )
)

print(summary_df)

write.csv(summary_df, file.path(OUT_DIR, "results", "DEG_summary.csv"), 
          row.names = FALSE)

# ============================================================================
# SECTION 6: Save R Objects
# ============================================================================
cat("\n=== Saving R Objects ===\n")

save(dds_mcap, res_mcap_AvD, res_mcap_BvD, res_mcap_CvD, res_mcap_season,
     file = file.path(OUT_DIR, "Mcap_DESeq2_objects.RData"))
cat("Saved: Mcap_DESeq2_objects.RData\n")

# ============================================================================
# Final Summary
# ============================================================================
cat("\n============================================\n")
cat("DESeq2 Analysis Complete!\n")
cat("============================================\n\n")

cat("Output directory:", OUT_DIR, "\n\n")

cat("Results files:\n")
list.files(file.path(OUT_DIR, "results"), full.names = FALSE)

cat("\nFigures:\n")
list.files(file.path(OUT_DIR, "figures"), full.names = FALSE)

cat("\nNormalized counts:\n")
list.files(file.path(OUT_DIR, "normalized_counts"), full.names = FALSE)

cat("\n\nTo load results in R:\n")
cat('  load("', file.path(OUT_DIR, "Mcap_DESeq2_objects.RData"), '")\n', sep = "")

cat("\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
