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
# ============================================================================

library(DESeq2)
library(tidyverse)
library(qvalue)

cat("============================================\n")
cat("DESeq2 Analysis - Within-Season B vs D\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n\n")

# --- Configuration ---
BASE_DIR <- "/home/darmstrong4/mc_rework"
COUNTS_DIR <- file.path(BASE_DIR, "06_featurecounts")
OUT_DIR <- file.path(BASE_DIR, "07_deseq2")

dir.create(file.path(OUT_DIR, "results"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figures"), showWarnings = FALSE)

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

run_season_analysis <- function(counts, sample_info, season_name, organism) {
    
    cat("--- ", organism, " - ", season_name, " ---\n\n", sep = "")
    
    # Subset to season
    season_samples <- rownames(sample_info)[sample_info$season == season_name]
    counts_season <- counts[, season_samples]
    info_season <- sample_info[season_samples, ]
    
    cat("Samples:", nrow(info_season), "\n")
    cat("  Treatment D (control):", sum(info_season$treatment == "D"), "\n")
    cat("  Treatment B (high CO2):", sum(info_season$treatment == "B"), "\n\n")
    
    # Design: just treatment
    dds <- DESeqDataSetFromMatrix(
        countData = counts_season,
        colData = info_season,
        design = ~ treatment
    )
    
    keep <- rowSums(counts(dds) >= 10) >= 3
    dds <- dds[keep, ]
    cat("Genes after filtering:", nrow(dds), "\n")
    
    cat("Running DESeq2...\n")
    dds <- DESeq(dds)
    cat("Complete.\n\n")
    
    # Results
    res <- results(dds, name = "treatment_B_vs_D", alpha = 0.05)
    res <- res[order(res$padj), ]
    
    cat("Treatment effect (B vs D):\n")
    cat("  p < 0.01:", sum(res$pvalue < 0.01, na.rm = TRUE), "\n")
    cat("  p < 0.05:", sum(res$pvalue < 0.05, na.rm = TRUE), "\n")
    cat("  padj < 0.05:", sum(res$padj < 0.05, na.rm = TRUE), "\n")
    cat("  padj < 0.10:", sum(res$padj < 0.10, na.rm = TRUE), "\n")
    
    if (sum(res$padj < 0.05, na.rm = TRUE) > 0) {
        cat("    Up in B:", sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE), "\n")
        cat("    Down in B:", sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE), "\n")
    }
    
    n_genes <- sum(!is.na(res$pvalue))
    cat("  Expected by chance: p<0.01 =", round(n_genes * 0.01),
        ", p<0.05 =", round(n_genes * 0.05), "\n\n")
    
    return(list(dds = dds, res = res, info = info_season))
}

# M. capitata
mcap_summer <- run_season_analysis(counts_mcap_BD, sample_info_BD, "Summer", "M. capitata")
mcap_winter <- run_season_analysis(counts_mcap_BD, sample_info_BD, "Winter", "M. capitata")

# D. trenchii
dtre_summer <- run_season_analysis(counts_dtre_BD, sample_info_BD, "Summer", "D. trenchii")
dtre_winter <- run_season_analysis(counts_dtre_BD, sample_info_BD, "Winter", "D. trenchii")

# ============================================================================
# SECTION 3: Save Results
# ============================================================================

cat("============================================\n")
cat("Saving Results\n")
cat("============================================\n\n")

# M. capitata
write.csv(as.data.frame(mcap_summer$res),
          file.path(OUT_DIR, "results", "Mcap_Summer_BvsD.csv"))
write.csv(as.data.frame(mcap_winter$res),
          file.path(OUT_DIR, "results", "Mcap_Winter_BvsD.csv"))

# D. trenchii
write.csv(as.data.frame(dtre_summer$res),
          file.path(OUT_DIR, "results", "Dtre_Summer_BvsD.csv"))
write.csv(as.data.frame(dtre_winter$res),
          file.path(OUT_DIR, "results", "Dtre_Winter_BvsD.csv"))

cat("Saved CSV results to: 07_deseq2/results/\n\n")

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

pi0_mcap_summer <- estimate_pi0(mcap_summer$res, "M. capitata Summer")
pi0_mcap_winter <- estimate_pi0(mcap_winter$res, "M. capitata Winter")
pi0_dtre_summer <- estimate_pi0(dtre_summer$res, "D. trenchii Summer")
pi0_dtre_winter <- estimate_pi0(dtre_winter$res, "D. trenchii Winter")

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
    p_01 = c(
        sum(mcap_summer$res$pvalue < 0.01, na.rm = TRUE),
        sum(mcap_winter$res$pvalue < 0.01, na.rm = TRUE),
        sum(dtre_summer$res$pvalue < 0.01, na.rm = TRUE),
        sum(dtre_winter$res$pvalue < 0.01, na.rm = TRUE)
    ),
    p_05 = c(
        sum(mcap_summer$res$pvalue < 0.05, na.rm = TRUE),
        sum(mcap_winter$res$pvalue < 0.05, na.rm = TRUE),
        sum(dtre_summer$res$pvalue < 0.05, na.rm = TRUE),
        sum(dtre_winter$res$pvalue < 0.05, na.rm = TRUE)
    ),
    padj_05 = c(
        sum(mcap_summer$res$padj < 0.05, na.rm = TRUE),
        sum(mcap_winter$res$padj < 0.05, na.rm = TRUE),
        sum(dtre_summer$res$padj < 0.05, na.rm = TRUE),
        sum(dtre_winter$res$padj < 0.05, na.rm = TRUE)
    ),
    padj_10 = c(
        sum(mcap_summer$res$padj < 0.10, na.rm = TRUE),
        sum(mcap_winter$res$padj < 0.10, na.rm = TRUE),
        sum(dtre_summer$res$padj < 0.10, na.rm = TRUE),
        sum(dtre_winter$res$padj < 0.10, na.rm = TRUE)
    )
)

print(summary_df, row.names = FALSE)

write.csv(summary_df,
          file.path(OUT_DIR, "results", "DEG_summary.csv"),
          row.names = FALSE)

# ============================================================================
# SECTION 6: P-value Histograms
# ============================================================================

cat("\n=== Creating P-value Histograms ===\n\n")

pdf(file.path(OUT_DIR, "figures", "pvalue_histograms.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

hist(mcap_summer$res$pvalue, breaks = 50, main = "M. capitata Summer - B vs D",
     xlab = "p-value", col = "steelblue", border = "white")
abline(h = sum(!is.na(mcap_summer$res$pvalue))/50, col = "red", lty = 2)

hist(mcap_winter$res$pvalue, breaks = 50, main = "M. capitata Winter - B vs D",
     xlab = "p-value", col = "steelblue", border = "white")
abline(h = sum(!is.na(mcap_winter$res$pvalue))/50, col = "red", lty = 2)

hist(dtre_summer$res$pvalue, breaks = 50, main = "D. trenchii Summer - B vs D",
     xlab = "p-value", col = "darkgreen", border = "white")
abline(h = sum(!is.na(dtre_summer$res$pvalue))/50, col = "red", lty = 2)

hist(dtre_winter$res$pvalue, breaks = 50, main = "D. trenchii Winter - B vs D",
     xlab = "p-value", col = "darkgreen", border = "white")
abline(h = sum(!is.na(dtre_winter$res$pvalue))/50, col = "red", lty = 2)

dev.off()
cat("Saved: figures/pvalue_histograms.pdf\n")

# ============================================================================
# SECTION 7: Volcano Plots
# ============================================================================

cat("\n=== Creating Volcano Plots ===\n\n")

make_volcano <- function(res, title, filename) {
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
    
    ggsave(file.path(OUT_DIR, "figures", filename), p, width = 8, height = 6)
}

make_volcano(mcap_summer$res, "M. capitata Summer - B vs D", "volcano_Mcap_Summer.pdf")
make_volcano(mcap_winter$res, "M. capitata Winter - B vs D", "volcano_Mcap_Winter.pdf")
make_volcano(dtre_summer$res, "D. trenchii Summer - B vs D", "volcano_Dtre_Summer.pdf")
make_volcano(dtre_winter$res, "D. trenchii Winter - B vs D", "volcano_Dtre_Winter.pdf")

cat("Saved volcano plots to: figures/\n")

# ============================================================================
# SECTION 8: Top DEGs
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

print_top_degs(mcap_summer$res, "M. capitata Summer")
print_top_degs(mcap_winter$res, "M. capitata Winter")
print_top_degs(dtre_summer$res, "D. trenchii Summer")
print_top_degs(dtre_winter$res, "D. trenchii Winter")

cat("Positive LFC = higher expression in B (high CO2/acidified)\n")
cat("Negative LFC = lower expression in B (high CO2/acidified)\n")

# ============================================================================
# SECTION 9: Save R Objects
# ============================================================================

save(mcap_summer, mcap_winter,
     file = file.path(OUT_DIR, "Mcap_DESeq2.RData"))

save(dtre_summer, dtre_winter,
     file = file.path(OUT_DIR, "Dtre_DESeq2.RData"))

cat("\n\nSaved R objects to: 07_deseq2/\n")

cat("\n============================================\n")
cat("Analysis Complete!\n")
cat("============================================\n")
