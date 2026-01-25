#!/usr/bin/env Rscript
# ==============================================================================
# GO_MWU Analysis: M. capitata Ocean Acidification Response
# Using signed -log10(pvalue) as the measure
# ==============================================================================

# Set working directory
setwd("/home/darmstrong4/mc_rework/10_GO_MWU")

# Load GO_MWU functions
source("gomwu.functions.R")

# ==============================================================================
# CREATE SIGNED P-VALUE INPUTS FROM DESEQ2 RESULTS
# ==============================================================================

cat("==============================================================================\n")
cat("Creating signed -log10(pvalue) inputs from DESeq2 results\n")
cat("==============================================================================\n\n")

# Function to create GO_MWU input from DESeq2 results
create_gomwu_input <- function(deseq_file, output_file) {
    
    # Read DESeq2 results
    deseq <- read.csv(deseq_file, row.names = 1)
    
    # Filter for genes with valid p-values and log2FC
    deseq_clean <- deseq[!is.na(deseq$pvalue) & !is.na(deseq$log2FoldChange), ]
    
    # Handle p-values of 0 (set to minimum non-zero p-value)
    min_pval <- min(deseq_clean$pvalue[deseq_clean$pvalue > 0])
    deseq_clean$pvalue[deseq_clean$pvalue == 0] <- min_pval
    
    # Create signed -log10(pvalue)
    # Positive values = upregulated, Negative values = downregulated
    signed_logP <- sign(deseq_clean$log2FoldChange) * -log10(deseq_clean$pvalue)
    
    # Create output dataframe
    gomwu_input <- data.frame(
        gene = rownames(deseq_clean),
        signed_logP = signed_logP
    )
    
    # Write to file
    write.csv(gomwu_input, output_file, row.names = FALSE, quote = FALSE)
    
    # Report stats
    cat("Input file:", basename(deseq_file), "\n")
    cat("  Total genes:", nrow(deseq), "\n")
    cat("  Genes with valid p-values:", nrow(deseq_clean), "\n")
    cat("  Signed -log10(p) range:", round(min(signed_logP), 2), "to", round(max(signed_logP), 2), "\n")
    cat("  Genes with |signed_logP| > 1.3 (p < 0.05):", sum(abs(signed_logP) > 1.3), "\n")
    cat("  Output:", output_file, "\n\n")
    
    return(gomwu_input)
}

# Create inputs for both seasons
summer_input <- create_gomwu_input(
    "/home/darmstrong4/mc_rework/07_deseq2/results/Mcap_Summer_BvsD.csv",
    "input/summer_signed_logP.csv"
)

winter_input <- create_gomwu_input(
    "/home/darmstrong4/mc_rework/07_deseq2/results/Mcap_Winter_BvsD.csv",
    "input/winter_signed_logP.csv"
)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Copy input files to working directory (GO_MWU expects files in current dir)
file.copy("input/go_annotations.tab", "go_annotations.tab", overwrite = TRUE)
file.copy("input/go.obo", "go.obo", overwrite = TRUE)
file.copy("input/summer_signed_logP.csv", "summer_signed_logP.csv", overwrite = TRUE)
file.copy("input/winter_signed_logP.csv", "winter_signed_logP.csv", overwrite = TRUE)

# Input files (now in working directory)
goAnnotations <- "go_annotations.tab"
goDatabase <- "go.obo"

# Seasons and GO divisions to analyze
seasons <- c("summer", "winter")
divisions <- c("BP", "MF", "CC")

# GO_MWU parameters
largest <- 0.1
smallest <- 10
clusterCutHeight <- 0.75

# absValue for signed -log10(pvalue):
# -log10(0.05) = 1.3  -> genes with p < 0.05 are "good candidates"
# -log10(0.01) = 2.0  -> genes with p < 0.01 are "good candidates"
absValue <- -log10(0.05)  # ~1.3

# FDR thresholds for plotting
level1 <- 0.05
level2 <- 0.01
level3 <- 0.001

# ==============================================================================
# RUN ANALYSIS FOR ALL COMBINATIONS
# ==============================================================================

results_summary <- data.frame()

for (season in seasons) {
    
    input_file <- paste0(season, "_signed_logP.csv")
    
    cat("\n")
    cat("##############################################################\n")
    cat("# Processing:", toupper(season), "\n")
    cat("##############################################################\n")
    
    for (div in divisions) {
        
        cat("\n=== ", season, " - ", div, " ===\n", sep="")
        
        prefix <- paste0(season, "_", div)
        
        tryCatch({
            # Run MWU statistics
            gomwuStats(input_file, goDatabase, goAnnotations, div,
                       perlPath = "perl",
                       largest = largest,
                       smallest = smallest,
                       clusterCutHeight = clusterCutHeight
            )
            
            # Generate plot
            pdf_file <- paste0("figures/", prefix, "_GO_MWU.pdf")
            pdf(pdf_file, width = 12, height = 10)
            
            results <- gomwuPlot(input_file, goAnnotations, div,
                                 absValue = absValue,
                                 level1 = level1,
                                 level2 = level2,
                                 level3 = level3,
                                 txtsize = 1.0,
                                 treeHeight = 0.5,
                                 colors = c("dodgerblue2", "firebrick1", "skyblue2", "lightcoral")
            )
            
            dev.off()
            cat("Saved plot:", pdf_file, "\n")
            
            # Move output files to output directory
            mwu_file <- list.files(pattern = paste0("^MWU_", div, "_"))
            if (length(mwu_file) > 0) {
                file.rename(mwu_file[1], paste0("output/", prefix, "_MWU_results.csv"))
            }
            
            dissim_file <- list.files(pattern = paste0("^dissim_", div, "_"))
            if (length(dissim_file) > 0) {
                file.rename(dissim_file[1], paste0("output/", prefix, "_dissimilarity.csv"))
            }
            
            # Record summary
            if (!is.null(results[[1]]) && nrow(results[[1]]) > 0) {
                n_sig <- nrow(results[[1]])
                n_up <- sum(results[[1]]$direction == 1)
                n_down <- sum(results[[1]]$direction == 0)
            } else {
                n_sig <- 0
                n_up <- 0
                n_down <- 0
            }
            
            results_summary <- rbind(results_summary, data.frame(
                Season = season,
                Division = div,
                Significant_Terms = n_sig,
                Upregulated = n_up,
                Downregulated = n_down
            ))
            
            cat("Significant GO terms:", n_sig, "(Up:", n_up, "Down:", n_down, ")\n")
            
        }, error = function(e) {
            cat("ERROR in", season, div, ":", conditionMessage(e), "\n")
            results_summary <<- rbind(results_summary, data.frame(
                Season = season,
                Division = div,
                Significant_Terms = NA,
                Upregulated = NA,
                Downregulated = NA
            ))
        })
    }
}

# Clean up intermediate files from working directory
int_files <- list.files(pattern = "^BP_|^MF_|^CC_|^dissim_|^MWU_|\\.tmp$")
if (length(int_files) > 0) file.remove(int_files)

# Clean up copied input files
file.remove(c("go_annotations.tab", "go.obo", "summer_signed_logP.csv", "winter_signed_logP.csv"))

# ==============================================================================
# SAVE SUMMARY
# ==============================================================================

cat("\n")
cat("##############################################################\n")
cat("# ANALYSIS COMPLETE\n")
cat("##############################################################\n\n")

cat("Method: Signed -log10(pvalue) from DESeq2 results\n")
cat("absValue threshold:", round(absValue, 3), "(corresponds to p < 0.05)\n\n")

print(results_summary)

write.csv(results_summary, "output/GO_MWU_summary.csv", row.names = FALSE)
cat("\nSummary saved to: output/GO_MWU_summary.csv\n")

cat("\n=== Session Info ===\n")
sessionInfo()
