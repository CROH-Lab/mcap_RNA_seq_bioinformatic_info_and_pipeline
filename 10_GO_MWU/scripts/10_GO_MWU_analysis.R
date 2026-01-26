#!/usr/bin/env Rscript
# ==============================================================================
# GO_MWU Analysis: M. capitata Ocean Acidification Response
# Using signed -log10(pvalue) as the measure
# Modular configuration per GO division
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
# MODULAR CONFIGURATION - ADJUST PARAMETERS PER DIVISION
# ==============================================================================

# Configuration list for each GO division
# Adjust these independently as needed
config <- list(
    BP = list(
        largest = 0.1,
        smallest = 10,
        clusterCutHeight = 0.9,
        level1 = 0.05,
        level2 = 0.01,
        level3 = 0.001,
        txtsize = 0.9,
        treeHeight = 0.9
    ),
    MF = list(
        largest = 0.1,
        smallest = 10,
        clusterCutHeight = 0.25,
        level1 = 0.05,
        level2 = 0.01,
        level3 = 0.001,
        txtsize = 1.0,
        treeHeight = 0.25
    ),
    CC = list(
        largest = 0.1,
        smallest = 10,
        clusterCutHeight = 0.75,
        level1 = 0.05,
        level2 = 0.01,
        level3 = 0.001,
        txtsize = 1.0,
        treeHeight = 0.75
    )
)

# Global settings
absValue <- -log10(0.05)  # ~1.3, genes with p < 0.05 are "good candidates"
seasons <- c("summer", "winter")
divisions <- c("BP", "MF", "CC")

# ==============================================================================
# COPY INPUT FILES TO WORKING DIRECTORY
# ==============================================================================

file.copy("input/go_annotations.tab", "go_annotations.tab", overwrite = TRUE)
file.copy("input/go.obo", "go.obo", overwrite = TRUE)
file.copy("input/summer_signed_logP.csv", "summer_signed_logP.csv", overwrite = TRUE)
file.copy("input/winter_signed_logP.csv", "winter_signed_logP.csv", overwrite = TRUE)

goAnnotations <- "go_annotations.tab"
goDatabase <- "go.obo"

# ==============================================================================
# RUN ANALYSIS FOR ALL COMBINATIONS
# ==============================================================================

results_summary <- data.frame()

# Store results for combined plots
plot_results <- list()

for (season in seasons) {
    
    input_file <- paste0(season, "_signed_logP.csv")
    plot_results[[season]] <- list()
    
    cat("\n")
    cat("##############################################################\n")
    cat("# Processing:", toupper(season), "\n")
    cat("##############################################################\n")
    
    for (div in divisions) {
        
        cat("\n=== ", season, " - ", div, " ===\n", sep="")
        cat("Parameters: smallest=", config[[div]]$smallest, 
            ", clusterCutHeight=", config[[div]]$clusterCutHeight,
            ", level1=", config[[div]]$level1, "\n", sep="")
        
        prefix <- paste0(season, "_", div)
        
        tryCatch({
            # Run MWU statistics with division-specific parameters
            gomwuStats(input_file, goDatabase, goAnnotations, div,
                       perlPath = "perl",
                       largest = config[[div]]$largest,
                       smallest = config[[div]]$smallest,
                       clusterCutHeight = config[[div]]$clusterCutHeight
            )
            
            # Generate individual plot
            pdf_file <- paste0("figures/", prefix, "_GO_MWU.pdf")
            pdf(pdf_file, width = 12, height = 10)
            
            results <- gomwuPlot(input_file, goAnnotations, div,
                                 absValue = absValue,
                                 level1 = config[[div]]$level1,
                                 level2 = config[[div]]$level2,
                                 level3 = config[[div]]$level3,
                                 txtsize = config[[div]]$txtsize,
                                 treeHeight = config[[div]]$treeHeight,
                                 colors = c("dodgerblue2", "firebrick1", "skyblue2", "lightcoral")
            )
            
            dev.off()
            cat("Saved individual plot:", pdf_file, "\n")
            
            # Store for combined plot
            plot_results[[season]][[div]] <- list(
                input_file = input_file,
                success = TRUE
            )
            
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
            plot_results[[season]][[div]] <<- list(success = FALSE)
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

# ==============================================================================
# CREATE COMBINED PLOTS (BP left, CC top-right, MF bottom-right)
# ==============================================================================

cat("\n")
cat("##############################################################\n")
cat("# Creating Combined Plots\n")
cat("##############################################################\n")

for (season in seasons) {
    
    cat("\nGenerating combined plot for", toupper(season), "...\n")
    
    input_file <- paste0(season, "_signed_logP.csv")
    combined_pdf <- paste0("figures/", season, "_combined_GO_MWU.pdf")
    
    tryCatch({
        # Create combined PDF with layout: BP (left), CC (top-right), MF (bottom-right)
        pdf(combined_pdf, width = 20, height = 14)
        
        # Set up layout matrix:
        # 1 1 2
        # 1 1 3
        layout(matrix(c(1, 1, 2,
                        1, 1, 3), nrow = 2, byrow = TRUE))
        
        # Plot 1: BP (large, left side)
        par(mar = c(4, 4, 4, 2))
        gomwuPlot(input_file, goAnnotations, "BP",
                  absValue = absValue,
                  level1 = config[["BP"]]$level1,
                  level2 = config[["BP"]]$level2,
                  level3 = config[["BP"]]$level3,
                  txtsize = config[["BP"]]$txtsize * 0.8,
                  treeHeight = config[["BP"]]$treeHeight,
                  colors = c("dodgerblue2", "firebrick1", "skyblue2", "lightcoral")
        )
        title(main = paste0(toupper(season), " - Biological Process (BP)"), 
              cex.main = 1.5, line = 1)
        
        # Plot 2: CC (top-right)
        par(mar = c(4, 4, 4, 2))
        gomwuPlot(input_file, goAnnotations, "CC",
                  absValue = absValue,
                  level1 = config[["CC"]]$level1,
                  level2 = config[["CC"]]$level2,
                  level3 = config[["CC"]]$level3,
                  txtsize = config[["CC"]]$txtsize * 0.8,
                  treeHeight = config[["CC"]]$treeHeight,
                  colors = c("dodgerblue2", "firebrick1", "skyblue2", "lightcoral")
        )
        title(main = paste0(toupper(season), " - Cellular Component (CC)"), 
              cex.main = 1.5, line = 1)
        
        # Plot 3: MF (bottom-right)
        par(mar = c(4, 4, 4, 2))
        gomwuPlot(input_file, goAnnotations, "MF",
                  absValue = absValue,
                  level1 = config[["MF"]]$level1,
                  level2 = config[["MF"]]$level2,
                  level3 = config[["MF"]]$level3,
                  txtsize = config[["MF"]]$txtsize * 0.8,
                  treeHeight = config[["MF"]]$treeHeight,
                  colors = c("dodgerblue2", "firebrick1", "skyblue2", "lightcoral")
        )
        title(main = paste0(toupper(season), " - Molecular Function (MF)"), 
              cex.main = 1.5, line = 1)
        
        dev.off()
        cat("Saved combined plot:", combined_pdf, "\n")
        
    }, error = function(e) {
        cat("ERROR creating combined plot for", season, ":", conditionMessage(e), "\n")
        if (dev.cur() > 1) dev.off()
    })
}

# ==============================================================================
# CLEANUP
# ==============================================================================

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

cat("Division-specific parameters:\n")
for (div in divisions) {
    cat("  ", div, ": smallest=", config[[div]]$smallest,
        ", clusterCutHeight=", config[[div]]$clusterCutHeight,
        ", level1=", config[[div]]$level1, "\n", sep="")
}

cat("\n")
print(results_summary)

write.csv(results_summary, "output/GO_MWU_summary.csv", row.names = FALSE)
cat("\nSummary saved to: output/GO_MWU_summary.csv\n")

cat("\n=== Output Files ===\n")
cat("Individual plots: figures/<season>_<division>_GO_MWU.pdf\n")
cat("Combined plots: figures/<season>_combined_GO_MWU.pdf\n")

cat("\n=== Session Info ===\n")
sessionInfo()
