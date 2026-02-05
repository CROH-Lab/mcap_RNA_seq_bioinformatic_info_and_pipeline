#!/usr/bin/env Rscript
# ==============================================================================
# GO_MWU Analysis: D. trenchii Symbiont Ocean Acidification Response
# Using signed -log10(pvalue) as the measure
# Modular configuration per GO division
# With robust representative GO extraction
# ==============================================================================

# Set working directory
setwd("/home/darmstrong4/mc_rework/11_symbiont_GO_MWU")

# Load GO_MWU functions
source("gomwu.functions.R")

# ==============================================================================
# CREATE SIGNED P-VALUE INPUTS FROM DESEQ2 RESULTS
# ==============================================================================

cat("==============================================================================\n")
cat("SYMBIONT GO_MWU Analysis: D. trenchii\n")
cat("Creating signed -log10(pvalue) inputs from DESeq2 results\n")
cat("==============================================================================\n\n")

# Function to create GO_MWU input from DESeq2 results
create_gomwu_input <- function(deseq_file, output_file) {
    deseq <- read.csv(deseq_file, row.names = 1)
    deseq_clean <- deseq[!is.na(deseq$pvalue) & !is.na(deseq$log2FoldChange), ]
    min_pval <- min(deseq_clean$pvalue[deseq_clean$pvalue > 0])
    deseq_clean$pvalue[deseq_clean$pvalue == 0] <- min_pval
    signed_logP <- sign(deseq_clean$log2FoldChange) * -log10(deseq_clean$pvalue)

    gomwu_input <- data.frame(
        gene = rownames(deseq_clean),
        signed_logP = signed_logP
    )
    write.csv(gomwu_input, output_file, row.names = FALSE, quote = FALSE)

    cat("Input file:", basename(deseq_file), "\n")
    cat("  Total genes:", nrow(deseq), "\n")
    cat("  Genes with valid p-values:", nrow(deseq_clean), "\n")
    cat("  Signed -log10(p) range:", round(min(signed_logP), 2), "to", round(max(signed_logP), 2), "\n")
    cat("  Genes with |signed_logP| > 1.3 (p < 0.05):", sum(abs(signed_logP) > 1.3), "\n\n")

    return(gomwu_input)
}

# Create inputs for both seasons - SYMBIONT (Dtre)
summer_input <- create_gomwu_input(
    "/home/darmstrong4/mc_rework/07_deseq2/results/Dtre_Summer_BvsD.csv",
    "input/symbiont_summer_signed_logP.csv"
)

winter_input <- create_gomwu_input(
    "/home/darmstrong4/mc_rework/07_deseq2/results/Dtre_Winter_BvsD.csv",
    "input/symbiont_winter_signed_logP.csv"
)

# ==============================================================================
# MODULAR CONFIGURATION
# ==============================================================================

config <- list(
    BP = list(
        largest = 0.1, smallest = 10, clusterCutHeight = 0.25,
        level1 = 0.05, level2 = 0.01, level3 = 0.001,
        txtsize = 1, treeHeight = 0.5
    ),
    MF = list(
        largest = 0.1, smallest = 10, clusterCutHeight = 0.25,
        level1 = 0.05, level2 = 0.01, level3 = 0.001,
        txtsize = 1.0, treeHeight = 0.5
    ),
    CC = list(
        largest = 0.1, smallest = 10, clusterCutHeight = 0.25,
        level1 = 0.05, level2 = 0.01, level3 = 0.001,
        txtsize = 1.0, treeHeight = 0.5
    )
)

absValue <- -log10(0.05)
seasons <- c("summer", "winter")
divisions <- c("BP", "MF", "CC")

# ==============================================================================
# COPY INPUT FILES TO WORKING DIRECTORY
# ==============================================================================

file.copy("input/symbiont_go_annotations.tab", "symbiont_go_annotations.tab", overwrite = TRUE)
file.copy("input/go.obo", "go.obo", overwrite = TRUE)
file.copy("input/symbiont_summer_signed_logP.csv", "symbiont_summer_signed_logP.csv", overwrite = TRUE)
file.copy("input/symbiont_winter_signed_logP.csv", "symbiont_winter_signed_logP.csv", overwrite = TRUE)

goAnnotations <- "symbiont_go_annotations.tab"
goDatabase <- "go.obo"

# ==============================================================================
# RUN ANALYSIS FOR ALL COMBINATIONS
# ==============================================================================

results_summary <- data.frame()

for (season in seasons) {
    input_file <- paste0("symbiont_", season, "_signed_logP.csv")

    cat("\n##############################################################\n")
    cat("# Processing SYMBIONT:", toupper(season), "\n")
    cat("##############################################################\n")

    for (div in divisions) {
        cat("\n=== Symbiont ", season, " - ", div, " ===\n", sep="")
        prefix <- paste0("symbiont_", season, "_", div)

        tryCatch({
            # Run MWU statistics
            gomwuStats(input_file, goDatabase, goAnnotations, div,
                       perlPath = "perl",
                       largest = config[[div]]$largest,
                       smallest = config[[div]]$smallest,
                       clusterCutHeight = config[[div]]$clusterCutHeight)

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
                                 colors = c("dodgerblue2", "firebrick1", "skyblue2", "lightcoral"))
            dev.off()
            cat("Saved plot:", pdf_file, "\n")

            # ==================================================================
            # MOVE OUTPUT FILES - Using pattern that works for host analysis
            # ==================================================================
            mwu_file <- list.files(pattern = paste0("^MWU_", div, "_"))
            if (length(mwu_file) > 0) {
                file.rename(mwu_file[1], paste0("output/", prefix, "_MWU_results.csv"))
                cat("Moved MWU results to output/\n")
            } else {
                cat("Warning: MWU file not found with pattern: ^MWU_", div, "_\n", sep="")
            }

            dissim_file <- list.files(pattern = paste0("^dissim_", div, "_"))
            if (length(dissim_file) > 0) {
                file.rename(dissim_file[1], paste0("output/", prefix, "_dissimilarity.csv"))
                cat("Moved dissimilarity to output/\n")
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
                Organism = "Symbiont",
                Season = season,
                Division = div,
                Significant_Terms = n_sig,
                Upregulated = n_up,
                Downregulated = n_down
            ))

        }, error = function(e) {
            cat("ERROR in symbiont", season, div, ":", conditionMessage(e), "\n")
            results_summary <<- rbind(results_summary, data.frame(
                Organism = "Symbiont", Season = season, Division = div,
                Significant_Terms = NA, Upregulated = NA, Downregulated = NA
            ))
        })
    }
}

# ==============================================================================
# CREATE COMBINED PLOTS
# ==============================================================================

cat("\n##############################################################\n")
cat("# Creating Combined Plots - SYMBIONT\n")
cat("##############################################################\n")

for (season in seasons) {
    cat("\nGenerating combined plot for SYMBIONT", toupper(season), "...\n")
    input_file <- paste0("symbiont_", season, "_signed_logP.csv")
    combined_pdf <- paste0("figures/symbiont_", season, "_combined_GO_MWU.pdf")

    tryCatch({
        pdf(combined_pdf, width = 20, height = 14)
        layout(matrix(c(1, 1, 2, 1, 1, 3), nrow = 2, byrow = TRUE))

        par(mar = c(4, 4, 4, 2))
        gomwuPlot(input_file, goAnnotations, "BP",
                  absValue = absValue, level1 = config[["BP"]]$level1,
                  level2 = config[["BP"]]$level2, level3 = config[["BP"]]$level3,
                  txtsize = config[["BP"]]$txtsize * 0.8,
                  treeHeight = config[["BP"]]$treeHeight,
                  colors = c("dodgerblue2", "firebrick1", "skyblue2", "lightcoral"))
        title(main = paste0("SYMBIONT ", toupper(season), " - BP"), cex.main = 1.5)

        par(mar = c(4, 4, 4, 2))
        gomwuPlot(input_file, goAnnotations, "CC",
                  absValue = absValue, level1 = config[["CC"]]$level1,
                  level2 = config[["CC"]]$level2, level3 = config[["CC"]]$level3,
                  txtsize = config[["CC"]]$txtsize * 0.8,
                  treeHeight = config[["CC"]]$treeHeight,
                  colors = c("dodgerblue2", "firebrick1", "skyblue2", "lightcoral"))
        title(main = paste0("SYMBIONT ", toupper(season), " - CC"), cex.main = 1.5)

        par(mar = c(4, 4, 4, 2))
        gomwuPlot(input_file, goAnnotations, "MF",
                  absValue = absValue, level1 = config[["MF"]]$level1,
                  level2 = config[["MF"]]$level2, level3 = config[["MF"]]$level3,
                  txtsize = config[["MF"]]$txtsize * 0.8,
                  treeHeight = config[["MF"]]$treeHeight,
                  colors = c("dodgerblue2", "firebrick1", "skyblue2", "lightcoral"))
        title(main = paste0("SYMBIONT ", toupper(season), " - MF"), cex.main = 1.5)

        dev.off()
        cat("Saved combined plot:", combined_pdf, "\n")
    }, error = function(e) {
        cat("ERROR creating combined plot:", conditionMessage(e), "\n")
        if (dev.cur() > 1) dev.off()
    })
}

# ==============================================================================
# CLEANUP
# ==============================================================================

int_files <- list.files(pattern = "^BP_|^MF_|^CC_|^dissim_|^MWU_|^cl_|\\.tmp$")
if (length(int_files) > 0) {
    cat("\nCleaning up", length(int_files), "intermediate files...\n")
    file.remove(int_files)
}

input_copies <- c("symbiont_go_annotations.tab", "go.obo",
                  "symbiont_summer_signed_logP.csv", "symbiont_winter_signed_logP.csv")
existing_copies <- input_copies[file.exists(input_copies)]
if (length(existing_copies) > 0) file.remove(existing_copies)

# ==============================================================================
# SAVE SUMMARY
# ==============================================================================

cat("\n##############################################################\n")
cat("# SYMBIONT GO_MWU ANALYSIS COMPLETE\n")
cat("##############################################################\n\n")

print(results_summary)
write.csv(results_summary, "output/symbiont_GO_MWU_summary.csv", row.names = FALSE)

cat("\n=== Output Files ===\n")
cat("Plots: figures/symbiont_<season>_<division>_GO_MWU.pdf\n")
cat("Combined: figures/symbiont_<season>_combined_GO_MWU.pdf\n")
cat("Results: output/symbiont_<season>_<division>_MWU_results.csv\n")

cat("\n=== Session Info ===\n")
sessionInfo()
