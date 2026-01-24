#!/usr/bin/env Rscript
# ==============================================================================
# GO_MWU Analysis: M. capitata Ocean Acidification Response
# ==============================================================================

# Set working directory
setwd("/home/darmstrong4/mc_rework/10_GO_MWU")

# Load GO_MWU functions
source("gomwu.functions.R")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Copy input files to working directory (GO_MWU expects files in current dir)
file.copy("input/go_annotations.tab", "go_annotations.tab", overwrite = TRUE)
file.copy("input/go.obo", "go.obo", overwrite = TRUE)
file.copy("input/summer_log2fc.csv", "summer_log2fc.csv", overwrite = TRUE)
file.copy("input/winter_log2fc.csv", "winter_log2fc.csv", overwrite = TRUE)

# Input files (now in working directory)
goAnnotations <- "go_annotations.tab"
goDatabase <- "go.obo"

# Seasons and GO divisions to analyze
seasons <- c("summer", "winter")
divisions <- c("BP", "MF", "CC")

# GO_MWU parameters
largest <- 0.1
smallest <- 5
clusterCutHeight <- 0.25
absValue <- 1

# FDR thresholds for plotting
level1 <- 0.1
level2 <- 0.05
level3 <- 0.01

# ==============================================================================
# RUN ANALYSIS FOR ALL COMBINATIONS
# ==============================================================================

results_summary <- data.frame()

for (season in seasons) {
    
    input_file <- paste0(season, "_log2fc.csv")  # No path prefix now
    
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
file.remove(c("go_annotations.tab", "go.obo", "summer_log2fc.csv", "winter_log2fc.csv"))

# ==============================================================================
# SAVE SUMMARY
# ==============================================================================

cat("\n")
cat("##############################################################\n")
cat("# ANALYSIS COMPLETE\n")
cat("##############################################################\n\n")

print(results_summary)

write.csv(results_summary, "output/GO_MWU_summary.csv", row.names = FALSE)
cat("\nSummary saved to: output/GO_MWU_summary.csv\n")

cat("\n=== Session Info ===\n")
sessionInfo()
