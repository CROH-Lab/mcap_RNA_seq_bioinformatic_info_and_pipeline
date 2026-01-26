#!/usr/bin/env Rscript
# ==============================================================================
# GO_MWU Analysis: M. capitata Ocean Acidification Response
# Using signed -log10(pvalue) as the measure
# Modular configuration per GO division
# With representative GO extraction
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

config <- list(
    BP = list(
        largest = 0.1,
        smallest = 10,
        clusterCutHeight = 0.5,
        level1 = 0.05,
        level2 = 0.01,
        level3 = 0.001,
        txtsize = 0.9,
        treeHeight = 0.5
    ),
    MF = list(
        largest = 0.1,
        smallest = 10,
        clusterCutHeight = 0.5,
        level1 = 0.05,
        level2 = 0.01,
        level3 = 0.001,
        txtsize = 1.0,
        treeHeight = 0.5
    ),
    CC = list(
        largest = 0.1,
        smallest = 10,
        clusterCutHeight = 0.5,
        level1 = 0.05,
        level2 = 0.01,
        level3 = 0.001,
        txtsize = 1.0,
        treeHeight = 0.5
    )
)

# Representative GO extraction parameters
rep_GO_config <- list(
    pcut = 0.05,   # adjusted p-value cutoff for representative GOs
    hcut = 0.9     # height to cut tree for independent groups
)

absValue <- -log10(0.05)
seasons <- c("summer", "winter")
divisions <- c("BP", "MF", "CC")

# ==============================================================================
# FUNCTION: Extract Representative GOs
# ==============================================================================

extract_representative_GOs <- function(results, pcut = 0.01, hcut = 0.9) {
    # results is the output from gomwuPlot()
    # results[[1]] = data frame of significant GO terms
    # results[[2]] = hierarchical clustering tree
    
    if (is.null(results[[1]]) || nrow(results[[1]]) == 0) {
        return(NULL)
    }
    
    # Cut tree to get independent groups
    ct <- cutree(results[[2]], h = hcut)
    
    rep_GOs <- data.frame()
    
    for (ci in unique(ct)) {
        rn <- names(ct)[ct == ci]
        
        # Remove obsolete terms
        obs <- grep("obsolete", rn)
        if (length(obs) > 0) { rn <- rn[-obs] }
        if (length(rn) == 0) { next }
        
        # Get results for this cluster
        rr <- results[[1]][rn, , drop = FALSE]
        
        # Find best (lowest p-value) term
        bestrr <- rr[which(rr$pval == min(rr$pval)), , drop = FALSE]
        best <- 1
        
        # If tie, choose the one with most genes
        if (nrow(bestrr) > 1) {
            nns <- sub(" .+", "", row.names(bestrr))
            fr <- c()
            for (i in 1:length(nns)) { 
                fr <- c(fr, eval(parse(text = nns[i]))) 
            }
            best <- which(fr == max(fr))[1]
        }
        
        # Add if passes p-value cutoff
        if (bestrr$pval[best] <= pcut) {
            term_name <- sub("\\d+\\/\\d+ ", "", row.names(bestrr)[best])
            rep_GOs <- rbind(rep_GOs, data.frame(
                cluster = ci,
                term = term_name,
                pval = bestrr$pval[best],
                p.adj = bestrr$p.adj[best],
                direction = ifelse(bestrr$direction[best] == 1, "up", "down"),
                nseqs = bestrr$nseqs[best],
                stringsAsFactors = FALSE
            ))
        }
    }
    
    # Sort by p-value
    if (nrow(rep_GOs) > 0) {
        rep_GOs <- rep_GOs[order(rep_GOs$pval), ]
        rownames(rep_GOs) <- NULL
    }
    
    return(rep_GOs)
}

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
all_rep_GOs <- list()

for (season in seasons) {
    
    input_file <- paste0(season, "_signed_logP.csv")
    all_rep_GOs[[season]] <- list()
    
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
            cat("Saved plot:", pdf_file, "\n")
            
            # Extract representative GOs
            rep_GOs <- extract_representative_GOs(results, 
                                                   pcut = rep_GO_config$pcut, 
                                                   hcut = rep_GO_config$hcut)
            
            if (!is.null(rep_GOs) && nrow(rep_GOs) > 0) {
                rep_GOs$season <- season
                rep_GOs$division <- div
                all_rep_GOs[[season]][[div]] <- rep_GOs
                
                # Save representative GOs
                write.csv(rep_GOs, paste0("output/", prefix, "_representative_GOs.csv"), 
                          row.names = FALSE)
                cat("Representative GO terms:", nrow(rep_GOs), "\n")
            }
            
            # Move output files
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
            
            n_rep <- ifelse(!is.null(rep_GOs), nrow(rep_GOs), 0)
            
            results_summary <- rbind(results_summary, data.frame(
                Season = season,
                Division = div,
                Significant_Terms = n_sig,
                Representative_Terms = n_rep,
                Upregulated = n_up,
                Downregulated = n_down
            ))
            
        }, error = function(e) {
            cat("ERROR in", season, div, ":", conditionMessage(e), "\n")
            results_summary <<- rbind(results_summary, data.frame(
                Season = season,
                Division = div,
                Significant_Terms = NA,
                Representative_Terms = NA,
                Upregulated = NA,
                Downregulated = NA
            ))
        })
    }
}

# ==============================================================================
# CREATE COMBINED PLOTS
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
        pdf(combined_pdf, width = 20, height = 14)
        layout(matrix(c(1, 1, 2,
                        1, 1, 3), nrow = 2, byrow = TRUE))
        
        # BP (left)
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
        title(main = paste0(toupper(season), " - Biological Process (BP)"), cex.main = 1.5)
        
        # CC (top-right)
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
        title(main = paste0(toupper(season), " - Cellular Component (CC)"), cex.main = 1.5)
        
        # MF (bottom-right)
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
        title(main = paste0(toupper(season), " - Molecular Function (MF)"), cex.main = 1.5)
        
        dev.off()
        cat("Saved combined plot:", combined_pdf, "\n")
        
    }, error = function(e) {
        cat("ERROR creating combined plot:", conditionMessage(e), "\n")
        if (dev.cur() > 1) dev.off()
    })
}

# ==============================================================================
# CREATE REPRESENTATIVE GO SUMMARY TABLE
# ==============================================================================

cat("\n")
cat("##############################################################\n")
cat("# Creating Representative GO Summary\n")
cat("##############################################################\n")

# Combine all representative GOs into one table
all_rep_combined <- do.call(rbind, lapply(seasons, function(s) {
    do.call(rbind, all_rep_GOs[[s]])
}))

if (!is.null(all_rep_combined) && nrow(all_rep_combined) > 0) {
    # Reorder columns
    all_rep_combined <- all_rep_combined[, c("season", "division", "term", "direction", 
                                              "pval", "p.adj", "nseqs", "cluster")]
    
    write.csv(all_rep_combined, "output/all_representative_GOs.csv", row.names = FALSE)
    cat("\nCombined representative GOs saved to: output/all_representative_GOs.csv\n")
    
    # Print summary by season and direction
    cat("\n=== Representative GO Summary ===\n\n")
    for (s in seasons) {
        cat(toupper(s), ":\n")
        for (d in divisions) {
            subset_data <- all_rep_combined[all_rep_combined$season == s & 
                                            all_rep_combined$division == d, ]
            if (nrow(subset_data) > 0) {
                n_up <- sum(subset_data$direction == "up")
                n_down <- sum(subset_data$direction == "down")
                cat("  ", d, ": ", nrow(subset_data), " representative terms (", 
                    n_up, " up, ", n_down, " down)\n", sep="")
            }
        }
        cat("\n")
    }
    
    # Print top representative terms for each
    cat("=== Top Representative GO Terms ===\n")
    for (s in seasons) {
        cat("\n", toupper(s), ":\n", sep="")
        for (d in divisions) {
            subset_data <- all_rep_combined[all_rep_combined$season == s & 
                                            all_rep_combined$division == d, ]
            if (nrow(subset_data) > 0) {
                cat("  ", d, ":\n", sep="")
                top_terms <- head(subset_data[order(subset_data$pval), ], 5)
                for (i in 1:nrow(top_terms)) {
                    cat("    ", top_terms$direction[i], ": ", top_terms$term[i], 
                        " (p.adj=", format(top_terms$p.adj[i], digits=3), ")\n", sep="")
                }
            }
        }
    }
}

# ==============================================================================
# CLEANUP
# ==============================================================================

int_files <- list.files(pattern = "^BP_|^MF_|^CC_|^dissim_|^MWU_|\\.tmp$")
if (length(int_files) > 0) file.remove(int_files)
file.remove(c("go_annotations.tab", "go.obo", "summer_signed_logP.csv", "winter_signed_logP.csv"))

# ==============================================================================
# SAVE SUMMARY
# ==============================================================================

cat("\n")
cat("##############################################################\n")
cat("# ANALYSIS COMPLETE\n")
cat("##############################################################\n\n")

print(results_summary)

write.csv(results_summary, "output/GO_MWU_summary.csv", row.names = FALSE)

cat("\n=== Output Files ===\n")
cat("Individual plots: figures/<season>_<division>_GO_MWU.pdf\n")
cat("Combined plots: figures/<season>_combined_GO_MWU.pdf\n")
cat("Full results: output/<season>_<division>_MWU_results.csv\n")
cat("Representative GOs: output/<season>_<division>_representative_GOs.csv\n")
cat("All representative GOs: output/all_representative_GOs.csv\n")

cat("\n=== Session Info ===\n")
sessionInfo()
