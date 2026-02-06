#!/usr/bin/env Rscript
# =============================================================================
# Diagnostic: Check file formats for representative GO dendrogram plotting
# =============================================================================

cat("=== Diagnostic: Checking file formats ===\n\n")

# Test files
rep_file <- "/home/darmstrong4/mc_rework/10_GO_MWU/output/summer_BP_representative_GOs.csv"
mwu_file <- "/home/darmstrong4/mc_rework/10_GO_MWU/output/summer_BP_MWU_results.csv"
dissim_file <- "/home/darmstrong4/mc_rework/10_GO_MWU/output/summer_BP_dissimilarity.csv"

# 1. Check representative GOs
cat("1. REPRESENTATIVE GOs FILE\n")
cat("   Path:", rep_file, "\n")
if (file.exists(rep_file)) {
    rep_gos <- read.csv(rep_file)
    cat("   Columns:", paste(names(rep_gos), collapse = ", "), "\n")
    cat("   Rows:", nrow(rep_gos), "\n")
    cat("   Sample terms:\n")
    print(head(rep_gos[, c("term", "pval", "direction")], 3))
} else {
    cat("   FILE NOT FOUND\n")
}

# 2. Check MWU results
cat("\n2. MWU RESULTS FILE\n")
cat("   Path:", mwu_file, "\n")
if (file.exists(mwu_file)) {
    # MWU results are space/tab-separated, not CSV!
    mwu <- read.table(mwu_file, header = TRUE, stringsAsFactors = FALSE,
                      quote = "\"", fill = TRUE)
    cat("   Columns:", paste(names(mwu), collapse = ", "), "\n")
    cat("   Rows:", nrow(mwu), "\n")
    cat("   Sample:\n")
    print(head(mwu[, c("term", "name", "delta.rank", "pval")], 3))
} else {
    cat("   FILE NOT FOUND\n")
}

# 3. Check dissimilarity matrix
cat("\n3. DISSIMILARITY MATRIX FILE\n")
cat("   Path:", dissim_file, "\n")
if (file.exists(dissim_file)) {
    # Read first few lines to check format
    first_lines <- readLines(dissim_file, n = 3)
    cat("   First line (truncated):", substr(first_lines[1], 1, 100), "...\n")
    cat("   Detected separator: TAB\n")
    
    # Read WITHOUT row.names first to handle potential duplicates
    diss <- read.table(dissim_file, sep = "\t", header = TRUE, 
                       check.names = FALSE, nrows = 5)
    
    cat("   Matrix dimensions:", nrow(diss), "x", ncol(diss), "(showing first 5 rows)\n")
    cat("   First few col names:", paste(head(colnames(diss), 3), collapse = ", "), "\n")
    
    # Check for compound GO IDs (semicolon-separated)
    compound_ids <- sum(grepl(";", colnames(diss)))
    cat("   Compound GO IDs (with ';'):", compound_ids, "\n")
    
    # Check if column names look like GO IDs
    if (any(grepl("^GO:", colnames(diss)))) {
        cat("   Column names appear to be GO IDs ✓\n")
    } else {
        cat("   WARNING: Column names don't look like GO IDs\n")
        cat("   Sample col names:", paste(head(colnames(diss)), collapse = ", "), "\n")
    }
} else {
    cat("   FILE NOT FOUND\n")
}

# 4. Check term matching
cat("\n4. TERM MATCHING TEST\n")
if (exists("rep_gos") && exists("mwu")) {
    # Try to match first representative term to MWU results
    test_term <- rep_gos$term[1]
    cat("   Test term from rep_gos:", test_term, "\n")
    
    match_idx <- which(mwu$name == test_term)
    if (length(match_idx) > 0) {
        cat("   FOUND in MWU results ✓\n")
        cat("   GO ID:", mwu$term[match_idx[1]], "\n")
        cat("   delta.rank:", mwu$delta.rank[match_idx[1]], "\n")
    } else {
        cat("   NOT FOUND in MWU results - trying fuzzy match\n")
        # Try partial match
        partial <- grep(test_term, mwu$name, ignore.case = TRUE)
        if (length(partial) > 0) {
            cat("   Partial matches found:\n")
            print(mwu$name[partial[1:min(3, length(partial))]])
        }
    }
    
    # Check if GO ID is in dissimilarity matrix
    if (exists("diss") && length(match_idx) > 0) {
        go_id <- mwu$term[match_idx[1]]
        cat("\n   Checking GO ID in dissimilarity matrix:", go_id, "\n")
        
        # Check against column names (GO IDs are in columns, not row names)
        # Also check for compound IDs containing our target
        direct_match <- go_id %in% colnames(diss)
        compound_match <- any(grepl(go_id, colnames(diss), fixed = TRUE))
        
        if (direct_match) {
            cat("   FOUND as direct match in dissimilarity matrix ✓\n")
        } else if (compound_match) {
            cat("   FOUND as part of compound GO ID ✓\n")
            matching_cols <- colnames(diss)[grepl(go_id, colnames(diss), fixed = TRUE)]
            cat("   Compound ID:", matching_cols[1], "\n")
        } else {
            cat("   NOT FOUND in first 5 rows - loading more to check...\n")
            full_diss <- read.table(dissim_file, sep = "\t", header = TRUE,
                                    check.names = FALSE, nrows = 100)
            if (go_id %in% colnames(full_diss)) {
                cat("   FOUND in larger sample ✓\n")
            } else if (any(grepl(go_id, colnames(full_diss), fixed = TRUE))) {
                cat("   FOUND as compound ID in larger sample ✓\n")
            } else {
                cat("   NOT FOUND\n")
            }
        }
    }
}

cat("\n=== Diagnostic Complete ===\n")
