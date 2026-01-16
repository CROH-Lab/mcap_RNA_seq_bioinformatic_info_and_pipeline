# ============================================================================
# Functional Enrichment Analysis - M. capitata DEGs
# ============================================================================
#
# Simplified version using Fisher's exact test (no clusterProfiler dependency)
#
# Input:
#   - DESeq2 results from 07_deseq2/
#   - eggNOG annotations from 04_reference/
#
# Output:
#   - GO enrichment results (BP, MF, CC)
#   - KEGG pathway enrichment results
#   - Publication-ready figures
#
# ============================================================================

# --- Load Libraries ---
cat("============================================\n")
cat("Functional Enrichment Analysis - M. capitata\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n\n")

library(tidyverse)

cat("Packages loaded successfully.\n\n")

# --- Configuration ---
BASE_DIR <- "/home/darmstrong4/mc_rework"
DESEQ_DIR <- file.path(BASE_DIR, "07_deseq2")
REF_DIR <- file.path(BASE_DIR, "04_reference")
OUT_DIR <- file.path(BASE_DIR, "08_functional_enrichment")

# Create output directories
dir.create(file.path(OUT_DIR, "results"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)

# Enrichment parameters
PADJ_CUTOFF <- 0.05      # DEG significance threshold
ENRICH_PVAL <- 0.05      # Enrichment p-value cutoff
MIN_GENES <- 3           # Minimum genes per term

# ============================================================================
# SECTION 1: Parse eggNOG Annotations
# ============================================================================

cat("=== SECTION 1: Parsing eggNOG Annotations ===\n\n")

# Read eggNOG results
eggnog_file <- file.path(REF_DIR, "Mcap_eggnog_results.txt")
cat("Reading:", eggnog_file, "\n")

eggnog <- read_tsv(eggnog_file, comment = "##", show_col_types = FALSE) %>%
    rename(gene_id = `#query`)

cat("Total genes with eggNOG annotations:", nrow(eggnog), "\n\n")

# --- Create Gene-to-GO mapping ---
cat("Creating gene-to-GO mapping...\n")

gene2go <- eggnog %>%
    filter(!is.na(GOs) & GOs != "-" & GOs != "") %>%
    select(gene_id, GOs) %>%
    separate_rows(GOs, sep = ",") %>%
    rename(GO = GOs) %>%
    filter(GO != "")

cat("  Genes with GO annotations:", n_distinct(gene2go$gene_id), "\n")
cat("  Total gene-GO associations:", nrow(gene2go), "\n")
cat("  Unique GO terms:", n_distinct(gene2go$GO), "\n\n")

# --- Create Gene-to-KEGG KO mapping ---
cat("Creating gene-to-KEGG mapping...\n")

gene2kegg <- eggnog %>%
    filter(!is.na(KEGG_ko) & KEGG_ko != "-" & KEGG_ko != "") %>%
    select(gene_id, KEGG_ko) %>%
    separate_rows(KEGG_ko, sep = ",") %>%
    mutate(KEGG_ko = gsub("ko:", "", KEGG_ko)) %>%
    rename(KEGG = KEGG_ko) %>%
    filter(KEGG != "")

cat("  Genes with KEGG annotations:", n_distinct(gene2kegg$gene_id), "\n")
cat("  Total gene-KEGG associations:", nrow(gene2kegg), "\n")
cat("  Unique KEGG orthologs:", n_distinct(gene2kegg$KEGG), "\n\n")

# --- Create Gene-to-Pathway mapping ---
cat("Creating gene-to-Pathway mapping...\n")

gene2pathway <- eggnog %>%
    filter(!is.na(KEGG_Pathway) & KEGG_Pathway != "-" & KEGG_Pathway != "") %>%
    select(gene_id, KEGG_Pathway) %>%
    separate_rows(KEGG_Pathway, sep = ",") %>%
    filter(grepl("^ko", KEGG_Pathway)) %>%
    rename(Pathway = KEGG_Pathway) %>%
    filter(Pathway != "")

cat("  Genes with Pathway annotations:", n_distinct(gene2pathway$gene_id), "\n")
cat("  Total gene-Pathway associations:", nrow(gene2pathway), "\n")
cat("  Unique KEGG pathways:", n_distinct(gene2pathway$Pathway), "\n\n")

# --- Store gene descriptions ---
gene2name <- eggnog %>%
    select(gene_id, Description, Preferred_name) %>%
    mutate(
        Description = ifelse(is.na(Description) | Description == "-", NA, Description),
        Preferred_name = ifelse(is.na(Preferred_name) | Preferred_name == "-", NA, Preferred_name)
    )

# ============================================================================
# SECTION 2: Load DESeq2 Results
# ============================================================================

cat("=== SECTION 2: Loading DESeq2 Results ===\n\n")

# Load M. capitata DESeq2 results
load(file.path(DESEQ_DIR, "Mcap_DESeq2_complete.RData"))

# Extract DEG lists
extract_degs <- function(res, padj_cut = 0.05) {
    res_df <- as.data.frame(res) %>%
        rownames_to_column("gene_id") %>%
        filter(!is.na(padj))
    
    list(
        all_genes = res_df$gene_id,
        all_tested = res_df %>% filter(!is.na(pvalue)) %>% pull(gene_id),
        degs = res_df %>% filter(padj < padj_cut) %>% pull(gene_id),
        degs_up = res_df %>% filter(padj < padj_cut & log2FoldChange > 0) %>% pull(gene_id),
        degs_down = res_df %>% filter(padj < padj_cut & log2FoldChange < 0) %>% pull(gene_id),
        res_df = res_df
    )
}

mcap_summer_degs <- extract_degs(mcap_summer$res_05, PADJ_CUTOFF)
mcap_winter_degs <- extract_degs(mcap_winter$res_05, PADJ_CUTOFF)

cat("M. capitata Summer:\n")
cat("  Total tested:", length(mcap_summer_degs$all_tested), "\n")
cat("  DEGs (padj <", PADJ_CUTOFF, "):", length(mcap_summer_degs$degs), "\n")
cat("    Up in high CO2:", length(mcap_summer_degs$degs_up), "\n")
cat("    Down in high CO2:", length(mcap_summer_degs$degs_down), "\n\n")

cat("M. capitata Winter:\n")
cat("  Total tested:", length(mcap_winter_degs$all_tested), "\n")
cat("  DEGs (padj <", PADJ_CUTOFF, "):", length(mcap_winter_degs$degs), "\n")
cat("    Up in high CO2:", length(mcap_winter_degs$degs_up), "\n")
cat("    Down in high CO2:", length(mcap_winter_degs$degs_down), "\n\n")

# ============================================================================
# SECTION 3: Fisher's Exact Test Enrichment Function
# ============================================================================

cat("=== SECTION 3: Running Enrichment Analysis ===\n\n")

# Fisher's exact test enrichment
run_enrichment <- function(degs, background, term2gene, label, min_genes = 3) {
    
    cat("Running enrichment for:", label, "\n")
    
    # Filter to annotated genes only
    degs_annotated <- degs[degs %in% term2gene$gene_id]
    background_annotated <- background[background %in% term2gene$gene_id]
    
    cat("  DEGs with annotations:", length(degs_annotated), "\n")
    cat("  Background with annotations:", length(background_annotated), "\n")
    
    if (length(degs_annotated) < min_genes) {
        cat("  WARNING: Too few annotated DEGs\n\n")
        return(NULL)
    }
    
    # Get all terms
    all_terms <- unique(term2gene[[1]])  # Second column is term
    term_col <- names(term2gene)[1]
    
    # Calculate enrichment for each term
    results <- data.frame()
    
    for (term in all_terms) {
        # Genes in this term
        genes_in_term <- term2gene %>% 
            filter(!!sym(term_col) == term) %>% 
            pull(gene_id) %>% 
            unique()
        
        # Filter to background
        genes_in_term <- genes_in_term[genes_in_term %in% background_annotated]
        
        if (length(genes_in_term) < min_genes) next
        
        # 2x2 contingency table
        # DEG & in term | DEG & not in term
        # Not DEG & in term | Not DEG & not in term
        
        deg_in_term <- sum(degs_annotated %in% genes_in_term)
        deg_not_in_term <- length(degs_annotated) - deg_in_term
        notdeg_in_term <- sum(genes_in_term %in% background_annotated) - deg_in_term
        notdeg_not_in_term <- length(background_annotated) - deg_in_term - deg_not_in_term - notdeg_in_term
        
        if (deg_in_term < min_genes) next
        
        # Fisher's exact test (one-sided, greater)
        mat <- matrix(c(deg_in_term, deg_not_in_term, notdeg_in_term, notdeg_not_in_term), 
                      nrow = 2, byrow = TRUE)
        
        test <- fisher.test(mat, alternative = "greater")
        
        # Gene ratio
        gene_ratio <- deg_in_term / length(degs_annotated)
        bg_ratio <- length(genes_in_term) / length(background_annotated)
        
        # Get DEG gene IDs in this term
        deg_genes <- degs_annotated[degs_annotated %in% genes_in_term]
        
        results <- rbind(results, data.frame(
            ID = term,
            Description = term,
            GeneRatio = paste0(deg_in_term, "/", length(degs_annotated)),
            BgRatio = paste0(length(genes_in_term), "/", length(background_annotated)),
            pvalue = test$p.value,
            OddsRatio = test$estimate,
            Count = deg_in_term,
            TermSize = length(genes_in_term),
            GeneIDs = paste(deg_genes, collapse = "/"),
            stringsAsFactors = FALSE
        ))
    }
    
    if (nrow(results) == 0) {
        cat("  No terms passed filters\n\n")
        return(NULL)
    }
    
    # Multiple testing correction
    results$p.adjust <- p.adjust(results$pvalue, method = "BH")
    
    # Sort by p-value
    results <- results %>% arrange(pvalue)
    
    n_sig <- sum(results$p.adjust < ENRICH_PVAL)
    cat("  Enriched terms (padj <", ENRICH_PVAL, "):", n_sig, "\n\n")
    
    return(results)
}

# ============================================================================
# SECTION 4: Run GO Enrichment
# ============================================================================

cat("=== SECTION 4: GO Enrichment ===\n\n")

# Prepare term2gene format (term, gene)
go_term2gene <- gene2go %>% select(GO, gene_id)

# Summer
go_summer_all <- run_enrichment(mcap_summer_degs$degs, mcap_summer_degs$all_tested, 
                                 go_term2gene, "Summer - All DEGs")
go_summer_up <- run_enrichment(mcap_summer_degs$degs_up, mcap_summer_degs$all_tested, 
                                go_term2gene, "Summer - Upregulated")
go_summer_down <- run_enrichment(mcap_summer_degs$degs_down, mcap_summer_degs$all_tested, 
                                  go_term2gene, "Summer - Downregulated")

# Winter
go_winter_all <- run_enrichment(mcap_winter_degs$degs, mcap_winter_degs$all_tested, 
                                 go_term2gene, "Winter - All DEGs")
go_winter_up <- run_enrichment(mcap_winter_degs$degs_up, mcap_winter_degs$all_tested, 
                                go_term2gene, "Winter - Upregulated")
go_winter_down <- run_enrichment(mcap_winter_degs$degs_down, mcap_winter_degs$all_tested, 
                                  go_term2gene, "Winter - Downregulated")

# ============================================================================
# SECTION 5: Run KEGG Pathway Enrichment
# ============================================================================

cat("=== SECTION 5: KEGG Pathway Enrichment ===\n\n")

# Prepare term2gene format
pathway_term2gene <- gene2pathway %>% select(Pathway, gene_id)

# Summer
kegg_summer_all <- run_enrichment(mcap_summer_degs$degs, mcap_summer_degs$all_tested, 
                                   pathway_term2gene, "Summer - All DEGs")
kegg_summer_up <- run_enrichment(mcap_summer_degs$degs_up, mcap_summer_degs$all_tested, 
                                  pathway_term2gene, "Summer - Upregulated")
kegg_summer_down <- run_enrichment(mcap_summer_degs$degs_down, mcap_summer_degs$all_tested, 
                                    pathway_term2gene, "Summer - Downregulated")

# Winter
kegg_winter_all <- run_enrichment(mcap_winter_degs$degs, mcap_winter_degs$all_tested, 
                                   pathway_term2gene, "Winter - All DEGs")
kegg_winter_up <- run_enrichment(mcap_winter_degs$degs_up, mcap_winter_degs$all_tested, 
                                  pathway_term2gene, "Winter - Upregulated")
kegg_winter_down <- run_enrichment(mcap_winter_degs$degs_down, mcap_winter_degs$all_tested, 
                                    pathway_term2gene, "Winter - Downregulated")

# ============================================================================
# SECTION 6: Add GO Term Names (from GO.db if available)
# ============================================================================

cat("=== SECTION 6: Adding GO Term Descriptions ===\n\n")

add_go_descriptions <- function(results) {
    if (is.null(results)) return(NULL)
    
    # Try to use GO.db
    if (requireNamespace("GO.db", quietly = TRUE)) {
        library(GO.db)
        
        go_terms <- results$ID
        term_info <- AnnotationDbi::select(GO.db, keys = go_terms, 
                                           columns = c("TERM", "ONTOLOGY"), 
                                           keytype = "GOID")
        
        results <- results %>%
            left_join(term_info, by = c("ID" = "GOID")) %>%
            mutate(Description = ifelse(!is.na(TERM), TERM, ID))
    }
    
    return(results)
}

go_summer_all <- add_go_descriptions(go_summer_all)
go_summer_up <- add_go_descriptions(go_summer_up)
go_summer_down <- add_go_descriptions(go_summer_down)
go_winter_all <- add_go_descriptions(go_winter_all)
go_winter_up <- add_go_descriptions(go_winter_up)
go_winter_down <- add_go_descriptions(go_winter_down)

cat("GO descriptions added.\n\n")

# ============================================================================
# SECTION 7: Save Results Tables
# ============================================================================

cat("=== SECTION 7: Saving Results Tables ===\n\n")

save_results <- function(results, filename) {
    if (!is.null(results) && nrow(results) > 0) {
        sig_results <- results %>% filter(p.adjust < ENRICH_PVAL)
        if (nrow(sig_results) > 0) {
            write_csv(sig_results, file.path(OUT_DIR, "tables", filename))
            cat("Saved:", filename, "(", nrow(sig_results), "terms )\n")
            return(sig_results)
        }
    }
    cat("No significant results for:", filename, "\n")
    return(NULL)
}

# GO results
go_summer_all_sig <- save_results(go_summer_all, "GO_Summer_AllDEGs.csv")
go_summer_up_sig <- save_results(go_summer_up, "GO_Summer_Upregulated.csv")
go_summer_down_sig <- save_results(go_summer_down, "GO_Summer_Downregulated.csv")
go_winter_all_sig <- save_results(go_winter_all, "GO_Winter_AllDEGs.csv")
go_winter_up_sig <- save_results(go_winter_up, "GO_Winter_Upregulated.csv")
go_winter_down_sig <- save_results(go_winter_down, "GO_Winter_Downregulated.csv")

# KEGG results
kegg_summer_all_sig <- save_results(kegg_summer_all, "KEGG_Summer_AllDEGs.csv")
kegg_summer_up_sig <- save_results(kegg_summer_up, "KEGG_Summer_Upregulated.csv")
kegg_summer_down_sig <- save_results(kegg_summer_down, "KEGG_Summer_Downregulated.csv")
kegg_winter_all_sig <- save_results(kegg_winter_all, "KEGG_Winter_AllDEGs.csv")
kegg_winter_up_sig <- save_results(kegg_winter_up, "KEGG_Winter_Upregulated.csv")
kegg_winter_down_sig <- save_results(kegg_winter_down, "KEGG_Winter_Downregulated.csv")

# Save full results (all terms, not just significant)
cat("\nSaving full results (all tested terms)...\n")
if (!is.null(go_summer_all)) write_csv(go_summer_all, file.path(OUT_DIR, "results", "GO_Summer_AllDEGs_full.csv"))
if (!is.null(go_summer_up)) write_csv(go_summer_up, file.path(OUT_DIR, "results", "GO_Summer_Up_full.csv"))
if (!is.null(go_summer_down)) write_csv(go_summer_down, file.path(OUT_DIR, "results", "GO_Summer_Down_full.csv"))
if (!is.null(go_winter_all)) write_csv(go_winter_all, file.path(OUT_DIR, "results", "GO_Winter_AllDEGs_full.csv"))
if (!is.null(go_winter_up)) write_csv(go_winter_up, file.path(OUT_DIR, "results", "GO_Winter_Up_full.csv"))
if (!is.null(go_winter_down)) write_csv(go_winter_down, file.path(OUT_DIR, "results", "GO_Winter_Down_full.csv"))

if (!is.null(kegg_summer_all)) write_csv(kegg_summer_all, file.path(OUT_DIR, "results", "KEGG_Summer_AllDEGs_full.csv"))
if (!is.null(kegg_winter_all)) write_csv(kegg_winter_all, file.path(OUT_DIR, "results", "KEGG_Winter_AllDEGs_full.csv"))

cat("\n")

# ============================================================================
# SECTION 8: Generate Figures
# ============================================================================

cat("=== SECTION 8: Generating Figures ===\n\n")

# Dotplot function
make_dotplot <- function(results, title, filename, n_show = 20) {
    if (is.null(results) || nrow(results) == 0) {
        cat("No results for:", title, "\n")
        return(NULL)
    }
    
    sig_results <- results %>% filter(p.adjust < ENRICH_PVAL)
    
    if (nrow(sig_results) == 0) {
        cat("No significant results for:", title, "\n")
        return(NULL)
    }
    
    # Prepare data
    plot_data <- sig_results %>%
        head(n_show) %>%
        mutate(
            GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
            Description = ifelse(nchar(Description) > 50, 
                                 paste0(substr(Description, 1, 47), "..."), 
                                 Description),
            Description = factor(Description, levels = rev(Description))
        )
    
    p <- ggplot(plot_data, aes(x = GeneRatio_num, y = Description)) +
        geom_point(aes(size = Count, color = -log10(p.adjust))) +
        scale_color_gradient(low = "blue", high = "red", name = "-log10(padj)") +
        scale_size_continuous(name = "Gene Count", range = c(2, 8)) +
        labs(
            title = title,
            x = "Gene Ratio",
            y = NULL
        ) +
        theme_bw(base_size = 10) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
            axis.text.y = element_text(size = 8)
        )
    
    ggsave(file.path(OUT_DIR, "figures", filename), p, 
           width = 8, height = min(4 + nrow(plot_data) * 0.25, 12))
    cat("Saved:", filename, "\n")
    
    return(p)
}

# Barplot function
make_barplot <- function(results, title, filename, n_show = 15) {
    if (is.null(results) || nrow(results) == 0) {
        cat("No results for:", title, "\n")
        return(NULL)
    }
    
    sig_results <- results %>% filter(p.adjust < ENRICH_PVAL)
    
    if (nrow(sig_results) == 0) {
        cat("No significant results for:", title, "\n")
        return(NULL)
    }
    
    plot_data <- sig_results %>%
        head(n_show) %>%
        mutate(
            Description = ifelse(nchar(Description) > 40, 
                                 paste0(substr(Description, 1, 37), "..."), 
                                 Description),
            Description = factor(Description, levels = rev(Description))
        )
    
    p <- ggplot(plot_data, aes(x = Count, y = Description, fill = -log10(p.adjust))) +
        geom_bar(stat = "identity") +
        scale_fill_gradient(low = "steelblue", high = "darkred", name = "-log10(padj)") +
        labs(
            title = title,
            x = "Gene Count",
            y = NULL
        ) +
        theme_bw(base_size = 10) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
            axis.text.y = element_text(size = 8)
        )
    
    ggsave(file.path(OUT_DIR, "figures", filename), p,
           width = 8, height = min(4 + nrow(plot_data) * 0.3, 10))
    cat("Saved:", filename, "\n")
    
    return(p)
}

# Generate GO plots
make_dotplot(go_summer_all, "GO Enrichment - Summer All DEGs", "GO_dotplot_Summer_All.pdf")
make_dotplot(go_summer_up, "GO Enrichment - Summer Upregulated", "GO_dotplot_Summer_Up.pdf")
make_dotplot(go_summer_down, "GO Enrichment - Summer Downregulated", "GO_dotplot_Summer_Down.pdf")
make_dotplot(go_winter_all, "GO Enrichment - Winter All DEGs", "GO_dotplot_Winter_All.pdf")
make_dotplot(go_winter_up, "GO Enrichment - Winter Upregulated", "GO_dotplot_Winter_Up.pdf")
make_dotplot(go_winter_down, "GO Enrichment - Winter Downregulated", "GO_dotplot_Winter_Down.pdf")

# Generate KEGG plots
make_barplot(kegg_summer_all, "KEGG Pathways - Summer All DEGs", "KEGG_barplot_Summer_All.pdf")
make_barplot(kegg_summer_up, "KEGG Pathways - Summer Upregulated", "KEGG_barplot_Summer_Up.pdf")
make_barplot(kegg_summer_down, "KEGG Pathways - Summer Downregulated", "KEGG_barplot_Summer_Down.pdf")
make_barplot(kegg_winter_all, "KEGG Pathways - Winter All DEGs", "KEGG_barplot_Winter_All.pdf")
make_barplot(kegg_winter_up, "KEGG Pathways - Winter Upregulated", "KEGG_barplot_Winter_Up.pdf")
make_barplot(kegg_winter_down, "KEGG Pathways - Winter Downregulated", "KEGG_barplot_Winter_Down.pdf")

cat("\n")

# ============================================================================
# SECTION 9: Summary Statistics
# ============================================================================

cat("=== SECTION 9: Summary Statistics ===\n\n")

count_enriched <- function(results, label) {
    if (is.null(results)) {
        return(data.frame(Analysis = label, Enriched_Terms = 0, Genes_in_Terms = 0))
    }
    
    sig <- results %>% filter(p.adjust < ENRICH_PVAL)
    n_terms <- nrow(sig)
    
    if (n_terms > 0) {
        n_genes <- length(unique(unlist(strsplit(sig$GeneIDs, "/"))))
    } else {
        n_genes <- 0
    }
    
    data.frame(Analysis = label, Enriched_Terms = n_terms, Genes_in_Terms = n_genes)
}

summary_go <- bind_rows(
    count_enriched(go_summer_all, "Summer - All DEGs"),
    count_enriched(go_summer_up, "Summer - Upregulated"),
    count_enriched(go_summer_down, "Summer - Downregulated"),
    count_enriched(go_winter_all, "Winter - All DEGs"),
    count_enriched(go_winter_up, "Winter - Upregulated"),
    count_enriched(go_winter_down, "Winter - Downregulated")
)

summary_kegg <- bind_rows(
    count_enriched(kegg_summer_all, "Summer - All DEGs"),
    count_enriched(kegg_summer_up, "Summer - Upregulated"),
    count_enriched(kegg_summer_down, "Summer - Downregulated"),
    count_enriched(kegg_winter_all, "Winter - All DEGs"),
    count_enriched(kegg_winter_up, "Winter - Upregulated"),
    count_enriched(kegg_winter_down, "Winter - Downregulated")
)

cat("GO Enrichment Summary:\n")
print(summary_go, row.names = FALSE)
cat("\n")

cat("KEGG Pathway Enrichment Summary:\n")
print(summary_kegg, row.names = FALSE)
cat("\n")

# Save summaries
write_csv(summary_go, file.path(OUT_DIR, "tables", "GO_enrichment_summary.csv"))
write_csv(summary_kegg, file.path(OUT_DIR, "tables", "KEGG_enrichment_summary.csv"))

# ============================================================================
# SECTION 10: Top Enriched Terms Preview
# ============================================================================

cat("=== SECTION 10: Top Enriched Terms ===\n\n")

print_top_terms <- function(results, label, n = 10) {
    if (is.null(results)) {
        cat(label, ": No results\n\n")
        return()
    }
    
    sig <- results %>% filter(p.adjust < ENRICH_PVAL)
    
    if (nrow(sig) == 0) {
        cat(label, ": No significant terms\n\n")
        return()
    }
    
    cat(label, " (", nrow(sig), " enriched terms):\n", sep = "")
    
    top <- sig %>%
        head(n) %>%
        select(ID, Description, Count, pvalue, p.adjust)
    
    print(top, row.names = FALSE)
    cat("\n")
}

print_top_terms(go_summer_all, "GO - Summer All DEGs")
print_top_terms(go_winter_all, "GO - Winter All DEGs")
print_top_terms(kegg_summer_all, "KEGG - Summer All DEGs")
print_top_terms(kegg_winter_all, "KEGG - Winter All DEGs")

# ============================================================================
# SECTION 11: Save R Objects
# ============================================================================

cat("=== SECTION 11: Saving R Objects ===\n\n")

enrichment_results <- list(
    # GO results
    go_summer_all = go_summer_all,
    go_summer_up = go_summer_up,
    go_summer_down = go_summer_down,
    go_winter_all = go_winter_all,
    go_winter_up = go_winter_up,
    go_winter_down = go_winter_down,
    # KEGG results
    kegg_summer_all = kegg_summer_all,
    kegg_summer_up = kegg_summer_up,
    kegg_summer_down = kegg_summer_down,
    kegg_winter_all = kegg_winter_all,
    kegg_winter_up = kegg_winter_up,
    kegg_winter_down = kegg_winter_down,
    # Annotation mappings
    gene2go = gene2go,
    gene2kegg = gene2kegg,
    gene2pathway = gene2pathway,
    gene2name = gene2name,
    # DEG lists
    mcap_summer_degs = mcap_summer_degs,
    mcap_winter_degs = mcap_winter_degs
)

save(enrichment_results, file = file.path(OUT_DIR, "Mcap_enrichment_results.RData"))
cat("Saved: Mcap_enrichment_results.RData\n\n")

# ============================================================================
# SECTION 12: Session Info
# ============================================================================

cat("=== Session Info ===\n\n")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "session_info.txt"))
cat("Saved: session_info.txt\n")

cat("\n============================================\n")
cat("Analysis Complete!\n")
cat("============================================\n\n")

cat("Output directory:", OUT_DIR, "\n\n")

cat("Tables (significant terms only):\n")
cat("  tables/GO_*.csv           - GO enrichment results\n")
cat("  tables/KEGG_*.csv         - KEGG pathway results\n")
cat("  tables/*_summary.csv      - Summary statistics\n\n")

cat("Full results (all tested terms):\n")
cat("  results/*_full.csv        - Complete results for filtering\n\n")

cat("Figures:\n")
cat("  figures/GO_dotplot_*.pdf  - GO enrichment dotplots\n")
cat("  figures/KEGG_barplot_*.pdf - KEGG pathway barplots\n\n")

cat("R Objects:\n")
cat("  Mcap_enrichment_results.RData\n")
