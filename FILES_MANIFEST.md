# Files Manifest

## What's in this Repository

### Included (version controlled)

- `README.md` - Project documentation
- `FILES_MANIFEST.md` - This file
- `scripts/` - All pipeline scripts
- `sample_info.txt` - Sample metadata
- `sample_info_with_chemistry.txt` - Sample metadata with carbonate chemistry

- `05_star_align/` - Alignment logs only (not BAMs)
  - `*/Log.final.out` - STAR alignment statistics
  - `alignment_summary.tsv` - Combined alignment stats
  - `organism_distribution.tsv` - Per-sample organism proportions

- `06_featurecounts/` - Count matrices
  - `Mcap_counts.txt` - Host counts (no multi-mappers)
  - `Cgor_counts.txt` - Clade C counts (with multi-mappers)
  - `Dtre_counts.txt` - Clade D counts (with multi-mappers)
  - `organism_proportions.txt` - Summary statistics
  - `*_assignment_summary.txt` - Assignment rates per organism

- `07_deseq2/` - Differential expression analysis
  - `results/Mcap_Summer_BvsD.csv` - Host summer DEG results
  - `results/Mcap_Winter_BvsD.csv` - Host winter DEG results
  - `results/Dtre_Summer_BvsD.csv` - Symbiont summer DEG results
  - `results/Dtre_Winter_BvsD.csv` - Symbiont winter DEG results
  - `results/DEG_summary.csv` - Summary statistics
  - `objects/*.rds` - DESeq2 and VST R objects
  - `figures/*.pdf` - QC plots (PCA, dispersion, p-value histograms)

- `08_host_deg_annotation/` - Host DEG functional annotation
  - `sequences/all_degs.txt` - List of 972 host DEG IDs
  - `sequences/degs_cds.fa` - CDS sequences for DEGs
  - `results/sprot_annotations_full.tsv` - SwissProt annotations
  - `results/trembl_annotations_full.tsv` - TrEMBL annotations
  - `results/all_annotations_full.tsv` - Combined annotations (905 genes)
  - `results/unannotated_degs.txt` - DEGs with no BLAST hits (67 genes)

- `09_symbiont_deg_annotation/` - Symbiont DEG functional annotation
  - `sequences/all_degs.txt` - List of 281 symbiont DEG IDs
  - `sequences/degs_cds.fa` - CDS sequences for DEGs
  - `results/sprot_annotations_full.tsv` - SwissProt annotations
  - `results/trembl_annotations_full.tsv` - TrEMBL annotations
  - `results/all_annotations_full.tsv` - Combined annotations (163 genes)
  - `results/unannotated_degs.txt` - DEGs with no BLAST hits (118 genes)

- `10_GO_MWU/` - Host GO enrichment analysis
  - `input/go_annotations.tab` - Gene-to-GO mapping
  - `input/*_signed_logP.csv` - GO_MWU input files
  - `output/*_MWU_results.csv` - Full enrichment results
  - `output/GO_MWU_summary.csv` - Summary statistics
  - `figures/*_GO_MWU.pdf` - Enrichment plots
  - `scripts/10_GO_MWU_analysis.R` - Analysis script
  - `gomwu.functions.R` - GO_MWU R functions
  - `gomwu_a.pl`, `gomwu_b.pl` - GO_MWU Perl scripts

- `11_symbiont_GO_MWU/` - Symbiont GO enrichment analysis
  - `input/symbiont_go_annotations.tab` - Gene-to-GO mapping
  - `input/*_signed_logP.csv` - GO_MWU input files
  - `output/*_MWU_results.csv` - Full enrichment results
  - `output/symbiont_GO_MWU_summary.csv` - Summary statistics
  - `figures/*_GO_MWU.pdf` - Enrichment plots
  - `scripts/11_symbiont_GO_MWU_analysis.R` - Analysis script

- `12_publication_figures/` - Publication-ready visualizations
  - `data/*.csv` - Summary data files
  - `figures/Fig1_ridgeline_*.pdf/png` - Expression distributions
  - `figures/Fig2*_venn.pdf` - Venn diagrams
  - `figures/Fig3*_heatmap_*.pdf` - Expression heatmaps
  - `figures/Fig4_*_sankey_*.html` - Interactive Sankey diagrams
  - `figures/Fig4D-E_bubble_*.pdf/png` - GO enrichment bubble plots
  - `figures/Fig5*_volcano_*.pdf/png` - Volcano plots
  - `figures/Fig6*_pca_*.pdf/png` - PCA plots
  - `scripts/12_publication_figures_*.R` - Figure generation scripts

- `logs/` - SLURM job logs

### Excluded (too large for GitHub)

- `00_raw_fastq/` - Raw FASTQ files (~400GB)
- `01_qc_raw/` - FastQC reports on raw reads
- `02_trimmed/` - Trimmed FASTQ files (~400GB)
- `03_qc_trimmed/` - FastQC reports on trimmed reads
- `04_reference/` - Reference genomes and indices (~60GB)
  - `holobiont_genome.fa` - Combined genome
  - `star_index_holobiont/` - STAR index (50GB)
- `05_star_align/*/*.bam` - BAM files (~150GB)
- `06_featurecounts/*_counts_raw.txt` - Raw featureCounts output with metadata
- `10_GO_MWU/input/go.obo` - GO database (31MB)
- `10_GO_MWU/output/*_dissimilarity.csv` - Large clustering matrices
- `11_symbiont_GO_MWU/input/go.obo` - GO database
- `11_symbiont_GO_MWU/output/*_dissimilarity.csv` - Large clustering matrices

### Data Availability

Large files are stored on TAMU HPC at:
`/home/darmstrong4/mc_rework/`

Reference genomes can be downloaded from:
- M. capitata: GigaDB 102268 (Stephens et al. 2022)
- C. goreaui: UQ eSpace (Chen et al. 2022)
- D. trenchii: UQ eSpace (Dougan et al. 2022)

GO database:
- go.obo: http://purl.obolibrary.org/obo/go.obo (release 2025-10-10)

### File Size Summary

| Directory | Included | Excluded |
|-----------|----------|----------|
| 00_raw_fastq | - | ~400 GB |
| 02_trimmed | - | ~400 GB |
| 04_reference | - | ~60 GB |
| 05_star_align | Logs (~1 MB) | BAMs (~150 GB) |
| 06_featurecounts | Counts (~5 MB) | Raw output (~50 MB) |
| 07_deseq2 | All (~20 MB) | - |
| 08_host_deg_annotation | All (~15 MB) | - |
| 09_symbiont_deg_annotation | All (~5 MB) | - |
| 10_GO_MWU | Results (~2 MB) | Dissimilarity (~15 MB) |
| 11_symbiont_GO_MWU | Results (~2 MB) | Dissimilarity (~160 MB) |
| 12_publication_figures | All (~50 MB) | - |
