# Files Manifest

Current state after subsetting to treatments B (OA) and D (control), 12 samples.

## What's in this Repository

### Included (version controlled)

- `README.md` — Project documentation
- `FILES_MANIFEST.md` — This file
- `scripts/` — All pipeline scripts
- `sample_info.txt` — Sample metadata (12 B/D samples)

- `05_star_align/` — Alignment logs only (not BAMs)
  - `*/Log.final.out` — STAR alignment statistics
  - `alignment_summary.tsv` — Combined alignment stats
  - `organism_distribution.tsv` — Per-sample organism proportions

- `06_featurecounts/` — Count matrices
  - `Mcap_counts.txt` — Host counts (no multi-mappers)
  - `Cgor_counts.txt` — Clade C counts (with multi-mappers)
  - `Dtre_counts.txt` — Clade D counts (with multi-mappers)
  - `organism_proportions.txt` — Summary statistics
  - `*_assignment_summary.txt` — Assignment rates per organism

- `07_deseq2/` — Differential expression analysis
  - `results/Mcap_Summer_BvsD.csv` — Host summer DEG results
  - `results/Mcap_Winter_BvsD.csv` — Host winter DEG results
  - `results/Dtre_Summer_BvsD.csv` — Symbiont summer DEG results
  - `results/Dtre_Winter_BvsD.csv` — Symbiont winter DEG results
  - `results/DEG_summary.csv` — Summary statistics
  - `objects/*.rds` — DESeq2 and VST R objects
  - `figures/*.pdf` — QC plots (PCA, dispersion, p-value histograms)

- `08_host_deg_annotation/` — Host DEG functional annotation
  - `sequences/all_degs.txt` — List of 972 host DEG IDs
  - `sequences/degs_cds.fa` — CDS sequences for DEGs
  - `results/sprot_annotations_full.tsv` — SwissProt annotations
  - `results/trembl_annotations_full.tsv` — TrEMBL annotations
  - `results/all_annotations_full.tsv` — Combined annotations (905 genes)
  - `results/unannotated_degs.txt` — DEGs with no BLAST hits (67 genes)

- `09_symbiont_deg_annotation/` — Symbiont DEG functional annotation
  - `sequences/all_degs.txt` — List of 281 symbiont DEG IDs
  - `sequences/degs_cds.fa` — CDS sequences for DEGs
  - `results/sprot_annotations_full.tsv` — SwissProt annotations
  - `results/trembl_annotations_full.tsv` — TrEMBL annotations
  - `results/all_annotations_full.tsv` — Combined annotations (163 genes)
  - `results/unannotated_degs.txt` — DEGs with no BLAST hits (118 genes)

- `10_host_GO_MWU/` — Host GO enrichment analysis
  - `input/go_annotations.tab` — Gene-to-GO mapping
  - `input/*_signed_logP.csv` — GO_MWU input files
  - `output/*_MWU_results.csv` — Full enrichment results
  - `output/*_representative_GOs.csv` — Clustered representative terms
  - `output/GO_MWU_summary.csv` — Summary statistics
  - `figures/*_GO_MWU.pdf` — Enrichment plots
  - `scripts/10_host_GO_MWU_analysis.R` — Analysis script
  - `gomwu.functions.R` — GO_MWU R functions
  - `gomwu_a.pl`, `gomwu_b.pl` — GO_MWU Perl scripts

- `11_symbiont_GO_MWU/` — Symbiont GO enrichment analysis
  - `input/symbiont_go_annotations.tab` — Gene-to-GO mapping
  - `input/*_signed_logP.csv` — GO_MWU input files
  - `output/symbiont_*_MWU_results.csv` — Full enrichment results
  - `output/symbiont_*_representative_GOs.csv` — Clustered representative terms
  - `output/symbiont_GO_MWU_summary.csv` — Summary statistics
  - `figures/symbiont_*_GO_MWU.pdf` — Enrichment plots
  - `scripts/11_symbiont_GO_MWU_analysis.R` — Analysis script

- `12_publication_figures/` — Publication-ready figures and tables
  - `figures/Fig2_combined_venn_ridgeline.pdf/png` — Venn + ridgeline plots
  - `figures/Fig3_heatmap_combined.pdf/png` — Combined DEG heatmaps
  - `figures/Fig4_pca_holobiont_combined.pdf/png` — Holobiont PCA
  - `figures/Fig5_Summer_Representative_GO_Dendrogram.pdf` — Summer GO dendrogram
  - `figures/Fig6_Winter_Representative_GO_Dendrogram.pdf` — Winter GO dendrogram
  - `figures/Fig7_Shared_GO_Chord_Diagram.pdf` — Shared GO chord diagram
  - `figures/FigS1_pca_individual_combined.pdf/png` — Individual PCA grid
  - `tables/Table_S_*_Shared_GO_Terms.html/docx/csv` — Supplementary GO tables
  - `data/DEG_summary_statistics.csv` — DEG counts summary
  - `scripts/12_publication_figures.R` — Unified figure generation script

- `13_cal_chem/` — Carbonate chemistry and calcification analysis
  - `input/sw_chem.csv` — Seawater chemistry data
  - `input/calc_data.csv` — Calcification measurements
  - `output/carbonate_chemistry_table.docx` — Chemistry summary table
  - `output/carbonate_chemistry_summary.csv` — Chemistry summary CSV
  - `output/pH_anova_results.txt` — pH ANOVA results
  - `output/calcification_anova_results.txt` — Calcification ANOVA results
  - `figures/calcification_barplot.pdf/png` — Fig 1: Calcification barplot
  - `scripts/13_cal_chem_analysis.R` — Analysis script

- `logs/` — SLURM job output and error logs

### Excluded (too large for GitHub, in .gitignore)

- `00_raw_fastq/` — Raw FASTQ files (~400 GB)
- `01_qc_raw/` — FastQC reports on raw reads
- `02_trimmed/` — Trimmed FASTQ files (~400 GB)
- `03_qc_trimmed/` — FastQC reports on trimmed reads
- `04_reference/` — Reference genomes and STAR index (~60 GB)
- `05_star_align/*/*.bam` — BAM files (~150 GB)
- `06_featurecounts/*_counts_raw.txt` — Raw featureCounts output with metadata
- `10_host_GO_MWU/input/go.obo` — GO database (31 MB)
- `10_host_GO_MWU/output/*_dissimilarity.csv` — Large clustering matrices
- `11_symbiont_GO_MWU/input/go.obo` — GO database
- `11_symbiont_GO_MWU/output/*_dissimilarity.csv` — Large clustering matrices

### Data Availability

Large files stored on TAMU HPC: `/home/darmstrong4/mc_rework/`

Reference genomes:
- *M. capitata*: GigaDB 102268 (Stephens et al. 2022)
- *C. goreaui*: UQ eSpace (Chen et al. 2022)
- *D. trenchii*: UQ eSpace (Dougan et al. 2022)

GO database: http://purl.obolibrary.org/obo/go.obo (release 2025-10-10)
