# Files Manifest

## What's in this Repository

### Included (version controlled)
- `README.md` - Project documentation
- `scripts/` - All pipeline scripts
- `metadata/` - Sample information
- `06_featurecounts/` - Count matrices (except raw counts)
  - `gene_counts.txt` - Integer counts for DESeq2
  - `gene_counts_fractional.txt` - Fractional counts
  - `gene_counts_Mcap.txt` - Host-only counts
  - `gene_counts_Cgor.txt` - Clade C counts
  - `gene_counts_Dtre.txt` - Clade D counts
  - `organism_proportions.txt` - Summary statistics
  - `assignment_summary.txt` - Assignment rates
- `05_star_align/` - Alignment logs only (not BAMs)
- `.gitignore` - Git ignore rules

### Excluded (too large for GitHub)
- `00_raw_fastq/` - Raw FASTQ files (~400GB)
- `02_trimmed/` - Trimmed FASTQ files (~400GB)
- `04_reference/` - Reference genomes and indices (~60GB)
- `05_star_align/*/*.bam` - BAM files (~150GB)
- `06_featurecounts/gene_counts_raw.txt` - Raw featureCounts output

### Data Availability
Large files are stored on TAMU HPC at:
`/home/darmstrong4/mc_rework/`

Reference genomes can be downloaded from:
- M. capitata: GigaDB 102268 (Stephens et al. 2022)
- C. goreaui: UQ eSpace (Chen et al. 2022)
- D. trenchii: UQ eSpace (Dougan et al. 2022)
