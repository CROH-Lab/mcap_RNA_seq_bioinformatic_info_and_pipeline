# M. capitata Holobiont RNA-seq Pipeline

## Overview

RNA-seq analysis of *Montipora capitata* coral holobiont response to seawater acidification. This pipeline uses genome-based alignment with STAR to quantify gene expression for both the coral host and its algal symbionts (Symbiodiniaceae) simultaneously.

## Sequencing Quality Summary

| Metric | Value |
|--------|-------|
| **Samples** | 24 |
| **Total Reads** | 1.37 billion |
| **Total Yield** | 412.5 Gb |
| **Read Length** | 150 bp PE |
| **Mean Quality Score** | Q38.9 |
| **% Bases ≥ Q30** | 94.6% |
| **GC Content** | 45-50% |

### Per-Sample Raw Read Statistics

| Metric | Min | Max | Mean |
|--------|-----|-----|------|
| Reads (M) | 44.9 | 80.9 | 57.3 |
| Quality Score | Q38.84 | Q38.96 | Q38.92 |
| % ≥ Q30 | 94.2% | 94.8% | 94.6% |
| Duplication | 14.4% | 19.5% | 16.3% |

## Trimming Summary (Cutadapt)

High-quality reads were retained after adapter removal and quality filtering.

| Metric | Value |
|--------|-------|
| **Read Retention** | 97.5% - 98.4% |
| **Adapter Detection** | 7.9% - 9.7% |
| **Quality Trimmed** | 0.3% of bases |
| **Too Short (<50bp)** | ~2.2% |

### Per-Sample Trimmed Read Counts

| Sample | Trimmed Reads | Retention |
|--------|---------------|-----------|
| 1AS | 48,178,543 | 97.8% |
| 1AW | 63,516,651 | 97.9% |
| 1BS | 53,749,496 | 97.5% |
| 1BW | 55,594,613 | 97.7% |
| 1CS | 56,353,424 | 97.9% |
| 1CW | 49,333,570 | 98.2% |
| 1DS | 49,208,380 | 98.0% |
| 1DW | 58,708,738 | 97.9% |
| 2AS | 55,551,688 | 98.4% |
| 2AW | 59,451,422 | 98.0% |
| 2BS | 50,120,093 | 97.6% |
| 2BW | 43,974,226 | 98.0% |
| 2CS | 54,206,111 | 97.6% |
| 2CW | 61,035,727 | 98.1% |
| 2DS | 54,211,740 | 97.6% |
| 2DW | 55,691,279 | 97.9% |
| 3AS | 55,993,112 | 97.9% |
| 3AW | 55,773,556 | 98.0% |
| 3BS | 53,177,277 | 98.0% |
| 3BW | 54,895,066 | 97.6% |
| 3CS | 57,196,082 | 97.8% |
| 3CW | 64,117,197 | 98.0% |
| 3DS | 56,598,638 | 97.9% |
| 3DW | 78,933,530 | 97.6% |

### Trimming Parameters

```bash
cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \  # TruSeq R1 adapter
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \  # TruSeq R2 adapter
    -a "A{10}" \                             # Poly-A tail removal
    -A "T{10}" \                             # Poly-T tail removal
    -q 20,20 \                               # Quality threshold
    --minimum-length 50                      # Minimum read length
```

## STAR Alignment Summary

### Holobiont Genome Reference

Reads were aligned to a combined holobiont genome containing host and symbiont sequences:

| Organism | Prefix | Scaffolds | Role |
|----------|--------|-----------|------|
| *Montipora capitata* | Mcap_ | 1,699 | Host coral |
| *Cladocopium goreaui* | Cgor_ | 6,843 | Symbiont (Clade C) |
| *Durusdinium trenchii* | Dtre_ | 44,682 | Symbiont (Clade D) |
| **Total** | | **53,224** | Combined holobiont |

### STAR Index Statistics

| Metric | Value |
|--------|-------|
| **Index Size** | 50 GB |
| **Genome Size** | ~16.1 Gb |
| **Splice Junctions** | 1,663,303 |
| **Exons Indexed** | 1,991,029 |
| **sjdbOverhang** | 149 bp |

### Alignment Statistics

| Metric | Mean | Range |
|--------|------|-------|
| **Uniquely Mapped** | 68.0% | 61.0% - 74.3% |
| **Multi-Mapped** | 18.8% | 12.9% - 23.1% |
| **Total Mapped** | 86.8% | 84.3% - 89.0% |
| **Unmapped (too short)** | 12.9% | 10.9% - 15.6% |

### Per-Sample Alignment Results

| Sample | Total Reads | Unique % | Multi % | Total Mapped % |
|--------|-------------|----------|---------|----------------|
| 1AS | 48,178,543 | 60.96% | 23.11% | 84.26% |
| 1AW | 63,516,651 | 74.27% | 13.46% | 87.99% |
| 1BS | 53,749,496 | 74.08% | 14.60% | 88.91% |
| 1BW | 55,594,613 | 70.21% | 17.67% | 88.03% |
| 1CS | 56,353,424 | 68.13% | 18.55% | 86.83% |
| 1CW | 49,333,570 | 62.99% | 22.27% | 85.41% |
| 1DS | 49,208,380 | 62.73% | 22.47% | 85.33% |
| 1DW | 58,708,738 | 66.45% | 19.91% | 86.53% |
| 2AS | 55,551,688 | 71.73% | 12.89% | 85.10% |
| 2AW | 59,451,422 | 64.75% | 21.39% | 86.25% |
| 2BS | 50,120,093 | 70.84% | 17.29% | 88.29% |
| 2BW | 43,974,226 | 65.27% | 20.81% | 86.21% |
| 2CS | 54,206,111 | 70.49% | 17.36% | 88.01% |
| 2CW | 61,035,727 | 71.61% | 15.18% | 87.15% |
| 2DS | 54,211,740 | 70.05% | 17.72% | 87.96% |
| 2DW | 55,691,279 | 67.36% | 19.53% | 87.03% |
| 3AS | 55,993,112 | 65.18% | 21.49% | 86.83% |
| 3AW | 55,773,556 | 65.68% | 20.69% | 86.49% |
| 3BS | 53,177,277 | 63.98% | 22.14% | 86.25% |
| 3BW | 54,895,066 | 70.07% | 17.42% | 87.68% |
| 3CS | 57,196,082 | 68.37% | 18.88% | 87.40% |
| 3CW | 64,117,197 | 66.97% | 20.19% | 87.29% |
| 3DS | 56,598,638 | 67.98% | 19.53% | 87.72% |
| 3DW | 78,933,530 | 72.44% | 16.01% | 88.58% |

### STAR Alignment Parameters

```bash
STAR --runMode alignReads \
    --runThreadN 64 \
    --genomeDir star_index_holobiont \
    --readFilesIn ${R1} ${R2} \
    --readFilesCommand zcat \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 1

# BAM sorting performed separately with samtools for stability
samtools sort -@ 64 -m 4G -o sorted.bam unsorted.bam
```

### Symbiont Community Profiles

Three distinct profiles were identified:

| Profile | Samples (n) | Characteristics |
|---------|-------------|-----------------|
| **Cgor-dominant** | 3 | High Clade C (*Cladocopium*), >25% |
| **Mixed** | 2 | Detectable both Cgor & Dtre |
| **Dtre-dominant** | 19 | Almost exclusively Clade D (*Durusdinium*) |

### Per-Sample Organism Distribution

| Sample | Host (Mcap) % | Clade C (Cgor) % | Clade D (Dtre) % | Profile |
|--------|---------------|------------------|------------------|---------|
| 1AS | 36.7 | 4.8 | 46.2 | Mixed |
| 1AW | 53.5 | 26.0 | 10.6 | Cgor-dominant |
| 1BS | 63.5 | 8.8 | 18.5 | Mixed |
| 1BW | 59.5 | 0.0 | 30.6 | Dtre-dominant |
| 1CS | 54.7 | 0.0 | 34.6 | Dtre-dominant |
| 1CW | 43.3 | 0.1 | 45.1 | Dtre-dominant |
| 1DS | 42.3 | 0.0 | 46.0 | Dtre-dominant |
| 1DW | 51.1 | 0.0 | 38.0 | Dtre-dominant |
| 2AS | 26.8 | 61.8 | 0.0 | Cgor-dominant |
| 2AW | 46.2 | 0.0 | 42.8 | Dtre-dominant |
| 2BS | 60.4 | 0.0 | 29.9 | Dtre-dominant |
| 2BW | 48.6 | 0.0 | 40.3 | Dtre-dominant |
| 2CS | 60.1 | 0.0 | 30.0 | Dtre-dominant |
| 2CW | 43.3 | 33.4 | 13.1 | Cgor-dominant |
| 2DS | 59.6 | 0.0 | 30.5 | Dtre-dominant |
| 2DW | 52.9 | 0.1 | 36.6 | Dtre-dominant |
| 3AS | 47.4 | 0.0 | 42.1 | Dtre-dominant |
| 3AW | 48.7 | 0.0 | 40.5 | Dtre-dominant |
| 3BS | 44.0 | 0.0 | 45.0 | Dtre-dominant |
| 3BW | 59.5 | 0.0 | 30.3 | Dtre-dominant |
| 3CS | 54.8 | 0.0 | 34.9 | Dtre-dominant |
| 3CW | 50.8 | 0.1 | 38.9 | Dtre-dominant |
| 3DS | 54.0 | 0.6 | 35.4 | Dtre-dominant |
| 3DW | 64.5 | 0.0 | 26.0 | Dtre-dominant |

### Summary

1. **D. trenchii dominance**: 79% of samples (19/24) are dominated by Clade D symbionts with no mapping to Clade C
2. **Sample 2AS**: 61.8% Clade C with no Clade D
3. **Symbiont composition may contribute to mapping variation**: Samples with more *D. trenchii* show higher multi-mapping rates (fragmented genome? 44,682 scaffolds)

## featureCounts Gene Quantification

### Quantification Strategy

Gene counts were generated **separately for each organism** with organism-specific multi-mapper handling:

| Organism | Multi-mapper Strategy | Rationale |
|----------|----------------------|-----------|
| *M. capitata* | **Discard** (default) | High-quality chromosome-level assembly (99.2% BUSCO), minimal true gene duplication |
| *C. goreaui* | **Count all** (`-M`) | Potential fragmented assembly (6,843 scaffolds), some genes split across contigs |
| *D. trenchii* | **Count all** (`-M`) | Whole genome duplication (WGD) creates real paralogs; discarding multi-mappers would lose 30-40% of signal |

**Reference**: Deschamps-Francoeur G et al. (2020) Handling multi-mapped reads in RNA-seq. Comput Struct Biotechnol J 18:1569-1576

### Matrix Dimensions

| Organism | Genes | Samples | Multi-mappers |
|----------|-------|---------|---------------|
| M. capitata | 54,384 | 24 | Discarded |
| C. goreaui | 45,322 | 24 | Counted |
| D. trenchii | 62,182 | 24 | Counted |

### Assignment Statistics by Organism

**M. capitata (Host) - No Multi-mappers**

| Metric | Mean | Range |
|--------|------|-------|
| Assignment Rate | 36.5% | 17.9% - 45.5% |
| Assigned Reads | 17.0M | 8.6M - 30.2M |

*Note: Assignment rate reflects host proportion of holobiont reads. Lower rates indicate higher symbiont content.*

**C. goreaui (Clade C) - With Multi-mappers**

| Sample Type | n | Assignment Rate | Assigned Reads |
|-------------|---|-----------------|----------------|
| Cgor-dominant | 3 | 18.6% - 44.3% | 15.1M - 33.9M |
| Mixed | 2 | 3.5% - 6.4% | 2.2M - 4.3M |
| Dtre-dominant | 19 | ~0% | <15K |

**D. trenchii (Clade D) - With Multi-mappers**

| Sample Type | n | Assignment Rate | Assigned Reads |
|-------------|---|-----------------|----------------|
| Dtre-dominant | 19 | 14.6% - 26.5% | 10.7M - 18.1M |
| Low Dtre | 5 | 0% - 10.4% | 2.9K - 7.0M |

### Symbiont Community Structure

| Profile | Samples | Description |
|---------|---------|-------------|
| **Cgor-dominant** | 2AS, 2CW, 1AW | >15% C. goreaui reads |
| **Mixed** | 1AS, 1BS | Detectable both symbionts |
| **Dtre-dominant** | 19 samples | Almost exclusively D. trenchii |

### Per-Sample Organism Proportions

| Sample | Mcap | Cgor | Dtre | Profile |
|--------|------|------|------|---------|
| 1AS | 34.0% | 7.6% | 58.4% | Mixed |
| 1AW | 50.6% | 37.3% | 12.1% | Cgor-dominant |
| 1BS | 62.7% | 14.1% | 23.1% | Mixed |
| 1BW | 61.8% | 0.0% | 38.2% | Dtre-dominant |
| 1CS | 57.3% | 0.0% | 42.7% | Dtre-dominant |
| 1CW | 44.6% | 0.0% | 55.4% | Dtre-dominant |
| 1DS | 40.6% | 0.0% | 59.4% | Dtre-dominant |
| 1DW | 53.2% | 0.0% | 46.8% | Dtre-dominant |
| 2AS | 20.3% | 79.7% | 0.0% | Cgor-dominant |
| 2AW | 46.6% | 0.0% | 53.3% | Dtre-dominant |
| 2BS | 61.7% | 0.0% | 38.3% | Dtre-dominant |
| 2BW | 50.2% | 0.0% | 49.8% | Dtre-dominant |
| 2CS | 60.5% | 0.0% | 39.5% | Dtre-dominant |
| 2CW | 39.5% | 45.9% | 14.6% | Cgor-dominant |
| 2DS | 61.2% | 0.0% | 38.8% | Dtre-dominant |
| 2DW | 54.2% | 0.0% | 45.7% | Dtre-dominant |
| 3AS | 47.3% | 0.0% | 52.7% | Dtre-dominant |
| 3AW | 49.6% | 0.0% | 50.4% | Dtre-dominant |
| 3BS | 43.7% | 0.0% | 56.3% | Dtre-dominant |
| 3BW | 61.0% | 0.0% | 39.0% | Dtre-dominant |
| 3CS | 55.5% | 0.0% | 44.5% | Dtre-dominant |
| 3CW | 52.5% | 0.0% | 47.5% | Dtre-dominant |
| 3DS | 55.2% | 1.0% | 43.8% | Dtre-dominant |
| 3DW | 67.8% | 0.0% | 32.1% | Dtre-dominant |

### Genes with Detected Expression

| Organism | Genes with Counts | Total Genes | Detection Rate |
|----------|-------------------|-------------|----------------|
| M. capitata | 50,031 | 54,384 | 92.0% |
| C. goreaui | 34,145 | 45,322 | 75.3% |
| D. trenchii | 42,997 | 62,182 | 69.1% |

### featureCounts Parameters

**M. capitata (Host)**
```bash
featureCounts -T 16 -p --countReadPairs \
    -t exon -g Parent \
    -a Mcap_genes_prefixed.gff3 \
    -o Mcap_counts_raw.txt *.bam
```

**C. goreaui and D. trenchii (Symbionts)**
```bash
featureCounts -T 16 -p --countReadPairs \
    -M \                           # Count multi-mapping reads
    -t exon -g Parent \
    -a {organism}_genes_prefixed.gff3 \
    -o {organism}_counts_raw.txt *.bam
```

### Output Files

| File | Description |
|------|-------------|
| `Mcap_counts.txt` | Host counts (no multi-mappers) |
| `Cgor_counts.txt` | Clade C counts (with multi-mappers) |
| `Dtre_counts.txt` | Clade D counts (with multi-mappers) |
| `*_counts_raw.txt` | Raw featureCounts output with metadata |
| `*_assignment_summary.txt` | Per-organism assignment statistics |
| `organism_proportions.txt` | Per-sample organism breakdown |

### Notes for DESeq2 Analysis

1. **M. capitata**: All 24 samples have sufficient coverage (8.6M - 30.2M reads)
2. **C. goreaui**: Only 3-5 samples have sufficient reads for DE analysis
3. **D. trenchii**: 21 samples have good coverage; exclude 2AS (2,899 reads)

## Experimental Design

- **Species**: *Montipora capitata* (Hawaiian rice coral)
- **Treatments**: 4 acdification treatment levels (A, B, C, D 'control')
- **Seasons**: Summer (S) and Winter (W)
- **Replicates**: 3 biological replicates per treatment × season
- **Sequencing**: Illumina NovaSeq, 150 bp paired-end

## DESeq2 Differential Expression Analysis

### Analysis Overview

Differential gene expression analysis was performed using DESeq2 to identify genes responsive to seawater acidification (OA) in both the coral host (*M. capitata*) and its dominant symbiont (*D. trenchii*). Analysis was conducted **within each season** to isolate treatment effects from strong seasonal signals.

### Treatment Subsetting Rationale

Analysis focused on **treatments B (high CO₂) and D (control)** only, excluding treatments A and C for the following reasons:

#### Ocean Acidification Contrast

Treatments A and C involved total alkalinity (TA) manipulation in addition to CO₂ changes, our goal was to decouple DIC and H+ for a more indepth analysis at calcification. This was achieved in the larger project by which the RNA subsamples were collected from but fialed to occur within the limited timeframe for the RNA-seq experiment. Treatment B (elevated CO₂, ambient TA) vs Treatment D (ambient CO₂, ambient TA) provide a clean test of ocean acidification effects and remains as an impactful scope for downstream functional analysis.

#### Symbiont Community Consistency

Three samples were identified as *Cladocopium goreaui* (Cgor)-dominant rather than *Durusdinium trenchii*-dominant:

| Sample | Treatment | Cgor % | Status |
|--------|-----------|--------|--------|
| 2AS | A | 79.7% | Excluded (Treatment A) |
| 2CW | C | 45.9% | Excluded (Treatment C) |
| 1AW | A | 37.3% | Excluded (Treatment A) |

All Cgor-dominant samples were in treatments A and C (not in the ocean acidification samples). By subsetting to treatments B and D, these samples were automatically excluded, ensuring consistent symbiont community composition across all analyzed samples to avoid confounding results with host gene expression responses.


### Experimental Design (B + D Subset)

| Factor | Levels | Description |
|--------|--------|-------------|
| Treatment | B, D | High CO₂ vs Control |
| Season | Summer, Winter | Seasonal replication |
| Replicates | 3 per group | Biological replicates |
| **Total samples** | **12** | 6 per treatment |
```
         Summer  Winter
    D       3       3      (Control)
    B       3       3      (High CO₂)
```

### Chemistry Conditions (B vs D)

| Treatment | Season | H⁺ (nmol kg-1) | pH | DIC (µmol kg-1) |
|-----------|--------|--------------|-----|---------------|
| D (Control) | Summer | 8.82 | 8.05 | 1860 |
| D (Control) | Winter | 7.95 | 8.10 | 1869 |
| B (High CO₂) | Summer | 16.19 | 7.79 | 1992 |
| B (High CO₂) | Winter | 14.99 | 7.82 | 2010 |

**H⁺ contrast**: ~7-8 nmol/kg difference (~0.25 pH units) between treatments.

### DESeq2 Model

Within-season analysis using treatment as a categorical factor:
```r
design = ~ treatment   # B vs D within each season
```

- **Reference level**: Treatment D (control)
- **Contrast**: Treatment B vs D (log2 fold change)

### Results Summary

| Organism | Season | Genes Tested | DEGs (padj < 0.05) | Up in B | Down in B |
|----------|--------|--------------|-------------------|---------|-----------|
| **M. capitata** | **Summer** | 23,294 | **772** | 405 | 367 |
| **M. capitata** | **Winter** | 23,062 | **215** | 132 | 83 |
| D. trenchii | Summer | 35,008 | 2 | 1 | 1 |
| **D. trenchii** | **Winter** | 35,147 | **279** | 200 | 79 |

### Signal Detection Statistics

| Organism | Season | p < 0.01 | Expected | p < 0.05 | Expected |
|----------|--------|----------|----------|----------|----------|
| M. capitata | Summer | 1,234 | 233 | 2,343 | 1,165 |
| M. capitata | Winter | 623 | 231 | 1,579 | 1,153 |
| D. trenchii | Summer | 253 | 350 | 1,088 | 1,750 |
| D. trenchii | Winter | 1,409 | 351 | 3,772 | 1,757 |


### Summary

#### Host Response (*M. capitata*)
- **Strong summer response**: 772 DEGs with balanced regulation (405 up, 367 down)
- **Moderate winter response**: 215 DEGs with slight up-regulation (132 up, 83 down)
- **Seasonal amplification**: OA effects ~3.6× stronger in summer, possibly due to compounding thermal stress
- **Large effect sizes**: Top DEGs show |LFC| > 5

#### Symbiont Response (*D. trenchii*)
- **Minimal summer response**: Only 2 DEGs — symbiont appears buffered from OA during warm season
- **Strong winter response**: 279 DEGs (estimated 3,963 true DEGs based on Pi0)
- **Predominantly up-regulated**: 200 up vs 79 down under acidification
- **Opposite seasonal pattern from host**: Symbiont most responsive when host is least responsive

### DESeq2 Parameters
```r
# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = sample_info,
    design = ~ treatment
)

# Filter low-count genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# Run DESeq2 (Wald test)
dds <- DESeq(dds)

# Extract results
results(dds, name = "treatment_B_vs_D", alpha = 0.05)
```

### Output Files

| File | Description |
|------|-------------|
| `results/Mcap_Summer_BvsD.csv` | Host summer DEG results |
| `results/Mcap_Winter_BvsD.csv` | Host winter DEG results |
| `results/Dtre_Summer_BvsD.csv` | Symbiont summer DEG results |
| `results/Dtre_Winter_BvsD.csv` | Symbiont winter DEG results |
| `results/DEG_summary.csv` | Summary statistics |
| `figures/pvalue_histograms.pdf` | P-value distributions |
| `figures/volcano_*.pdf` | Volcano plots per analysis |
| `Mcap_DESeq2.RData` | R objects for downstream analysis |
| `Dtre_DESeq2.RData` | R objects for downstream analysis |


## Directory Structure

```
mc_rework/
├── 00_raw_fastq/           # Symlinks to original FASTQ files
├── 01_qc_raw/              # FastQC/MultiQC on raw reads
├── 02_trimmed/             # Cutadapt-trimmed reads
├── 03_qc_trimmed/          # FastQC/MultiQC on trimmed reads
├── 04_reference/           # Holobiont genome & STAR index
│   ├── holobiont_genome.fa           # Combined genome (Mcap + Cgor + Dtre)
│   ├── holobiont_genes.gff3          # Combined annotations
│   ├── star_index_holobiont/         # STAR genome index (50 GB)
│   ├── Mcap_transcript_to_gene.tsv   # Mcap transcript-gene mapping
│   └── organism_distribution.tsv     # Per-sample organism stats
├── 05_star_align/          # STAR alignment output (per sample)
│   ├── {sample}/
│   │   ├── {sample}_Aligned.sortedByCoord.out.bam
│   │   ├── {sample}_Aligned.sortedByCoord.out.bam.bai
│   │   ├── {sample}_Log.final.out
│   │   └── {sample}_Log.out
│   ├── alignment_summary.tsv
│   └── organism_distribution.tsv
├── 06_featurecounts/       # Gene-level counts (pending)
├── 07_deseq2/              # DESeq2 analysis (pending)
├── logs/                   # SLURM job logs
├── scripts/                # Pipeline scripts
└── metadata/               # Sample information
```


## Reference Sources

### Host: *Montipora capitata*
- **Assembly**: HIv3 (chromosome-level, 781 Mb, 14 chromosomes + 1,685 scaffolds)
- **Scaffolds**: 1,699
- **Source**: Stephens et al. 2022 (GigaDB 102268)
- **BUSCO**: 99.2% complete
- **Citation**: Stephens TG et al. GigaScience. 2022;11:giac098

### Symbiont: *Cladocopium goreaui* (Clade C)
- **Strain**: SCF055
- **Scaffolds**: 6,843
- **Source**: Chen et al. 2022 (UQ eSpace)
- **BUSCO**: 82.4% complete
- **Citation**: Chen Y et al. Microorganisms. 2022;10(8):1662

### Symbiont: *Durusdinium trenchii* (Clade D)
- **Strain**: SCF082
- **Scaffolds**: 44,682
- **Source**: Dougan et al. 2022 (UQ eSpace)
- **Citation**: Dougan KE et al. Science Advances. 2024;10:eadn2218

## Conda Environment

The `mcap_rnaseq` environment includes:
- **QC**: fastqc, multiqc
- **Trimming**: cutadapt
- **Alignment**: star (2.7.11b), samtools
- **Counting**: subread (featureCounts)
- **Analysis**: R, DESeq2, tximport
- **Visualization**: pheatmap, EnhancedVolcano

## Troubleshooting

### STAR memory requirements
The holobiont STAR index requires ~50GB RAM to load. Use exclusive node access for stability:
```bash
#SBATCH --exclusive
#SBATCH --mem=320G
```

### Multi-mapping variation
Variation in unique mapping rates (61-74%) reflects biological differences in symbiont composition, not technical issues. Samples with more *D. trenchii* show higher multi-mapping due to its potentially fragmented genome assembly.

## References

- Drury C et al. (2022) Intrapopulation adaptive variance supports thermal tolerance in a reef-building coral. Communications Biology 5:486
- Love MI et al. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15:550
- Dobin A et al. (2013) STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29:15-21

## Contact

David Armstrong  
CROH Lab, Texas A&M University-Corpus Christi  
darmstrong4@islander.tamucc.edu
