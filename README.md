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

### Why featureCounts Instead of STAR Counting?

While STAR has built-in read counting (`--quantMode GeneCounts`), featureCounts was chosen for this holobiont analysis because:

| Feature | STAR Counting | featureCounts |
|---------|---------------|---------------|
| Multi-mapper handling | Discard or count all | **Fractional counting** (1/n per location) |
| GFF3 support | Limited (prefers GTF) | **Native GFF3 support** |
| Attribute flexibility | gene_id only | **Any attribute** (we used `Parent`) |

With 18.8% multi-mapping reads (due to *D. trenchii*'s fragmented genome), fractional counting recovered ~20% more assigned reads compared to discarding multi-mappers.

### Quantification Summary

| Metric | Value |
|--------|-------|
| **Total Features** | 161,888 |
| **Mcap Transcripts** | 54,384 |
| **Cgor Genes** | 45,322 |
| **Dtre Genes** | 62,182 |
| **Mean Assignment Rate** | 51.5% |
| **Assignment Range** | 48.0% - 58.3% |

### featureCounts Parameters
```bash
featureCounts \
    -T 16 \                    # Threads
    -p \                       # Paired-end mode
    --countReadPairs \         # Count read pairs (fragments)
    -M \                       # Count multi-mapping reads
    --fraction \               # Fractional counting (1/n per location)
    -t exon \                  # Feature type
    -g Parent \                # Group by Parent attribute (GFF3)
    -a holobiont_genes.gff3 \  # Annotation file
    -o gene_counts_raw.txt \   # Output file
    *.bam                      # Input BAM files
```

### Per-Sample Assignment Statistics

| Sample | Assigned | NoFeatures | Ambiguous | Unmapped | Rate |
|--------|----------|------------|-----------|----------|------|
| 1AS | 31,148,831 | 21,306,639 | 3,247,789 | 7,673,841 | 49.1% |
| 1AW | 45,335,945 | 25,223,762 | 2,910,269 | 7,794,751 | 55.8% |
| 1BS | 34,645,649 | 24,427,279 | 2,144,208 | 6,085,946 | 51.5% |
| 1BW | 36,206,909 | 24,106,263 | 2,594,029 | 6,737,239 | 52.0% |
| 1CS | 37,020,783 | 23,584,788 | 2,925,364 | 7,507,502 | 52.1% |
| 1CW | 32,368,661 | 21,309,434 | 2,986,602 | 7,272,925 | 50.6% |
| 1DS | 30,618,218 | 22,773,030 | 3,066,892 | 7,281,418 | 48.0% |
| 1DW | 38,035,924 | 25,645,097 | 3,125,121 | 8,008,415 | 50.8% |
| 2AS | 44,631,300 | 19,842,296 | 3,492,749 | 8,542,628 | 58.3% |
| 2AW | 38,049,781 | 26,361,815 | 3,488,080 | 8,239,307 | 50.0% |
| 2BS | 31,947,007 | 22,588,136 | 2,330,219 | 5,946,296 | 50.9% |
| 2BW | 28,605,974 | 19,129,262 | 2,464,752 | 6,121,201 | 50.8% |
| 2CS | 33,497,267 | 25,564,737 | 2,470,362 | 6,583,859 | 49.2% |
| 2CW | 45,627,254 | 24,065,449 | 3,466,083 | 8,064,204 | 56.2% |
| 2DS | 34,723,505 | 24,380,928 | 2,530,382 | 6,629,449 | 50.9% |
| 2DW | 35,901,420 | 24,776,177 | 2,851,139 | 7,302,194 | 50.7% |
| 3AS | 36,198,770 | 25,337,564 | 3,250,764 | 7,463,313 | 50.1% |
| 3AW | 36,334,992 | 24,029,361 | 3,214,474 | 7,599,981 | 51.0% |
| 3BS | 34,470,100 | 23,469,047 | 3,285,333 | 7,382,960 | 50.2% |
| 3BW | 34,831,546 | 24,693,139 | 2,538,496 | 6,868,446 | 50.5% |
| 3CS | 36,274,411 | 25,959,223 | 2,903,512 | 7,293,066 | 50.1% |
| 3CW | 42,849,489 | 27,022,045 | 3,446,366 | 8,228,912 | 52.5% |
| 3DS | 37,089,706 | 25,394,063 | 2,975,559 | 7,069,212 | 51.1% |
| 3DW | 51,244,118 | 34,066,965 | 3,230,998 | 9,120,165 | 52.5% |

### Organism-Specific Gene Counts

| Sample | Mcap | Cgor | Dtre | Total | %Mcap | %Cgor | %Dtre |
|--------|------|------|------|-------|-------|-------|-------|
| 1AS | 10,675,410 | 1,553,766 | 11,789,148 | 24,018,324 | 44.4 | 6.5 | 49.1 |
| 1AW | 22,372,859 | 10,573,624 | 3,424,078 | 36,370,561 | 61.5 | 29.1 | 9.4 |
| 1BS | 20,657,552 | 3,022,219 | 4,911,062 | 28,590,833 | 72.3 | 10.6 | 17.2 |
| 1BW | 21,318,437 | 2,197 | 8,459,339 | 29,779,973 | 71.6 | 0.0 | 28.4 |
| 1CS | 20,395,171 | 2,472 | 9,873,888 | 30,271,531 | 67.4 | 0.0 | 32.6 |
| 1CW | 14,226,534 | 3,135 | 11,395,876 | 25,625,545 | 55.5 | 0.0 | 44.5 |
| 1DS | 12,367,355 | 3,180 | 11,675,765 | 24,046,300 | 51.4 | 0.0 | 48.6 |
| 1DW | 19,612,480 | 3,314 | 11,106,059 | 30,721,853 | 63.8 | 0.0 | 36.2 |
| 2AS | 9,444,842 | 23,690,411 | 2,016 | 33,137,269 | 28.5 | 71.5 | 0.0 |
| 2AW | 17,444,003 | 3,336 | 12,874,168 | 30,321,507 | 57.5 | 0.0 | 42.5 |
| 2BS | 18,740,704 | 1,974 | 7,522,035 | 26,264,713 | 71.4 | 0.0 | 28.6 |
| 2BW | 13,996,172 | 2,429 | 8,922,941 | 22,921,542 | 61.1 | 0.0 | 38.9 |
| 2CS | 19,348,380 | 2,641 | 8,155,507 | 27,506,528 | 70.3 | 0.0 | 29.6 |
| 2CW | 17,990,413 | 12,974,883 | 4,249,479 | 35,214,775 | 51.1 | 36.8 | 12.1 |
| 2DS | 20,230,661 | 2,216 | 8,227,462 | 28,460,339 | 71.1 | 0.0 | 28.9 |
| 2DW | 18,816,864 | 2,942 | 10,196,366 | 29,016,172 | 64.8 | 0.0 | 35.1 |
| 3AS | 16,781,942 | 3,378 | 12,030,023 | 28,815,343 | 58.2 | 0.0 | 41.7 |
| 3AW | 17,600,824 | 3,122 | 11,540,746 | 29,144,692 | 60.4 | 0.0 | 39.6 |
| 3BS | 14,864,575 | 3,554 | 12,372,262 | 27,240,391 | 54.6 | 0.0 | 45.4 |
| 3BW | 20,263,141 | 4,949 | 8,356,873 | 28,624,963 | 70.8 | 0.0 | 29.2 |
| 3CS | 19,397,476 | 2,725 | 9,989,128 | 29,389,329 | 66.0 | 0.0 | 34.0 |
| 3CW | 21,827,910 | 3,438 | 12,653,031 | 34,484,379 | 63.3 | 0.0 | 36.7 |
| 3DS | 19,749,923 | 217,022 | 10,094,553 | 30,061,498 | 65.7 | 0.7 | 33.6 |
| 3DW | 32,770,872 | 2,934 | 9,982,608 | 42,756,414 | 76.6 | 0.0 | 23.3 |

### Output Files

| File | Description |
|------|-------------|
| `gene_counts.txt` | Integer counts for DESeq2 (all organisms) |
| `gene_counts_fractional.txt` | Original fractional counts |
| `gene_counts_Mcap.txt` | Host coral counts only |
| `gene_counts_Cgor.txt` | Clade C symbiont counts only |
| `gene_counts_Dtre.txt` | Clade D symbiont counts only |
| `organism_proportions.txt` | Per-sample organism breakdown |
| `assignment_summary.txt` | Per-sample assignment statistics |

### Notes on Assignment Rates

The ~51% mean assignment rate is typical for holobiont RNA-seq and reflects:

1. **Intergenic/intronic reads** (~35% NoFeatures): Reads mapping outside annotated exons
2. **Unmapped reads** (~12%): Reads not aligning to the holobiont reference
3. **Ambiguous reads** (~5%): Reads overlapping multiple features

Higher assignment rates in samples 2AS (58.3%) and 1AW (55.8%) correlate with Cgor-dominant profiles, as the *C. goreaui* genome is potentially better assembled than *D. trenchii*.

## Experimental Design

- **Species**: *Montipora capitata* (Hawaiian rice coral)
- **Treatments**: 4 acdification treatment levels (A, B, C, D 'control')
- **Seasons**: Summer (S) and Winter (W)
- **Replicates**: 3 biological replicates per treatment × season
- **Sequencing**: Illumina NovaSeq, 150 bp paired-end


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
