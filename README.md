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
Variation in unique mapping rates (61-74%) reflects biological differences in symbiont composition, not technical issues. Samples with more *D. trenchii* show higher multi-mapping due to its fragmented genome assembly.

## References

- Drury C et al. (2022) Intrapopulation adaptive variance supports thermal tolerance in a reef-building coral. Communications Biology 5:486
- Love MI et al. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15:550
- Dobin A et al. (2013) STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29:15-21

## Contact

David Armstrong  
CROH Lab, Texas A&M University-Corpus Christi  
darmstrong4@islander.tamucc.edu
