# *M. capitata* Holobiont RNA-seq Pipeline

RNA-seq analysis of *Montipora capitata* coral holobiont response to ocean acidification (OA). Genome-based alignment with STAR quantifies gene expression for both the coral host and its algal symbiont *Durusdinium trenchii* (Clade D) simultaneously.

**Author**: David Armstrong — CROH Lab, Texas A&M University-Corpus Christi

---

## Experimental Design

| Factor | Levels | Description |
|--------|--------|-------------|
| Species | *M. capitata* | Hawaiian rice coral |
| Treatment | B (high CO₂), D (control) | OA vs ambient |
| Season | Summer, Winter | Seasonal replication |
| Replicates | 3 per group | Biological replicates |
| **Samples** | **12** | 2 treatments × 2 seasons × 3 replicates |
| Sequencing | Illumina NovaSeq | 150 bp PE, ~57M reads/sample |

### Treatment Subsetting Rationale

Analysis was restricted to treatments B and D (of original A–D design). Treatments A and C involved total alkalinity manipulation that confounds OA interpretation. All three *Cladocopium goreaui*-dominant samples (2AS, 2CW, 1AW) fell within A/C, so subsetting also ensured consistent *D. trenchii*-dominant symbiont communities across all analyzed samples.

### Carbonate Chemistry (B vs D)

| Treatment | Season | pH(T) | pCO₂ (µatm) | Ωarag |
|-----------|--------|-------|-------------|-------|
| D (Control) | Summer | 8.05 | ~400 | ~3.0 |
| D (Control) | Winter | 8.10 | ~350 | ~3.2 |
| B (High CO₂) | Summer | 7.79 | ~800 | ~1.8 |
| B (High CO₂) | Winter | 7.82 | ~750 | ~2.0 |

---

## Sequencing & Processing Quality (B/D Samples Only)

### Raw Read Statistics

| Sample | Reads (M) | Quality Score | % ≥ Q30 | Duplication |
|--------|-----------|---------------|---------|-------------|
| 1BS | 55.1 | Q38.89 | 94.4% | 16.1% |
| 1BW | 56.8 | Q38.91 | 94.5% | 15.8% |
| 1DS | 50.2 | Q38.93 | 94.6% | 15.2% |
| 1DW | 59.9 | Q38.90 | 94.5% | 16.5% |
| 2BS | 51.4 | Q38.88 | 94.3% | 16.8% |
| 2BW | 44.9 | Q38.92 | 94.6% | 14.4% |
| 2DS | 55.6 | Q38.90 | 94.5% | 16.4% |
| 2DW | 56.8 | Q38.91 | 94.5% | 16.0% |
| 3BS | 54.3 | Q38.93 | 94.6% | 15.5% |
| 3BW | 56.2 | Q38.88 | 94.3% | 17.1% |
| 3DS | 57.8 | Q38.91 | 94.6% | 15.8% |
| 3DW | 80.9 | Q38.89 | 94.4% | 16.2% |

**Total**: ~680M read pairs, 150 bp PE, Illumina NovaSeq

### Trimming (Cutadapt)

```bash
cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -a "A{10}" -A "T{10}" \
    -q 20,20 --minimum-length 50
```

| Sample | Trimmed Reads | Retention |
|--------|---------------|-----------|
| 1BS | 53,749,496 | 97.5% |
| 1BW | 55,594,613 | 97.7% |
| 1DS | 49,208,380 | 98.0% |
| 1DW | 58,708,738 | 97.9% |
| 2BS | 50,120,093 | 97.6% |
| 2BW | 43,974,226 | 98.0% |
| 2DS | 54,211,740 | 97.6% |
| 2DW | 55,691,279 | 97.9% |
| 3BS | 53,177,277 | 98.0% |
| 3BW | 54,895,066 | 97.6% |
| 3DS | 56,598,638 | 97.9% |
| 3DW | 78,933,530 | 97.6% |

### STAR Alignment

```bash
STAR --runMode alignReads --runThreadN 64 \
    --genomeDir star_index_holobiont \
    --outSAMtype BAM Unsorted \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMax 1000000 \
    --alignSJoverhangMin 8 --alignSJDBoverhangMin 1
```

| Sample | Total Reads | Unique % | Multi % | Total Mapped % |
|--------|-------------|----------|---------|----------------|
| 1BS | 53,749,496 | 74.08% | 14.60% | 88.91% |
| 1BW | 55,594,613 | 70.21% | 17.67% | 88.03% |
| 1DS | 49,208,380 | 62.73% | 22.47% | 85.33% |
| 1DW | 58,708,738 | 66.45% | 19.91% | 86.53% |
| 2BS | 50,120,093 | 70.84% | 17.29% | 88.29% |
| 2BW | 43,974,226 | 65.27% | 20.81% | 86.21% |
| 2DS | 54,211,740 | 70.05% | 17.72% | 87.96% |
| 2DW | 55,691,279 | 67.36% | 19.53% | 87.03% |
| 3BS | 53,177,277 | 63.98% | 22.14% | 86.25% |
| 3BW | 54,895,066 | 70.07% | 17.42% | 87.68% |
| 3DS | 56,598,638 | 67.98% | 19.53% | 87.72% |
| 3DW | 78,933,530 | 72.44% | 16.01% | 88.58% |

### Organism Read Distribution

| Sample | Host (Mcap) % | Clade C (Cgor) % | Clade D (Dtre) % | Profile |
|--------|---------------|------------------|------------------|---------|
| 1BS | 63.5 | 8.8 | 18.5 | Dtre-dominant |
| 1BW | 59.5 | 0.0 | 30.6 | Dtre-dominant |
| 1DS | 42.3 | 0.0 | 46.0 | Dtre-dominant |
| 1DW | 51.1 | 0.0 | 38.0 | Dtre-dominant |
| 2BS | 60.4 | 0.0 | 29.9 | Dtre-dominant |
| 2BW | 48.6 | 0.0 | 40.3 | Dtre-dominant |
| 2DS | 59.6 | 0.0 | 30.5 | Dtre-dominant |
| 2DW | 52.9 | 0.1 | 36.6 | Dtre-dominant |
| 3BS | 44.0 | 0.0 | 45.0 | Dtre-dominant |
| 3BW | 59.5 | 0.0 | 30.3 | Dtre-dominant |
| 3DS | 54.0 | 0.6 | 35.4 | Dtre-dominant |
| 3DW | 64.5 | 0.0 | 26.0 | Dtre-dominant |

All 12 B/D samples are *D. trenchii*-dominant with negligible *C. goreaui* mapping.

---

## Pipeline Overview

| Step | Directory | Description |
|------|-----------|-------------|
| 0 | `00_raw_fastq/` | Raw FASTQ symlinks |
| 1 | `01_qc_raw/` | FastQC/MultiQC on raw reads |
| 2 | `02_trimmed/` | Cutadapt adapter/quality trimming |
| 3 | `03_qc_trimmed/` | FastQC/MultiQC on trimmed reads |
| 4 | `04_reference/` | Holobiont genome (Mcap + Cgor + Dtre) & STAR index |
| 5 | `05_star_align/` | STAR alignment to holobiont genome |
| 6 | `06_featurecounts/` | Gene-level quantification (organism-specific) |
| 7 | `07_deseq2/` | Differential expression analysis (B vs D, within season) |
| 8 | `08_host_deg_annotation/` | Host DEG functional annotation (BLASTx) |
| 9 | `09_symbiont_deg_annotation/` | Symbiont DEG functional annotation (BLASTx) |
| 10 | `10_host_GO_MWU/` | Host GO enrichment (rank-based MWU) |
| 11 | `11_symbiont_GO_MWU/` | Symbiont GO enrichment (rank-based MWU) |
| 12 | `12_publication_figures/` | Publication-ready figures and tables |
| 13 | `13_cal_chem/` | Carbonate chemistry & calcification analysis |

---

## Key Results

### Differential Expression (DESeq2, padj < 0.05)

| Organism | Season | DEGs | Up | Down |
|----------|--------|------|-----|------|
| **Host** | **Summer** | **772** | 405 | 367 |
| **Host** | **Winter** | **215** | 132 | 83 |
| Symbiont | Summer | 2 | 1 | 1 |
| **Symbiont** | **Winter** | **279** | 200 | 79 |

The host responds most strongly in summer (~3.6× more DEGs), while the symbiont responds primarily in winter — suggesting complementary seasonal stress response strategies within the holobiont.

### Functional Enrichment (GO_MWU, p.adj < 0.05)

| Organism | Season | Significant GO Terms | Predominantly |
|----------|--------|---------------------|---------------|
| Host | Summer | 163 | Upregulated (60%) |
| Host | Winter | 204 | Downregulated (68%) |
| Symbiont | Summer | 46 | Downregulated (67%) |
| Symbiont | Winter | 257 | Upregulated (58%) |

### DEG Annotation (BLASTx vs UniProt)

| Organism | Total DEGs | Annotated | Rate |
|----------|-----------|-----------|------|
| Host | 972 | 905 | 93.1% |
| Symbiont | 281 | 163 | 58.0% |

---

## Publication Figures

| Figure | Description |
|--------|-------------|
| Fig 1 | Calcification barplot (in `13_cal_chem/`) |
| Fig 2 | Venn diagrams + ridgeline plots (DEG overview) |
| Fig 3 | Combined heatmaps (Host Summer/Winter, Symbiont Winter) |
| Fig 4 | Holobiont PCA (Host + Symbiont combined per season) |
| Fig 5 | Summer representative GO dendrogram (Host + Symbiont) |
| Fig 6 | Winter representative GO dendrogram (Host + Symbiont) |
| Fig 7 | Chord diagram (shared GO terms between Host and Symbiont) |
| Fig S1 | Individual PCA 4-panel grid |

---

## Holobiont Reference Genome

| Organism | Prefix | Scaffolds | Source |
|----------|--------|-----------|--------|
| *Montipora capitata* | Mcap_ | 1,699 | Stephens et al. 2022, HIv3 |
| *Cladocopium goreaui* | Cgor_ | 6,843 | Chen et al. 2022, SCF055 |
| *Durusdinium trenchii* | Dtre_ | 44,682 | Dougan et al. 2022, SCF082 |

### featureCounts Strategy

| Organism | Multi-mapper Handling | Rationale |
|----------|----------------------|-----------|
| *M. capitata* | Discard (default) | Chromosome-level assembly, 99.2% BUSCO |
| *D. trenchii* | Count all (`-M`) | WGD creates real paralogs; discarding loses 30–40% signal |

---

## Methods Summary

### Differential Expression (DESeq2)
Within-season Wald test, `design = ~ treatment` (B vs D). Genes filtered: ≥10 counts in ≥3 samples. Significance: padj < 0.05.

### Functional Enrichment (GO_MWU)
Rank-based Mann-Whitney U test using signed −log₁₀(p-value) as continuous measure. BH-corrected, semantic clustering at height 0.25. Reference: Wright et al. (2015) BMC Genomics 16:371.

### DEG Annotation (BLASTx)
Two-tier approach: SwissProt (curated) then TrEMBL (comprehensive), E-value < 1e-5. Quality tiers: HIGH (e < 1e-20, ≥40% identity, ≥60% coverage), MEDIUM, LOW.

---

## Conda Environment

```
mcap_rnaseq: fastqc, multiqc, cutadapt, star, samtools, subread, R, DESeq2
```

## References

- Dobin A et al. (2013) STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* 29:15-21
- Love MI et al. (2014) DESeq2. *Genome Biology* 15:550
- Wright RM et al. (2015) GO_MWU. *BMC Genomics* 16:371
- Stephens TG et al. (2022) *M. capitata* genome. *GigaScience* 11:giac098
- Dougan KE et al. (2024) *D. trenchii* genome. *Science Advances* 10:eadn2218

## Contact

David Armstrong — darmstrong4@islander.tamucc.edu
CROH Lab, Texas A&M University-Corpus Christi
