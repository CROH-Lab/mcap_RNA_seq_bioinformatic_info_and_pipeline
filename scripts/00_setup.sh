#!/bin/bash
# ============================================================================
# M. capitata RNA-seq Holobiont Pipeline - Initial Setup
# Run this INTERACTIVELY (not as SLURM job) to set up the environment
#
# Pipeline Overview:
#   00 - Setup (this script)
#   01 - Raw read QC (FastQC/MultiQC)
#   02 - Trimming (Cutadapt)
#   03 - Trimmed read QC (FastQC/MultiQC)
#   04 - STAR index (holobiont genome)
#   05 - STAR alignment
#   06 - featureCounts
#   07 - DESeq2 (R script)
# ============================================================================

set -euo pipefail

# --- Configuration ---
BASE_DIR="/home/darmstrong4/mc_rework"
RAW_FASTQ_SOURCE="/home/darmstrong4/Montipora_capitata_ocean_acidifcation_RNA_seq/00_fastq"
SCRIPTS_DIR="${BASE_DIR}/scripts"

echo "=============================================="
echo "M. capitata Holobiont RNA-seq Pipeline Setup"
echo "=============================================="

# --- Create Directory Structure ---
echo ""
echo "=== Creating directory structure ==="
mkdir -p ${BASE_DIR}/{00_raw_fastq,01_qc_raw,02_trimmed,03_qc_trimmed,04_reference,05_star_align,06_counts,07_deseq2,logs,scripts}

cd ${BASE_DIR}
echo "Working directory: $(pwd)"

# --- Copy scripts ---
echo ""
echo "=== Copying pipeline scripts ==="
SCRIPT_SOURCE="$(dirname "$0")"
cp ${SCRIPT_SOURCE}/*.sh ${SCRIPTS_DIR}/ 2>/dev/null || echo "Scripts already in place or copy manually"
cp ${SCRIPT_SOURCE}/*.yml ${SCRIPTS_DIR}/ 2>/dev/null || echo "YAML file already in place or copy manually"
cp ${SCRIPT_SOURCE}/*.R ${SCRIPTS_DIR}/ 2>/dev/null || echo "R scripts already in place or copy manually"
chmod +x ${SCRIPTS_DIR}/*.sh

# --- Create or update conda environment ---
echo ""
echo "=== Setting up Conda environment ==="
echo "This may take 15-30 minutes..."

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found. Please load anaconda/miniconda first."
    echo "Try: module load anaconda3 or similar"
    exit 1
fi

# Create environment if it doesn't exist
if conda env list | grep -q "mcap_rnaseq"; then
    echo "Environment 'mcap_rnaseq' already exists."
    echo "To update: conda env update -n mcap_rnaseq -f ${SCRIPTS_DIR}/mcap_rnaseq_env.yml"
    echo "To recreate: conda env remove -n mcap_rnaseq && conda env create -f ${SCRIPTS_DIR}/mcap_rnaseq_env.yml"
else
    echo "Creating conda environment 'mcap_rnaseq'..."
    conda env create -f ${SCRIPTS_DIR}/mcap_rnaseq_env.yml
fi

# --- Link raw fastq files ---
echo ""
echo "=== Linking raw FASTQ files ==="
if [[ -d ${RAW_FASTQ_SOURCE} ]]; then
    for f in ${RAW_FASTQ_SOURCE}/*.fastq.gz; do
        if [[ -f "$f" ]]; then
            ln -sf "$f" 00_raw_fastq/ 2>/dev/null || true
        fi
    done
    echo "Linked $(ls 00_raw_fastq/*.fastq.gz 2>/dev/null | wc -l) FASTQ files"
else
    echo "WARNING: Source directory not found: ${RAW_FASTQ_SOURCE}"
    echo "Please link FASTQ files manually to ${BASE_DIR}/00_raw_fastq/"
fi

# --- Create sample info file ---
echo ""
echo "=== Creating sample info file ==="
if [[ ! -f ${BASE_DIR}/sample_info.txt ]]; then
    cat > ${BASE_DIR}/sample_info.txt << 'EOF'
sample	genotype	treatment	season
1AS	1	A	Summer
1AW	1	A	Winter
1BS	1	B	Summer
1BW	1	B	Winter
1CS	1	C	Summer
1CW	1	C	Winter
1DS	1	D	Summer
1DW	1	D	Winter
2AS	2	A	Summer
2AW	2	A	Winter
2BS	2	B	Summer
2BW	2	B	Winter
2CS	2	C	Summer
2CW	2	C	Winter
2DS	2	D	Summer
2DW	2	D	Winter
3AS	3	A	Summer
3AW	3	A	Winter
3BS	3	B	Summer
3BW	3	B	Winter
3CS	3	C	Summer
3CW	3	C	Winter
3DS	3	D	Summer
3DW	3	D	Winter
EOF
    echo "Created sample_info.txt"
else
    echo "sample_info.txt already exists"
fi

# --- Verify setup ---
echo ""
echo "=============================================="
echo "Setup Complete!"
echo "=============================================="
echo ""
echo "Directory structure:"
ls -la ${BASE_DIR}/
echo ""
echo "FASTQ files linked:"
ls 00_raw_fastq/*.fastq.gz 2>/dev/null | head -5
echo "... ($(ls 00_raw_fastq/*.fastq.gz 2>/dev/null | wc -l) total files)"
echo ""
echo "Reference files needed in 04_reference/:"
echo "  - Montipora_capitata_HIv3.assembly.fasta.gz  (host genome)"
echo "  - Mcap.genes.gff3                            (host annotation)"
echo "  - Cladocopium_goreaui_genome_fa.gz           (Clade C genome)"
echo "  - Cladocopium_goreaui_genes_gff3.gz          (Clade C annotation)"
echo "  - Dtrenchii_SCF082_ASSEMBLY_fasta.gz         (Clade D genome)"
echo "  - Dtrenchii_SCF082_ANNOT_gff.gz              (Clade D annotation)"
echo ""
echo "Pipeline steps:"
echo "  1. conda activate mcap_rnaseq"
echo "  2. sbatch scripts/01_fastqc_raw.sh        # Raw QC"
echo "  3. sbatch scripts/02_trim.sh              # Trimming"
echo "  4. sbatch scripts/03_fastqc_trimmed.sh    # Trimmed QC"
echo "  5. sbatch scripts/04_star_index.sh        # Build STAR index"
echo "  6. sbatch scripts/05_star_align.sh        # Align reads"
echo "  7. sbatch scripts/06_featurecounts.sh     # Count reads"
echo "  8. Rscript scripts/07_deseq2.R            # Differential expression"
echo ""
