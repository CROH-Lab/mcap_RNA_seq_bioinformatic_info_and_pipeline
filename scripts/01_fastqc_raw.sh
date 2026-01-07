#!/bin/bash
#SBATCH -J mcap_fastqc_raw           # Job name
#SBATCH -o logs/fastqc_raw_%j.out    # Standard output
#SBATCH -e logs/fastqc_raw_%j.err    # Standard error
#SBATCH --partition=normal           # Use the normal partition
#SBATCH --nodes=1                    # Use a single node
#SBATCH --cpus-per-task=8            # Allocate 8 CPUs
#SBATCH --time=0-04:00:00            # Runtime (4 hours)
#SBATCH --mail-user=darmstrong4@islander.tamucc.edu
#SBATCH --mail-type=BEGIN,END,FAIL   # Notifications for job status

# ============================================================================
# M. capitata RNA-seq Rework Pipeline - Step 1: FastQC on Raw Reads
# ============================================================================

# Don't use -u flag (causes issues with conda activation)
set -eo pipefail

# --- Configuration ---
BASE_DIR="/home/darmstrong4/mc_rework"
THREADS=8

cd ${BASE_DIR}

# --- Activate conda environment ---
# Use direct conda initialization instead of sourcing bashrc
eval "$(conda shell.bash hook)"
conda activate mcap_rnaseq

echo "=== FastQC on Raw Reads ==="
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo "FastQC version: $(fastqc --version)"

# --- Run FastQC ---
mkdir -p 01_qc_raw

fastqc \
    -t ${THREADS} \
    -o 01_qc_raw \
    00_raw_fastq/*.fastq.gz

# --- Generate MultiQC report ---
echo ""
echo "=== Generating MultiQC Report ==="
multiqc 01_qc_raw \
    -o 01_qc_raw \
    -n multiqc_raw_reads \
    --title "M. capitata Raw Reads Quality Report"

echo ""
echo "=== FastQC Complete ==="
echo "End time: $(date)"
echo "Reports saved to: ${BASE_DIR}/01_qc_raw/"
echo "MultiQC report: ${BASE_DIR}/01_qc_raw/multiqc_raw_reads.html"
