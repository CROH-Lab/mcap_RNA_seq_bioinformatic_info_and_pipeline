#!/bin/bash
#SBATCH -J mcap_fastqc_trim          # Job name
#SBATCH -o logs/fastqc_trim_%j.out   # Standard output
#SBATCH -e logs/fastqc_trim_%j.err   # Standard error
#SBATCH --partition=normal           # Use the normal partition
#SBATCH --nodes=1                    # Use a single node
#SBATCH --cpus-per-task=8            # Allocate 8 CPUs
#SBATCH --time=0-04:00:00            # Runtime (4 hours)
#SBATCH --mail-user=darmstrong4@islander.tamucc.edu
#SBATCH --mail-type=BEGIN,END,FAIL   # Notifications for job status

# ============================================================================
# M. capitata RNA-seq Rework Pipeline - Step 3: FastQC on Trimmed Reads
# ============================================================================

set -eo pipefail

# --- Configuration ---
BASE_DIR="/home/darmstrong4/mc_rework"
THREADS=8

cd ${BASE_DIR}

# --- Activate conda environment ---
eval "$(conda shell.bash hook)"
conda activate mcap_rnaseq

echo "=== FastQC on Trimmed Reads ==="
echo "Start time: $(date)"

# --- Run FastQC on trimmed reads ---
mkdir -p 03_qc_trimmed

# Count total files to process
TOTAL_FILES=$(ls 02_trimmed/*_trimmed.fastq.gz 2>/dev/null | wc -l)
echo "Processing ${TOTAL_FILES} trimmed FASTQ files..."

fastqc \
    -t ${THREADS} \
    -o 03_qc_trimmed \
    02_trimmed/*_trimmed.fastq.gz

echo "FastQC complete: $(date)"

# --- Generate comparison MultiQC report ---
echo ""
echo "=== Generating Comparative MultiQC Report ==="
multiqc \
    01_qc_raw \
    02_trimmed \
    03_qc_trimmed \
    -o 03_qc_trimmed \
    -n multiqc_trimming_comparison \
    --title "M. capitata: Raw vs Trimmed Reads Comparison" \
    --force

echo "MultiQC complete: $(date)"

# --- Calculate read retention statistics ---
# Use FastQC data files instead of re-reading FASTQs (MUCH faster)
echo ""
echo "=== Calculating Read Retention Statistics ==="
echo -e "Sample\tRaw_Reads\tTrimmed_Reads\tRetention_Pct" > 03_qc_trimmed/read_retention_stats.txt

for SAMPLE in 1AS 1AW 1BS 1BW 1CS 1CW 1DS 1DW 2AS 2AW 2BS 2BW 2CS 2CW 2DS 2DW 3AS 3AW 3BS 3BW 3CS 3CW 3DS 3DW; do
    echo "Processing ${SAMPLE}..."
    
    # Get raw read count from FastQC data (R1 only, since PE)
    RAW_ZIP="01_qc_raw/${SAMPLE}_R1_001_fastqc.zip"
    if [[ -f "${RAW_ZIP}" ]]; then
        RAW_READS=$(unzip -p ${RAW_ZIP} */fastqc_data.txt | grep "Total Sequences" | cut -f2)
    else
        RAW_READS="NA"
    fi
    
    # Get trimmed read count from FastQC data
    TRIM_ZIP="03_qc_trimmed/${SAMPLE}_R1_trimmed_fastqc.zip"
    if [[ -f "${TRIM_ZIP}" ]]; then
        TRIM_READS=$(unzip -p ${TRIM_ZIP} */fastqc_data.txt | grep "Total Sequences" | cut -f2)
    else
        TRIM_READS="NA"
    fi
    
    # Calculate retention
    if [[ "${RAW_READS}" != "NA" ]] && [[ "${TRIM_READS}" != "NA" ]] && [[ ${RAW_READS} -gt 0 ]]; then
        RETENTION=$(echo "scale=2; ${TRIM_READS} * 100 / ${RAW_READS}" | bc)
    else
        RETENTION="NA"
    fi
    
    echo -e "${SAMPLE}\t${RAW_READS}\t${TRIM_READS}\t${RETENTION}%"
done >> 03_qc_trimmed/read_retention_stats.txt

echo ""
echo "=== Read Retention Summary ==="
column -t 03_qc_trimmed/read_retention_stats.txt

# Calculate overall statistics
echo ""
echo "=== Overall Statistics ==="
awk -F'\t' 'NR>1 && $2!="NA" && $3!="NA" {
    raw += $2
    trim += $3
    n++
}
END {
    if (n > 0) {
        printf "Total samples: %d\n", n
        printf "Total raw reads: %.2f M\n", raw/1000000
        printf "Total trimmed reads: %.2f M\n", trim/1000000
        printf "Overall retention: %.2f%%\n", trim*100/raw
    }
}' 03_qc_trimmed/read_retention_stats.txt

echo ""
echo "=== Post-Trimming QC Complete ==="
echo "End time: $(date)"
echo ""
echo "Output files:"
echo "  FastQC reports: ${BASE_DIR}/03_qc_trimmed/"
echo "  Comparison report: ${BASE_DIR}/03_qc_trimmed/multiqc_trimming_comparison.html"
echo "  Retention stats: ${BASE_DIR}/03_qc_trimmed/read_retention_stats.txt"
