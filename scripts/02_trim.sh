#!/bin/bash
#SBATCH -J mcap_trim                 # Job name
#SBATCH -o logs/trim_%A_%a.out       # Standard output (array-aware)
#SBATCH -e logs/trim_%A_%a.err       # Standard error (array-aware)
#SBATCH --partition=normal           # Use the normal partition
#SBATCH --nodes=1                    # Use a single node per task
#SBATCH --cpus-per-task=8            # Allocate 8 CPUs
#SBATCH --array=1-24                 # Job array for 24 samples
#SBATCH --time=0-06:00:00            # Runtime (6 hours)
#SBATCH --mail-user=darmstrong4@islander.tamucc.edu
#SBATCH --mail-type=BEGIN,END,FAIL   # Notifications for job status

# ============================================================================
# M. capitata RNA-seq Rework Pipeline - Step 2: Adapter Trimming
# Following Drury et al. 2022: Cutadapt to remove poly-A tails and low quality bases
# ============================================================================

set -eo pipefail

# --- Configuration ---
BASE_DIR="/home/darmstrong4/mc_rework"
THREADS=8

cd ${BASE_DIR}

# --- Activate conda environment ---
eval "$(conda shell.bash hook)"
conda activate mcap_rnaseq

# --- Sample array ---
SAMPLES=(1AS 1AW 1BS 1BW 1CS 1CW 1DS 1DW 2AS 2AW 2BS 2BW 2CS 2CW 2DS 2DW 3AS 3AW 3BS 3BW 3CS 3CW 3DS 3DW)
SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "=== Trimming Sample: ${SAMPLE} ==="
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Start time: $(date)"
echo "Cutadapt version: $(cutadapt --version)"

# --- Input/Output files ---
R1_IN="00_raw_fastq/${SAMPLE}_R1_001.fastq.gz"
R2_IN="00_raw_fastq/${SAMPLE}_R2_001.fastq.gz"

mkdir -p 02_trimmed

R1_OUT="02_trimmed/${SAMPLE}_R1_trimmed.fastq.gz"
R2_OUT="02_trimmed/${SAMPLE}_R2_trimmed.fastq.gz"
REPORT="02_trimmed/${SAMPLE}_cutadapt_report.txt"

# --- Verify input files exist ---
if [[ ! -f ${R1_IN} ]] || [[ ! -f ${R2_IN} ]]; then
    echo "ERROR: Input files not found for sample ${SAMPLE}"
    echo "Expected: ${R1_IN} and ${R2_IN}"
    exit 1
fi

# --- Run Cutadapt ---
# Following Drury et al. 2022:
#   - Remove Illumina TruSeq adapters
#   - Remove poly-A tails
#   - Trim low quality bases (Q20)
#   - Discard reads shorter than 50bp

echo "Running Cutadapt..."
cutadapt \
    -j ${THREADS} \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -a "A{10}" \
    -A "T{10}" \
    -q 20,20 \
    --minimum-length 50 \
    --trim-n \
    --max-n 0.1 \
    -o ${R1_OUT} \
    -p ${R2_OUT} \
    ${R1_IN} \
    ${R2_IN} \
    > ${REPORT} 2>&1

# --- Extract key statistics ---
echo ""
echo "=== Cutadapt Summary ==="
grep -E "Total read pairs|Read 1 with adapter|Read 2 with adapter|Pairs written|Total basepairs|Quality-trimmed" ${REPORT} || cat ${REPORT}

# --- Append to summary file (with lock to prevent race condition) ---
{
    flock -x 200
    echo -e "${SAMPLE}\t$(grep "Pairs written" ${REPORT} | grep -oP '[\d,]+' | head -1)" >> 02_trimmed/trimming_summary.txt
} 200>02_trimmed/.summary.lock

echo ""
echo "=== Sample ${SAMPLE} Complete ==="
echo "End time: $(date)"
echo "Output: ${R1_OUT}, ${R2_OUT}"
ls -lh ${R1_OUT} ${R2_OUT}
