#!/bin/bash
#SBATCH -J star_align                # Job name
#SBATCH -o logs/star_align_%A_%a.out # Standard output (array-aware)
#SBATCH -e logs/star_align_%A_%a.err # Standard error (array-aware)
#SBATCH --partition=normal           # Normal partition (350GB available)
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --mem=320G
#SBATCH --array=1-24                 # 24 samples
#SBATCH --time=12:00:00
#SBATCH --mail-user=darmstrong4@islander.tamucc.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# ============================================================================
# M. capitata Holobiont RNA-seq Pipeline - Step 5: STAR Alignment
#
# Aligns trimmed reads to the holobiont genome using STAR.
#
# Strategy:
#   - Uses normal partition with highmem QOS for 120GB RAM
#   - Limits concurrent jobs to 2 to avoid memory contention
#   - STAR outputs UNSORTED BAM (avoids STAR's internal sorter issues)
#   - samtools sorts separately (more robust and memory-efficient)
#   - Skips samples that have already completed successfully
#
# Input:  Trimmed FASTQ files (02_trimmed/)
# Output: Sorted BAM files (05_star_align/)
#
# Expected mapping rates for holobiont:
#   - Uniquely mapped: 60-75%
#   - Multi-mapped: 15-25%
#   - Unmapped: 10-20%
# ============================================================================

set -eo pipefail

# --- Configuration ---
BASE_DIR="/home/darmstrong4/mc_rework"
TRIMMED_DIR="${BASE_DIR}/02_trimmed"
STAR_INDEX="${BASE_DIR}/04_reference/star_index_holobiont"
OUTPUT_DIR="${BASE_DIR}/05_star_align"
THREADS=${SLURM_CPUS_PER_TASK:-21}

cd ${BASE_DIR}

# --- Activate conda environment ---
eval "$(conda shell.bash hook)"
conda activate mcap_rnaseq

# --- Sample array (24 samples) ---
SAMPLES=(
    1AS 1AW 1BS 1BW 1CS 1CW 1DS 1DW
    2AS 2AW 2BS 2BW 2CS 2CW 2DS 2DW
    3AS 3AW 3BS 3BW 3CS 3CW 3DS 3DW
)
SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}

# --- Output paths ---
OUTDIR="${OUTPUT_DIR}/${SAMPLE}"
FINAL_BAM="${OUTDIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
UNSORTED_BAM="${OUTDIR}/${SAMPLE}_Aligned.out.bam"
LOG_FINAL="${OUTDIR}/${SAMPLE}_Log.final.out"

echo "=============================================="
echo "STAR Alignment: ${SAMPLE}"
echo "=============================================="
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Start time: $(date)"
echo "Node: $(hostname)"
echo "STAR version: $(STAR --version)"
echo "samtools version: $(samtools --version | head -1)"
echo "Threads: ${THREADS}"
echo "Memory requested: 96GB"
echo ""

# --- Check if already completed ---
if [[ -f "${FINAL_BAM}" ]]; then
    SIZE=$(stat -c%s "${FINAL_BAM}" 2>/dev/null || echo "0")
    if [[ ${SIZE} -gt 100000000 ]]; then
        echo "=============================================="
        echo "SKIPPING: ${SAMPLE} already completed"
        echo "BAM file: ${FINAL_BAM}"
        echo "Size: $(ls -lh ${FINAL_BAM} | awk '{print $5}')"
        echo "=============================================="
        
        # Still print stats if available
        if [[ -f "${LOG_FINAL}" ]]; then
            echo ""
            echo "Previous alignment statistics:"
            grep -E "Uniquely mapped reads %" ${LOG_FINAL} || true
        fi
        exit 0
    else
        echo "Found incomplete BAM ($(ls -lh ${FINAL_BAM} 2>/dev/null | awk '{print $5}' || echo '0')), will re-run..."
    fi
fi

# --- Input files ---
R1="${TRIMMED_DIR}/${SAMPLE}_R1_trimmed.fastq.gz"
R2="${TRIMMED_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"

# --- Verify inputs exist ---
if [[ ! -f "${R1}" ]] || [[ ! -f "${R2}" ]]; then
    echo "ERROR: Trimmed FASTQ files not found for ${SAMPLE}"
    echo "  R1: ${R1}"
    echo "  R2: ${R2}"
    exit 1
fi

if [[ ! -d "${STAR_INDEX}" ]]; then
    echo "ERROR: STAR index not found at ${STAR_INDEX}"
    exit 1
fi

echo "Input files:"
echo "  R1: ${R1} ($(ls -lh ${R1} | awk '{print $5}'))"
echo "  R2: ${R2} ($(ls -lh ${R2} | awk '{print $5}'))"
echo "  Index: ${STAR_INDEX}"
echo ""

# --- Clean up any previous failed attempt ---
rm -rf ${OUTDIR}
mkdir -p ${OUTDIR}

# ============================================================================
# STEP 1: STAR Alignment (output UNSORTED BAM)
# ============================================================================
echo "=== Step 1: STAR Alignment ==="
echo "Outputting unsorted BAM (sorting handled by samtools)"
echo ""

STAR --runMode alignReads \
    --runThreadN ${THREADS} \
    --genomeDir ${STAR_INDEX} \
    --readFilesIn ${R1} ${R2} \
    --readFilesCommand zcat \
    --outFileNamePrefix ${OUTDIR}/${SAMPLE}_ \
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

STAR_EXIT=$?

if [[ ${STAR_EXIT} -ne 0 ]]; then
    echo "ERROR: STAR alignment failed for ${SAMPLE} (exit code: ${STAR_EXIT})"
    exit 1
fi

echo ""
echo "STAR alignment complete!"
echo "Unsorted BAM: $(ls -lh ${UNSORTED_BAM} | awk '{print $5}')"
echo ""

# ============================================================================
# STEP 2: Sort BAM with samtools
# ============================================================================
echo "=== Step 2: Sorting BAM with samtools ==="

# Use 4GB per thread for sorting (16 threads × 4GB = 64GB max)
samtools sort \
    -@ ${THREADS} \
    -m 4G \
    -T ${OUTDIR}/${SAMPLE}_sort_tmp \
    -o ${FINAL_BAM} \
    ${UNSORTED_BAM}

SORT_EXIT=$?

if [[ ${SORT_EXIT} -ne 0 ]]; then
    echo "ERROR: samtools sort failed for ${SAMPLE} (exit code: ${SORT_EXIT})"
    exit 1
fi

echo "Sorting complete!"
echo "Sorted BAM: $(ls -lh ${FINAL_BAM} | awk '{print $5}')"
echo ""

# ============================================================================
# STEP 3: Index BAM
# ============================================================================
echo "=== Step 3: Indexing BAM ==="

samtools index -@ ${THREADS} ${FINAL_BAM}

echo "BAM index created: ${FINAL_BAM}.bai"
echo ""

# ============================================================================
# STEP 4: Clean up
# ============================================================================
echo "=== Step 4: Cleaning up ==="

# Remove unsorted BAM to save space
rm -f ${UNSORTED_BAM}
echo "Removed unsorted BAM"

# Remove STAR temp directory
rm -rf ${OUTDIR}/${SAMPLE}__STARtmp
echo "Removed STAR temp directory"
echo ""

# ============================================================================
# STEP 5: Extract and display statistics
# ============================================================================
echo "=== Alignment Statistics ==="

if [[ -f "${LOG_FINAL}" ]]; then
    echo "Sample: ${SAMPLE}"
    echo ""
    
    # Display key metrics
    grep -E "Number of input reads|Uniquely mapped reads|% of reads mapped|Mapping speed" ${LOG_FINAL}
    echo ""
    
    # Extract metrics for summary
    TOTAL_READS=$(grep "Number of input reads" ${LOG_FINAL} | awk '{print $NF}')
    UNIQUE_READS=$(grep "Uniquely mapped reads number" ${LOG_FINAL} | awk '{print $NF}')
    UNIQUE_PCT=$(grep "Uniquely mapped reads %" ${LOG_FINAL} | awk '{print $NF}')
    MULTI_PCT=$(grep "% of reads mapped to multiple loci" ${LOG_FINAL} | awk '{print $NF}')
    UNMAPPED_SHORT=$(grep "% of reads unmapped: too short" ${LOG_FINAL} | awk '{print $NF}')
    UNMAPPED_OTHER=$(grep "% of reads unmapped: other" ${LOG_FINAL} | awk '{print $NF}')
    
    echo "Summary:"
    echo "  Total reads: ${TOTAL_READS}"
    echo "  Uniquely mapped: ${UNIQUE_READS} (${UNIQUE_PCT})"
    echo "  Multi-mapped: ${MULTI_PCT}"
    echo "  Unmapped (too short): ${UNMAPPED_SHORT}"
    echo "  Unmapped (other): ${UNMAPPED_OTHER}"
    
    # Append to summary file (with file locking for concurrent jobs)
    SUMMARY_FILE="${OUTPUT_DIR}/alignment_summary.tsv"
    {
        flock -x 200
        if [[ ! -f "${SUMMARY_FILE}" ]]; then
            echo -e "Sample\tTotal_Reads\tUnique_Mapped\tUnique_Pct\tMulti_Pct\tUnmapped_Short\tUnmapped_Other" > ${SUMMARY_FILE}
        fi
        # Remove any existing entry for this sample (in case of re-run)
        grep -v "^${SAMPLE}	" ${SUMMARY_FILE} > ${SUMMARY_FILE}.tmp || true
        mv ${SUMMARY_FILE}.tmp ${SUMMARY_FILE}
        echo -e "${SAMPLE}\t${TOTAL_READS}\t${UNIQUE_READS}\t${UNIQUE_PCT}\t${MULTI_PCT}\t${UNMAPPED_SHORT}\t${UNMAPPED_OTHER}" >> ${SUMMARY_FILE}
    } 200>${OUTPUT_DIR}/.alignment.lock
    
else
    echo "WARNING: Log file not found at ${LOG_FINAL}"
fi

# ============================================================================
# Final output summary
# ============================================================================
echo ""
echo "=== Output Files ==="
ls -lh ${OUTDIR}/

echo ""
echo "=============================================="
echo "Sample ${SAMPLE} Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Final BAM: ${FINAL_BAM}"
echo "BAM Index: ${FINAL_BAM}.bai"
echo "Log file: ${LOG_FINAL}"
