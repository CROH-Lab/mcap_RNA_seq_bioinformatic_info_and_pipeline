#!/bin/bash
#SBATCH --job-name=featurecounts
#SBATCH --output=logs/featurecounts_%j.out
#SBATCH --error=logs/featurecounts_%j.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --mail-user=darmstrong4@islander.tamucc.edu
#SBATCH --mail-type=END,FAIL

# ============================================================================
# M. capitata Holobiont RNA-seq Pipeline - Step 6: featureCounts
# ============================================================================

set -e

echo "=============================================="
echo "featureCounts - Read Counting"
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "=============================================="

# Configuration
BASEDIR="/home/darmstrong4/mc_rework"
BAMDIR="${BASEDIR}/05_star_align"
REFDIR="${BASEDIR}/04_reference"
OUTDIR="${BASEDIR}/06_featurecounts"
GFFFILE="${REFDIR}/holobiont_genes.gff3"
THREADS=16

# Setup
cd "${BASEDIR}"
mkdir -p "${OUTDIR}"

# Activate conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate mcap_rnaseq

echo "featureCounts version: $(featureCounts -v 2>&1 | head -1)"
echo "Working directory: $(pwd)"
echo ""

# Check inputs
if [[ ! -f "${GFFFILE}" ]]; then
    echo "ERROR: GFF file not found: ${GFFFILE}"
    exit 1
fi

# Get BAM files
BAMFILES=$(ls ${BAMDIR}/*/[0-9]*_Aligned.sortedByCoord.out.bam | sort)
NUMBAMS=$(echo "${BAMFILES}" | wc -l)
echo "Found ${NUMBAMS} BAM files"

# ============================================================================
# STEP 1: Run featureCounts (skip if already done)
# ============================================================================
RAWCOUNTS="${OUTDIR}/gene_counts_raw.txt"

if [[ -f "${RAWCOUNTS}" ]]; then
    echo ""
    echo "=== Raw counts file exists, skipping featureCounts ==="
    echo "Delete ${RAWCOUNTS} to re-run"
else
    echo ""
    echo "=== Running featureCounts (fractional multi-mapper counting) ==="

    featureCounts \
        -T ${THREADS} \
        -p \
        --countReadPairs \
        -M \
        --fraction \
        -t exon \
        -g Parent \
        -a "${GFFFILE}" \
        -o "${RAWCOUNTS}" \
        ${BAMFILES}

    echo "featureCounts complete!"
fi

# ============================================================================
# STEP 2: Create clean count matrices
# ============================================================================
echo ""
echo "=== Creating clean count matrices ==="

# Clean fractional counts
tail -n +2 "${RAWCOUNTS}" | cut -f1,7- > "${OUTDIR}/gene_counts_fractional_full.txt"

# Simplify column headers (remove paths)
head -1 "${OUTDIR}/gene_counts_fractional_full.txt" | \
    sed 's|/home/darmstrong4/mc_rework/05_star_align/||g' | \
    sed 's|/[^/]*_Aligned.sortedByCoord.out.bam||g' > "${OUTDIR}/gene_counts_fractional.txt"

tail -n +2 "${OUTDIR}/gene_counts_fractional_full.txt" >> "${OUTDIR}/gene_counts_fractional.txt"
rm "${OUTDIR}/gene_counts_fractional_full.txt"

# Create integer version for DESeq2
head -1 "${OUTDIR}/gene_counts_fractional.txt" > "${OUTDIR}/gene_counts.txt"
tail -n +2 "${OUTDIR}/gene_counts_fractional.txt" | \
    awk 'BEGIN {FS=OFS="\t"} {
        printf "%s", $1
        for (i=2; i<=NF; i++) printf "\t%d", int($i + 0.5)
        printf "\n"
    }' >> "${OUTDIR}/gene_counts.txt"

echo "Created: gene_counts.txt (integer counts for DESeq2)"
echo "Created: gene_counts_fractional.txt (original fractional counts)"

# ============================================================================
# STEP 3: Create organism-specific matrices
# ============================================================================
echo ""
echo "=== Creating organism-specific matrices ==="

# Mcap (host) - IDs start with "Montipora_capitata"
head -1 "${OUTDIR}/gene_counts.txt" > "${OUTDIR}/gene_counts_Mcap.txt"
grep "^Montipora_capitata" "${OUTDIR}/gene_counts.txt" >> "${OUTDIR}/gene_counts_Mcap.txt"
MCAP=$(tail -n +2 "${OUTDIR}/gene_counts_Mcap.txt" | wc -l)
echo "Mcap: ${MCAP} transcripts"

# Cgor (Clade C) - IDs start with "evm.model.scf"
head -1 "${OUTDIR}/gene_counts.txt" > "${OUTDIR}/gene_counts_Cgor.txt"
grep "^evm\.model\.scf" "${OUTDIR}/gene_counts.txt" >> "${OUTDIR}/gene_counts_Cgor.txt"
CGOR=$(tail -n +2 "${OUTDIR}/gene_counts_Cgor.txt" | wc -l)
echo "Cgor: ${CGOR} genes"

# Dtre (Clade D) - IDs start with "SCF082"
head -1 "${OUTDIR}/gene_counts.txt" > "${OUTDIR}/gene_counts_Dtre.txt"
grep "^SCF082" "${OUTDIR}/gene_counts.txt" >> "${OUTDIR}/gene_counts_Dtre.txt"
DTRE=$(tail -n +2 "${OUTDIR}/gene_counts_Dtre.txt" | wc -l)
echo "Dtre: ${DTRE} genes"

TOTAL=$((MCAP + CGOR + DTRE))
echo "Total: ${TOTAL} features"

# ============================================================================
# STEP 4: Calculate statistics
# ============================================================================
echo ""
echo "=== Per-sample statistics ==="

# Get sample names
HEADER=$(head -1 "${OUTDIR}/gene_counts.txt")
SAMPLES=$(echo "${HEADER}" | cut -f2-)

echo -e "Sample\tMcap\tCgor\tDtre\tTotal\tPct_Mcap\tPct_Cgor\tPct_Dtre" > "${OUTDIR}/organism_proportions.txt"

NUMCOLS=$(head -1 "${OUTDIR}/gene_counts.txt" | awk '{print NF}')

for i in $(seq 2 ${NUMCOLS}); do
    SAMPLE=$(head -1 "${OUTDIR}/gene_counts.txt" | cut -f${i})

    MCAP_SUM=$(tail -n +2 "${OUTDIR}/gene_counts_Mcap.txt" | cut -f${i} | awk '{s+=$1} END {print s}')
    CGOR_SUM=$(tail -n +2 "${OUTDIR}/gene_counts_Cgor.txt" | cut -f${i} | awk '{s+=$1} END {print s}')
    DTRE_SUM=$(tail -n +2 "${OUTDIR}/gene_counts_Dtre.txt" | cut -f${i} | awk '{s+=$1} END {print s}')

    TOTAL_SUM=$((MCAP_SUM + CGOR_SUM + DTRE_SUM))

    if [[ ${TOTAL_SUM} -gt 0 ]]; then
        PCT_MCAP=$(awk "BEGIN {printf \"%.1f\", ${MCAP_SUM}*100/${TOTAL_SUM}}")
        PCT_CGOR=$(awk "BEGIN {printf \"%.1f\", ${CGOR_SUM}*100/${TOTAL_SUM}}")
        PCT_DTRE=$(awk "BEGIN {printf \"%.1f\", ${DTRE_SUM}*100/${TOTAL_SUM}}")
    else
        PCT_MCAP=0; PCT_CGOR=0; PCT_DTRE=0
    fi

    echo -e "${SAMPLE}\t${MCAP_SUM}\t${CGOR_SUM}\t${DTRE_SUM}\t${TOTAL_SUM}\t${PCT_MCAP}\t${PCT_CGOR}\t${PCT_DTRE}" >> "${OUTDIR}/organism_proportions.txt"
done

echo ""
echo "Organism proportions:"
column -t "${OUTDIR}/organism_proportions.txt"

# ============================================================================
# STEP 5: Assignment rate summary
# ============================================================================
echo ""
echo "=== Assignment rates from featureCounts ==="

SUMMARY="${RAWCOUNTS}.summary"
if [[ -f "${SUMMARY}" ]]; then
    echo -e "Sample\tAssigned\tNoFeatures\tAmbiguous\tUnmapped\tRate" > "${OUTDIR}/assignment_summary.txt"

    for i in $(seq 2 ${NUMCOLS}); do
        SAMPLE=$(head -1 "${OUTDIR}/gene_counts.txt" | cut -f${i})
        ASSIGNED=$(sed -n '2p' "${SUMMARY}" | cut -f${i})
        NOFEAT=$(grep "^Unassigned_NoFeatures" "${SUMMARY}" | cut -f${i})
        AMBIG=$(grep "^Unassigned_Ambiguity" "${SUMMARY}" | cut -f${i})
        UNMAP=$(grep "^Unassigned_Unmapped" "${SUMMARY}" | cut -f${i})

        TOTAL=$((ASSIGNED + NOFEAT + AMBIG + UNMAP))
        RATE=$(awk "BEGIN {printf \"%.1f\", ${ASSIGNED}*100/${TOTAL}}")

        echo -e "${SAMPLE}\t${ASSIGNED}\t${NOFEAT}\t${AMBIG}\t${UNMAP}\t${RATE}%" >> "${OUTDIR}/assignment_summary.txt"
    done

    echo ""
    column -t "${OUTDIR}/assignment_summary.txt"
fi

# ============================================================================
# Final Summary
# ============================================================================
echo ""
echo "=============================================="
echo "featureCounts Processing Complete!"
echo "=============================================="
echo ""
echo "Output files:"
ls -lh "${OUTDIR}/"
echo ""
echo "Count matrices for DESeq2:"
echo "  All (integer):  ${OUTDIR}/gene_counts.txt"
echo "  All (fraction): ${OUTDIR}/gene_counts_fractional.txt"
echo "  Host only:      ${OUTDIR}/gene_counts_Mcap.txt"
echo "  Clade C only:   ${OUTDIR}/gene_counts_Cgor.txt"
echo "  Clade D only:   ${OUTDIR}/gene_counts_Dtre.txt"
echo ""
echo "Matrix dimensions:"
echo "  Features: $(tail -n +2 ${OUTDIR}/gene_counts.txt | wc -l)"
echo "  Samples:  $((NUMCOLS - 1))"
echo ""
echo "End time: $(date)"
