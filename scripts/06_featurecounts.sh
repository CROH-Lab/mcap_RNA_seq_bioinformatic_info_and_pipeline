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
#SBATCH --mail-type=BEGIN,END,FAIL

# ============================================================================
# M. capitata Holobiont RNA-seq Pipeline - Step 6: featureCounts
# ============================================================================
#
# PURPOSE:
#   Quantify gene expression from aligned RNA-seq reads using featureCounts.
#   Runs SEPARATELY for each organism with organism-specific multi-mapper
#   handling strategies.
#
# MULTI-MAPPER STRATEGY:
#   - M. capitata (host): DISCARD multi-mappers
#     → High-quality chromosome-level assembly (99.2% BUSCO)
#     → Minimal true gene duplication
#     → Multi-mappers likely represent alignment artifacts
#
#   - C. goreaui (Clade C symbiont): COUNT ALL multi-mappers (-M)
#     → Fragmented assembly (6,843 scaffolds)
#     → Some genes split across scaffolds
#
#   - D. trenchii (Clade D symbiont): COUNT ALL multi-mappers (-M)
#     → Whole genome duplication (WGD) creates real paralogs
#     → Highly fragmented assembly (44,682 scaffolds)
#     → Discarding multi-mappers would lose 30-40% of signal
#     → DESeq2 size factor normalization handles count inflation
#
# INPUT:
#   - BAM files from STAR alignment (05_star_align/*/sample_Aligned.sortedByCoord.out.bam)
#   - Organism-specific GFF3 files (04_reference/*_genes_prefixed.gff3)
#
# OUTPUT:
#   - Mcap_counts.txt: M. capitata gene counts (no multi-mappers)
#   - Cgor_counts.txt: C. goreaui gene counts (with multi-mappers)
#   - Dtre_counts.txt: D. trenchii gene counts (with multi-mappers)
#   - *_counts_raw.txt: Raw featureCounts output with metadata columns
#   - *_counts_raw.txt.summary: Assignment statistics per organism
#
# USAGE:
#   sbatch scripts/06_featurecounts.sh
#
# AUTHOR: David Armstrong
# DATE: January 2025
# ============================================================================

set -e

echo "=============================================="
echo "featureCounts - Organism-Specific Read Counting"
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "=============================================="

# ============================================================================
# Configuration
# ============================================================================
BASEDIR="/home/darmstrong4/mc_rework"
BAMDIR="${BASEDIR}/05_star_align"
REFDIR="${BASEDIR}/04_reference"
OUTDIR="${BASEDIR}/06_featurecounts"
THREADS=16

# Organism-specific GFF files
GFF_MCAP="${REFDIR}/Mcap_genes_prefixed.gff3"
GFF_CGOR="${REFDIR}/Cgor_genes_prefixed.gff3"
GFF_DTRE="${REFDIR}/Dtre_genes_prefixed.gff3"

# ============================================================================
# Setup
# ============================================================================
cd "${BASEDIR}"
mkdir -p "${OUTDIR}"

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate mcap_rnaseq

echo ""
echo "featureCounts version: $(featureCounts -v 2>&1 | head -1)"
echo "Working directory: $(pwd)"
echo ""

# ============================================================================
# Verify input files
# ============================================================================
echo "=== Checking input files ==="

# Check GFF files exist
for gff in "${GFF_MCAP}" "${GFF_CGOR}" "${GFF_DTRE}"; do
    if [[ ! -f "${gff}" ]]; then
        echo "ERROR: GFF file not found: ${gff}"
        exit 1
    fi
    echo "  Found: $(basename ${gff})"
done

# Get BAM files
BAMFILES=$(ls ${BAMDIR}/*/[0-9]*_Aligned.sortedByCoord.out.bam 2>/dev/null | sort)
NUMBAMS=$(echo "${BAMFILES}" | wc -l)

if [[ ${NUMBAMS} -eq 0 ]]; then
    echo "ERROR: No BAM files found in ${BAMDIR}"
    exit 1
fi

echo ""
echo "Found ${NUMBAMS} BAM files:"
echo "${BAMFILES}" | xargs -n1 basename | head -5
echo "..."
echo ""

# ============================================================================
# STEP 1: M. capitata (Host) - NO multi-mappers
# ============================================================================
echo "=============================================="
echo "STEP 1: M. capitata (Host)"
echo "  Strategy: DISCARD multi-mappers (default)"
echo "  Rationale: High-quality assembly, minimal duplication"
echo "=============================================="

MCAP_RAW="${OUTDIR}/Mcap_counts_raw.txt"

featureCounts \
    -T ${THREADS} \
    -p \
    --countReadPairs \
    -t exon \
    -g Parent \
    -a "${GFF_MCAP}" \
    -o "${MCAP_RAW}" \
    ${BAMFILES}

echo ""
echo "M. capitata featureCounts complete"
echo ""

# ============================================================================
# STEP 2: C. goreaui (Clade C) - WITH multi-mappers
# ============================================================================
echo "=============================================="
echo "STEP 2: C. goreaui (Clade C Symbiont)"
echo "  Strategy: COUNT ALL multi-mappers (-M)"
echo "  Rationale: Fragmented assembly, preserve signal"
echo "=============================================="

CGOR_RAW="${OUTDIR}/Cgor_counts_raw.txt"

featureCounts \
    -T ${THREADS} \
    -p \
    --countReadPairs \
    -M \
    -t exon \
    -g Parent \
    -a "${GFF_CGOR}" \
    -o "${CGOR_RAW}" \
    ${BAMFILES}

echo ""
echo "C. goreaui featureCounts complete"
echo ""

# ============================================================================
# STEP 3: D. trenchii (Clade D) - WITH multi-mappers
# ============================================================================
echo "=============================================="
echo "STEP 3: D. trenchii (Clade D Symbiont)"
echo "  Strategy: COUNT ALL multi-mappers (-M)"
echo "  Rationale: Whole genome duplication creates real paralogs"
echo "=============================================="

DTRE_RAW="${OUTDIR}/Dtre_counts_raw.txt"

featureCounts \
    -T ${THREADS} \
    -p \
    --countReadPairs \
    -M \
    -t exon \
    -g Parent \
    -a "${GFF_DTRE}" \
    -o "${DTRE_RAW}" \
    ${BAMFILES}

echo ""
echo "D. trenchii featureCounts complete"
echo ""

# ============================================================================
# STEP 4: Create clean count matrices for DESeq2
# ============================================================================
echo "=============================================="
echo "STEP 4: Creating clean count matrices"
echo "=============================================="

# Function to clean featureCounts output
# Removes metadata columns (Chr, Start, End, Strand, Length)
# Simplifies sample names (removes path and suffix)
clean_counts() {
    local input=$1
    local output=$2
    local organism=$3
    
    echo "  Processing ${organism}..."
    
    # Extract header and simplify sample names
    # featureCounts output: Geneid | Chr | Start | End | Strand | Length | Sample1 | Sample2 | ...
    # We want: Geneid | Sample1 | Sample2 | ...
    
    # Get header, extract Geneid + sample columns (7+), simplify paths
    head -2 "${input}" | tail -1 | cut -f1,7- | \
        sed 's|/home/darmstrong4/mc_rework/05_star_align/||g' | \
        sed 's|/[^/]*_Aligned.sortedByCoord.out.bam||g' > "${output}"
    
    # Get count data (skip first 2 header lines)
    tail -n +3 "${input}" | cut -f1,7- >> "${output}"
    
    # Report
    local ngenes=$(tail -n +2 "${output}" | wc -l)
    local nsamples=$(head -1 "${output}" | awk '{print NF-1}')
    echo "    ${organism}: ${ngenes} genes × ${nsamples} samples"
}

# Clean each organism's counts
clean_counts "${MCAP_RAW}" "${OUTDIR}/Mcap_counts.txt" "M. capitata"
clean_counts "${CGOR_RAW}" "${OUTDIR}/Cgor_counts.txt" "C. goreaui"
clean_counts "${DTRE_RAW}" "${OUTDIR}/Dtre_counts.txt" "D. trenchii"

echo ""

# ============================================================================
# STEP 5: Generate assignment statistics
# ============================================================================
echo "=============================================="
echo "STEP 5: Assignment Statistics"
echo "=============================================="

# Function to parse featureCounts summary
parse_summary() {
    local summary=$1
    local organism=$2
    local outfile=$3
    
    echo ""
    echo "=== ${organism} Assignment Summary ==="
    
    # Get sample names from main output
    local counts_file="${summary%.summary}"
    local header=$(head -2 "${counts_file}" | tail -1)
    
    # Create summary table
    echo -e "Sample\tAssigned\tMultiMapping\tNoFeatures\tAmbiguous\tTotal\tAssign_Rate\tMulti_Rate" > "${outfile}"
    
    # Get number of samples
    local ncols=$(head -1 "${summary}" | awk '{print NF}')
    
    for i in $(seq 2 ${ncols}); do
        # Extract sample name from header
        local sample=$(head -2 "${counts_file}" | tail -1 | cut -f$i | \
            sed 's|/home/darmstrong4/mc_rework/05_star_align/||g' | \
            sed 's|/[^/]*_Aligned.sortedByCoord.out.bam||g')
        
        # Extract counts from summary file
        local assigned=$(grep "^Assigned" "${summary}" | cut -f$i)
        local multi=$(grep "^Unassigned_MultiMapping" "${summary}" | cut -f$i || echo "0")
        local nofeat=$(grep "^Unassigned_NoFeatures" "${summary}" | cut -f$i)
        local ambig=$(grep "^Unassigned_Ambiguity" "${summary}" | cut -f$i)
        local unmap=$(grep "^Unassigned_Unmapped" "${summary}" | cut -f$i)
        
        # Handle missing values
        [[ -z "$multi" ]] && multi=0
        [[ -z "$ambig" ]] && ambig=0
        
        # Calculate total and rates
        local total=$((assigned + multi + nofeat + ambig + unmap))
        local assign_rate=$(awk "BEGIN {printf \"%.1f\", ${assigned}*100/${total}}")
        local multi_rate=$(awk "BEGIN {printf \"%.1f\", ${multi}*100/${total}}")
        
        echo -e "${sample}\t${assigned}\t${multi}\t${nofeat}\t${ambig}\t${total}\t${assign_rate}%\t${multi_rate}%"
    done >> "${outfile}"
    
    # Display summary
    column -t "${outfile}"
    
    # Calculate averages
    local avg_assign=$(tail -n +2 "${outfile}" | awk -F'\t' '{sum+=$7} END {printf "%.1f", sum/NR}')
    local avg_multi=$(tail -n +2 "${outfile}" | awk -F'\t' '{sum+=$8} END {printf "%.1f", sum/NR}')
    echo ""
    echo "  Average assignment rate: ${avg_assign}%"
    echo "  Average multi-mapping rate: ${avg_multi}%"
}

# Parse each organism's summary
parse_summary "${MCAP_RAW}.summary" "M. capitata (no multi-mappers)" "${OUTDIR}/Mcap_assignment_summary.txt"
parse_summary "${CGOR_RAW}.summary" "C. goreaui (with multi-mappers)" "${OUTDIR}/Cgor_assignment_summary.txt"
parse_summary "${DTRE_RAW}.summary" "D. trenchii (with multi-mappers)" "${OUTDIR}/Dtre_assignment_summary.txt"

# ============================================================================
# STEP 6: Calculate organism proportions per sample
# ============================================================================
echo ""
echo "=============================================="
echo "STEP 6: Organism Proportions per Sample"
echo "=============================================="

echo -e "Sample\tMcap_Assigned\tCgor_Assigned\tDtre_Assigned\tTotal\tPct_Mcap\tPct_Cgor\tPct_Dtre" > "${OUTDIR}/organism_proportions.txt"

# Get sample list from Mcap counts
SAMPLES=$(head -1 "${OUTDIR}/Mcap_counts.txt" | cut -f2-)

# Get number of samples
NUMSAMPLES=$(head -1 "${OUTDIR}/Mcap_counts.txt" | awk '{print NF-1}')

for i in $(seq 2 $((NUMSAMPLES + 1))); do
    SAMPLE=$(head -1 "${OUTDIR}/Mcap_counts.txt" | cut -f${i})
    
    # Sum counts for each organism
    MCAP_SUM=$(tail -n +2 "${OUTDIR}/Mcap_counts.txt" | cut -f${i} | awk '{s+=$1} END {print s}')
    CGOR_SUM=$(tail -n +2 "${OUTDIR}/Cgor_counts.txt" | cut -f${i} | awk '{s+=$1} END {print s}')
    DTRE_SUM=$(tail -n +2 "${OUTDIR}/Dtre_counts.txt" | cut -f${i} | awk '{s+=$1} END {print s}')
    
    TOTAL=$((MCAP_SUM + CGOR_SUM + DTRE_SUM))
    
    if [[ ${TOTAL} -gt 0 ]]; then
        PCT_MCAP=$(awk "BEGIN {printf \"%.1f\", ${MCAP_SUM}*100/${TOTAL}}")
        PCT_CGOR=$(awk "BEGIN {printf \"%.1f\", ${CGOR_SUM}*100/${TOTAL}}")
        PCT_DTRE=$(awk "BEGIN {printf \"%.1f\", ${DTRE_SUM}*100/${TOTAL}}")
    else
        PCT_MCAP=0; PCT_CGOR=0; PCT_DTRE=0
    fi
    
    echo -e "${SAMPLE}\t${MCAP_SUM}\t${CGOR_SUM}\t${DTRE_SUM}\t${TOTAL}\t${PCT_MCAP}\t${PCT_CGOR}\t${PCT_DTRE}"
done >> "${OUTDIR}/organism_proportions.txt"

echo "Organism proportions per sample:"
column -t "${OUTDIR}/organism_proportions.txt"

# ============================================================================
# Final Summary
# ============================================================================
echo ""
echo "=============================================="
echo "featureCounts Processing Complete!"
echo "=============================================="
echo ""
echo "Output directory: ${OUTDIR}"
echo ""
echo "Count matrices for DESeq2:"
echo "  M. capitata: ${OUTDIR}/Mcap_counts.txt (no multi-mappers)"
echo "  C. goreaui:  ${OUTDIR}/Cgor_counts.txt (with multi-mappers)"
echo "  D. trenchii: ${OUTDIR}/Dtre_counts.txt (with multi-mappers)"
echo ""
echo "Matrix dimensions:"
for org in Mcap Cgor Dtre; do
    ngenes=$(tail -n +2 "${OUTDIR}/${org}_counts.txt" | wc -l)
    nsamples=$(head -1 "${OUTDIR}/${org}_counts.txt" | awk '{print NF-1}')
    echo "  ${org}: ${ngenes} genes × ${nsamples} samples"
done
echo ""
echo "Assignment statistics:"
ls -lh "${OUTDIR}"/*_assignment_summary.txt
echo ""
echo "NOTE: Symbiont counts include multi-mappers (-M flag)"
echo "      DESeq2 size factors will normalize for inflated library sizes"
echo ""
echo "End time: $(date)"
