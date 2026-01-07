#!/bin/bash
# ============================================================================
# M. capitata Holobiont RNA-seq Pipeline - Step 4 QC: STAR Index Verification
#
# This script performs comprehensive quality control checks on the STAR index
# and associated reference files to ensure everything is ready for alignment.
#
# Checks performed:
#   1. Critical index files present and sized appropriately
#   2. Genome parameters verification
#   3. GFF annotation integrity
#   4. Scaffold consistency between genome, GFF, and index
#   5. Organism distribution (host + symbionts)
#   6. Sample entries from each organism
#   7. Index build log review
#
# Usage: bash scripts/04_qc_star_index.sh
# ============================================================================

set -eo pipefail

# --- Configuration ---
BASE_DIR="/home/darmstrong4/mc_rework"
REF_DIR="${BASE_DIR}/04_reference"
INDEX_DIR="${REF_DIR}/star_index_holobiont"
GFF_FILE="${REF_DIR}/holobiont_genes.gff3"
GENOME_FILE="${REF_DIR}/holobiont_genome.fa"

cd ${REF_DIR}

echo "=============================================="
echo "STAR INDEX QC CHECKS"
echo "Date: $(date)"
echo "=============================================="

# Track pass/fail status
QC_PASS=true

# ============================================
# 1. CRITICAL INDEX FILES
# ============================================
echo ""
echo "=== 1. CRITICAL INDEX FILES ==="
echo ""

CRITICAL_FILES="SA SAindex Genome genomeParameters.txt chrName.txt chrLength.txt sjdbList.out.tab exonInfo.tab"
for f in ${CRITICAL_FILES}; do
    if [[ -f "${INDEX_DIR}/$f" ]]; then
        SIZE=$(ls -lh "${INDEX_DIR}/$f" | awk '{print $5}')
        echo "  ✓ $f ($SIZE)"
    else
        echo "  ✗ $f MISSING!"
        QC_PASS=false
    fi
done

# Check minimum file sizes for critical files
SA_SIZE=$(stat -c%s "${INDEX_DIR}/SA" 2>/dev/null || echo "0")
GENOME_SIZE=$(stat -c%s "${INDEX_DIR}/Genome" 2>/dev/null || echo "0")
SAINDEX_SIZE=$(stat -c%s "${INDEX_DIR}/SAindex" 2>/dev/null || echo "0")

echo ""
echo "Size validation:"
if [[ ${SA_SIZE} -gt 1000000000 ]]; then
    echo "  ✓ SA file size OK (>1GB)"
else
    echo "  ✗ SA file too small - index may be corrupt"
    QC_PASS=false
fi

if [[ ${GENOME_SIZE} -gt 1000000000 ]]; then
    echo "  ✓ Genome file size OK (>1GB)"
else
    echo "  ✗ Genome file too small - index may be corrupt"
    QC_PASS=false
fi

if [[ ${SAINDEX_SIZE} -gt 100000000 ]]; then
    echo "  ✓ SAindex file size OK (>100MB)"
else
    echo "  ✗ SAindex file too small - index may be corrupt"
    QC_PASS=false
fi

# ============================================
# 2. GENOME PARAMETERS
# ============================================
echo ""
echo "=== 2. GENOME PARAMETERS ==="
echo ""
cat ${INDEX_DIR}/genomeParameters.txt

# Verify key parameters
echo ""
echo "Parameter validation:"
if grep -q "sjdbOverhang 149" ${INDEX_DIR}/genomeParameters.txt; then
    echo "  ✓ sjdbOverhang = 149 (correct for 150bp reads)"
else
    echo "  ⚠ sjdbOverhang may not be optimal for your read length"
fi

if grep -q "sjdbGTFtagExonParentTranscript Parent" ${INDEX_DIR}/genomeParameters.txt; then
    echo "  ✓ GFF3 Parent attribute configured correctly"
else
    echo "  ⚠ GFF3 Parent attribute setting not found"
fi

# ============================================
# 3. GFF INTEGRITY CHECK
# ============================================
echo ""
echo "=== 3. GFF INTEGRITY CHECK ==="
echo ""

# Count lines by type
TOTAL_LINES=$(wc -l < ${GFF_FILE})
HEADER_LINES=$(grep -c '^##' ${GFF_FILE} || echo "0")
COMMENT_LINES=$(grep -c '^#' ${GFF_FILE} || echo "0")
VALID_FEATURES=$(awk -F'\t' 'NF >= 8' ${GFF_FILE} | wc -l)

echo "Line counts:"
echo "  Total lines: ${TOTAL_LINES}"
echo "  Header lines (##): ${HEADER_LINES}"
echo "  All comment lines (#): ${COMMENT_LINES}"
echo "  Valid feature lines (8+ columns): ${VALID_FEATURES}"

# Check for malformed lines
MALFORMED=$(awk -F'\t' '!/^#/ && NF > 0 && NF < 8' ${GFF_FILE} | wc -l)
echo ""
if [[ ${MALFORMED} -gt 0 ]]; then
    echo "  ⚠ Found ${MALFORMED} lines with <8 columns (excluding headers)"
    echo "  First 3 examples:"
    awk -F'\t' '!/^#/ && NF > 0 && NF < 8 {print "    "$0}' ${GFF_FILE} | head -3
else
    echo "  ✓ No malformed feature lines detected"
fi

echo ""
echo "Feature types in GFF:"
cut -f3 ${GFF_FILE} | grep -v "^#" | sort | uniq -c | sort -rn | head -10

# Check Parent attributes
EXONS_TOTAL=$(awk -F'\t' '$3 == "exon"' ${GFF_FILE} | wc -l)
EXONS_WITH_PARENT=$(grep -c 'exon.*Parent=' ${GFF_FILE} || echo "0")
EXONS_NO_PARENT=$((EXONS_TOTAL - EXONS_WITH_PARENT))

echo ""
echo "Exon Parent attribute check:"
echo "  Total exons: ${EXONS_TOTAL}"
echo "  Exons with Parent: ${EXONS_WITH_PARENT}"
if [[ ${EXONS_NO_PARENT} -gt 0 ]]; then
    echo "  ✗ Exons without Parent: ${EXONS_NO_PARENT}"
    QC_PASS=false
else
    echo "  ✓ All exons have Parent attributes"
fi

# ============================================
# 4. SCAFFOLD CONSISTENCY CHECK
# ============================================
echo ""
echo "=== 4. SCAFFOLD CONSISTENCY CHECK ==="
echo ""

# Get scaffolds from each source
cut -f1 ${GFF_FILE} | grep -v "^#" | sort -u > /tmp/gff_scaffolds.txt
grep "^>" ${GENOME_FILE} | sed 's/>//' | sort -u > /tmp/genome_scaffolds.txt
sort ${INDEX_DIR}/chrName.txt > /tmp/index_scaffolds.txt

GFF_SCAFFOLDS=$(wc -l < /tmp/gff_scaffolds.txt)
GENOME_SCAFFOLDS=$(wc -l < /tmp/genome_scaffolds.txt)
INDEX_SCAFFOLDS=$(wc -l < /tmp/index_scaffolds.txt)

echo "Scaffold counts:"
echo "  In genome FASTA: ${GENOME_SCAFFOLDS}"
echo "  In GFF annotations: ${GFF_SCAFFOLDS}"
echo "  In STAR index: ${INDEX_SCAFFOLDS}"

# Check genome vs index match
if [[ ${GENOME_SCAFFOLDS} -eq ${INDEX_SCAFFOLDS} ]]; then
    echo "  ✓ Genome and index scaffold counts match"
else
    echo "  ✗ Mismatch between genome (${GENOME_SCAFFOLDS}) and index (${INDEX_SCAFFOLDS})"
    QC_PASS=false
fi

# Check for scaffolds in GFF but not in index (problematic)
GFF_NOT_IN_INDEX=$(comm -23 /tmp/gff_scaffolds.txt /tmp/index_scaffolds.txt | wc -l)
echo ""
echo "Scaffolds in GFF but NOT in index:"
if [[ ${GFF_NOT_IN_INDEX} -gt 0 ]]; then
    echo "  ✗ Found ${GFF_NOT_IN_INDEX} scaffolds in GFF missing from index!"
    echo "  Examples:"
    comm -23 /tmp/gff_scaffolds.txt /tmp/index_scaffolds.txt | head -5 | sed 's/^/    /'
    QC_PASS=false
else
    echo "  ✓ All GFF scaffolds present in index (0 missing)"
fi

# Check for scaffolds in index but not in GFF (okay - unannotated)
INDEX_NOT_IN_GFF=$(comm -13 /tmp/gff_scaffolds.txt /tmp/index_scaffolds.txt | wc -l)
echo ""
echo "Scaffolds in index but NOT in GFF (unannotated - OK):"
echo "  Count: ${INDEX_NOT_IN_GFF}"

# ============================================
# 5. ORGANISM DISTRIBUTION
# ============================================
echo ""
echo "=== 5. ORGANISM DISTRIBUTION ==="
echo ""

# Scaffolds per organism
MCAP_SCAFFOLDS=$(grep -c '^Mcap_' /tmp/index_scaffolds.txt)
CGOR_SCAFFOLDS=$(grep -c '^Cgor_' /tmp/index_scaffolds.txt)
DTRE_SCAFFOLDS=$(grep -c '^Dtre_' /tmp/index_scaffolds.txt)

echo "Scaffolds per organism:"
echo "  Mcap (M. capitata - host):     ${MCAP_SCAFFOLDS}"
echo "  Cgor (C. goreaui - Clade C):   ${CGOR_SCAFFOLDS}"
echo "  Dtre (D. trenchii - Clade D):  ${DTRE_SCAFFOLDS}"

# Exons per organism
MCAP_EXONS=$(awk -F'\t' '$1 ~ /^Mcap_/ && $3 == "exon"' ${GFF_FILE} | wc -l)
CGOR_EXONS=$(awk -F'\t' '$1 ~ /^Cgor_/ && $3 == "exon"' ${GFF_FILE} | wc -l)
DTRE_EXONS=$(awk -F'\t' '$1 ~ /^Dtre_/ && $3 == "exon"' ${GFF_FILE} | wc -l)

echo ""
echo "Exons per organism:"
echo "  Mcap: ${MCAP_EXONS}"
echo "  Cgor: ${CGOR_EXONS}"
echo "  Dtre: ${DTRE_EXONS}"

# Splice junctions per organism
MCAP_SJ=$(grep -c '^Mcap_' ${INDEX_DIR}/sjdbList.out.tab)
CGOR_SJ=$(grep -c '^Cgor_' ${INDEX_DIR}/sjdbList.out.tab)
DTRE_SJ=$(grep -c '^Dtre_' ${INDEX_DIR}/sjdbList.out.tab)
TOTAL_SJ=$((MCAP_SJ + CGOR_SJ + DTRE_SJ))

echo ""
echo "Splice junctions per organism:"
echo "  Mcap: ${MCAP_SJ}"
echo "  Cgor: ${CGOR_SJ}"
echo "  Dtre: ${DTRE_SJ}"
echo "  Total: ${TOTAL_SJ}"

# Genes per organism
MCAP_GENES=$(awk -F'\t' '$1 ~ /^Mcap_/ && $3 == "gene"' ${GFF_FILE} | wc -l)
CGOR_GENES=$(awk -F'\t' '$1 ~ /^Cgor_/ && $3 == "gene"' ${GFF_FILE} | wc -l)
DTRE_GENES=$(awk -F'\t' '$1 ~ /^Dtre_/ && $3 == "gene"' ${GFF_FILE} | wc -l)

echo ""
echo "Genes per organism:"
echo "  Mcap: ${MCAP_GENES}"
echo "  Cgor: ${CGOR_GENES}"
echo "  Dtre: ${DTRE_GENES}"

# Verify all organisms are represented
echo ""
if [[ ${MCAP_SCAFFOLDS} -gt 0 && ${CGOR_SCAFFOLDS} -gt 0 && ${DTRE_SCAFFOLDS} -gt 0 ]]; then
    echo "  ✓ All three organisms represented in index"
else
    echo "  ✗ Missing organism(s) in index!"
    QC_PASS=false
fi

# ============================================
# 6. SAMPLE ENTRIES CHECK
# ============================================
echo ""
echo "=== 6. SAMPLE ENTRIES FROM EACH ORGANISM ==="
echo ""

echo "Sample Mcap exon:"
awk -F'\t' '$1 ~ /^Mcap_/ && $3 == "exon"' ${GFF_FILE} | head -1
echo ""

echo "Sample Cgor exon:"
awk -F'\t' '$1 ~ /^Cgor_/ && $3 == "exon"' ${GFF_FILE} | head -1
echo ""

echo "Sample Dtre exon:"
awk -F'\t' '$1 ~ /^Dtre_/ && $3 == "exon"' ${GFF_FILE} | head -1
echo ""

echo "Sample splice junctions from index:"
echo "  Mcap: $(grep '^Mcap_' ${INDEX_DIR}/sjdbList.out.tab | head -1)"
echo "  Cgor: $(grep '^Cgor_' ${INDEX_DIR}/sjdbList.out.tab | head -1)"
echo "  Dtre: $(grep '^Dtre_' ${INDEX_DIR}/sjdbList.out.tab | head -1)"

# ============================================
# 7. INDEX LOG CHECK
# ============================================
echo ""
echo "=== 7. INDEX BUILD LOG CHECK ==="
echo ""

LOG_FILE="${INDEX_DIR}/Log.out"

if [[ -f ${LOG_FILE} ]]; then
    # Check for fatal errors
    FATAL_ERRORS=$(grep -ci "fatal\|error" ${LOG_FILE} || echo "0")
    # Exclude "SOLUTION" lines which contain "error" but aren't errors
    REAL_ERRORS=$(grep -i "fatal\|error" ${LOG_FILE} | grep -cv "SOLUTION" || echo "0")
    
    if [[ ${REAL_ERRORS} -gt 0 ]]; then
        echo "  ✗ Found ${REAL_ERRORS} error(s) in log:"
        grep -i "fatal\|error" ${LOG_FILE} | grep -v "SOLUTION" | head -5
        QC_PASS=false
    else
        echo "  ✓ No fatal errors in build log"
    fi
    
    # Count warnings
    WARNING_COUNT=$(grep -c "WARNING" ${LOG_FILE} || echo "0")
    echo "  Warning count: ${WARNING_COUNT}"
    
    # Long repeat warnings are expected
    LONG_REPEAT=$(grep -c "long repeat" ${LOG_FILE} || echo "0")
    if [[ ${LONG_REPEAT} -gt 0 ]]; then
        echo "  Long repeat warnings: ${LONG_REPEAT} (expected for complex genomes)"
    fi
    
    # Check for successful completion
    if grep -q "finished successfully" ${LOG_FILE}; then
        echo "  ✓ Index build completed successfully"
    else
        echo "  ✗ 'finished successfully' message not found in log"
        QC_PASS=false
    fi
else
    echo "  ⚠ Log file not found: ${LOG_FILE}"
fi

# ============================================
# CLEANUP
# ============================================
rm -f /tmp/gff_scaffolds.txt /tmp/genome_scaffolds.txt /tmp/index_scaffolds.txt

# ============================================
# FINAL SUMMARY
# ============================================
echo ""
echo "=============================================="
echo "QC SUMMARY"
echo "=============================================="
echo ""
echo "Index directory: ${INDEX_DIR}"
echo "Index size: $(du -sh ${INDEX_DIR}/ | cut -f1)"
echo ""
echo "Holobiont composition:"
echo "  Total scaffolds: ${INDEX_SCAFFOLDS}"
echo "  Total splice junctions: ${TOTAL_SJ}"
echo "  Total exons: ${EXONS_TOTAL}"
echo ""
echo "Organism breakdown:"
printf "  %-35s %10s %10s %10s\n" "" "Mcap" "Cgor" "Dtre"
printf "  %-35s %10s %10s %10s\n" "Scaffolds:" "${MCAP_SCAFFOLDS}" "${CGOR_SCAFFOLDS}" "${DTRE_SCAFFOLDS}"
printf "  %-35s %10s %10s %10s\n" "Genes:" "${MCAP_GENES}" "${CGOR_GENES}" "${DTRE_GENES}"
printf "  %-35s %10s %10s %10s\n" "Exons:" "${MCAP_EXONS}" "${CGOR_EXONS}" "${DTRE_EXONS}"
printf "  %-35s %10s %10s %10s\n" "Splice junctions:" "${MCAP_SJ}" "${CGOR_SJ}" "${DTRE_SJ}"
echo ""

if [[ "${QC_PASS}" == "true" ]]; then
    echo "=============================================="
    echo "✓ ALL QC CHECKS PASSED"
    echo "=============================================="
    echo ""
    echo "Index is ready for alignment!"
    echo "Next step: sbatch scripts/05_star_align.sh"
    exit 0
else
    echo "=============================================="
    echo "✗ SOME QC CHECKS FAILED"
    echo "=============================================="
    echo ""
    echo "Review the issues above before proceeding."
    exit 1
fi
