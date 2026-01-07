#!/bin/bash
#SBATCH -J star_index                # Job name
#SBATCH -o logs/star_index_%j.out    # Standard output
#SBATCH -e logs/star_index_%j.err    # Standard error
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G                   # STAR index needs lots of RAM
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=darmstrong4@islander.tamucc.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# ============================================================================
# M. capitata Holobiont RNA-seq Pipeline - Step 4: STAR Index
#
# Creates a combined holobiont reference genome and STAR index containing:
#   - Host: Montipora capitata (prefix: Mcap_)
#   - Symbiont: Cladocopium goreaui (prefix: Cgor_)
#   - Symbiont: Durusdinium trenchii (prefix: Dtre_)
#
# Input:  Individual genome FASTAs and GFF annotations
# Output: Combined holobiont genome, GFF, and STAR index
#
# CHANGES FROM ORIGINAL:
#   1. Added --limitSjdbInsertNsj 2000000 (holobiont has 1.66M splice junctions)
#   2. Added GFF filtering to remove malformed lines (protein seqs, comments)
#   3. Added index verification step to confirm critical files exist
#   4. Force regeneration of combined GFF with proper filtering
# ============================================================================

set -eo pipefail

echo "=============================================="
echo "STAR Index Build - Holobiont Reference"
echo "Date: $(date)"
echo "=============================================="

# --- Configuration ---
BASE_DIR="/home/darmstrong4/mc_rework"
REF_DIR="${BASE_DIR}/04_reference"
THREADS=${SLURM_CPUS_PER_TASK:-16}

cd ${REF_DIR}

# --- Activate conda environment ---
eval "$(conda shell.bash hook)"
conda activate mcap_rnaseq

echo "STAR version: $(STAR --version)"
echo "Threads: ${THREADS}"
echo ""

# ============================================================================
# STEP 1: Decompress and Prepare Genome Files
# ============================================================================
echo "=== STEP 1: Preparing Genome Files ==="

# --- M. capitata (Host) ---
echo "Processing M. capitata genome..."
if [[ ! -f Mcap_genome.fa ]]; then
    gunzip -c Montipora_capitata_HIv3.assembly.fasta.gz > Mcap_genome_raw.fa
    # Add Mcap_ prefix to chromosome/scaffold names
    sed 's/^>/>Mcap_/' Mcap_genome_raw.fa > Mcap_genome.fa
    rm Mcap_genome_raw.fa
fi
MCAP_SCAFFOLDS=$(grep -c "^>" Mcap_genome.fa)
MCAP_SIZE=$(grep -v "^>" Mcap_genome.fa | tr -d '\n' | wc -c)
echo "  Scaffolds: ${MCAP_SCAFFOLDS}"
echo "  Size: $(numfmt --to=iec ${MCAP_SIZE})"

# --- C. goreaui (Clade C Symbiont) ---
echo "Processing C. goreaui genome..."
if [[ ! -f Cgor_genome.fa ]]; then
    gunzip -c Cladocopium_goreaui_genome_fa.gz > Cgor_genome_raw.fa
    # Add Cgor_ prefix to chromosome/scaffold names
    sed 's/^>/>Cgor_/' Cgor_genome_raw.fa > Cgor_genome.fa
    rm Cgor_genome_raw.fa
fi
CGOR_SCAFFOLDS=$(grep -c "^>" Cgor_genome.fa)
CGOR_SIZE=$(grep -v "^>" Cgor_genome.fa | tr -d '\n' | wc -c)
echo "  Scaffolds: ${CGOR_SCAFFOLDS}"
echo "  Size: $(numfmt --to=iec ${CGOR_SIZE})"

# --- D. trenchii (Clade D Symbiont) ---
echo "Processing D. trenchii genome..."
if [[ ! -f Dtre_genome.fa ]]; then
    gunzip -c Dtrenchii_SCF082_ASSEMBLY_fasta.gz > Dtre_genome_raw.fa
    # Add Dtre_ prefix to chromosome/scaffold names
    sed 's/^>/>Dtre_/' Dtre_genome_raw.fa > Dtre_genome.fa
    rm Dtre_genome_raw.fa
fi
DTRE_SCAFFOLDS=$(grep -c "^>" Dtre_genome.fa)
DTRE_SIZE=$(grep -v "^>" Dtre_genome.fa | tr -d '\n' | wc -c)
echo "  Scaffolds: ${DTRE_SCAFFOLDS}"
echo "  Size: $(numfmt --to=iec ${DTRE_SIZE})"

# --- Combine into holobiont genome ---
echo ""
echo "Creating combined holobiont genome..."
cat Mcap_genome.fa Cgor_genome.fa Dtre_genome.fa > holobiont_genome.fa

TOTAL_SCAFFOLDS=$(grep -c "^>" holobiont_genome.fa)
TOTAL_SIZE=$(grep -v "^>" holobiont_genome.fa | tr -d '\n' | wc -c)
echo "Combined holobiont genome:"
echo "  Total scaffolds: ${TOTAL_SCAFFOLDS}"
echo "  Total size: $(numfmt --to=iec ${TOTAL_SIZE})"
echo ""

# ============================================================================
# STEP 2: Prepare GFF Annotation Files
# ============================================================================
# CHANGE: Added filtering to remove malformed lines from GFF files
# The Cladocopium GFF contains embedded protein sequences and comment lines
# that were causing issues when concatenated. Now we filter to keep only
# valid GFF feature lines (gene, mRNA, exon, CDS, UTRs) plus proper headers.
# ============================================================================
echo "=== STEP 2: Preparing GFF Annotations ==="

# Define valid GFF3 feature types to keep
VALID_FEATURES="gene|mRNA|transcript|exon|CDS|five_prime_UTR|three_prime_UTR|start_codon|stop_codon"

# --- M. capitata GFF ---
echo "Processing M. capitata GFF..."
# Force regeneration to ensure clean output
rm -f Mcap_genes_prefixed.gff3
# Add prefix and filter to valid features only
awk -F'\t' -v prefix="Mcap_" -v features="$VALID_FEATURES" '
    BEGIN {OFS="\t"}
    /^##gff-version/ {print; next}
    /^##sequence-region/ {$0 = prefix $0; print; next}
    /^#/ {next}  # Skip other comments
    NF >= 8 && $3 ~ "^("features")$" {$1 = prefix $1; print}
' Mcap.genes.gff3 > Mcap_genes_prefixed.gff3

MCAP_EXONS=$(grep -c $'\texon\t' Mcap_genes_prefixed.gff3)
echo "  Exons: ${MCAP_EXONS}"

# --- C. goreaui GFF ---
# CHANGE: This GFF has embedded protein sequences that need filtering
echo "Processing C. goreaui GFF..."
rm -f Cgor_genes_prefixed.gff3
gunzip -c Cladocopium_goreaui_genes_gff3.gz | \
awk -F'\t' -v prefix="Cgor_" -v features="$VALID_FEATURES" '
    BEGIN {OFS="\t"}
    /^##gff-version/ {print; next}
    /^##sequence-region/ {$0 = prefix $0; print; next}
    /^#/ {next}  # Skip other comments (including #PROT lines)
    NF >= 8 && $3 ~ "^("features")$" {$1 = prefix $1; print}
' > Cgor_genes_prefixed.gff3

CGOR_EXONS=$(grep -c $'\texon\t' Cgor_genes_prefixed.gff3)
echo "  Exons: ${CGOR_EXONS}"

# --- D. trenchii GFF ---
echo "Processing D. trenchii GFF..."
rm -f Dtre_genes_prefixed.gff3
gunzip -c Dtrenchii_SCF082_ANNOT_gff.gz | \
awk -F'\t' -v prefix="Dtre_" -v features="$VALID_FEATURES" '
    BEGIN {OFS="\t"}
    /^##gff-version/ {print; next}
    /^##sequence-region/ {$0 = prefix $0; print; next}
    /^#/ {next}  # Skip other comments
    NF >= 8 && $3 ~ "^("features")$" {$1 = prefix $1; print}
' > Dtre_genes_prefixed.gff3

DTRE_EXONS=$(grep -c $'\texon\t' Dtre_genes_prefixed.gff3)
echo "  Exons: ${DTRE_EXONS}"

# --- Combine GFF files ---
echo ""
echo "Creating combined holobiont GFF..."
# Remove old combined GFF if exists
rm -f holobiont_genes.gff3
cat Mcap_genes_prefixed.gff3 Cgor_genes_prefixed.gff3 Dtre_genes_prefixed.gff3 > holobiont_genes.gff3

# Count features in combined file
TOTAL_GENES=$(grep -c $'\tgene\t' holobiont_genes.gff3 || echo "0")
TOTAL_MRNA=$(grep -c $'\tmRNA\t' holobiont_genes.gff3 || echo "0")
TOTAL_EXONS=$(grep -c $'\texon\t' holobiont_genes.gff3 || echo "0")
TOTAL_CDS=$(grep -c $'\tCDS\t' holobiont_genes.gff3 || echo "0")

echo "Combined holobiont annotation:"
echo "  Total genes: ${TOTAL_GENES}"
echo "  Total mRNA: ${TOTAL_MRNA}"
echo "  Total exons: ${TOTAL_EXONS}"
echo "  Total CDS: ${TOTAL_CDS}"

# Verify GFF is clean
echo ""
echo "Verifying GFF integrity..."
MALFORMED=$(awk -F'\t' '!/^#/ && NF < 8' holobiont_genes.gff3 | wc -l)
if [[ ${MALFORMED} -gt 0 ]]; then
    echo "WARNING: Found ${MALFORMED} malformed lines in GFF!"
else
    echo "  ✓ No malformed lines detected"
fi

# Verify all exons have Parent attributes (required for STAR)
EXONS_NO_PARENT=$(awk -F'\t' '$3 == "exon" && $9 !~ /Parent=/' holobiont_genes.gff3 | wc -l)
if [[ ${EXONS_NO_PARENT} -gt 0 ]]; then
    echo "WARNING: Found ${EXONS_NO_PARENT} exons without Parent attribute!"
else
    echo "  ✓ All exons have Parent attributes"
fi
echo ""

# ============================================================================
# STEP 3: Build STAR Index
# ============================================================================
echo "=== STEP 3: Building STAR Index ==="
echo "This may take 1-3 hours..."
echo ""

# Remove old index if exists
rm -rf star_index_holobiont

# Create index directory
mkdir -p star_index_holobiont

# Calculate genomeSAindexNbases based on genome size
# Formula: min(14, log2(GenomeLength)/2 - 1)
GENOME_LENGTH=${TOTAL_SIZE}
SA_INDEX_NBASES=$(python3 -c "import math; print(min(14, int(math.log2(${GENOME_LENGTH})/2 - 1)))")
echo "Calculated --genomeSAindexNbases: ${SA_INDEX_NBASES}"

# CHANGE: Added --limitSjdbInsertNsj 2000000
# The holobiont has ~1.66 million splice junctions from the combined annotations.
# Default limit is 1,000,000 which causes the index build to fail.
echo "Using --limitSjdbInsertNsj 2000000 (holobiont has ~1.66M splice junctions)"
echo ""

# Build STAR index
STAR --runMode genomeGenerate \
    --runThreadN ${THREADS} \
    --limitGenomeGenerateRAM 42886059616 \
    --genomeDir star_index_holobiont \
    --genomeFastaFiles holobiont_genome.fa \
    --sjdbGTFfile holobiont_genes.gff3 \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbOverhang 149 \
    --genomeSAindexNbases ${SA_INDEX_NBASES} \
    --limitSjdbInsertNsj 2000000

STAR_EXIT=$?

if [[ ${STAR_EXIT} -ne 0 ]]; then
    echo "ERROR: STAR index build failed with exit code ${STAR_EXIT}!"
    exit 1
fi

echo ""
echo "STAR command completed. Verifying index..."

# ============================================================================
# STEP 3b: Verify Index Files
# ============================================================================
# CHANGE: Added verification step to confirm all critical files exist
# Previous run appeared successful but was missing SA, SAindex, Genome files
# ============================================================================
echo ""
echo "=== Verifying STAR Index ==="

CRITICAL_FILES="SA SAindex Genome genomeParameters.txt chrName.txt chrLength.txt sjdbList.out.tab"
INDEX_OK=true

for f in ${CRITICAL_FILES}; do
    if [[ -f "star_index_holobiont/$f" ]]; then
        SIZE=$(ls -lh "star_index_holobiont/$f" | awk '{print $5}')
        echo "  ✓ $f ($SIZE)"
    else
        echo "  ✗ $f MISSING!"
        INDEX_OK=false
    fi
done

if [[ "${INDEX_OK}" != "true" ]]; then
    echo ""
    echo "ERROR: STAR index is incomplete! Critical files are missing."
    echo "Check logs/star_index_*.err for error messages."
    exit 1
fi

# Additional checks
echo ""
echo "Index statistics:"
INDEXED_SCAFFOLDS=$(wc -l < star_index_holobiont/chrName.txt)
INDEXED_SJ=$(wc -l < star_index_holobiont/sjdbList.out.tab)
echo "  Scaffolds indexed: ${INDEXED_SCAFFOLDS}"
echo "  Splice junctions: ${INDEXED_SJ}"

# Verify scaffold counts match
if [[ ${INDEXED_SCAFFOLDS} -ne ${TOTAL_SCAFFOLDS} ]]; then
    echo "WARNING: Scaffold count mismatch! Expected ${TOTAL_SCAFFOLDS}, got ${INDEXED_SCAFFOLDS}"
fi

echo ""
echo "STAR index created successfully!"
echo ""
ls -lh star_index_holobiont/

# ============================================================================
# STEP 4: Create Reference Info File
# ============================================================================
echo ""
echo "=== STEP 4: Creating Reference Documentation ==="

cat > reference_info.txt << EOF
M. capitata Holobiont Reference - STAR Pipeline
================================================

Created: $(date)

GENOME STATISTICS:
------------------
Organism                    Scaffolds       Size
M. capitata (Mcap_)         ${MCAP_SCAFFOLDS}           $(numfmt --to=iec ${MCAP_SIZE})
C. goreaui (Cgor_)          ${CGOR_SCAFFOLDS}           $(numfmt --to=iec ${CGOR_SIZE})
D. trenchii (Dtre_)         ${DTRE_SCAFFOLDS}           $(numfmt --to=iec ${DTRE_SIZE})
-------------------------------------------------
TOTAL                       ${TOTAL_SCAFFOLDS}           $(numfmt --to=iec ${TOTAL_SIZE})

ANNOTATION STATISTICS:
----------------------
Total genes: ${TOTAL_GENES}
Total mRNA: ${TOTAL_MRNA}
Total exons: ${TOTAL_EXONS}
Total CDS: ${TOTAL_CDS}
Splice junctions: ${INDEXED_SJ}

OUTPUT FILES:
-------------
holobiont_genome.fa         Combined genome FASTA (all organisms)
holobiont_genes.gff3        Combined gene annotations (filtered)
star_index_holobiont/       STAR index directory

INDIVIDUAL PREFIXED FILES:
--------------------------
Mcap_genome.fa              M. capitata genome (prefixed)
Mcap_genes_prefixed.gff3    M. capitata annotations (prefixed, filtered)
Cgor_genome.fa              C. goreaui genome (prefixed)
Cgor_genes_prefixed.gff3    C. goreaui annotations (prefixed, filtered)
Dtre_genome.fa              D. trenchii genome (prefixed)
Dtre_genes_prefixed.gff3    D. trenchii annotations (prefixed, filtered)

SOURCE DATA:
------------
1. HOST: Montipora capitata
   Genome: Stephens et al. 2022 (GigaDB 102268)
   Assembly: HIv3 (chromosome-level)
   Citation: Stephens TG et al. GigaScience. 2022;11:giac098

2. SYMBIONT: Cladocopium goreaui (Clade C)
   Strain: SCF055
   Source: Chen et al. 2022 (UQ eSpace)
   Citation: Chen Y et al. Microorganisms. 2022;10(8):1662

3. SYMBIONT: Durusdinium trenchii (Clade D)
   Strain: SCF082
   Source: Dougan et al. 2022 (UQ eSpace)
   Citation: Dougan KE et al. Science Advances. 2024;10:eadn2218

STAR INDEX PARAMETERS:
----------------------
--sjdbOverhang 149                    (read length - 1)
--genomeSAindexNbases ${SA_INDEX_NBASES}
--sjdbGTFtagExonParentTranscript Parent   (for GFF3 format)
--limitSjdbInsertNsj 2000000          (increased for holobiont)

USAGE:
------
# Alignment
STAR --genomeDir star_index_holobiont --readFilesIn R1.fq.gz R2.fq.gz ...

# Counting (use -t exon -g Parent for GFF3)
featureCounts -a holobiont_genes.gff3 -t exon -g Parent -o counts.txt *.bam
EOF

echo "Reference info saved to reference_info.txt"

# ============================================================================
# STEP 5: Cleanup (Optional)
# ============================================================================
echo ""
echo "=== Cleanup ==="
echo "Keeping individual genome files for potential future use."
echo "To save space, you can remove them with:"
echo "  rm Mcap_genome.fa Cgor_genome.fa Dtre_genome.fa"
echo ""

# ============================================================================
# Final Summary
# ============================================================================
echo "=============================================="
echo "STAR Index Build Complete!"
echo "=============================================="
echo ""
echo "Holobiont reference created:"
echo "  Genome: ${REF_DIR}/holobiont_genome.fa"
echo "  Annotation: ${REF_DIR}/holobiont_genes.gff3"
echo "  STAR index: ${REF_DIR}/star_index_holobiont/"
echo ""
echo "Index verification:"
echo "  ✓ All critical files present"
echo "  ✓ ${INDEXED_SCAFFOLDS} scaffolds indexed"
echo "  ✓ ${INDEXED_SJ} splice junctions loaded"
echo ""
echo "Index size:"
du -sh star_index_holobiont/
echo ""
echo "Next step: sbatch scripts/05_star_align.sh"
echo ""
echo "End time: $(date)"
