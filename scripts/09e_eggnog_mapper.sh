#!/bin/bash
#SBATCH --job-name=eggnog_Dtre
#SBATCH --output=/home/darmstrong4/mc_rework/logs/eggnog_Dtre_%j.out
#SBATCH --error=/home/darmstrong4/mc_rework/logs/eggnog_Dtre_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-user=darmstrong4@islander.tamucc.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Load conda environment
source ~/.bashrc
conda activate mcap_rnaseq

# Directories
WORK_DIR="/home/darmstrong4/mc_rework/09_symbiont_deg_annotation"
DB_DIR="/home/darmstrong4/mc_rework/09_symbiont_deg_annotation/databases/eggnog"
QUERY="${WORK_DIR}/sequences/expressed_proteins.faa"
OUT_PREFIX="Dtre_expressed"

# CORRECT base URL (the bug uses eggnogdb.embl.de, but it should be eggnog5.embl.de)
BASE_URL="http://eggnog5.embl.de/download/emapperdb-5.0.2"

echo "============================================"
echo "eggNOG-mapper Annotation Pipeline"
echo "============================================"
echo "Query: ${QUERY}"
echo "Database dir: ${DB_DIR}"
echo "Output prefix: ${OUT_PREFIX}"
echo "Threads: ${SLURM_CPUS_PER_TASK}"
echo "Start time: $(date)"
echo "============================================"

# Count input sequences
echo "Input sequences: $(grep -c '>' ${QUERY})"
echo ""

# ============================================
# Step 1: Clean up failed downloads
# ============================================
echo "Cleaning up any failed/empty downloads..."
rm -f ${DB_DIR}/*.gz 2>/dev/null
rm -f ${DB_DIR}/eggnog.db 2>/dev/null
rm -f ${DB_DIR}/eggnog_proteins.dmnd 2>/dev/null
mkdir -p ${DB_DIR}

# ============================================
# Step 2: Manual database download (correct URLs)
# ============================================
echo "============================================"
echo "Downloading eggNOG databases (manual method)"
echo "Using correct URL: ${BASE_URL}"
echo "============================================"

cd ${DB_DIR}

# Download annotation database
if [ ! -f "${DB_DIR}/eggnog.db" ]; then
    echo "Downloading eggnog.db (~7GB compressed)..."
    wget -c ${BASE_URL}/eggnog.db.gz
    if [ $? -eq 0 ] && [ -s eggnog.db.gz ]; then
        echo "Decompressing eggnog.db.gz..."
        gunzip -f eggnog.db.gz
    else
        echo "ERROR: Failed to download eggnog.db.gz"
        exit 1
    fi
else
    echo "eggnog.db already exists"
fi

# Download taxa database
if [ ! -f "${DB_DIR}/eggnog.taxa.db" ]; then
    echo "Downloading eggnog.taxa.tar.gz..."
    wget -c ${BASE_URL}/eggnog.taxa.tar.gz
    if [ $? -eq 0 ] && [ -s eggnog.taxa.tar.gz ]; then
        echo "Extracting eggnog.taxa.tar.gz..."
        tar -zxf eggnog.taxa.tar.gz
        rm eggnog.taxa.tar.gz
    else
        echo "ERROR: Failed to download eggnog.taxa.tar.gz"
        exit 1
    fi
else
    echo "eggnog.taxa.db already exists"
fi

# Download diamond database
if [ ! -f "${DB_DIR}/eggnog_proteins.dmnd" ]; then
    echo "Downloading eggnog_proteins.dmnd.gz (~4GB)..."
    wget -c ${BASE_URL}/eggnog_proteins.dmnd.gz
    if [ $? -eq 0 ] && [ -s eggnog_proteins.dmnd.gz ]; then
        echo "Decompressing eggnog_proteins.dmnd.gz..."
        gunzip -f eggnog_proteins.dmnd.gz
    else
        echo "ERROR: Failed to download eggnog_proteins.dmnd.gz"
        exit 1
    fi
else
    echo "eggnog_proteins.dmnd already exists"
fi

# Verify downloads
echo ""
echo "Verifying database files..."
ls -lh ${DB_DIR}/

if [ ! -f "${DB_DIR}/eggnog.db" ] || [ ! -f "${DB_DIR}/eggnog_proteins.dmnd" ]; then
    echo "ERROR: Database files missing after download"
    exit 1
fi

echo "Database download complete!"
echo ""

# ============================================
# Step 3: Run eggNOG-mapper
# ============================================
echo "============================================"
echo "Running eggNOG-mapper"
echo "============================================"
echo "Taxonomic scope: Eukaryota"
echo "E-value threshold: 1e-5"
echo ""

cd ${WORK_DIR}

emapper.py \
    -i ${QUERY} \
    -o ${OUT_PREFIX} \
    --output_dir ${WORK_DIR}/results \
    --data_dir ${DB_DIR} \
    --cpu ${SLURM_CPUS_PER_TASK} \
    --tax_scope Eukaryota \
    --go_evidence non-electronic \
    --target_orthologs all \
    --seed_ortholog_evalue 1e-5 \
    --override

echo ""
echo "============================================"
echo "eggNOG-mapper complete"
echo "============================================"
echo "Output files:"
ls -lh ${WORK_DIR}/results/${OUT_PREFIX}*
echo ""
echo "Annotation summary:"
ANNOT_FILE="${WORK_DIR}/results/${OUT_PREFIX}.emapper.annotations"
if [ -f "${ANNOT_FILE}" ]; then
    echo "Total annotated: $(tail -n +6 ${ANNOT_FILE} | wc -l)"
    echo "With GO terms: $(tail -n +6 ${ANNOT_FILE} | awk -F'\t' '$10 != "-"' | wc -l)"
    echo "With KEGG KO: $(tail -n +6 ${ANNOT_FILE} | awk -F'\t' '$12 != "-"' | wc -l)"
    echo "With KEGG pathway: $(tail -n +6 ${ANNOT_FILE} | awk -F'\t' '$13 != "-"' | wc -l)"
fi
echo ""
echo "End time: $(date)"
echo "============================================"
