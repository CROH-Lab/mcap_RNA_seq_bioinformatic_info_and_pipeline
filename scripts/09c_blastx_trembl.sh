#!/bin/bash
#SBATCH --job-name=dtre_blastx_trembl
#SBATCH --output=/home/darmstrong4/mc_rework/logs/dtre_blastx_trembl_%j.out
#SBATCH --error=/home/darmstrong4/mc_rework/logs/dtre_blastx_trembl_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-user=darmstrong4@islander.tamucc.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Load conda environment
source ~/.bashrc
conda activate mcap_rnaseq

# Directories
WORK_DIR="/home/darmstrong4/mc_rework/09_symbiont_deg_annotation"
DB="${WORK_DIR}/databases/trembl_db"
QUERY="${WORK_DIR}/sequences/degs_no_sprot_hit.fa"
OUT="${WORK_DIR}/results/degs_vs_trembl.tsv"

# Check if input exists
if [ ! -f ${QUERY} ]; then
    echo "ERROR: No-hit file not found: ${QUERY}"
    echo "Run 09b_parse_sprot.sh first to generate this file"
    exit 1
fi

echo "============================================"
echo "D. trenchii BLASTx vs TrEMBL"
echo "============================================"
echo "Query: ${QUERY}"
echo "Database: ${DB}"
echo "Threads: ${SLURM_CPUS_PER_TASK}"
echo "Max target seqs: 5"
echo "E-value threshold: 1e-5"
echo "Start time: $(date)"
echo "============================================"

# Count input sequences
echo "Input sequences: $(grep -c '>' ${QUERY})"

# Run BLASTx
blastx -query ${QUERY} \
    -db ${DB} \
    -out ${OUT} \
    -evalue 1e-5 \
    -max_target_seqs 5 \
    -num_threads ${SLURM_CPUS_PER_TASK} \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs stitle"

echo ""
echo "============================================"
echo "BLASTx complete"
echo "============================================"
echo "Output: ${OUT}"
echo "Total hits: $(wc -l < ${OUT})"
echo "Unique queries with hits: $(cut -f1 ${OUT} | sort -u | wc -l)"
echo "End time: $(date)"
echo "============================================"
