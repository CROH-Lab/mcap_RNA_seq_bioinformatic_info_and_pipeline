#!/bin/bash
#SBATCH --job-name=h_GO_MWU
#SBATCH --output=/home/darmstrong4/mc_rework/logs/GO_MWU_%j.out
#SBATCH --error=/home/darmstrong4/mc_rework/logs/GO_MWU_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-user=darmstrong4@islander.tamucc.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Load environment
source ~/.bashrc
conda activate mcap_rnaseq

# Run analysis
cd /home/darmstrong4/mc_rework/10_GO_MWU
Rscript scripts/10_GO_MWU_analysis.R

echo ""
echo "=== Output files ==="
ls -lh output/
ls -lh figures/
