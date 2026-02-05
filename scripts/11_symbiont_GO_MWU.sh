#!/bin/bash
#SBATCH --job-name=s_GO_MWU
#SBATCH --output=logs/symbiont_GO_MWU_%j.out
#SBATCH --error=logs/symbiont_GO_MWU_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=normal
#SBATCH --mail-user=darmstrong4@islander.tamucc.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate mcap_rnaseq

# Navigate to working directory
cd /home/darmstrong4/mc_rework/11_symbiont_GO_MWU

# Run the analysis
echo "Starting SYMBIONT GO_MWU analysis..."
echo "Date: $(date)"
echo ""

Rscript scripts/11_symbiont_GO_MWU_analysis.R

echo ""
echo "=== Output files ===" 
ls -lh output/
ls -lh figures/

echo ""
echo "Analysis complete: $(date)"
