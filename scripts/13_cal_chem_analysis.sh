#!/bin/bash
#SBATCH --job-name=cal_chem
#SBATCH --output=logs/cal_chem_%j.out
#SBATCH --error=logs/cal_chem_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=normal
#SBATCH --mail-user=darmstrong4@islander.tamucc.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate mcap_rnaseq

# Navigate to working directory
cd /home/darmstrong4/mc_rework/13_cal_chem

# Run the analysis
echo "Starting publication figure generation..."
echo "Date: $(date)"
echo ""

Rscript scripts/13_cal_chem_analysis.R

echo ""
echo "=== Output files ===" 
echo "Figures:"
ls -lh figures/
echo ""
echo "Output:"
ls -lh output/

echo ""
echo "Figure generation complete: $(date)"
