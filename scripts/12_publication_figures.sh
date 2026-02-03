#!/bin/bash
#SBATCH --job-name=pub_figures
#SBATCH --output=logs/pub_figures_%j.out
#SBATCH --error=logs/pub_figures_%j.err
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
cd /home/darmstrong4/mc_rework/12_publication_figures

# Run the analysis
echo "Starting publication figure generation..."
echo "Date: $(date)"
echo ""

Rscript scripts/12_publication_figures.R

echo ""
echo "=== Output files ===" 
echo "Figures:"
ls -lh figures/
echo ""
echo "Data:"
ls -lh data/

echo ""
echo "Figure generation complete: $(date)"
