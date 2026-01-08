#!/bin/bash
#SBATCH -J star_div                # Job name
#SBATCH -o logs/star_diversity.out # Standard output (array-aware)
#SBATCH -e logs/star_diveristy.err # Standard error (array-aware)
#SBATCH --partition=normal           # Normal partition (350GB available)
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=12:00:00
#SBATCH --mail-user=darmstrong4@islander.tamucc.edu
#SBATCH --mail-type=BEGIN,END,FAIL

echo "=============================================="
echo "SYMBIONT COMPOSITION ACROSS ALL SAMPLES"
echo "=============================================="
echo ""
echo "This will take a few minutes..."
echo ""

# Create output file
echo -e "Sample\tMcap_reads\tCgor_reads\tDtre_reads\tUnmapped\tMcap%\tCgor%\tDtre%\tSymbiont_Ratio" > 05_star_align/organism_distribution.tsv

for SAMPLE in 1AS 1AW 1BS 1BW 1CS 1CW 1DS 1DW 2AS 2AW 2BS 2BW 2CS 2CW 2DS 2DW 3AS 3AW 3BS 3BW 3CS 3CW 3DS 3DW; do
    echo "Processing ${SAMPLE}..."
    BAM="05_star_align/${SAMPLE}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    
    # Count reads per organism
    COUNTS=$(samtools view ${BAM} | cut -f3 | cut -d'_' -f1 | sort | uniq -c | awk '{print $2":"$1}' | tr '\n' '\t')
    
    # Extract individual counts
    MCAP=$(samtools view ${BAM} | cut -f3 | grep "^Mcap_" | wc -l)
    CGOR=$(samtools view ${BAM} | cut -f3 | grep "^Cgor_" | wc -l)
    DTRE=$(samtools view ${BAM} | cut -f3 | grep "^Dtre_" | wc -l)
    UNMAPPED=$(samtools view ${BAM} | cut -f3 | grep "^\*" | wc -l)
    
    # Calculate percentages
    TOTAL=$((MCAP + CGOR + DTRE + UNMAPPED))
    MCAP_PCT=$(echo "scale=1; $MCAP * 100 / $TOTAL" | bc)
    CGOR_PCT=$(echo "scale=1; $CGOR * 100 / $TOTAL" | bc)
    DTRE_PCT=$(echo "scale=1; $DTRE * 100 / $TOTAL" | bc)
    
    # Symbiont ratio (Cgor:Dtre)
    if [[ $DTRE -gt 0 ]]; then
        SYM_RATIO=$(echo "scale=2; $CGOR / $DTRE" | bc)
    else
        SYM_RATIO="NA"
    fi
    
    echo -e "${SAMPLE}\t${MCAP}\t${CGOR}\t${DTRE}\t${UNMAPPED}\t${MCAP_PCT}\t${CGOR_PCT}\t${DTRE_PCT}\t${SYM_RATIO}" >> 05_star_align/organism_distribution.tsv
done

echo ""
echo "=== ORGANISM DISTRIBUTION SUMMARY ==="
column -t -s $'\t' 05_star_align/organism_distribution.tsv

echo ""
echo "=== INTERPRETATION GUIDE ==="
echo "Cgor:Dtre Ratio > 1 = More Clade C (typical healthy)"
echo "Cgor:Dtre Ratio < 1 = More Clade D (stress-associated)"
echo ""
echo "Results saved to: 05_star_align/organism_distribution.tsv"
