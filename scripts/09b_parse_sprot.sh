#!/bin/bash
# Parse SwissProt results and identify no-hit genes for TrEMBL

cd /home/darmstrong4/mc_rework/09_symbiont_deg_annotation

echo "=== Processing SwissProt Results ==="

# Check if results exist
if [ ! -f results/degs_vs_sprot.tsv ]; then
    echo "ERROR: SwissProt results not found"
    exit 1
fi

# Count results
echo "Total hits: $(wc -l < results/degs_vs_sprot.tsv)"
echo "Unique genes with hits: $(cut -f1 results/degs_vs_sprot.tsv | sort -u | wc -l)"

# Create annotation table - keep best hit per gene
echo -e "gene_id\tuniprot_acc\tuniprot_entry\tpident\talign_length\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tqcovs\tquality_tier\tdescription" > results/sprot_annotations_full.tsv

awk -F'\t' '!seen[$1]++ {
    # Extract UniProt accession from sseqid (sp|P12345|NAME_SPECIES)
    split($2, acc, "|")
    uniprot = acc[2]
    entry = acc[3]
    
    # Assign quality tier
    if ($11 < 1e-20 && $3 >= 40 && $15 >= 60) tier="HIGH"
    else if ($11 < 1e-10 && $3 >= 30 && $15 >= 50) tier="MEDIUM"
    else tier="LOW"
    
    # Extract description (remove OS=... onwards)
    desc = $16
    sub(/ OS=.*/, "", desc)
    
    print $1"\t"uniprot"\t"entry"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"tier"\t"desc
}' results/degs_vs_sprot.tsv >> results/sprot_annotations_full.tsv

echo ""
echo "=== SwissProt Annotation Summary ==="
echo "Total annotated genes: $(tail -n +2 results/sprot_annotations_full.tsv | wc -l)"
echo ""
echo "Quality tier breakdown:"
tail -n +2 results/sprot_annotations_full.tsv | cut -f17 | sort | uniq -c

# Find genes without SwissProt hits
echo ""
echo "=== Identifying genes without SwissProt hits ==="
comm -23 <(sort sequences/all_degs.txt) \
         <(cut -f1 results/degs_vs_sprot.tsv | sort -u) > sequences/no_sprot_hit_ids.txt

NO_HIT=$(wc -l < sequences/no_sprot_hit_ids.txt)
echo "Genes without SwissProt hits: $NO_HIT"

# Extract sequences for TrEMBL search
if [ $NO_HIT -gt 0 ]; then
    awk 'NR==FNR {ids[$1]; next} 
         /^>/ {id = $1; gsub(/^>/, "", id); keep = (id in ids)} 
         keep' sequences/no_sprot_hit_ids.txt sequences/degs_cds.fa > sequences/degs_no_sprot_hit.fa
    echo "Extracted $(grep -c '>' sequences/degs_no_sprot_hit.fa) sequences for TrEMBL search"
    echo "Ready to run: sbatch scripts/09c_blastx_trembl.sh"
else
    echo "All genes have SwissProt hits - no TrEMBL search needed!"
fi
