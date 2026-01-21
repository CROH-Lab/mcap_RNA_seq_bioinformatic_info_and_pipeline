#!/bin/bash
# Parse TrEMBL results identically to SwissProt

cd /home/darmstrong4/mc_rework/08_host_deg_annotation

# Check if TrEMBL results exist
if [ ! -f results/degs_vs_trembl.tsv ]; then
    echo "ERROR: TrEMBL results not found"
    exit 1
fi

echo "=== Processing TrEMBL Results ==="

# Count results
echo "Total hits: $(wc -l < results/degs_vs_trembl.tsv)"
echo "Unique genes with hits: $(cut -f1 results/degs_vs_trembl.tsv | sort -u | wc -l)"

# Create annotation table with all columns and quality tier
echo -e "gene_id\tuniprot_acc\tuniprot_entry\tpident\talign_length\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tqcovs\tquality_tier\tdescription\tsource" > results/trembl_annotations_full.tsv

awk -F'\t' '!seen[$1]++ {
    # Extract UniProt accession from sseqid (tr|A0A0G2K1R2|A0A0G2K1R2_DANRE)
    split($2, acc, "|")
    db = acc[1]
    uniprot = acc[2]
    entry = acc[3]
    
    # Assign quality tier
    if ($11 < 1e-20 && $3 >= 40 && $15 >= 60) tier="HIGH"
    else if ($11 < 1e-10 && $3 >= 30 && $15 >= 50) tier="MEDIUM"
    else tier="LOW"
    
    # Extract short description (remove OS=... onwards)
    desc = $16
    sub(/ OS=.*/, "", desc)
    
    # Print all columns plus source
    print $1"\t"uniprot"\t"entry"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"tier"\t"desc"\tTrEMBL"
}' results/degs_vs_trembl.tsv >> results/trembl_annotations_full.tsv

echo ""
echo "=== TrEMBL Annotation Summary ==="
echo "Total annotated genes: $(tail -n +2 results/trembl_annotations_full.tsv | wc -l)"
echo ""
echo "Quality tier breakdown:"
tail -n +2 results/trembl_annotations_full.tsv | cut -f17 | sort | uniq -c
echo ""
echo "Preview:"
head -6 results/trembl_annotations_full.tsv | column -t -s$'\t'
