#!/bin/bash
# Parse TrEMBL results and merge all annotations

cd /home/darmstrong4/mc_rework/09_symbiont_deg_annotation

echo "=== Processing TrEMBL Results ==="

# Check if TrEMBL results exist
if [ -f results/degs_vs_trembl.tsv ] && [ -s results/degs_vs_trembl.tsv ]; then
    echo "Total hits: $(wc -l < results/degs_vs_trembl.tsv)"
    echo "Unique genes with hits: $(cut -f1 results/degs_vs_trembl.tsv | sort -u | wc -l)"
    
    # Create TrEMBL annotation table
    echo -e "gene_id\tuniprot_acc\tuniprot_entry\tpident\talign_length\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tqcovs\tquality_tier\tdescription\tsource" > results/trembl_annotations_full.tsv
    
    awk -F'\t' '!seen[$1]++ {
        split($2, acc, "|")
        uniprot = acc[2]
        entry = acc[3]
        
        if ($11 < 1e-20 && $3 >= 40 && $15 >= 60) tier="HIGH"
        else if ($11 < 1e-10 && $3 >= 30 && $15 >= 50) tier="MEDIUM"
        else tier="LOW"
        
        desc = $16
        sub(/ OS=.*/, "", desc)
        
        print $1"\t"uniprot"\t"entry"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"tier"\t"desc"\tTrEMBL"
    }' results/degs_vs_trembl.tsv >> results/trembl_annotations_full.tsv
    
    echo ""
    echo "TrEMBL annotations: $(tail -n +2 results/trembl_annotations_full.tsv | wc -l)"
else
    echo "No TrEMBL results found (or empty file)"
fi

echo ""
echo "=== Merging All Annotations ==="

# Create combined file
echo -e "gene_id\tuniprot_acc\tuniprot_entry\tpident\talign_length\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tqcovs\tquality_tier\tdescription\tsource" > results/all_annotations_full.tsv

# Add SwissProt results
tail -n +2 results/sprot_annotations_full.tsv | awk -F'\t' '{print $0"\tSwissProt"}' >> results/all_annotations_full.tsv

# Add TrEMBL results if they exist
if [ -f results/trembl_annotations_full.tsv ]; then
    tail -n +2 results/trembl_annotations_full.tsv >> results/all_annotations_full.tsv
fi

# Summary
echo ""
echo "=== Combined Annotation Summary ==="
echo "Total annotated DEGs: $(tail -n +2 results/all_annotations_full.tsv | wc -l)"
echo ""
echo "By source:"
tail -n +2 results/all_annotations_full.tsv | cut -f19 | sort | uniq -c
echo ""
echo "By quality tier:"
tail -n +2 results/all_annotations_full.tsv | cut -f17 | sort | uniq -c

# Count unannotated
echo ""
echo "=== Annotation Coverage ==="
TOTAL_DEGS=$(wc -l < sequences/all_degs.txt)
ANNOTATED=$(tail -n +2 results/all_annotations_full.tsv | wc -l)
UNANNOTATED=$((TOTAL_DEGS - ANNOTATED))
echo "Total DEGs: $TOTAL_DEGS"
echo "Annotated: $ANNOTATED ($(awk "BEGIN {printf \"%.1f\", $ANNOTATED/$TOTAL_DEGS*100}")%)"
echo "Unannotated: $UNANNOTATED ($(awk "BEGIN {printf \"%.1f\", $UNANNOTATED/$TOTAL_DEGS*100}")%)"

# Create unannotated list
comm -23 <(sort sequences/all_degs.txt) \
         <(tail -n +2 results/all_annotations_full.tsv | cut -f1 | sort) > results/unannotated_degs.txt
echo "Unannotated gene list saved to: results/unannotated_degs.txt"
