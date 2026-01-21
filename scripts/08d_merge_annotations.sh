#!/bin/bash
# Merge SwissProt and TrEMBL annotations

cd /home/darmstrong4/mc_rework/09_blast_annotation

echo "=== Merging SwissProt and TrEMBL Annotations ==="

# Add source column to SwissProt annotations
echo -e "gene_id\tuniprot_acc\tuniprot_entry\tpident\talign_length\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tqcovs\tquality_tier\tdescription\tsource" > results/all_annotations_full.tsv

# Add SwissProt results (with source column)
tail -n +2 results/sprot_annotations_full.tsv | awk -F'\t' '{print $0"\tSwissProt"}' >> results/all_annotations_full.tsv

# Add TrEMBL results (already has source column)
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
echo ""
echo "By source and quality:"
tail -n +2 results/all_annotations_full.tsv | awk -F'\t' '{print $19"\t"$17}' | sort | uniq -c

# Count genes still without annotation
echo ""
echo "=== Annotation Coverage ==="
TOTAL_DEGS=$(wc -l < sequences/all_degs.txt)
ANNOTATED=$(tail -n +2 results/all_annotations_full.tsv | wc -l)
UNANNOTATED=$((TOTAL_DEGS - ANNOTATED))
echo "Total DEGs: $TOTAL_DEGS"
echo "Annotated: $ANNOTATED ($(awk "BEGIN {printf \"%.1f\", $ANNOTATED/$TOTAL_DEGS*100}")%)"
echo "Unannotated: $UNANNOTATED ($(awk "BEGIN {printf \"%.1f\", $UNANNOTATED/$TOTAL_DEGS*100}")%)"

# Create list of unannotated genes
comm -23 <(sort sequences/all_degs.txt) \
         <(tail -n +2 results/all_annotations_full.tsv | cut -f1 | sort) > results/unannotated_degs.txt
echo ""
echo "Unannotated gene list saved to: results/unannotated_degs.txt"
