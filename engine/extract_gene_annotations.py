#!/usr/bin/env python3
"""
Extract gene annotations from GFF3 file (fast, no R dependency).
Usage: python3 extract_gene_annotations.py <gff3_file> [output_csv]
"""
import gzip, re, csv, sys

gff_file = sys.argv[1]
output = sys.argv[2] if len(sys.argv) > 2 else "gene_annotations.csv"

# Determine if gzipped
open_func = gzip.open if gff_file.endswith('.gz') else open

gene_count = 0
with open_func(gff_file, 'rt') as f, open(output, 'w', newline='') as out:
    writer = csv.writer(out)
    writer.writerow(['gene_id', 'name', 'product', 'go_terms', 'note'])

    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9 or parts[2] != 'gene':
            continue

        attrs = parts[8]
        # Extract key fields
        id_match = re.search(r'ID=([^;]+)', attrs)
        name_match = re.search(r'Name=([^;]+)', attrs)
        prod_match = re.search(r'product=([^;]+)', attrs)
        go_match = re.search(r'Ontology_term=([^;]+)', attrs)
        note_match = re.search(r'Note=([^;]+)', attrs)

        gene_id = id_match.group(1) if id_match else ''
        name = name_match.group(1) if name_match else ''
        product = prod_match.group(1) if prod_match else ''
        go_terms = go_match.group(1) if go_match else ''
        note = note_match.group(1) if note_match else ''

        # Use Note as product if product field is empty
        if not product and note:
            product = note

        writer.writerow([gene_id, name, product, go_terms, note])
        gene_count += 1

print(f"Extracted {gene_count} genes -> {output}")
