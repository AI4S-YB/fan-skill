# PAV Detection (PanGene / PanGenie)

**Goal:** Detect presence/absence variation (PAV) of genes across a plant pan-genome
**Best for:** Identifying genes that differ between accessions — key for understanding trait variation

## Prerequisites

- Reference genome and annotation (GFF3/GTF)
- Other genome assemblies (FASTA) OR resequencing reads (FASTQ)
- Python 3.6+, BLAST/diamond for sequence similarity

## PanGene Method

### Gene-Sequence-Based PAV

```bash
# Extract CDS sequences from reference
gffread reference.gff3 -g reference.fa -x reference.cds.fa

# Map CDS to each genome using minimap2
for genome in genome*.fa; do
  minimap2 -t 8 -x asm5 "$genome" reference.cds.fa > "${genome%.fa}.paf"
done

# Parse PAF to generate PAV matrix
python3 << 'EOF'
import pandas as pd

def parse_paf_to_pav(paf_file, min_coverage=0.5, min_identity=0.8):
    """Parse minimap2 PAF to gene presence/absence"""
    pav = {}
    with open(paf_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            gene = parts[0]
            qlen = int(parts[1])
            matches = int(parts[9])
            alen = int(parts[10])
            coverage = alen / qlen
            identity = matches / alen if alen > 0 else 0
            if coverage >= min_coverage and identity >= min_identity:
                pav[gene] = 1
    return pav

# Build PAV matrix across all genomes
all_genes = set()
genome_pavs = {}
for paf in paf_files:
    gname = paf.replace('.paf', '')
    genome_pavs[gname] = parse_paf_to_pav(paf)

# Create matrix
pav_df = pd.DataFrame(genome_pavs).fillna(0).astype(int)
pav_df.to_csv("pav_matrix.csv")
EOF
```

## PanGenie Method

### K-mer-Based PAV (no assembly needed)

```bash
# Build k-mer index from reference annotation
pangenie build -r reference.fa -a reference.gff3 \
  -o pangenie_index/

# Genotype PAV in resequencing samples
pangenie genotype -i pangenie_index/ \
  -s samples.txt \
  -t 16 \
  -o pav_genotypes.vcf

# Convert VCF to presence/absence matrix
pangenie matrix -i pav_genotypes.vcf \
  -o pav_matrix.csv
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| min_coverage (PanGene) | 0.5 | Gene must cover >50% length to count as present |
| min_identity (PanGene) | 0.8 | Alignment identity threshold |
| k-mer size (PanGenie) | 21-31 | Shorter: more sensitive; longer: more specific |
| min k-mer count (PanGenie) | 5 | Filter noise in low-coverage samples |

## Plant-Specific Notes

- Tandemly duplicated genes: PAV detection may miss one copy while detecting another
- Gene family PAV: for large gene families, detect at family level rather than individual gene level
- For wheat and other polyploids: separate subgenomes for PAV calling
- TE-related PAV: may detect TE as "gene" — filter transposon-annotated regions

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| All genes "present" in all genomes | min_coverage too low | Raise to 0.7-0.8 |
| Too many "absent" genes | Poor genome assembly contiguity | Use PanGenie k-mer approach instead |
| Inconsistent PAV between replicates | Mapping parameters inconsistent | Standardize minimap2 parameters across runs |
| PanGenie: high false negative | Low sequencing depth | Increase to >15x coverage |
