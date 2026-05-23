# Assembly Preprocessing for Annotation

**Goal:** Prepare genome assembly for annotation — quality check, contig filtering, and format conversion
**Best for:** All genome annotation projects, especially draft assemblies

## Prerequisites
- seqkit 2.0+
- samtools 1.15+
- Python 3 with biopython

## Assembly Quality Assessment

### Contig Statistics

```bash
seqkit stats -a genome.fasta
```

Check these metrics:
- **Total size**: Expected genome size for the species?
- **N50**: Scaffold N50 > 100Kb for reasonable gene prediction
- **Number of contigs**: Fewer is better; >10,000 is fragmented
- **GC content**: Plant genomes typically 35-45% (rice 43%, maize 47%, wheat 45%)

### Identify Foreign Contamination

```bash
# BLAST a random sample of contigs against nt
seqkit sample -n 1000 genome.fasta > sample_contigs.fasta
blastn -query sample_contigs.fasta -db nt -outfmt '6 qseqid stitle' -num_threads 8 | head -50
```

### Mitochondrial and Chloroplast Contigs

```bash
# Run BLAST against organelle genomes
makeblastdb -in plant_organelle_db.fasta -dbtype nucl
blastn -query genome.fasta -db plant_organelle_db.fasta \
  -outfmt '6 qseqid sseqid pident length' -out organelle_hits.txt

# Extract contigs with strong organelle hits (>90% identity, >500bp)
awk '$3 > 90 && $4 > 500 {print $1}' organelle_hits.txt | sort -u > organelle_contigs.txt
```

## Contig Filtering

### Remove Short Contigs

```bash
# Remove contigs shorter than 500bp (too short for gene prediction)
seqkit seq -m 500 genome.fasta > genome_filtered.fasta
```

### Remove Organelle Contigs (Optional)

```bash
seqkit grep -v -f organelle_contigs.txt genome_filtered.fasta > genome_nuclear.fasta
```

## Genome Sorting and Renaming

### Sort by Size (Helps BUSCO and RepeatMasker)

```bash
seqkit sort -l -r genome.fasta > genome_sorted.fasta
```

### Rename Contigs with Simple Names

```bash
# Long contig names cause issues with some tools
seqkit replace -p ".*" -r "contig_{nr}" genome.fasta > genome_renamed.fasta
```

## Key Parameters

| Parameter | Purpose |
|-----------|---------|
| -m 500 | Minimum contig length for filtering (default for annotation) |
| -M | Maximum contig length |
| -l | Sort by sequence length |
| -r | Reverse sort order (largest first) |
| -p | Pattern for sequence ID replacement |

## Plant-Specific Considerations

- Organelle genomes are often assembled alongside nuclear DNA; filter them out to avoid spurious gene predictions
- For polyploid genomes, check for haplotype-resolved assemblies; individual haplotypes should be annotated separately
- Chloroplast genomes may be present at 100-1000x coverage in WGS data; verify with coverage analysis
- Plant nuclear genomes often contain large NUPTs (nuclear plastid DNA insertions); these are legitimate targets for annotation

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Fasta index exceeds memory" | Assembly file too large | Use `samtools faidx` for indexed access |
| Contaminant sequences | Bacterial/vector contamination | Screen with VecScreen or BLAST against UniVec |
| Ambiguous bases (N) | Assembly gaps | Count with `seqtk comp genome.fasta`; if >5%, re-assemble |
| Duplicate contig names | Merging multiple assemblies | Rename with unique prefixes before merging |
