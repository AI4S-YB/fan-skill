# Repeat Masking

**Goal:** Mask repetitive elements (transposons, LTRs, simple repeats) in plant genomes before gene prediction
**Best for:** All plant genome annotation projects

## Prerequisites
- RepeatMasker 4.1+
- RepeatModeler 2.0+ (for de novo repeat library)
- RMBlast or other search engine
- Plant repeat database (e.g., MIPS REdat, Repbase plants)

## Standard Workflow

### 1. Build Species-Specific Repeat Library

```bash
# De novo repeat identification
BuildDatabase -name species_db -engine ncbi genome.fasta
RepeatModeler -database species_db -engine ncbi -pa 16 -LTRStruct

# Combine with known plant repeats
cat species_db-families.fa plant_repeats.lib > combined_repeats.lib
```

### 2. Run RepeatMasker

```bash
# Mask using combined library
RepeatMasker -pa 16 \
  -lib combined_repeats.lib \
  -xsmall \
  -gff \
  -dir repeatmasker_output \
  genome.fasta
```

### 3. Hard vs Soft Masking

```bash
# Hard masking (replace repeats with N) — use for gene prediction
RepeatMasker -pa 16 -lib combined_repeats.lib genome.fasta

# Soft masking (lowercase repeats) — use for BLAST/repeat analysis
RepeatMasker -pa 16 -lib combined_repeats.lib -xsmall genome.fasta
```

## Key Parameters

| Parameter | Purpose |
|-----------|---------|
| -pa | Number of parallel threads |
| -lib | Custom repeat library (preferred for plants) |
| -species | Use built-in RepeatMasker species library |
| -xsmall | Soft mask (lowercase) instead of N |
| -gff | Output GFF annotation of repeats |
| -s | Slow/sensitive mode (use for final annotation) |
| -nolow | Skip low-complexity/simple repeats |

## Plant-Specific Considerations

- **LTR retrotransposons** are the dominant repeat class in most plant genomes
- Use `-LTRStruct` in RepeatModeler to leverage LTR structural features
- For large genomes (>1Gb), consider EDTA or LTR_FINDER_parallel before RepeatMasker
- Monocot vs dicot: TE families differ substantially; use taxonomy-appropriate databases
- Polyploid genomes: expect duplicated TE families across subgenomes

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "No sequences in database" | Empty or misformatted repeat library | Validate with `grep -c ">" combined_repeats.lib` |
| "Out of memory" | Large genome with -s flag | Reduce threads or mask by chromosome |
| "RMBlast not found" | Search engine not in PATH | Set `RMBLAST_DIR` or use `-engine ncbi` |
| "RepeatModeler timeout" | Large genome assembly | Run on subset (longest 100 scaffolds) first |
| Low masking percentage | Repeat library missing plant-specific TEs | Add MIPS/Repbase plant sequences to library |
