# Functional Annotation

**Goal:** Assign functional terms (GO, KEGG, Pfam, InterPro) to predicted proteins
**Best for:** Post gene-prediction annotation of plant proteomes

## Prerequisites
- Predicted protein sequences (FASTA)
- EggNOG-mapper 2.1+
- InterProScan 5.0+
- Sufficient disk space for databases (50GB+)

## EggNOG-mapper

### Download Plant Databases

```bash
# Download EggNOG database (Viridiplantae level recommended for plants)
download_eggnog_data.py --taxa Viridiplantae
```

### Run Annotation

```bash
emapper.py \
  -i predicted_proteins.fasta \
  --output eggnog_results \
  --tax_scope Viridiplantae \
  --cpu 16 \
  --sensmode diamond \
  --go_evidence non-electronic \
  --pfam_realign realign
```

### Output Files

| File | Content |
|------|---------|
| eggnog_results.emapper.annotations | Full annotation table |
| eggnog_results.emapper.seed_orthologs | Best seed orthologs |
| eggnog_results.emapper.hits | All database hits |

## InterProScan

### Run Annotation

```bash
interproscan.sh \
  -i predicted_proteins.fasta \
  -f TSV,GFF3,XML,JSON \
  -goterms \
  -pa \
  -cpu 16 \
  -o interpro_results \
  -dp
```

### Key Analyses Included

| Analysis | What It Detects |
|----------|-----------------|
| Pfam | Protein domain families |
| Gene3D | Structural domain assignments |
| SUPERFAMILY | SCOP superfamily assignments |
| CDD | Conserved domains from NCBI |
| TIGRFAMs | Prokaryotic protein families (limited plant utility) |
| PANTHER | Protein family/subfamily classifications |
| SignalP | Signal peptide prediction |
| TMHMM | Transmembrane helix prediction |
| Phobius | Combined signal peptide + transmembrane |
| Coils | Coiled-coil region prediction |

## Merging Results

```python
# Merge EggNOG and InterProScan annotations
import pandas as pd

eggnog = pd.read_csv('eggnog_results.emapper.annotations', sep='\t',
                      comment='#', skiprows=3)
interpro = pd.read_csv('interpro_results.tsv', sep='\t',
                        names=['protein', 'md5', 'length', 'analysis',
                               'db_id', 'description', 'start', 'end',
                               'score', 'status', 'date', 'ipr_id', 'ipr_desc',
                               'go_terms', 'pathways'])

# Generate summarized annotation per protein
functional_summary = eggnog[['#query', 'GOs', 'KEGG_ko', 'Preferred_name', 'Description']]
# Combine with InterPro GO terms
```

## Key Parameters

| Parameter | Purpose |
|-----------|---------|
| --tax_scope | Limit search to specific taxonomic group |
| --sensmode | Sensitivity mode (faster: diamond / more sensitive: diamond-sensitive) |
| --go_evidence | Filter GO terms by evidence code |
| -goterms | Enable GO term mapping in InterProScan |
| -f | Output formats (TSV recommended for parsing) |
| -dp | Disable pre-calculated match lookup (use for first run) |

## Plant-Specific Considerations

- Use `--tax_scope Viridiplantae` (green plants) or narrower for EggNOG-mapper
- PANTHER 18.0+ includes better plant pathway coverage
- For crop species, supplement with species-specific Pfam families if available
- Plant-specific enzyme families (e.g., cytochrome P450, glycosyltransferases) need proper taxonomic context

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Database not found" | EggNOG DB not downloaded | Run `download_eggnog_data.py` first |
| "Out of memory" | Large proteome | Split input into chunks (5000 proteins each) |
| InterProScan slow | Full suite on large proteome | Use `-appl Pfam,Gene3D,SUPERFAMILY` to limit analyses |
| "No hits found" | Wrong taxonomic scope | Use broader `--tax_scope` or check protein quality |
