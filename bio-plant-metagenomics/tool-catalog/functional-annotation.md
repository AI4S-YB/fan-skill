# Functional Annotation (Prokka / DRAM / EggNOG-mapper)

**Goal:** Annotate MAGs and metagenomic contigs with functional categories, metabolic pathways
**Best for:** Understanding the functional potential of plant-associated microbiomes

## Prerequisites

- Prokka (https://github.com/tseemann/prokka)
- DRAM (https://github.com/WrightonLabCSU/DRAM)
- EggNOG-mapper (https://github.com/eggnogdb/eggnog-mapper)
- Download required databases before running

## Prokka

### Genome Annotation (for complete/near-complete MAGs)

```bash
# Annotate a single MAG
prokka --outdir prokka_mag1 \
  --prefix mag1 \
  --kingdom Bacteria \
  --genus Pseudomonas \
  --cpus 8 \
  --rfam \
  mag1.fa

# Batch annotate multiple MAGs
for mag in quality_bins/*.fa; do
  name=$(basename $mag .fa)
  prokka --outdir prokka_${name} \
    --prefix ${name} \
    --kingdom Bacteria \
    --cpus 4 \
    ${mag}
done

# Extract all protein sequences for downstream analysis
cat prokka_*/*.faa > all_proteins.faa
```

## DRAM

### Distilled and Refined Annotation of Metabolism (MAG-optimized)

```bash
# Step 1: Annotate
DRAM.py annotate \
  -i 'quality_bins/*.fa' \
  -o dram_annotation \
  --threads 32 \
  --min_contig_size 2500

# Step 2: Distill
DRAM.py distill \
  -i dram_annotation/annotations.tsv \
  -o dram_distilled \
  --genome_stats mag_quality_report.tsv

# Extract plant-relevant metabolism
grep -E 'Nitrogen|Phosphorus|Sulfur|Iron|Plant_hormone' \
  dram_distilled/metabolism_summary.tsv
```

## EggNOG-mapper

### Functional Annotation for Large Metagenomes

```bash
# Download EggNOG database (run once)
download_eggnog_data.py

# Annotate proteins
emapper.py \
  -i all_proteins.faa \
  -o eggnog_output \
  --cpu 32 \
  --itype proteins \
  --tax_scope Bacteria \
  --dbmem

# Parse results for specific pathways
python3 << 'EOF'
import pandas as pd

df = pd.read_csv("eggnog_output.emapper.annotations",
                 sep="\t", skiprows=3, comment="#")

# Filter plant-relevant functions
plant_kegg = ["nitrogen", "phosphorus", "sulfur",
              "tryptophan", "siderophore", "antibiotic"]
for pathway in plant_kegg:
    hits = df[df["KEGG_Pathway"].str.contains(pathway, na=False)]
    print(f"{pathway}: {len(hits)} genes")
EOF
```

## Key Parameters

| Tool | Parameter | Recommended | Rationale |
|------|-----------|-------------|-----------|
| Prokka | `--kingdom` | Bacteria | Set domain for correct genetic code |
| Prokka | `--rfam` | enabled | Infer non-coding RNA |
| DRAM | `--min_contig_size` | 2500 | Remove short/plant contigs |
| EggNOG | `--dbmem` | enabled | Load DB into memory (faster) |
| EggNOG | `--tax_scope` | Bacteria | Restrict to bacterial orthologs |

## Plant-Specific Notes

- **Plant hormone pathways**: Search for tryptophan-dependent IAA synthesis (ipdC, iaaM, iaaH), cytokinin (ipt), ACC deaminase (acdS)
- **Nutrient cycling**: N-fixation (nifHDK), P-solubilization (pqq, gcd), S-oxidation (sox)
- **Iron acquisition**: Siderophore biosynthesis (sid, ent, pch), iron transport (feo, fhu)
- **Plant cell wall degradation**: Cellulases, xylanases, pectinases — differentiate pathogen vs saprophyte

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| Prokka: "contig too short" | MAG too fragmented | Filter contigs > 500bp and try again |
| DRAM: memory error | Too many proteins | Split input or increase memory |
| EggNOG: "no hits" | Wrong domain or too-divergent proteins | Try broader `--tax_scope` or diamond blastx |
| Missing plant-relevant functions | Incomplete MAG or wrong compartment | Expect fewer functions in low-quality MAGs |
