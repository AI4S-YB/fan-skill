---
name: bio-plant-variant-calling
description: >
  Plant variant discovery from resequencing data — FASTQ to VCF.
  Covers BWA-MEM/minimap2 alignment, GATK/bcftools/DeepVariant
  calling, hard filtering/VQSR, and VCF QC. Dual-mode B+C.
tool_type: mixed
primary_tool: GATK4
workflow: true
---

# Plant Variant Calling

SNP/InDel discovery from resequencing data for plant genomes.

## Decision Modes

Set in `params.yaml`:

```yaml
decision_mode: hybrid   # rule | expert | hybrid (recommended)
```

| Mode | Behavior |
|------|----------|
| **rule** | Strictly follow `decision-matrix.yaml` conditions -> methods |
| **expert** | Read `analyst-notebook.md` and reason like a plant genomicist |
| **hybrid** | Matrix first -> fall back to expert when no rule matches |

## Before Analysis

1. Profile the data: `bash bio-plant-infra/scripts/inspect_data.sh <input_path>`
2. Check the environment: `bash bio-plant-infra/scripts/check_env.sh`
3. If species detected, read `bio-plant-infra/references/species-cheatsheet.md`
4. Verify reference genome: index must include `.fa`, `.fai`, and `.dict`

## Analysis Flow

### Step 1: Alignment
Choose aligner based on sequencing platform:
- Illumina short reads -> BWA-MEM
- PacBio/ONT long reads -> minimap2
See `decision-matrix.yaml` > `alignment`.

### Step 2: BAM Processing
Sort, mark duplicates, index BAM files.
Optionally run BQSR if known sites are available.
See `tool-catalog/bwa-mem.md` for commands.

### Step 3: Variant Calling
Choose caller based on sample count, ploidy, and platform.
See `decision-matrix.yaml` > `variant_caller`.
- GATK HaplotypeCaller (diploid, >= 10 samples) -> `tool-catalog/gatk-haplotype.md`
- bcftools mpileup (small cohorts, any ploidy) -> `tool-catalog/bcftools.md`
- DeepVariant (long reads, maximum accuracy) -> `tool-catalog/deepvariant.md`

### Step 4: Variant Filtering
Hard filtering or VQSR based on sample count and known sites.
See `decision-matrix.yaml` > `variant_filtering`.
All parameters in `tool-catalog/vcf-qc.md`.

### Step 5: VCF QC and Visualization
Validate Ti/Tv ratio, missing rate, heterozygosity.
See `tool-catalog/visualization.md`.

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| After alignment | Mapping rate | > 90% (diploid), > 70% (polyploid) |
| After calling | Ti/Tv ratio | 2.0-2.5 |
| After filtering | SNP retention | > 80% |
| After QC | Heterozygosity | < 5% (inbred), 10-30% (outcross) |
| After QC | Missing rate per sample | < 20% |

## References

- `bio-plant-infra/references/species-cheatsheet.md` — species-specific guidance
- `bio-plant-infra/references/plant-databases.md` — genome & annotation sources
- `bio-plant-infra/references/qc-thresholds.yaml` — QC thresholds
- `references/variant-plant-special.md` — plant variant calling specifics
