# SV Genotyping (vg giraffe)

**Goal:** Genotype structural variants in new samples by mapping reads to a pan-genome graph
**Best for:** Large-scale SV genotyping across populations when a pan-genome graph exists

## Prerequisites

- vg toolkit (https://github.com/vgteam/vg)
- Pan-genome graph (GFA or vg format)
- Sequencing reads (FASTQ) for samples to genotype
- vg autoindex for graph indexing

## Basic Usage

```bash
# Step 1: Build graph index
vg autoindex --workflow giraffe \
  -g pan_graph.gfa \
  -p pan_index \
  -t 16

# Step 2: Map reads to graph
vg giraffe -Z pan_index.giraffe.gbz \
  -G pan_index.giraffe.gbwt \
  -d pan_index.dist \
  -f sample_R1.fq.gz -f sample_R2.fq.gz \
  -t 16 \
  > mapped.gam

# Step 3: Call variants from mapped reads
vg pack -x pan_index.xg \
  -g mapped.gam \
  -o sample.pack \
  -t 16

vg call pan_index.xg \
  -k sample.pack \
  -a \
  > sv_calls.vcf

# Step 4: Genotype SVs across multiple samples
vg autoindex --workflow giraffe \
  -g pan_graph.gfa -p graph_index

for sample in sample1 sample2 sample3; do
  vg giraffe -Z graph_index.giraffe.gbz \
    -f ${sample}_R1.fq.gz -f ${sample}_R2.fq.gz \
    -t 16 > ${sample}.gam
done

# Joint genotyping
vg pack -x graph_index.xg \
  -g sample1.gam -g sample2.gam -g sample3.gam \
  -o joint.pack

vg call graph_index.xg \
  -k joint.pack \
  -a > joint_svs.vcf
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| `-t` | 16-32 | Thread count |
| `-a` (vg call) | enabled | Output all alleles, not just variants |
| `--min-total-depth` | 5 | Minimum read depth for calling |
| `--min-alt-fraction` | 0.15 | Minimum alternate allele fraction |

## Plant-Specific Notes

- For polyploids: expect higher depth per locus due to homeolog mapping
- Large plant genomes: increase graph index building memory to 128-256 GB
- Repetitive regions: vg may produce ambiguous mappings — filter calls with low QUAL
- For population-scale genotyping: batch 50-100 samples per joint calling run

## Alternative: GraphAligner (Long Reads)

```bash
# Better for long reads and large SVs
GraphAligner -g pan_graph.gfa \
  -f longreads.fastq \
  -a aligned.gaf \
  -t 16 \
  -x vg
```

## Alternative: Paragraph (Small SVs)

```bash
# Specialized for small SV genotyping
multigrmpy.py -i sample.bam \
  -g graph_manifest.json \
  -o sv_genotypes.vcf
```

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "gbwt index too large" | Graph too big | Split by chromosome or use pruned graph |
| Low mapping rate | Graph doesn't represent sample diversity | Include more genomes in the graph |
| Excessive memory in vg pack | Too many samples in one pack run | Batch samples, 20-50 at a time |
| "fragment too long" warning | Reads exceed graph segment size | Increase segment length in graph construction |
