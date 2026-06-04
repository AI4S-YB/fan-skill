# Bacterial and Fungal Genome Consultation Guide

## Quick Decision Framework

### Step 1: Identify Organism Type

| Characteristic | Bacteria | Fungi |
|---------------|----------|-------|
| Expected size | 2-12 Mb | 10-100 Mb |
| Introns | Rare | Common |
| Gene count | 2,000-8,000 | 5,000-15,000 |
| BUSCO lineage | bacteria_odb10 | fungi_odb10 |

### Step 2: Choose Assembly Strategy

| Data Available | Recommended Tool |
|---------------|------------------|
| Illumina + Long-read | Unicycler (hybrid) |
| Long-read only | Flye + Polishing |
| Illumina only | SPAdes |

### Step 3: Select Annotation Pipeline

| Organism | Recommended Tool |
|----------|-----------------|
| Bacteria | Prokka |
| Fungi + RNA-seq | BRAKER3 |
| Fungi, no RNA-seq | AUGUSTUS (fungal training) |

## Common Consultation Scenarios

### Q1: "I have bacterial Illumina and ONT data"

**Recommended workflow:**
1. Quality control with FastQC
2. Hybrid assembly with Unicycler
3. Assess with QUAST and BUSCO
4. Annotate with Prokka
5. Run antiSMASH for BGCs
6. Run CARD RGI for resistance genes

### Q2: "How do I know if assembly is complete?"

**Check for:**
- Single circular chromosome (most bacteria)
- No gaps (N characters)
- BUSCO complete > 95%
- Circularization confirmed by overlap

### Q3: "What secondary metabolites does my strain produce?"

**Use antiSMASH:**
- Predicts BGC types (PKS, NRPS, terpenes, etc.)
- Compares to known clusters
- For fungi, use antiSMASH fungal mode

### Q4: "Is my strain antibiotic resistant?"

**Run CARD RGI:**
- Identifies resistance genes
- Predicts resistance phenotype
- Reports confidence levels

### Q5: "How do I compare multiple strains?"

**Comparative workflow:**
1. Annotate all genomes uniformly
2. Build pan-genome with Roary
3. Calculate ANI with FastANI
4. Build core-genome phylogeny
5. Visualize synteny with Mauve

## Parameter Selection Guide

### SPAdes Parameters

| Parameter | Typical Value | Notes |
|-----------|---------------|-------|
| --careful | Yes | Reduces mismatches |
| --cov-cutoff | auto | Automatic coverage cutoff |
| -k | 21,33,55,77 | Multiple k-mer sizes |

### Prokka Parameters

| Parameter | Typical Value | Notes |
|-----------|---------------|-------|
| --kingdom | Bacteria/Archaea | Required |
| --genus | Species genus | Improves annotation |
| --evalue | 1e-06 | BLAST e-value cutoff |

### antiSMASH Parameters

| Parameter | Bacteria | Fungi |
|-----------|----------|-------|
| --taxon | bacteria | fungi |
| --cb-knownclusters | Yes | Yes |
| --smcog-trees | Yes | No |

## Risk Assessment

### High Risk
- Multiple strains/species mixed
- High contamination (>10%)
- Very low coverage (<20x)
- No reference for quality check

### Medium Risk
- Moderate contamination (5-10%)
- Single platform only
- Distant reference genome

### Low Risk
- Pure culture
- Hybrid data available
- Close reference genome

## Result Interpretation

### Assembly Quality

| N50 | Interpretation |
|-----|----------------|
| > chromosome size | Excellent (circular) |
| > 100 kb | Good |
| < 50 kb | Poor (may need reassembly) |

### BUSCO Results

| Complete | Fragmented | Missing | Interpretation |
|----------|------------|---------|----------------|
| >95% | <5% | <5% | Excellent |
| 90-95% | 5-10% | <5% | Good |
| <90% | >10% | >5% | Incomplete |

### antiSMASH BGCs

| Match Type | Interpretation |
|------------|----------------|
| 100% to known | Known compound |
| 70-99% | Variant of known |
| <70% | Novel cluster |

## When to Escalate

- Metagenome-assembled genome (MAG)
- Polyploid fungal genome
- Plasmid-rich bacteria
- Extreme GC content (>70% or <30%)
- Symbiont with reduced genome

## 数据状态标签规范 (C4强制)

在生成方案的数据画像部分，必须显式标注数据状态：
- **数据状态**: FULL — [说明哪些数据完整可用，如"细菌/真菌基因组测序数据(Illumina+长读长混合)、纯培养物"]
- **数据状态**: PARTIAL — [说明哪些数据缺失，如"仅有Illumina短读长数据；或存在污染(>10%)"]
- **数据状态**: EMPTY — [说明数据不可用原因，如"菌株尚未分离培养"]

违例判定: 仅列出文件名/大小但无显式FULL/PARTIAL/EMPTY状态标签 → C4=0分
