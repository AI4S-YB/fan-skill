# BSA Analysis Notebook

## Overview

This notebook provides step-by-step guidance for Bulk Segregant Analysis (BSA), including MutMap, QTL-seq, and BSR-seq methods.

---

## Phase 0: Data Preparation

### 0.1 Verify Sample Information

```bash
# Check FASTQ files
ls -la *.fastq.gz

# Expected: 2+ pools (e.g., mutant_pool_R1.fastq.gz, wild_pool_R1.fastq.gz)
# Minimum: 20x coverage per pool for initial analysis
# Recommended: 30-50x coverage per pool
```

### 0.2 Reference Genome Preparation

```bash
# Index reference genome
bwa index reference.fasta
samtools faidx reference.fasta

# Create sequence dictionary (for GATK)
gatk CreateSequenceDictionary -R reference.fasta
```

---

## Phase 1: Read Alignment

### 1.1 BWA-MEM Alignment (DNA-seq)

```bash
# Align each pool
# Pool A (mutant/extreme phenotype 1)
bwa mem -t 16 -M -R '@RG\tID:poolA\tSM:poolA\tPL:ILLUMINA' \
    reference.fasta poolA_R1.fastq.gz poolA_R2.fastq.gz | \
    samtools sort -@ 8 -o poolA.sorted.bam

# Pool B (wild-type/extreme phenotype 2)
bwa mem -t 16 -M -R '@RG\tID:poolB\tSM:poolB\tPL:ILLUMINA' \
    reference.fasta poolB_R1.fastq.gz poolB_R2.fastq.gz | \
    samtools sort -@ 8 -o poolB.sorted.bam

# Index BAM files
samtools index poolA.sorted.bam
samtools index poolB.sorted.bam
```

### 1.2 STAR Alignment (RNA-seq for BSR-seq)

```bash
# Generate genome index
STAR --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir star_index/ \
    --genomeFastaFiles reference.fasta \
    --sjdbGTFfile annotation.gtf \
    --sjdbOverhang 149

# Align reads (two-pass mode recommended)
# First pass
STAR --runThreadN 16 \
    --genomeDir star_index/ \
    --readFilesIn poolA_R1.fastq.gz poolA_R2.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix poolA_ \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic

# Repeat for poolB
```

---

## Phase 2: Variant Calling

### 2.1 GATK HaplotypeCaller (Recommended)

```bash
# Call variants jointly for all pools
gatk HaplotypeCaller \
    -R reference.fasta \
    -I poolA.sorted.bam \
    -I poolB.sorted.bam \
    -O raw_variants.vcf \
    --emit-ref-confidence GVCF \
    --ploidy 2

# For BSR-seq RNA data, add:
# --dont-use-soft-clipped-bases
# --standard-min-confidence-threshold-for-calling 20
```

### 2.2 bcftools mpileup (Alternative)

```bash
# Joint variant calling
bcftools mpileup -f reference.fasta \
    poolA.sorted.bam poolB.sorted.bam | \
    bcftools call -mv -Oz -o raw_variants.vcf.gz

# Index VCF
bcftools index raw_variants.vcf.gz
```

### 2.3 Variant Filtering

```bash
# GATK hard filtering
gatk VariantFiltration \
    -R reference.fasta \
    -V raw_variants.vcf \
    -O filtered_variants.vcf \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "SOR > 3.0" --filter-name "SOR3"

# Keep only PASS variants
bcftools view -f PASS filtered_variants.vcf -Oz -o final_variants.vcf.gz
```

---

## Phase 3: SNP-index Calculation

### 3.1 MutMap SNP-index

```python
# Python script for MutMap SNP-index calculation
import pysam
import numpy as np

def calculate_mutmap_snpindex(vcf_file, pool_name, min_dp=3, max_dp=100):
    """
    Calculate SNP-index for MutMap analysis.
    SNP-index = proportion of alternate (mutant) allele in mutant pool.
    For causal mutations, SNP-index should approach 1.0.
    """
    results = []
    vcf = pysam.VariantFile(vcf_file)
    
    for rec in vcf:
        if rec.filter.keys() != {'PASS'}:
            continue
        
        # Get depth and allele counts for the mutant pool
        sample = rec.samples[pool_name]
        dp = sample['DP']
        
        if dp < min_dp or dp > max_dp:
            continue
        
        # Count alternate alleles
        alt_count = sum(sample.get('AD', [0, 0])[1:])  # Handle multi-allelic
        snp_index = alt_count / dp if dp > 0 else 0
        
        results.append({
            'chrom': rec.chrom,
            'pos': rec.pos,
            'snp_index': snp_index,
            'depth': dp
        })
    
    return results

# Usage
# snp_indices = calculate_mutmap_snpindex('final_variants.vcf.gz', 'poolA')
```

### 3.2 QTL-seq ΔSNP-index

```python
def calculate_delta_snpindex(vcf_file, poolA_name, poolB_name, min_dp=3):
    """
    Calculate ΔSNP-index for QTL-seq analysis.
    ΔSNP-index = SNP-index_poolA - SNP-index_poolB
    At QTL positions, ΔSNP-index deviates significantly from 0.
    """
    results = []
    vcf = pysam.VariantFile(vcf_file)
    
    for rec in vcf:
        if rec.filter.keys() != {'PASS'}:
            continue
        
        # Get allele counts for both pools
        sampleA = rec.samples[poolA_name]
        sampleB = rec.samples[poolB_name]
        
        dpA = sampleA['DP']
        dpB = sampleB['DP']
        
        if dpA < min_dp or dpB < min_dp:
            continue
        
        # Calculate SNP-index for each pool
        # Assuming parent A allele is reference
        alt_countA = sum(sampleA.get('AD', [0, 0])[1:])
        alt_countB = sum(sampleB.get('AD', [0, 0])[1:])
        
        snp_indexA = alt_countA / dpA if dpA > 0 else 0
        snp_indexB = alt_countB / dpB if dpB > 0 else 0
        
        delta_snp_index = snp_indexA - snp_indexB
        
        results.append({
            'chrom': rec.chrom,
            'pos': rec.pos,
            'snp_indexA': snp_indexA,
            'snp_indexB': snp_indexB,
            'delta_snp_index': delta_snp_index,
            'dpA': dpA,
            'dpB': dpB
        })
    
    return results
```

### 3.3 Sliding Window Analysis

```python
import pandas as pd

def sliding_window_average(snp_data, window_size=1000000, step_size=100000, min_snps=10):
    """
    Calculate sliding window average of SNP-index or ΔSNP-index.
    """
    df = pd.DataFrame(snp_data)
    results = []
    
    for chrom in df['chrom'].unique():
        chrom_data = df[df['chrom'] == chrom].sort_values('pos')
        positions = chrom_data['pos'].values
        
        for start in range(1, positions.max(), step_size):
            end = start + window_size
            window_snps = chrom_data[(chrom_data['pos'] >= start) & 
                                      (chrom_data['pos'] < end)]
            
            if len(window_snps) < min_snps:
                continue
            
            # Calculate average
            avg_snp_index = window_snps['snp_index'].mean() if 'snp_index' in window_snps else \
                           window_snps['delta_snp_index'].mean()
            
            results.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'midpoint': (start + end) // 2,
                'avg_snp_index': avg_snp_index,
                'n_snps': len(window_snps)
            })
    
    return pd.DataFrame(results)
```

---

## Phase 4: Confidence Intervals

### 4.1 Permutation Test

```python
import numpy as np
from scipy import stats

def permutation_confidence_interval(snp_data, n_permutations=1000, 
                                    confidence_level=0.99, window_size=1000000):
    """
    Calculate confidence intervals by permutation.
    Randomly shuffle SNP positions to generate null distribution.
    """
    df = pd.DataFrame(snp_data)
    
    # Observed window averages
    observed = sliding_window_average(df, window_size=window_size)
    
    # Permutation
    permuted_max = []
    permuted_min = []
    
    for i in range(n_permutations):
        # Shuffle positions within chromosomes
        shuffled = df.copy()
        for chrom in shuffled['chrom'].unique():
            mask = shuffled['chrom'] == chrom
            shuffled.loc[mask, 'pos'] = np.random.permutation(shuffled.loc[mask, 'pos'])
        
        # Calculate window averages for shuffled data
        perm_result = sliding_window_average(shuffled, window_size=window_size)
        permuted_max.append(perm_result['avg_snp_index'].max())
        permuted_min.append(perm_result['avg_snp_index'].min())
    
    # Calculate confidence thresholds
    alpha = (1 - confidence_level) / 2
    upper_threshold = np.percentile(permuted_max, (1 - alpha) * 100)
    lower_threshold = np.percentile(permuted_min, alpha * 100)
    
    return {
        'upper_ci': upper_threshold,
        'lower_ci': lower_threshold,
        'observed': observed
    }
```

### 4.2 Binomial Normal Approximation

```python
def binomial_ci_normal(n, p, confidence_level=0.95):
    """
    Calculate binomial confidence interval using normal approximation.
    For SNP-index, n = depth, p = SNP-index.
    """
    from scipy import stats
    
    z = stats.norm.ppf((1 + confidence_level) / 2)
    
    se = np.sqrt(p * (1 - p) / n)
    lower = max(0, p - z * se)
    upper = min(1, p + z * se)
    
    return lower, upper
```

---

## Phase 5: Candidate Region Identification

### 5.1 Identify Significant Peaks

```python
def identify_candidate_regions(window_data, upper_threshold, lower_threshold,
                               min_consecutive=3, max_gap=500000):
    """
    Identify candidate regions exceeding confidence thresholds.
    """
    candidates = []
    
    for chrom in window_data['chrom'].unique():
        chrom_data = window_data[window_data['chrom'] == chrom].sort_values('midpoint')
        
        # Find windows exceeding threshold
        if 'delta_snp_index' in chrom_data.columns:
            significant = chrom_data[
                (chrom_data['avg_snp_index'] > upper_threshold) |
                (chrom_data['avg_snp_index'] < lower_threshold)
            ]
        else:
            significant = chrom_data[chrom_data['avg_snp_index'] > upper_threshold]
        
        if len(significant) == 0:
            continue
        
        # Merge consecutive windows
        current_region = None
        for _, window in significant.iterrows():
            if current_region is None:
                current_region = {
                    'chrom': chrom,
                    'start': window['start'],
                    'end': window['end'],
                    'max_snp_index': window['avg_snp_index'],
                    'n_windows': 1
                }
            elif window['start'] - current_region['end'] <= max_gap:
                current_region['end'] = window['end']
                current_region['max_snp_index'] = max(
                    current_region['max_snp_index'], window['avg_snp_index']
                )
                current_region['n_windows'] += 1
            else:
                if current_region['n_windows'] >= min_consecutive:
                    candidates.append(current_region)
                current_region = {
                    'chrom': chrom,
                    'start': window['start'],
                    'end': window['end'],
                    'max_snp_index': window['avg_snp_index'],
                    'n_windows': 1
                }
        
        if current_region and current_region['n_windows'] >= min_consecutive:
            candidates.append(current_region)
    
    return pd.DataFrame(candidates)
```

---

## Phase 6: BSR-seq Specific Analysis

### 6.1 Euclidean Distance (ED)

```python
def calculate_ed(sampleA_depth, sampleB_depth, sampleA_alt, sampleB_alt):
    """
    Calculate Euclidean Distance for BSR-seq.
    Higher ED indicates greater differentiation between pools.
    """
    # Normalize by total depth
    freqA = sampleA_alt / sampleA_depth if sampleA_depth > 0 else 0
    freqB = sampleB_alt / sampleB_depth if sampleB_depth > 0 else 0
    
    # Euclidean distance
    ed = np.sqrt((freqA - freqB) ** 2)
    
    # Normalized ED (optional)
    ed_normalized = ed / np.sqrt(2)  # Max ED is sqrt(2)
    
    return ed, ed_normalized
```

### 6.2 G-statistic

```python
from scipy.stats import chi2

def calculate_g_statistic(sampleA_ref, sampleA_alt, sampleB_ref, sampleB_alt):
    """
    Calculate G-statistic for allele frequency difference.
    G = 2 * sum(O * ln(O/E))
    """
    # Observed counts
    obs = np.array([[sampleA_ref, sampleA_alt],
                    [sampleB_ref, sampleB_alt]])
    
    # Expected counts under null (equal frequencies)
    row_sums = obs.sum(axis=1)
    col_sums = obs.sum(axis=0)
    total = obs.sum()
    
    expected = np.outer(row_sums, col_sums) / total
    
    # G-statistic
    # Avoid log(0)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = obs / expected
        ratio[ratio == 0] = 1  # log(1) = 0
        g_stat = 2 * np.sum(obs * np.log(ratio))
    
    # P-value from chi-square distribution
    p_value = 1 - chi2.cdf(g_stat, df=1)
    
    return g_stat, p_value
```

---

## Phase 7: Visualization

### 7.1 Manhattan Plot

```python
import matplotlib.pyplot as plt

def plot_manhattan(window_data, threshold=None, output_file='manhattan.pdf'):
    """
    Create Manhattan plot of SNP-index or ΔSNP-index.
    """
    fig, ax = plt.subplots(figsize=(12, 4))
    
    chromosomes = sorted(window_data['chrom'].unique())
    colors = ['#1f77b4', '#ff7f0e'] * ((len(chromosomes) + 1) // 2)
    
    offset = 0
    x_ticks = []
    x_labels = []
    
    for i, chrom in enumerate(chromosomes):
        chrom_data = window_data[window_data['chrom'] == chrom]
        x_pos = chrom_data['midpoint'] + offset
        
        ax.scatter(x_pos, chrom_data['avg_snp_index'],
                   c=colors[i], s=5, alpha=0.7)
        
        x_ticks.append(offset + chrom_data['midpoint'].median())
        x_labels.append(chrom)
        offset += chrom_data['midpoint'].max() + 1000000
    
    if threshold:
        ax.axhline(y=threshold['upper'], color='red', linestyle='--', linewidth=1)
        ax.axhline(y=threshold['lower'], color='red', linestyle='--', linewidth=1)
    
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, rotation=45)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('SNP-index' if 'delta' not in window_data.columns else 'ΔSNP-index')
    ax.set_title('BSA Analysis')
    
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
```

### 7.2 Highlight Candidate Regions

```python
def plot_with_candidates(window_data, candidates, threshold, output_file='bsa_candidates.pdf'):
    """
    Plot SNP-index with candidate regions highlighted.
    """
    fig, ax = plt.subplots(figsize=(14, 5))
    
    # Plot SNP-index
    # ... (similar to manhattan plot)
    
    # Highlight candidate regions
    for _, region in candidates.iterrows():
        ax.axvspan(region['start'], region['end'], 
                   alpha=0.3, color='red',
                   label='Candidate region')
    
    plt.savefig(output_file)
    plt.close()
```

---

## Phase 8: Candidate Gene Extraction

### 8.1 Extract Genes in Candidate Regions

```python
import gffutils

def extract_candidate_genes(candidates, gff_file, upstream=5000, downstream=5000):
    """
    Extract genes within candidate regions from GFF file.
    """
    # Create temporary database
    db = gffutils.create_db(gff_file, ':memory:', 
                            force=True, keep_order=True)
    
    candidate_genes = []
    
    for _, region in candidates.iterrows():
        # Find genes overlapping the region
        region_start = region['start'] - upstream
        region_end = region['end'] + downstream
        
        for gene in db.region(seqid=region['chrom'],
                              start=region_start,
                              end=region_end,
                              featuretype='gene'):
            candidate_genes.append({
                'chrom': gene.chrom,
                'start': gene.start,
                'end': gene.end,
                'gene_id': gene.id,
                'strand': gene.strand,
                'region': f"{region['chrom']}:{region['start']}-{region['end']}"
            })
    
    return pd.DataFrame(candidate_genes)
```

### 8.2 Variant Annotation with SnpEff

```bash
# Annotate variants
snpEff -v plant_species final_variants.vcf > annotated_variants.vcf

# Extract HIGH and MODERATE impact variants in candidate regions
# Using bcftools
bcftools view -i 'INFO/ANN[0].IMPACT="HIGH" || INFO/ANN[0].IMPACT="MODERATE"' \
    -R candidate_regions.bed annotated_variants.vcf > candidate_variants.vcf
```

---

## FAQ

### Q1: SNP-index values are noisy

**Solutions:**
- Increase sliding window size
- Filter low-depth SNPs more stringently
- Check for sequencing quality issues
- Consider using larger pool size

### Q2: No clear peaks detected

**Possible causes:**
- Low heritability trait
- Multiple small-effect QTL
- Insufficient pool size
- Low sequencing depth

**Solutions:**
- Increase pool size
- Increase sequencing depth
- Use complementary mapping approach
- Fine-map with additional markers

### Q3: Multiple peaks across genome

**Possible causes:**
- Multiple QTL
- Population structure
- Linked loci

**Solutions:**
- Verify population is from single cross
- Check for seed contamination
- Consider epistatic interactions
