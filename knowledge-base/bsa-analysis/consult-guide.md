# BSA Analysis Consultation Guide

## Quick Decision Framework

### Step 1: Identify BSA Method

| Scenario | Method |
|----------|--------|
| Mutant mapping with mutant parent | MutMap |
| QTL mapping in biparental cross | QTL-seq |
| Only RNA-seq data available | BSR-seq |
| Multiple traits | QTL-seq (multiple pairs) |

### Step 2: Verify Data Quality

**Minimum Requirements:**
- [ ] Two pools with contrasting phenotypes
- [ ] Reference genome available
- [ ] Pool size ≥ 20 individuals
- [ ] Sequencing depth ≥ 20x per pool

**If any requirement missing →** Assess feasibility and recommend improvements.

### Step 3: Choose Analysis Pipeline

| Data Type | Recommended Pipeline |
|-----------|---------------------|
| DNA-seq, Illumina | BWA-MEM → GATK → QTL-seq |
| DNA-seq, quick analysis | BWA-MEM → bcftools → custom |
| RNA-seq | STAR → GATK → BSR-seq |
| Multiple platforms | Combine carefully, account for batch |

## Common Consultation Scenarios

### Q1: "I have an F2 population with a visible mutant phenotype"

**Recommended approach:**
1. If mutant parent available → **MutMap**
2. Cross mutant × wild-type
3. Select ~50 F2 mutants for pool
4. Sequence to 30-50x depth
5. Look for SNP-index > 0.9

### Q2: "I want to map QTL for a quantitative trait"

**Recommended approach:**
1. Use **QTL-seq**
2. Select extreme individuals (top/bottom 20%)
3. Minimum 30 individuals per pool
4. Sequence both pools to 30x
5. Calculate ΔSNP-index with 99% CI

### Q3: "Can I use RNA-seq for BSA?"

**Yes, with caveats:**
- Use **BSR-seq** method
- Only genes expressed in sampled tissue are informative
- Combine allele frequency with expression differences
- Good for identifying regulatory mutations

### Q4: "How do I handle polyploid species?"

**Special considerations:**
- Use subgenome-aware alignment
- Filter homoeologous regions carefully
- May need higher sequencing depth
- Consider using haplotype-specific markers

## Parameter Selection Guide

### Window Size

| Genome Size | Recombination Rate | Window Size |
|-------------|-------------------|-------------|
| < 500 Mb | High | 500 kb |
| 500 Mb - 1 Gb | Moderate | 1 Mb |
| > 1 Gb | Low | 2 Mb |

### Confidence Level

| Stage | Confidence Level |
|-------|-----------------|
| Initial scan | 95% |
| Publication | 99% |
| Multiple testing correction | 99.9% |

## Risk Assessment

### High Risk
- Pool size < 20
- Depth < 15x
- No biological replicates
- Complex polyploid

### Medium Risk
- Pool size 20-30
- Depth 15-25x
- Single pair of pools

### Low Risk
- Pool size > 50
- Depth > 30x
- Multiple independent pools

## Result Interpretation

### Good Signal
- Clear peak(s) exceeding CI
- Peak width consistent with recombination
- Single or few peaks

### Poor Signal
- No regions exceed CI
- Many small peaks across genome
- Signal matches genome-wide noise

### Troubleshooting Poor Results

| Issue | Possible Cause | Solution |
|-------|---------------|----------|
| No peaks | Low heritability | Increase pool size |
| Too many peaks | Population structure | Check population history |
| Broad peaks | Low recombination | Fine mapping needed |
| Inconsistent results | Batch effects | Re-sequence together |

## When to Escalate

- Complex polyploid genetics
- Multiple interacting QTL
- Epistatic effects expected
- Non-Mendelian inheritance
- Structural variants suspected

## 数据状态标签规范 (C4强制)

在生成方案的数据画像部分，必须显式标注数据状态：
- **数据状态**: FULL — [说明哪些数据完整可用，如"两个表型极端池的DNA-seq数据(各30x)、参考基因组完整"]
- **数据状态**: PARTIAL — [说明哪些数据缺失，如"仅有一个池的数据，缺乏极端表型池"]
- **数据状态**: EMPTY — [说明数据不可用原因，如"群体尚未构建"]

违例判定: 仅列出文件名/大小但无显式FULL/PARTIAL/EMPTY状态标签 → C4=0分
