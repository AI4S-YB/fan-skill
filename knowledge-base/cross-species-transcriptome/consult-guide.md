# Cross-Species Transcriptome Consultation Guide

## Quick Decision Framework

### Step 1: Verify Data Requirements

**Minimum Requirements:**
- [ ] Expression data from ≥2 species
- [ ] Ortholog mapping available or can be computed
- [ ] Sample metadata (tissue/condition labels)

**If any requirement missing →** Stop and request data before proceeding.

### Step 2: Choose Ortholog Mapping Strategy

| Scenario | Recommended Method |
|----------|-------------------|
| Species in Ensembl | biomaRt (quickest) |
| Protein sequences available | OrthoFinder (most accurate) |
| Need quick 1:1 mapping | Reciprocal best hit (BLAST) |
| Species in OrthoDB | Direct OrthoDB mapping |

### Step 3: Determine Normalization Strategy

| Expression Type | Recommended Normalization |
|----------------|--------------------------|
| Raw counts | TMM → log2(CPM+1) |
| TPM available | log2(TPM+1) directly |
| FPKM available | Convert to TPM or use quantile |
| Mixed platforms | Quantile normalization |

### Step 4: Batch Correction Decision

**ALWAYS apply batch correction when:**
- Data from different labs/sources
- Different library preparation protocols
- Different sequencing platforms

**Use ComBat when:** ≥3 samples per species
**Use MNN when:** Non-linear batch effects, high-dimensional data

### Step 5: Analysis Pathway Selection

```
Is tissue annotation available?
├── YES → Tissue-specificity analysis (tau/SPM)
│         └── Compare specificity across species
└── NO → Global expression comparison
          └── PCA/UMAP + correlation analysis
```

## Common Consultation Scenarios

### Q1: "I want to compare expression between human and mouse tissues"

**Recommended workflow:**
1. Download expression from GTEx (human) and mouse tissue atlas
2. Use biomaRt to get human-mouse 1:1 orthologs
3. Apply ComBat batch correction (species as batch)
4. Calculate tissue-specificity with tau index
5. Compare tau values for ortholog pairs

### Q2: "How do I know if expression is conserved?"

**Metrics to use:**
- Pearson/Spearman correlation of expression profiles
- Tissue-specificity conservation (tau difference)
- Expression divergence rate (phylogenetic contrast)

**Interpretation:**
- r > 0.8: Highly conserved expression
- r 0.5-0.8: Moderately conserved
- r < 0.5: Divergent expression

### Q3: "Can I compare expression levels directly across species?"

**NO** - Absolute expression levels are not comparable due to:
- Different gene lengths
- Different mRNA metabolism rates
- Platform/protocol differences

**Use relative comparisons instead:**
- Expression rank within species
- Z-score normalized expression
- Tissue-specificity indices

### Q4: "How many species can I compare simultaneously?"

| Species Count | Recommendation |
|--------------|----------------|
| 2-5 | Standard ComBat workflow |
| 6-15 | Use ComBat with careful QC |
| >15 | Consider phylogenetic methods (phylogenetic contrasts) |

## Risk Assessment Checklist

### High Risk (Results may be unreliable)
- [ ] No ortholog information available
- [ ] Mixed RNA-seq platforms without batch correction
- [ ] Single sample per species per tissue
- [ ] Unknown tissue homology

### Medium Risk (Proceed with caution)
- [ ] Limited ortholog coverage (<70% genes)
- [ ] Different annotation versions between species
- [ ] Batch effects not fully correctable

### Low Risk (Standard analysis)
- [ ] Well-annotated species with 1:1 orthologs
- [ ] Biological replicates available
- [ ] Comparable tissue sampling across species

## Output Interpretation Guide

### Tissue-Specificity (tau)
| tau Value | Interpretation |
|-----------|---------------|
| 0.0-0.3 | Ubiquitous expression |
| 0.3-0.7 | Moderately specific |
| 0.7-1.0 | Tissue-specific |

### Expression Conservation
| Correlation | Interpretation |
|-------------|---------------|
| >0.9 | Strong conservation |
| 0.7-0.9 | Moderate conservation |
| <0.7 | Expression divergence |

## When to Escalate

Escalate to domain expert when:
- Novel ortholog mapping required (no established resources)
- Complex polyploid species comparisons
- Cross-kingdom comparisons (plant vs animal)
- Non-model species without reference genomes

## 数据状态标签规范 (C4强制)

在生成方案的数据画像部分，必须显式标注数据状态：
- **数据状态**: FULL — [说明哪些数据完整可用，如"多物种表达数据、直系同源映射信息、样本metadata完整"]
- **数据状态**: PARTIAL — [说明哪些数据缺失，如"表达数据完整但缺乏直系同源映射；或组织注释不完整"]
- **数据状态**: EMPTY — [说明数据不可用原因，如"跨物种数据尚未收集"]

违例判定: 仅列出文件名/大小但无显式FULL/PARTIAL/EMPTY状态标签 → C4=0分
