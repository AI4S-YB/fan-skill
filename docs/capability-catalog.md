# Fan-Skill Analysis Capability Catalog / 分析能力目录

## By Research Goal / 按研究目标

### 🔍 Find Genes Controlling Traits / 找控制性状的基因

| Entry / 条目 | What it does / 做什么 | Typical input / 典型输入 |
|------|------|------|
| `gwas` | Genome-wide association — scan whole genome for trait-associated loci / 全基因组关联扫描 | VCF + phenotype data |
| `qtl-mapping` | Linkage-based QTL mapping for biparental populations / 双亲群体连锁QTL定位 | Genetic map + phenotype |
| `candidate-gene-association` | Association analysis within pre-selected gene regions / 候选基因区域关联分析 | VCF + candidate gene list + phenotype |
| `eqtl` | Expression QTL — genetic variants associated with gene expression / 遗传变异与基因表达的关联 | Genotype + expression matrix |
| `population` | Population structure — PCA, ADMIXTURE, Fst, phylogenetic tree / 群体遗传结构 | VCF/PLINK binary |

### 🧬 Understand Gene Function / 理解基因功能

| Entry / 条目 | What it does / 做什么 | Typical input / 典型输入 |
|------|------|------|
| `rnaseq` | Differential expression — which genes change under treatment / 差异表达分析 | Count matrix + sample metadata |
| `time-series` | Developmental/time-course expression patterns / 时间序列表达模式 | Expression matrix + time points |
| `grn` | Gene regulatory network inference / 基因调控网络推断 | Expression matrix + TF list |
| `small-rna` | miRNA prediction, target prediction, differential miRNA / 小RNA分析 | sRNA-seq reads + reference |
| `multi-omics` | Multi-omics integration (transcriptome + metabolome + proteome) / 多组学整合 | Multiple omics datasets |

### 🌾 Predict Breeding Value / 预测育种值

| Entry / 条目 | What it does / 做什么 | Typical input / 典型输入 |
|------|------|------|
| `genomic-selection` | Predict breeding values using genome-wide markers / 基因组选择预测育种值 | Genotype + phenotype |
| `hybrid-prediction` | Predict hybrid performance, GCA/SCA, heterotic groups / 杂种优势预测 | Parental genotypes + hybrid phenotypes |
| `phenotype` | Heritability, BLUP/BLUE, spatial correction, MET summary / 表型数据分析 | Phenotype + trial design |
| `enviromics` | Environmental + genomic integration, G×E modeling / 环境组+基因组整合 | Climate data + genotype + phenotype |

### 🛠️ Data Processing / 数据预处理

| Entry / 条目 | What it does / 做什么 | Typical input / 典型输入 |
|------|------|------|
| `variant-calling` | SNP/InDel discovery from resequencing / 变异检测 | FASTQ + reference genome |
| `genotype-imputation` | Fill missing genotypes using reference panels / 基因型填充 | VCF + reference panel |
| `genome-assembly` | De novo genome assembly (HiFi, ONT, hybrid) / 基因组从头组装 | Sequencing reads |
| `genome-annotation` | Gene prediction, functional annotation / 基因组注释 | Assembly + RNA-seq evidence |

### 🔬 Epigenomics & Chromatin / 表观组学

| Entry / 条目 | What it does / 做什么 | Typical input / 典型输入 |
|------|------|------|
| `chipseq` | ChIP-seq peak calling, differential binding / 蛋白-DNA结合分析 | ChIP + input FASTQ |
| `atacseq` | Open chromatin region identification / 染色质开放区域鉴定 | ATAC-seq FASTQ |
| `methylation` | DNA methylation analysis (WGBS/RRBS) / DNA甲基化分析 | Bisulfite-seq reads |

### 🦠 Microbiome & Metagenomics / 微生物组

| Entry / 条目 | What it does / 做什么 | Typical input / 典型输入 |
|------|------|------|
| `amplicon` | 16S/ITS amplicon analysis / 扩增子分析 | Amplicon FASTQ |
| `metagenomics` | Shotgun metagenomics — assembly, binning, MAGs / 宏基因组组装分箱 | Metagenomic reads |

### 🧪 Metabolomics & Proteomics / 代谢与蛋白组

| Entry / 条目 | What it does / 做什么 | Typical input / 典型输入 |
|------|------|------|
| `metabolomics` | LC-MS/GC-MS peak detection, differential metabolites / 代谢组差异分析 | Mass spectrometry data |
| `proteomics` | DDA/DIA protein quantification, differential proteins, PPI / 蛋白组差异分析 | Mass spectrometry data |

### 🧬 Evolution & Comparative Genomics / 进化与比较基因组

| Entry / 条目 | What it does / 做什么 | Typical input / 典型输入 |
|------|------|------|
| `comparative` | Synteny, Ks, gene families, selection pressure / 共线性、Ks、选择压力 | Genome + annotation |
| `pan-genome` | Core/variable genome, PAV detection / 泛基因组分析 | Multiple genome assemblies |

### ✂️ Genome Editing & Markers / 基因编辑与标记

| Entry / 条目 | What it does / 做什么 | Typical input / 典型输入 |
|------|------|------|
| `crispr` | sgRNA design, off-target prediction / CRISPR靶点设计 | Target gene sequence |
| `marker` | KASP/InDel/SSR marker development, parental recommendation / 分子标记开发 | GWAS peak / QTL interval |

### 📊 Visualization / 可视化

| Entry / 条目 | What it does / 做什么 | Typical input / 典型输入 |
|------|------|------|
| `visualization` | Publication-quality figures from analysis results / 发表级科研图表 | Analysis output files |

---

## Plant-Specific Capabilities / 植物特化能力

- **Polyploid-aware** / 多倍体感知: auto-detection, subgenome-specific analysis
- **Self/cross-pollinated distinction** / 自交/异交区分: different GWAS and GS strategies
- **12 crop species cheatsheet** / 12 作物速查: reference genomes, LD decay, breeding systems
- **Non-model species strategy** / 非模式物种策略: Mercator annotation, cross-species inference
- **Multi-environment trial (MET)** / 多环境试验: standard in plant breeding, built into phenotype + GWAS + GS

---

## Architecture / 架构

```
Each entry = 4 modules / 每个条目 = 4 模块:

consult-guide.md  →  "What to ask the user" / "问用户什么"
rules.yaml        →  "Which method to use" / "选什么方法" (C-layer + design gates)
notebook.md       →  "Why this method" / "为什么" (B-layer expert reasoning)
analysis-primer.md → "What the results mean" / "结果解读"
```

Fan-skill uses a B+C dual-mode decision architecture:
- **C-layer (rules.yaml)**: Deterministic — same data → same method. Rule engine with priority matching.
- **B-layer (notebook.md)**: Flexible — narrative expert reasoning for edge cases.
- **User control**: `decision_mode: rule | expert | hybrid`

Fan-skill 使用 B+C 双模式决策架构：
- **C 层 (rules.yaml)**: 确定性——同样的数据选同样的方法。优先级匹配规则引擎。
- **B 层 (notebook.md)**: 灵活性——叙事式专家推理处理规则盲区。
- **用户控制**: `decision_mode: rule | expert | hybrid`
