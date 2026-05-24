# Fan-Skill 分析能力目录

30 项分析条目，按研究目标分类组织。

## 寻找控制性状的基因

| 条目 | 功能 | 典型输入 |
|------|------|------|
| `gwas` | 全基因组关联分析——扫描性状相关位点 | VCF + 表型 |
| `qtl-mapping` | 基于连锁的 QTL 定位（双亲群体） | 遗传图谱 + 表型 |
| `candidate-gene-association` | 候选基因区域内的关联分析 | VCF + 基因列表 + 表型 |
| `eqtl` | 表达 QTL——与表达相关的变异 | 基因型 + 表达矩阵 |
| `population` | 群体结构——PCA、ADMIXTURE、Fst、系统发育 | VCF/PLINK |

## 理解基因功能

| 条目 | 功能 | 典型输入 |
|------|------|------|
| `rnaseq` | 差异表达分析 | 计数矩阵 + 元数据 |
| `time-series` | 发育/时间序列模式分析 | 表达量 + 时间点 |
| `grn` | 基因调控网络推断 | 表达量 + 转录因子列表 |
| `small-rna` | miRNA 预测、靶基因预测 | sRNA-seq 读数 + 参考基因组 |
| `multi-omics` | 多组学整合分析 | 多个组学数据集 |

## 预测育种值

| 条目 | 功能 | 典型输入 |
|------|------|------|
| `genomic-selection` | 预测育种值 | 基因型 + 表型 |
| `hybrid-prediction` | 杂种优势预测、GCA/SCA | 亲本基因型 + 杂种表型 |
| `phenotype` | 遗传力、BLUP/BLUE、MET 汇总 | 表型 + 试验设计 |
| `enviromics` | 环境 + 基因组整合分析 | 气候 + 基因型 + 表型 |

## 数据处理

| 条目 | 功能 | 典型输入 |
|------|------|------|
| `variant-calling` | SNP/InDel 检测 | FASTQ + 参考基因组 |
| `genotype-imputation` | 缺失基因型填补 | VCF + 参考面板 |
| `genome-assembly` | 从头组装（HiFi、ONT、混合） | 测序读数 |
| `genome-annotation` | 基因预测、功能注释 | 组装结果 + RNA-seq 证据 |

## 表观基因组与染色质

| 条目 | 功能 | 典型输入 |
|------|------|------|
| `chipseq` | ChIP-seq 峰检测、差异结合分析 | ChIP + input FASTQ |
| `atacseq` | 开放染色质区域鉴定 | ATAC-seq FASTQ |
| `methylation` | DNA 甲基化（WGBS/RRBS） | 亚硫酸氢盐测序读数 |

## 微生物组与宏基因组

| 条目 | 功能 | 典型输入 |
|------|------|------|
| `amplicon` | 16S/ITS 扩增子分析 | 扩增子 FASTQ |
| `metagenomics` | 宏基因组组装、分箱、MAG | 宏基因组读数 |

## 代谢组与蛋白质组

| 条目 | 功能 | 典型输入 |
|------|------|------|
| `metabolomics` | LC-MS/GC-MS、差异代谢物分析 | 质谱数据 |
| `proteomics` | DDA/DIA 定量、差异蛋白、PPI | 质谱数据 |

## 进化与比较基因组学

| 条目 | 功能 | 典型输入 |
|------|------|------|
| `comparative` | 共线性、Ks、基因家族、选择压力 | 基因组 + 注释 |
| `pan-genome` | 核心/可变基因组、PAV 检测 | 多个组装结果 |

## 基因组编辑与标记

| 条目 | 功能 | 典型输入 |
|------|------|------|
| `crispr` | sgRNA 设计、脱靶预测 | 目标基因序列 |
| `marker` | KASP/InDel/SSR 开发、亲本推荐 | GWAS 峰值 / QTL 区间 |

## 可视化

| 条目 | 功能 | 典型输入 |
|------|------|------|
| `visualization` | 发表级别的图表 | 分析输出 |

## 植物特有功能

- **多倍体感知**：自动检测、亚基因组特异性分析
- **自交/异交**：不同的 GWAS 和 GS 策略
- **12 种作物速查表**：参考基因组、LD 衰减、育种系统
- **非模式物种策略**：Mercator 注释、跨物种推断
- **多环境试验（MET）**：植物育种的标准配置
