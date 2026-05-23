# 植物蛋白质组学分析：分析笔记

## 在开始之前

植物蛋白质组学的核心挑战来自三个方面：样品复杂性（次生代谢物干扰）、数据库完整性（非模式物种注释不足）、以及翻译后修饰的多样性。先问自己：

### 你用了什么采集模式？

- **DDA (Data-Dependent Acquisition)**：传统模式，每个循环选 top-N 离子做 MS/MS → MaxQuant 分析
- **DIA (Data-Independent Acquisition)**：所有离子循环碎裂 → DIA-NN 或 Spectronaut 分析
- **TMT (Tandem Mass Tag)**：多重标记，可同时比较多达 18 个样品 → MaxQuant TMT 模式
- **PRM (Parallel Reaction Monitoring)**：靶向定量，验证候选蛋白

### 你的植物物种有好的蛋白数据库吗？

- **模式植物**（拟南芥、水稻、玉米、大豆、番茄、小麦）→ UniProt 参考蛋白质组完整，直接使用
- **非模式植物**→ 可能需要：
  - 从基因组注释预测蛋白质组
  - 用转录组（RNA-seq）构建样本特异性蛋白数据库
  - BLAST 到近缘物种的蛋白数据库

### 你想要什么样的结果？

- **全局蛋白组变化** → LFQ (Label-Free Quantification)
- **精确蛋白丰度变化** → TMT 更准（无 missing value 问题）
- **特定通路的磷酸化变化** → 磷酸化蛋白质组学
- **蛋白互作网络** → 差异蛋白 → STRING PPI → Cytoscape 网络分析

## 选择你的定量工具

### MaxQuant — DDA 标准

基于特征检测和 match-between-runs 的 label-free quantitation。

**适合**：DDA 数据，≥3 重复
**Plant 特殊优点**：match-between-runs 极大减少植物样本的 missing values（植物次生代谢物会干扰离子化）

### DIA-NN — DIA 标准

基于深度学习的 DIA 数据分析。

**适合**：DIA 数据，所有重复数
**Plant 优势**：产生极低 missing value 率的定量矩阵 → 对植物来说这是一大优势

### MaxQuant TMT 模式

**适合**：需要精确定量的实验，有 TMT 标记数据
**注意**：植物样本中 tryptophan 和 methionine 的氧化可能导致 ratio compression

## 预处理和数据清洗

### Missing Value 处理

植物蛋白质组学中 missing values 通常有 10-30%，来源：

1. **随机 missing (MAR)**：离子化抑制 → 用 imputation（如 KNN, MinProb）
2. **非随机 missing (MNAR)**：蛋白在某个条件下真的不存在 → 不应 impute！

**如何处理**：
- 先区分 MAR vs MNAR（检查 peptide intensity 分布）
- MAR：KNN imputation (k=10)
- MNAR：保留为 NA，或 MinDet imputation（左尾截断分布）

### 归一化

植物样本间蛋白量差异可能很大（如不同组织的 Rubisco 含量差 10 倍），需要仔细选择归一化：
- **Median normalization**：适合大多数植物实验
- **Quantile normalization**：适合样本蛋白分布差异大的情况
- **VSN (Variance Stabilizing Normalization)**：适合 low-intensity 蛋白多的情况

## 差异蛋白分析

### limma — 植物蛋白质组学默认选择

limma 的 empirical Bayes moderation 对蛋白质组学非常适用：

```r
library(limma)

# 模型矩阵
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)

# 拟合
fit <- lmFit(lfq_matrix, design)
contrast <- makeContrasts(Treatment - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# 结果
results <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
```

**关键阈值**：
- P-value adjusted < 0.05（FDR）
- |log2FC| > 1（但植物中 fold change 通常较小，0.58 也是合理的阈值）
- 与至少 2 个 unique peptides 匹配

### 植物特异性考虑

- 光合作用蛋白（特别是 Rubisco, LHC）通常丰度最高 → 在前处理时可能被 depletion 处理掉
- Rubisco 占了叶片总蛋白的 30-50% → 如果要看低丰度蛋白，必须做 depletion
- 次生代谢物（多酚、色素类）会干扰蛋白提取 → 需要 TCA/丙酮沉淀或 phenol 提取

## PPI 网络分析

### STRING DB 植物使用

STRING 数据库覆盖主要植物物种：
- 拟南芥 (Arabidopsis thaliana)
- 水稻 (Oryza sativa)
- 玉米 (Zea mays)
- 大豆 (Glycine max)
- 番茄 (Solanum lycopersicum)

**非模式植物怎么办？**
- 使用近缘模式物种的 STRING 网络（通过序列相似性映射）
- 或使用 co-expression 数据（如 ATTED-II for 拟南芥）构建功能网络

## PTM 分析

### 磷酸化蛋白组学

植物磷酸化是信号转导的核心。常见富集方法：
- TiO2 富集：适合植物，但会同时富集非磷酸化酸性肽段
- IMAC (Fe3+)：选择性更好，但对植物次生代谢物敏感
- 混合模式：TiO2 + 乳酸可显著改善特异性

### 分析工具

- **Motif-x**：从磷酸化位点周围提取 kinase 识别基序
- **kinase-substrate prediction**：预测哪些 kinase 磷酸化了这些位点
- **PhosphoSitePlus**：包含植物磷酸化数据（主要是拟南芥）

## 常见陷阱

### Rubisco 的"诅咒"

Rubisco (Ribulose-1,5-bisphosphate carboxylase/oxygenase) 是地球上最丰富的蛋白，在植物绿色组织中占总蛋白 30-50%。

- 如果不 depletion，低丰度蛋白几乎检测不到
- Depletion 方法：抗 Rubisco 抗体免疫亲和柱（拟南芥/水稻）或 PEG 沉淀法
- 如果没有 depletion，至少报告 Rubisco 占总谱图的百分比

### 多倍体的蛋白推断

同源基因编码的蛋白可能几乎相同（例如小麦 A/B/D 亚基因组的 Rubisco 活化酶）。
- MaxQuant 会将这些蛋白归入同一个 "protein group"
- 如果确实需要区分亚基因组来源，需要 unique peptides 验证
- 对于非模式多倍体，可能需要建立亚基因组特异性的蛋白数据库

### 数据库搜索空间

植物基因组大 → 蛋白数据库大 → 搜索空间大 → FDR 膨胀

**建议**：
- 添加常见 contaminant 数据库（cRAP: keratin, trypsin 等）
- 如果有 RNA-seq，用翻译后的蛋白序列缩小搜索空间
- 如果使用大数据库（>100K sequence），考虑用更严格的 FDR（1% at protein level）

### 光合作用相关的时间效应

- 采样时间对植物蛋白质组影响极大（光合作用蛋白有明显的昼夜节律）
- 必须在实验设计中记录采样时间
- 如果所有样本在不同时间采集 → 混淆昼夜节律和处理效应
