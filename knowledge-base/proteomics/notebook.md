# 植物蛋白质组学分析：分析笔记

## 在开始之前

植物蛋白质组学的核心挑战来自三个方面：样品复杂性（次生代谢物干扰）、数据库完整性（非模式物种注释不足）、以及翻译后修饰的多样性。先问自己：

### 你用了什么采集模式？

- **DDA (Data-Dependent Acquisition)**：传统模式，每个循环选 top-N 离子做 MS/MS → MaxQuant 分析
- **DIA (Data-Independent Acquisition)**：所有离子循环碎裂 → DIA-NN 或 Spectronaut 分析
- **TMT (Tandem Mass Tag)**：多重标记，可同时比较多达 18 个样品 → MaxQuant TMT 模式
- **PRM (Parallel Reaction Monitoring)**：靶向定量，验证候选蛋白

### DDA vs DIA 决策指南

选择 DDA 还是 DIA 取决于实验目标和资源：

| 特征 | DDA | DIA |
|------|-----|-----|
| 数据分析复杂度 | 低（MaxQuant 成熟） | 中等（需 DIA-NN 或 Spectronaut） |
| Missing values | 高（10-30%） | 极低（<5%） |
| 蛋白鉴定数量 | ~3000-4000（植物组织） | ~3500-5000（植物组织） |
| 重复性 | 中（stochastic precursor selection） | 高（所有离子都碎裂） |
| 数据分析深度 | 浅 — 标准分析即可 | 深 — 光谱库方法或 library-free 方法 |
| 植物适用性 | 经典选项，流程成熟 | 植物 matrix 中 missing value 少是一大优势 |

**植物 DDA 选择的理由**：
- 样品数量少（<20），重复数多（>4/组）
- 实验室已有 DDA 分析管道
- 非模式物种（DIA 的光谱库可能不完整）

**植物 DIA 选择的理由**：
- 样品数量多（>50），需要矩阵完整性好的定量结果
- 进行时间序列实验（missing values 会严重破坏趋势分析）
- 植物次生代谢物多 → DDA 的随机离子选择在干扰环境中更差

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

## Label-Free vs TMT 标记：成本效益分析

### Label-Free Quantification (MaxQuant LFQ)

**优点**：
- 样品制备简单，无需额外的标记试剂
- 每个样品独立上机 — 样品数量灵活
- 成本低（每个样品 ~$100-200）

**缺点**：
- Missing values 高（植物样品尤其明显）
- 运行间变异性大 → 需要更多的生物学重复
- 低丰度蛋白的定量不精确

**植物适用**：
- 叶片、根等次生代谢物较少的组织
- 大规模筛选实验（如数百个品系的蛋白组谱）
- 非模式物种（标记化学对非模式物种可能行为不同）

### TMT Labeling

**优点**：
- 无 missing values — 所有通道在同一 MS run 中定量
- 定量精度高（CV < 10% for medium-abundance proteins）
- 可多重分析（TMTpro 16-plex 允许 16 个样品同时定量）

**缺点**：
- 成本：TMT 试剂每标签 $50-100（16-plex ~$800-1600）
- Ratio compression：共分离前体的干扰信号导致 fold change 被压缩
- 需要更复杂的 MS 仪器（MS3 或 SPS-MS3 for accurate ratio）

**值得使用 TMT 的场景**：
- 精确定量至关重要的实验（如剂量-反应曲线）
- 样品数量恰好匹配标记通道数
- 需要检测微小的蛋白丰度变化（<1.5 fold）
- 实验中有珍贵的、材料量有限的样品（如激光显微切割的特定组织）

### 关键决策：成本阈值

如果研究预算有限：
- <12 样品且追求最高精度 → TMT 可能不值得（成本增幅 >50%）
- >30 样品 → LFQ 是唯一经济合理的选择
- 如果是验证性实验（已有候选蛋白）→ PRM 靶向定量最经济

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

## 植物特异性蛋白提取挑战

### 细胞壁的障碍

植物细胞壁是蛋白提取的第一道难关：
- 纤维素/半纤维素/果胶组成的刚性结构 → 需要物理破碎（液氮研磨、bead beating）
- **推荐方案**：液氮研磨 + 提取缓冲液（含有 0.5-1% SDS 或 urea/thiourea） → 超声破碎

### 次生代谢物的严重干扰

植物组织富含多酚类、色素、萜类等次生代谢物：
- 多酚氧化后与蛋白共价结合 → 蛋白沉淀、蛋白酶抑制
- 色素（叶绿素、花青素）干扰光谱检测和肽段离子化

**三大提取方法比较**：

| 方法 | 适用组织 | 优点 | 缺点 |
|------|---------|------|------|
| TCA/丙酮沉淀 | 叶片、幼苗 | 极好去除叶绿素和多酚 | 蛋白沉淀可能难以复溶 |
| 酚提取 (Phenol extraction) | 富含多酚的组织（果实、树皮） | 高效去除多酚和碳水化合物 | 步骤多，耗时长（~6-8小时） |
| SDS 提取 + FASP | 根、种子 | 蛋白回收率高，可溶性佳 | SDS 残留影响 MS |

**组织特异性推荐**：
- 绿色叶片：TCA/丙酮
- 果实（草莓、葡萄）：酚提取 + 醋酸铵/甲醇沉淀
- 种子（大豆、油菜）：SDS 提取 + FASP（高油脂组织需要先脱脂）
- 根：TCA/丙酮 或 SDS 提取（次生代谢物较少）

### 蛋白降解

植物组织中含有大量蛋白酶（液泡中的 cysteine proteases 在组织破碎后释放）：
- **缓解**：提取缓冲液中加入 protease inhibitor cocktail (含 PMSF, leupeptin, pepstatin)
- **关键步骤**：全程在 4°C 操作，液氮速冻组织后储存于 -80°C

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

**STRING DB 的植物覆盖局限性**：
- 非模式物种的 PPI 信息极其稀疏（仅通过直系同源转移，通常 < 20% 的查询蛋白有网络信息）
- 植物特有的蛋白-蛋白互作（如受体激酶-配体、E3 ligase-底物）在 STRING 中覆盖尤其差
- 亚细胞定位信息和共表达数据（植物的两个重要 PPI 证据来源）在非模式物种中缺失

**非模式植物怎么办？**
- 使用近缘模式物种的 STRING 网络（通过序列相似性映射）
- 或使用 co-expression 数据（如 ATTED-II for 拟南芥）构建功能网络
- **补充方案**：PlantReactome（植物特异性通路数据库）+ CORNET（植物共表达网络工具）

## PTM 分析

### 磷酸化蛋白组学

植物磷酸化是信号转导的核心。常见富集方法：
- TiO2 富集：适合植物，但会同时富集非磷酸化酸性肽段
- IMAC (Fe3+)：选择性更好，但对植物次生代谢物敏感
- 混合模式：TiO2 + 乳酸可显著改善特异性

### 泛素化 (Ubiquitination)

植物中泛素-蛋白酶体系统调控几乎所有发育过程：
- 富集：di-Gly 抗体（识别 trypsin 消化后留下的 K-ε-GG 残余）
- 重要性：E3 泛素连接酶是植物中最大的基因家族之一
- 挑战：di-Gly 抗体特异性有限（也富集 NEDD8 修饰）

### 乙酰化 (Acetylation)

植物蛋白乙酰化参与光合作用调控和代谢调控：
- 富集：anti-acetyl-lysine (Kac) 抗体
- 光合作用相关酶的乙酰化在昼夜转换中动态变化
- 注意：化学乙酰化（非酶促，在碱性条件下自发发生）可能造成背景

### 分析工具

- **Motif-x**：从磷酸化位点周围提取 kinase 识别基序
- **kinase-substrate prediction**：预测哪些 kinase 磷酸化了这些位点
- **PhosphoSitePlus**：包含植物磷酸化数据（主要是拟南芥）

## 细胞器蛋白质组学

### 叶绿体蛋白组

叶绿体编码 ~100 个蛋白，但 >3000 个核编码蛋白输入叶绿体：
- **富集方法**：Percoll 梯度离心分离完整叶绿体 → 分级分离（膜/基质/类囊体）
- **关键验证**：用已知的叶绿体 markers（LHCB1 类囊体膜, RBCL 基质, TOC75 外膜）验证富集效率
- **挑战**：分离过程中线粒体和过氧化物酶体的交叉污染 > 10%

### 线粒体蛋白组

植物线粒体蛋白组比动物更复杂（植物线粒体有特殊的代谢通路）：
- **富集方法**：差速离心 + Percoll 梯度
- **特殊功能**：光呼吸、alternative oxidase (AOX)、uncoupling proteins
- **标记蛋白**：VDAC1 (外膜), COXII (内膜), MDH (基质)

### 何时做细胞器蛋白组学？

- 感兴趣蛋白的丰度在全细胞裂解物中太低 → 富集可提高灵敏度
- 需要进行空间分辨的功能分析（如蛋白运输、亚细胞调控）
- 每次细胞器分离通常需要 2-5g 鲜组织作为起始材料

## 常见陷阱

### Rubisco 的"诅咒"

Rubisco (Ribulose-1,5-bisphosphate carboxylase/oxygenase) 是地球上最丰富的蛋白，在植物绿色组织中占总蛋白 30-50%。

- 如果不 depletion，低丰度蛋白几乎检测不到
- Depletion 方法：抗 Rubisco 抗体免疫亲和柱（拟南芥/水稻）或 PEG 沉淀法
- 如果没有 depletion，至少报告 Rubisco 占总谱图的百分比

### RuBisCO 以外的高丰度蛋白

除了 Rubisco，还有一组蛋白因其高丰度而掩盖低丰度蛋白的检测：
- LHCB (Light-Harvesting Chlorophyll a/b-binding) 蛋白家族 → 占总蛋白 10-20%
- 叶绿体 ATP synthase 亚基 → 占叶绿体蛋白 5-10%
- 贮藏蛋白（种子中）：大豆 glycinin/conglycinin → 占种子蛋白 40-60%

**对策**：对种子组织 → 用 hexane 脱脂 + 高盐缓冲液去除贮藏蛋白

### 蛋白酶降解的隐蔽性

即使冷冻保存，植物组织中的蛋白酶在提取缓冲液加入时立即被激活：
- 检测蛋白降解的线索：SDS-PAGE 上出现大量 <20 kDa 的条带，或蛋白鉴定中覆盖度不均匀
- **预防胜于治疗**：新鲜组织 -> 立即液氮速冻 -> -80°C -> 提取时加入 1mM PMSF + cocktail

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
