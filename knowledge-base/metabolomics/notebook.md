# 植物代谢组学：分析笔记

这份笔记帮助你像一个有经验的植物代谢组学家一样思考和做决策。

## 在开始之前：理解你的分析平台

### LC-MS vs GC-MS vs NMR

| 特性 | LC-MS | GC-MS | NMR |
|------|-------|-------|-----|
| 检测范围 | 广泛（极性到非极性） | 挥发性和可衍生化合物 | 所有含 H 的代谢物 |
| 灵敏度 | 极高（fg-pg） | 高（pg-ng） | 低（ug-mg） |
| 鉴定能力 | MS/MS 可鉴定 | EI 谱库成熟 | 结构鉴定强 |
| 定量能力 | 相对定量 | 相对/绝对定量 | 绝对定量 |
| 样品制备 | 简单（溶剂提取） | 需要衍生化 | 简单/无破坏性 |
| 通量 | 高 | 中 | 低 |
| 重现性 | 中（受基质效应影响） | 高 | 很高 |
| 植物适用性 | 次生代谢物（黄酮/生物碱/萜类） | 初级代谢物（糖/酸/氨基酸） | 初级代谢物 |

植物代谢组学最常用 LC-MS，因为植物次生代谢物（黄酮、生物碱、萜类、苯丙素等）种类极多，LC-MS 的覆盖范围最广。

### 离子模式选择

- **正离子模式（ESI+）**：检测生物碱、氨基酸、核苷、甾醇
- **负离子模式（ESI-）**：检测有机酸、黄酮、酚酸、糖类
- **建议两者都采集**。数据分析时可分别处理或合并。

---

## 植物代谢组学的特色挑战

### 代谢物的极端多样性

植物界有 20 万+ 种次生代谢物，而且大多数未包含在公共数据库中。同一个化合物可能：
- 以多种糖苷形式存在（如槲皮素-3-O-葡萄糖苷、槲皮素-3-O-芸香糖苷）
- 存在多种同分异构体（如绿原酸异构体）
- 在不同物种中以不同修饰形式出现

这不是 "注释失败"，这是植物代谢组的真实面貌。

### 基质效应

植物提取物中富含色素（叶绿素）、多酚和脂质，这些共流出物会显著影响离子化效率。表现：
- 同一化合物在不同样本中信号强度不同，即使是相同浓度
- 保留时间漂移（柱污染）
- 有些化合物只在一部分样本中能检测到

对策：
- 使用同位素内标（如果可用）
- 使用 QC 样本（等量混合所有样本）进行信号校正
- 在序列开始前运行 5-10 次 QC 以平衡色谱柱

### 组织特异性

不同植物组织的代谢物谱差异极大：
- **叶**：光合色素、黄酮、苯丙素
- **根**：生物碱、萜类、酚酸
- **花**：花青素、挥发性萜类、类胡萝卜素
- **果实**：糖、有机酸、类胡萝卜素、花青素
- **种子**：油脂、储存蛋白、植酸

---

## 选择你的峰检测策略

### XCMS：LC-MS 的首选

XCMS 是目前 LC-MS 数据最通用的峰检测和对齐工具：

**适用场景**：
- 高分辨 LC-MS（Q-TOF、Orbitrap）
- 非靶向代谢组学
- 需要检测大量特征（>1000 个）

**参数提示**：
- **ppm**：质量误差容忍度。Orbitrap 用 5-10ppm，Q-TOF 用 15-30ppm
- **peakwidth**：峰宽范围（秒）。HPLC 用 c(5, 20)，UPLC 用 c(2, 10)
- **snthresh**：信噪比阈值。非靶向用 3-5，靶向用 10
- **prefilter**：至少连续 3 个扫描点，强度 > 100

### MZmine：GC-MS 或需要交互式操作

MZmine 有图形界面，更适合：
- GC-MS 数据（EI 谱库搜索内置）
- 需要可视化检查峰积分
- 小批量数据的交互式分析
- 分子网络构建（GNPS 兼容）

---

## 差异代谢物：不仅仅是 p < 0.05

### Limma：最稳健的统计框架

limma（用于代谢组学的 limma-trend 或 voom）：
- 基于经验贝叶斯调节方差估计
- 对小样本量（n=3-5）表现好
- 适合有多个处理因素的实验设计

**检查清单**：
1. 确认数据已 log2 转换
2. 确认缺失值已合理填充（不是简单地填 0）
3. 确认实验设计矩阵正确指定
4. 报告 log2 fold change + adjusted p-value + 置信区间

### 多重假设检验校正

代谢组数据通常有 1000-10000 个特征，需要多重检验校正：
- **BH (Benjamini-Hochberg)**：常用，q < 0.05 或 q < 0.1
- **Bonferroni**：太保守，会漏掉大部分真实差异代谢物
- 植物代谢组学中，q < 0.1 是可接受的标准（因为效应往往中等且生物重复数有限）

---

## 代谢物鉴定：从置信度到化学结构

### MSI 置信度等级

| 等级 | 定义 | 证据 |
|------|------|------|
| **Level 1** | 确认鉴定 | 与标准品：RT + MS1 + MS/MS 全部匹配 |
| **Level 2** | 推定鉴定 | MS/MS 谱库匹配（无标准品验证） |
| **Level 3** | 初步候选 | 分子式 + 数据库候选（如基于 MS1 的 m/z匹配） |
| **Level 4** | 未知化合物 | 已知 m/z 和 RT，但无鉴定 |
| **Level 5** | 完全未知 | 仅检测到特征，无任何鉴定信息 |

大多数植物非靶向代谢组学数据处于 Level 2-3。诚实地报告置信度。

### SIRIUS + CSI:FingerID 的价值

当你采集了 MS/MS 数据时：
1. SIRIUS 计算候选分子式（基于同位素模式 + MS/MS 碎片树）
2. CSI:FingerID 将 MS/MS 谱与结构数据库比对，预测化学结构
3. CANOPUS 预测化合物类别（即使无数据库匹配）

这在植物代谢组学中特别有用，因为很多植物次生代谢物不在标准数据库中。

---

## 通路映射：植物特异的挑战

### PlantCyc vs KEGG Plant

- **PlantCyc**：专门为植物代谢通路设计。包括 600+ 植物特异通路（如黄酮生物合成、芥子油苷生物合成、苯丙素通路）。覆盖物种 500+。
- **KEGG Plant**：更通用的植物通路。适合与模型物种（拟南芥、水稻）比较。

**两者结合使用**效果最好。PlantCyc 找到植物特异通路，KEGG 提供跨物种比较。

### 通路富集分析注意事项

通路富集的前提是：你的代谢物已被鉴定。如果 80% 的特征是未知的（这在植物代谢组学中很常见），通路富集结果不可靠。报告中需要明确说明。

---

## 常见陷阱

### 盲信数据库匹配

m/z 匹配不等于鉴定。一个 m/z 可能对应数十种同分异构体。如果没有 RT 或 MS/MS 验证，只能说是"特征"或"候选化合物"，不能说是"代谢物 X"。

### 忽略缺失值模式

缺失（未检测到）可能有生物学意义：
- 随机缺失（检测限以下）→ 用极小值或 left-censored 方法填充
- 非随机缺失（某组中完全不存在）→ 可能是有意义的生物学信号
- 不要盲目用均值或 0 填充所有缺失值

### 批次效应忽略

代谢组数据极容易受批次影响：
- 样品制备日期、LC-MS 运行日期（批次效应）往往与生物因素混杂
- 随机化样本运行顺序（不要按组别顺序进样）
- 如果批次无法避免，用 ComBat 或 QC-RLSC 校正
- 报告中必须说明批次处理策略

### 植物色素干扰

叶绿素和花青素在 LC-MS 中产生强烈的背景信号：
- 用适当的 SPE（固相萃取）去除
- 注意色谱柱负载量，不要超载
- 使用合适的提取溶剂（如甲醇:水，尽量避免强色素溶解）

---

## LC-MS vs GC-MS vs NMR：决策指南

上述表格提供了技术特性比较。以下是实际选择流程：

### 选择决策树

1. **你的目标是什么？**
   - 想尽可能多地看到不同的代谢物（非靶向）→ LC-MS
   - 关注初级代谢物（糖/有机酸/氨基酸）的绝对定量 → GC-MS 或 NMR
   - 关注挥发性有机物（花香、果实香气）→ GC-MS（顶空进样）
   - 需要确定未知化合物的精确结构 → NMR + LC-MS/MS 组合

2. **你的样本类型是什么？**
   - 叶片/绿色组织：LC-MS（次生代谢物丰富），注意叶绿素干扰
   - 根/根茎：LC-MS 或 GC-MS（生物碱、酚酸）
   - 花/果实：GC-MS（挥发性成分）+ LC-MS（花青素、类胡萝卜素）
   - 种子：GC-MS（脂肪酸甲酯化）+ LC-MS（脂质组学）
   - 木质部/韧皮部汁液：LC-MS（低浓度、高极性化合物）

3. **你的预算和通量需求？**
   - 大样本量（> 200 样本）且预算有限 → LC-MS
   - 小样本量但对结构鉴定要求极高 → NMR + LC-MS/MS
   - 需要定量准确性 → GC-MS（有标准品的情况）或 NMR

### 植物代谢组学平台组合策略

植物代谢组学中，单一平台无法覆盖所有代谢物类别。推荐的组合方案：

| 项目类型 | 推荐平台组合 | 理由 |
|---------|------------|------|
| 常规非靶向 | LC-MS (ESI+/-) | 覆盖范围最广 |
| 综合非靶向 | LC-MS (ESI+/-) + GC-MS | 初级+次级代谢物全覆盖 |
| 深度鉴定 | LC-MS/MS + NMR | 定量+结构鉴定 |
| 大规模群体 | LC-MS (ESI+/-) | 高通量，成本可控 |
| 挥发性代谢物 | GC-MS (HS-SPME) | 专用于挥发物 |
| 脂质组 | LC-MS/MS (lipidomics模式) | 专用于膜脂和储存脂 |

---

## 植物代谢物鉴定的独特挑战

### 未知质谱谱图问题

植物次生代谢物种类远超动物。截至 2026 年：
- 公共质谱库（GNPS/MoNA/METLIN）覆盖的植物特异性代谢物不足 30%
- 许多植物物种特有化合物完全没有参考谱图
- 黄酮、萜类、生物碱的同分异构体数量极多

应对未知谱图的策略：
1. **分子网络**（GNPS feature-based molecular networking）：将相似 MS/MS 谱图的特征聚类，通过已知化合物的邻近节点推断未知化合物类别
2. **CANOPUS 预测**：即使无数据库匹配，SIRIUS CANOPUS 可以基于碎片树预测化合物的化学类别（如"黄酮苷"或"三萜皂苷"），为功能解读提供线索
3. **植物化学文献交叉引用**：结合已有的植物化学文献（如该属的已知化合物列表），缩小候选范围

### 异构体分离

同一个 m/z 可能对应十几种或更多的异构体（如绿原酸有 3-、4-、5-咖啡酰奎宁酸三种异构体，槲皮素糖苷有数十种）。在植物代谢组学中：

- **色谱分离至关重要**：UPLC 比 HPLC 能提供更好的异构体分离。离子淌度（ion mobility, IM）是区分同分异构体的有力工具
- **MS/MS 碎片差异**：即使是同分异构体，碎片谱通常有细微差异。不要自动将相似碎片谱的峰合并
- **保留时间和碎片谱联合建模**：使用 `MS-DIAL` 或 `MZmine 3` 的离子淌度+碎片谱+RT 联合鉴定，可大幅提高异构体区分能力

### 糖苷化修饰的多样性

植物次生代谢物最常见的是糖苷化修饰，这使鉴定极其复杂：
- 同一苷元（aglycone）可能连接不同糖基（葡萄糖、鼠李糖、半乳糖、阿拉伯糖、木糖）
- 同一糖基可能连接在不同位置（如 3-O-葡萄糖苷 vs 7-O-葡萄糖苷）
- 同一化合物可能有二糖或三糖修饰

**实用建议**：在没有标准品的情况下，不要声称确定了糖基的连接位置。报告时可使用 Level 2 或 Level 3 置信度标注。

---

## 通路映射：PlantCyc vs KEGG Plant vs MetaCyc

### MetaCyc 的角色

MetaCyc 是所有 Pathway Tools 数据库（包括 PlantCyc）的"母数据库"：
- **MetaCyc**：涵盖所有生物界（细菌、真菌、植物、动物）的通路。通路经过实验验证，质量最高。包含 3000+ 条通路。
- **PlantCyc**：从 MetaCyc 中选取植物相关通路，并添加了 600+ 条植物特异性通路。覆盖 500+ 种植物物种。
- **KEGG Plant**：独立的通路数据库。与 MetaCyc/PlantCyc 是竞争关系，通路定义和组织方式不同。

### 三者如何结合使用

| 分析目标 | 首选 | 次选 | 原因 |
|---------|------|------|------|
| 初级代谢通路 | MetaCyc | PlantCyc | 初级代谢通路在 MetaCyc 中更全面且经过实验验证 |
| 植物特异次生代谢 | PlantCyc | KEGG Plant | PlantCyc 的黄酮、萜类、生物碱通路更完整 |
| 跨物种代谢比较 | KEGG Plant | MetaCyc | KEGG 提供了一套通用的代谢通路参考框架 |
| 酶功能注释 | MetaCyc | PlantCyc | MetaCyc 的酶学注释（EC）更准确 |
| 比较植物与其他生物 | MetaCyc | KEGG | 跨界的通路比较只能用 MetaCyc |

**实际工作流程**：先用 PlantCyc 找植物特异性通路，用 KEGG Plant 做跨物种比较（如水稻 vs 玉米 vs 拟南芥），如果关注核心代谢通路（TCA 循环、糖酵解等），用 MetaCyc 作为参考。

---

## 昼夜节律和昼夜效应对植物代谢组的影响

植物代谢组在 24 小时内剧烈波动。这不是噪声，是有生物学意义的信号。

### 日间 vs 夜间差异

- **光合产物（蔗糖、葡萄糖、果糖）**：白天升高，夜间消耗。日间/夜间比值可达 3-10 倍
- **淀粉**：白天积累，夜间降解
- **氨基酸**：日间合成，夜间分解代谢增强。部分氨基酸（如脯氨酸）夜间可能升高（胁迫响应）
- **有机酸（苹果酸、柠檬酸）**：CAM 植物的昼夜波动尤为剧烈
- **次生代谢物**：黄酮和花青素通常在光照下积累；某些萜类和生物碱的昼夜规律则物种特异

### 实验设计中的时间控制

时间效应在代谢组学实验中极为重要：

1. **固定取样时间**：所有样本在同一时间窗口采集（如上午 9:00-11:00）。不同日期采集需标注时间。
2. **避免跨时间段比较**：不要在早上采对照组、下午采处理组 — 昼夜差异会与处理效应混杂。
3. **昼夜节律实验的特殊处理**：如果要研究昼夜节律本身，按 Zeitgeber Time (ZT) 设计取样时间点（如 ZT0, ZT3, ZT6, ZT9, ZT12, ZT15, ZT18, ZT21）
4. **光照/黑暗条件记录**：记录取样时的光照条件、光强、温度。温室和田间条件差异大。

### CAM 植物的特殊注意事项

景天科、仙人掌科和某些兰科植物采用 CAM（景天酸代谢）：
- 夜间：气孔开放，固定 CO2 为苹果酸（苹果酸含量夜间飙升 10-100 倍）
- 白天：气孔关闭，苹果酸脱羧释放 CO2 供光合作用
- 代谢组学取样必须严格控制在特定时间点，否则数据不可解释

---

## 常见陷阱补充

### 电离抑制效应

植物粗提物中含有大量容易离子化的化合物，导致电离抑制：

- **多酚**：在 ESI- 模式下竞争电离，抑制有机酸和黄酮的离子化
- **脂质和磷脂**：在 ESI+ 模式下抑制生物碱和氨基酸的检测
- **高浓度糖类**：虽然糖类本身在 ESI 条件下不易离子化，但高浓度基质会影响喷雾稳定性

**缓解措施**：
- SPE 净化（反相或混合模式）
- 稀释样本（牺牲低丰度代谢物的检测）
- 使用小内径色谱柱（提高分离度，减少共流出）
- 碰撞截面（CCS）预测过滤（离子淌度-MS）

### 加合物形成

相同化合物在不同条件下可能形成多种加合物，导致数据混乱：

- **[M+H]+**：最常见，应该用于定量
- **[M+Na]+** 和 **[M+K]+**：钠钾加合物，在生物样本中很常见（植物富含钾）
- **[M+NH4]+**：铵加合物（使用甲酸铵流动相时）
- **[M-H2O+H]+**：源内碎裂，脱水产物
- **[2M+H]+**：二聚体，在高浓度时出现
- **[M+Cl]-** 和 **[M+HCOO]-**：负离子模式下的加合物

**如何处理**：
1. XCMS 的 CAMERA 或 MZmine 的加合物检测模块可以自动识别和分组加合物
2. 检查总离子流图（TIC）中 m/z 差值是否为 21.98 (Na) 或 37.96 (K)
3. 在报告中说明加合物处理策略
4. 对于复杂的加合物模式，使用 RAMClustR 进行分组

---

## 多元统计分析

### PCA 无监督降维

PCA 是代谢组学分析的第一步，用于评估数据整体质量：

```R
library(FactoMineR)
library(factoextra)

# 加载峰表数据（行=样本，列=代谢物特征）
peak_table <- read.table("peak_table.txt", header=TRUE, row.names=1)

# PCA分析（去除缺失值）
peak_table_na <- na.omit(peak_table)
pca_result <- PCA(peak_table_na, scale.unit=TRUE, ncp=5, graph=FALSE)

# 可视化
fviz_pca_ind(pca_result,
             col.ind = sample_groups,
             palette = "jco",
             addEllipses = TRUE,
             ellipse.type = "confidence")

# 关键检查点：
# 1. QC样本是否紧密聚集（评估技术重复性）
# 2. 同组样本是否聚类
# 3. 是否存在明显的批次效应或离群样本
```

**PCA结果解读**：
- PC1和PC2解释的方差比例应报告
- QC样本聚集程度反映数据质量（理想情况下QC应形成紧密簇）
- 样本按组别分离是预期的生物学信号
- 如果按批次分离，需要批次校正

### OPLS-DA 两组判别分析

OPLS-DA（正交偏最小二乘判别分析）用于两组比较：

```R
library(ropls)

# 准备数据
X <- as.matrix(peak_table)  # 代谢物特征矩阵
Y <- factor(sample_groups)   # 分组信息

# OPLS-DA分析
oplsda_model <- opls(X, Y,
                     predI = 1,        # 预测成分数
                     orthoI = NA,      # 正交成分数（NA表示自动选择）
                     crossvalI = 7,    # 7折交叉验证
                     log10L = TRUE,    # log10转换p值
                     permI = 200)      # 置换检验次数

# 提取VIP值
vip_scores <- getVipVn(oplsda_model)
vip_df <- data.frame(feature = colnames(X), VIP = vip_scores)

# 筛选VIP > 1的特征
significant_features <- vip_df[vip_df$VIP > 1, ]

# S-plot可视化（相关性 vs 协方差）
plot(oplsda_model, typeVc = "x-loading")

# 置换检验图（验证模型显著性）
plot(oplsda_model, typeVc = "permutation")
```

**OPLS-DA模型评估指标**：
| 指标 | 含义 | 可接受阈值 |
|------|------|----------|
| R2X | X矩阵解释率 | > 0.3 |
| R2Y | Y变量解释率 | > 0.7 |
| Q2 | 预测能力 | > 0.5（优秀 > 0.9） |
| Q2_perm | 置换检验Q2 | 原始Q2显著高于置换值 |

**注意事项**：
- OPLS-DA易过拟合，必须进行置换检验
- 样本量 < 10时，模型不稳定，谨慎解读
- 报告时必须包含R2X, R2Y, Q2值

### PLS-DA 多组判别分析

PLS-DA用于三组及以上的分类问题：

```R
library(mixOmics)

# 准备数据
X <- as.matrix(peak_table)
Y <- factor(sample_groups)  # 3+ 组别

# PLS-DA分析
plsda_model <- plsda(X, Y, ncomp = 2, scale = TRUE)

# 交叉验证确定最佳成分数
set.seed(123)
tune_plsda <- tune(X, Y,
                   ncomp = 5,
                   logratio = "none",
                   test = "Mfold",
                   folds = 5,
                   dist = "max.dist",
                   progressBar = FALSE)

# 使用最佳成分数重新建模
optimal_ncomp <- tune_plsda$choice.ncomp
plsda_final <- plsda(X, Y, ncomp = optimal_ncomp)

# VIP值提取
vip_scores <- vip(plsda_final)

# 样本得分图
plotIndiv(plsda_final, ind.names = FALSE,
          group = Y, legend = TRUE,
          ellipse = TRUE, style = "ggplot2")

# 变量重要性图
plotVar(plsda_final, cutoff = 0.5)
```

### VIP筛选与单变量统计结合

```R
# 结合VIP和p值进行差异代谢物筛选
library(limma)

# 1. VIP筛选
vip_threshold <- 1.0
vip_selected <- colnames(X)[vip_scores[,1] > vip_threshold]

# 2. 单变量统计（limma）
design <- model.matrix(~ 0 + Y)
fit <- lmFit(t(X), design)
contrast <- makeContrasts(group2 - group1, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
pvals <- topTable(fit2, number = Inf)$P.Value
names(pvals) <- rownames(topTable(fit2, number = Inf))

# 3. 结合筛选
diff_metabolites <- intersect(
  vip_selected,
  names(pvals[pvals < 0.05])
)

# 报告筛选条件
cat("VIP阈值:", vip_threshold, "\n")
cat("p值阈值: 0.05\n")
cat("差异代谢物数量:", length(diff_metabolites), "\n")
```

---

## 机器学习分析

### 随机森林特征选择

随机森林适用于中等及以上规模数据（n >= 20）：

```R
library(randomForest)
library(caret)

# 准备数据
data <- data.frame(X, group = Y)

# 划分训练集和测试集（70/30）
set.seed(123)
train_index <- createDataPartition(Y, p = 0.7, list = FALSE)
train_data <- data[train_index, ]
test_data <- data[-train_index, ]

# 训练随机森林模型
set.seed(456)
rf_model <- randomForest(group ~ ., data = train_data,
                         ntree = 1000,
                         mtry = sqrt(ncol(X)) |> floor(),
                         importance = TRUE,
                         proximity = TRUE)

# 模型评估
rf_pred <- predict(rf_model, test_data)
confusion_matrix <- confusionMatrix(rf_pred, test_data$group)
print(confusion_matrix)

# OOB error
cat("OOB error rate:", rf_model$err.rate[nrow(rf_model$err.rate), "OOB"], "\n")

# 变量重要性
importance_df <- data.frame(
  feature = rownames(importance(rf_model)),
  MeanDecreaseAccuracy = importance(rf_model)[, "MeanDecreaseAccuracy"],
  MeanDecreaseGini = importance(rf_model)[, "MeanDecreaseGini"]
)

# 排序并展示Top 20
importance_df <- importance_df[order(-importance_df$MeanDecreaseGini), ]
head(importance_df, 20)

# 可视化
varImpPlot(rf_model, n.var = 30, main = "Variable Importance")
```

**随机森林参数调优**：

```R
# 使用caret进行mtry参数调优
tune_grid <- expand.grid(mtry = c(2, 5, 10, 15, 20))

ctrl <- trainControl(method = "cv",
                     number = 5,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

rf_tune <- train(group ~ ., data = train_data,
                 method = "rf",
                 tuneGrid = tune_grid,
                 trControl = ctrl,
                 metric = "ROC",
                 ntree = 1000)

# 最佳mtry值
print(rf_tune$bestTune)
```

### t-SNE 和 UMAP 降维可视化

```R
library(Rtsne)
library(umap)

# t-SNE分析
set.seed(789)
tsne_result <- Rtsne(X, dims = 2, perplexity = 30,
                     verbose = FALSE, max_iter = 1000)

# t-SNE可视化
tsne_df <- data.frame(
  tSNE1 = tsne_result$Y[, 1],
  tSNE2 = tsne_result$Y[, 2],
  group = Y
)

ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = group)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(title = "t-SNE Plot")

# UMAP分析
umap_result <- umap(X, n_neighbors = 15, min_dist = 0.1)

# UMAP可视化
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  group = Y
)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = group)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(title = "UMAP Plot")
```

**参数选择指南**：
| 参数 | 推荐值 | 说明 |
|------|--------|------|
| t-SNE perplexity | sqrt(n) 到 n/3 | 小样本用较小值 |
| UMAP n_neighbors | 5-50 | 影响局部/全局结构平衡 |
| UMAP min_dist | 0.001-0.5 | 影响点聚集程度 |

### ROC曲线与生物标志物验证

```R
library(pROC)

# 单代谢物ROC分析
roc_single <- roc(Y, X[, "feature_name"], levels = c("control", "case"))

# 计算最佳截断值
best_threshold <- coords(roc_single, "best", best.method = "youden",
                         ret = c("threshold", "sensitivity", "specificity"))

# AUC值及95%置信区间
auc_value <- auc(roc_single)
ci_value <- ci.auc(roc_single)

# ROC曲线绘制
ggroc(roc_single) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "gray") +
  annotate("text", x = 0.3, y = 0.2,
           label = paste("AUC =", round(auc_value, 3),
                        "\n95% CI:", round(ci_value[1], 3), "-", round(ci_value[3], 3))) +
  theme_minimal() +
  labs(title = "ROC Curve for Biomarker Evaluation")

# 多代谢物联合ROC（逻辑回归模型）
library(caret)

# 训练逻辑回归模型
set.seed(123)
train_control <- trainControl(method = "cv", number = 5,
                              classProbs = TRUE,
                              summaryFunction = twoClassSummary)

lr_model <- train(group ~ feature1 + feature2 + feature3,
                  data = train_data,
                  method = "glm",
                  family = "binomial",
                  trControl = train_control,
                  metric = "ROC")

print(lr_model)
```

---

## 可视化模块

### 火山图

```R
library(EnhancedVolcano)

# 准备差异分析结果
diff_results <- data.frame(
  feature = rownames(topTable(fit2, number = Inf)),
  log2FC = topTable(fit2, number = Inf)$logFC,
  pvalue = topTable(fit2, number = Inf)$P.Value,
  padj = topTable(fit2, number = Inf)$adj.P.Val
)

# 火山图绘制
EnhancedVolcano(diff_results,
                lab = diff_results$feature,
                x = "log2FC",
                y = "pvalue",
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~italic(P)),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 3.0,
                title = "Differential Metabolites",
                subtitle = "Volcano Plot",
                legendPosition = "right",
                selectLab = head(diff_results$feature[diff_results$pvalue < 0.001], 10))
```

### 热图聚类

```R
library(ComplexHeatmap)
library(circlize)

# 准备差异代谢物数据
diff_data <- X[diff_metabolites, ]

# Z-score标准化（行标准化）
diff_data_scaled <- t(scale(t(diff_data)))

# 颜色映射
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# 行注释（代谢物类别）
row_anno <- rowAnnotation(
  Class = metabolite_class,
  col = list(Class = c("Flavonoid" = "#E41A1C",
                       "Alkaloid" = "#377EB8",
                       "Terpenoid" = "#4DAF4A",
                       "Other" = "#984EA3"))
)

# 列注释（样本分组）
col_anno <- columnAnnotation(
  Group = sample_groups,
  col = list(Group = c("Control" = "#66C2A5",
                       "Treatment" = "#FC8D62"))
)

# 绘制热图
Heatmap(diff_data_scaled,
        name = "Z-score",
        col = col_fun,
        top_annotation = col_anno,
        left_annotation = row_anno,
        show_row_names = FALSE,
        show_column_names = FALSE,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "complete",
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 10))
```

**差异代谢物过多时的处理**：
- 按 p 值或 VIP 值排序，仅展示 Top 50-100
- 按代谢物类别分组展示
- 使用折线图或点图替代

---

## 差异代谢物分析完整流程

```R
# 完整分析流程示例
library(XCMS)
library(limma)
library(ropls)
library(pROC)
library(ggplot2)
library(ComplexHeatmap)

# Step 1: 数据预处理
# XCMS峰检测 → 峰对齐 → 缺失值填充

# Step 2: 数据过滤与标准化
# 去除低质量峰（在<50%样本中检测到）
# PQN或VSN标准化

# Step 3: PCA评估数据质量
# 检查QC聚集程度、批次效应

# Step 4: 差异分析
# 4a. 单变量统计（limma）
# 4b. 多变量分析（OPLS-DA/PLS-DA VIP）

# Step 5: 差异代谢物筛选
# 结合VIP、p值、fold change

# Step 6: 功能注释与通路分析
# SIRIUS+CSI:FingerID鉴定
# PlantCyc/KEGG通路富集

# Step 7: 可视化
# 火山图、热图、ROC曲线

# Step 8: 结果报告
# 差异代谢物表、统计指标、置信度标注
```
