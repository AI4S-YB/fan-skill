# 植物 CRISPR/Cas9 设计：分析笔记

## 在开始之前

植物基因编辑的核心挑战来自三个方面：基因组复杂性（多倍体、重复序列）、转化效率（不同物种差异巨大）、以及同源基因干扰。先问自己：

### 你想达成什么编辑目标？

- **基因敲除 (Knockout/KO)**：产生移码突变 → 设计 sgRNA 靶向早期外显子，PAM 后产生 DSB，NHEJ 修复引入 indel
- **基因敲入 (Knock-in/KI)**：插入外源序列 → 需要 HDR 供体模板 + 同源臂设计
- **碱基编辑 (Base Editing)**：C-to-T 或 A-to-G 单碱基替换 → 需要目标碱基位于编辑窗口内（通常为 protospacer position 4-8）
- **Prime Editing**：精确的短片段插入/删除/替换 → 需要 pegRNA（prime editing guide RNA）+ nicking sgRNA
- **多重编辑 (Multiplex)**：同时编辑多个基因 → 多个 sgRNA 共表达 + 可能需要 Cas12a（更简单的 crRNA array）

### 你的物种有什么资源？

- **主要作物**（水稻、玉米、拟南芥、大豆、番茄、小麦）→ CRISPOR 有预建参考基因组
- **非模式物种**→ CRISPR-P 允许上传自定义基因组，或使用 CHOPCHOP 命令行
- **多倍体**（小麦、棉花、油菜）→ 需要考虑同源基因共享 sgRNA 靶序列（同时编辑所有亚基因组拷贝）

### 你用什么 Cas 变体？

| Cas 变体 | PAM | 特点 | 植物应用 |
|----------|-----|------|---------|
| SpCas9 | NGG | 最常用 | 几乎所有植物物种 |
| Cas9-NG | NG | 放宽的 PAM | GC 含量高的植物基因组 |
| Cas12a (Cpf1) | TTTV | T-rich PAM, 交错切割 | AT-rich 区域 |
| Cas9-VQR | NGAN | 替代 PAM | 特定靶序列 |
| xCas9 | NG/GAA/GAT | 最宽 PAM | 靶序列 PAM 选择受限时 |

## 选择你的 sgRNA 设计工具

### CRISPOR — 主要作物首选

预建了主要植物物种的参考基因组，Web 界面易用。

**适合**：水稻、玉米、拟南芥、大豆、番茄、小麦
**优势**：综合评分（MIT specificity score + Doench efficiency score + 脱靶预测）

### CRISPR-P — 非模式植物

支持上传自定义基因组 FASTA。

**适合**：没有在 CRISPOR 中列出的植物物种
**优势**：支持自定义基因组、多倍体分析、基因家族考虑

### CHOPCHOP — 本地运行

命令行版可在本地/HPC 上运行。

**适合**：高通量设计（数百个基因）、自动化管道、自有基因组

## sgRNA 设计原则

### 最优 sgRNA 的特征

1. **位置**：CDS 前 1/3 处（产生早期的移码终止）
2. **特异性**：全基因组中独特的靶序列（避免同源基因同时被编辑）
3. **效率**：GC%=40-60%, 避免 4+ T 连续（Pol III 终止信号）
4. **PAM**：NGG（SpCas9）优先，检查是否有合适的 PAM 在基因组上下文中
5. **二级结构**：sgRNA 不应形成强发夹结构

### 植物特异性考虑

- **靶向所有同源拷贝** vs **靶向特定亚基因组拷贝**：需要在 sgRNA 选择时明确
- 如果目标是在多倍体中敲除所有亚基因组拷贝，选择保守区域的 sgRNA
- 如果要研究亚基因组特异性功能，在外显子序列不同的区域设计 sgRNA
- 小麦中 3 个亚基因组（A、B、D）→ 通常需要 1 个保守 sgRNA 或 3 个特异性 sgRNA

## 脱靶分析

### Cas-OFFinder

全基因组脱靶搜索工具，Cas-OFFinder 可以搜索最多 3 个错配 + 1 个 DNA/RNA bulge 的脱靶位点。

**注意**：
- 植物基因组中重复序列多 → 某些 sgRNA 可能有数百个脱靶位点 → 应该排除
- 启动子区域的脱靶可能影响其他基因的表达
- 基因间区域的脱靶通常无害

### 特异性评分

- **CFD (Cutting Frequency Determination)** score：评估脱靶位点被切割的可能性
- **MIT specificity score**：CRISPOR 中使用的综合评分
- 一般标准：CFD score > 80 为高特异性

## 编辑效率预测

### DeepSpCas9

基于深度学习预测 SpCas9 在特定 sgRNA 靶序列上的编辑效率。

输入：30bp 序列（4bp upstream + 20bp protospacer + 3bp PAM + 3bp downstream）
输出：预测的 indel 频率

### 影响编辑效率的因素

1. **染色质可及性**：开放染色质区域的 sgRNA 效率更高（植物的 DNase-seq 或 ATAC-seq 数据可以辅助判断）
2. **GC 含量**：40-60% 最优
3. **PAM 邻近序列**：位置 17-20（靠近 PAM 的 4bp）对效率影响最大
4. **sgRNA 的二级结构**：折叠后正确的 scaffold 结构对于 Cas9 loading 至关重要

## HDR 供体模板设计

### Knock-in 策略

植物中 HDR（Homology-Directed Repair）效率远低于 NHEJ。策略包括：

1. **同时抑制 NHEJ**：使用 NHEJ 抑制剂（如 Scr7）或 Ku70/Ku80 突变体背景
2. **增加供体模板量**：通过 geminivirus replicon 系统在植物体内大量扩增供体
3. **同源臂长度**：通常 400-800bp
4. **共表达 Cas9 + sgRNA + 供体**：使用双 T-DNA 或 polycistronic 系统

### 供体模板设计要点

- 左同源臂 + 插入序列 + 右同源臂（全长 800-2000bp）
- 供体 PAM 应该被突变掉（silent mutation）以防止 Cas9 切割供体
- 插入序列 < 2kb（HDR 效率随插入长度指数下降）
- 供体 DNA 可以通过 particle bombardment 或 Agrobacterium 递送

## 基因分型验证

### 编辑检测方法

1. **T7E1/Surveyor assay**：酶切错配检测，快速但半定量
2. **PCR + Sanger 测序**：金标准，可以确定确切的 indel 序列
3. **ICE (Inference of CRISPR Edits)** 分析：Sanger 测序色谱图分解，估计编辑效率和基因型比例
4. **扩增子 NGS**：多重编辑时使用，高灵敏度检测低频率编辑事件

### 植物特异性验证

- T0 代通常是嵌合体（多种基因型共存于同一植株）
- 需要 T1 或 T2 代才能获得纯合编辑植株
- 分离"转基因"（Cas9/sgRNA）和"编辑"（突变）：在 T2 代可以筛选 transgene-free 编辑事件
- 体细胞编辑事件不遗传到下一代 → T0 的基因分型不能完全预测 T1 的结果

## 常见陷阱

### 多倍体的 sgRNA 特异性

- 在一个亚基因组中"独特"的靶序列可能与另一个亚基因组中的序列只有 1-2 个错配 → 仍可能被切割
- **检查**：将 sgRNA 比对到所有亚基因组，查看错配分布
- **策略**：如果不想编辑所有拷贝，选择 2+ 错配位于 PAM 近端（seed region）的 sgRNA

### 基因冗余

- 同一基因家族的成员可能导致"敲除表型缺失"（另一个家族成员补偿了功能）
- **策略**：先检查基因家族（OrthoFinder），如果是大家族，考虑多重敲除多个成员

### 转化效率差异

不同植物物种的转化效率天差地别：

| 物种 | 转化方法 | 典型效率 |
|------|---------|---------|
| 水稻 | Agrobacterium | 50-90% |
| 玉米 | Agrobacterium | 5-20% |
| 小麦 | Particle bombardment | 1-10% |
| 大豆 | Agrobacterium (cotyledonary node) | 1-5% |
| 棉花 | Agrobacterium | 5-15% |
| 番茄 | Agrobacterium | 10-50% |

如果转化效率 <5%，可能需要在 sgRNA 设计阶段选择 3-5 条 sgRNA 以增加成功概率。
