# 通用排错指南

## 环境问题

### R 包缺失
```
Error in library(GAPIT3): there is no package called 'GAPIT3'
```
→ `devtools::install_github("jiabowang/GAPIT3")` 或使用 singularity 容器

### Conda 环境冲突
```
LibMamba UnsatisfiableError
```
→ 使用 `pixi` 代替 conda，或创建独立环境

### Singularity 镜像不可用
→ 检查 `test_config.yaml` 中 `container_image` 路径
→ 如镜像在 `/share/` 下，确认 `--bind` 参数

## 数据问题

### 染色体名不匹配
```
Warning: Chromosome name 'Chr01' not found in reference
```
→ 植物的染色体命名极不统一。常见模式：
- 水稻: Chr01..Chr12 或 chr01..chr12 或 1..12
- 玉米: chr1..chr10 或 1..10
- 小麦: chr1A..chr7D 或 1A..7D
→ 用 `inspect_data.sh` 的染色体名检测功能

### 编码问题 (中文 Windows 数据)
```
Error: line 1 contains embedded nulls
```
→ `iconv -f GBK -t UTF-8 input.csv > input_utf8.csv`

### 表型文件格式
常见陷阱：
- 表型值用中文标记（"高"/"低"而非数值）
- 缺失值用 "NA"/"."/"-" 不统一
- 品种名含特殊字符（空格、/、括号）
→ `inspect_data.sh` 会自动检测并报告

## 分析问题

### λ 值偏高 (GWAS)
- 加 PCA 协变量 (PC1-PC5)
- 加 K 矩阵 (亲缘关系)
- 检查是否存在极端样本 (>6 SD on any PC)
- 自交作物 λ>1.2 不一定是问题

### 零个显著位点 (GWAS)
可能原因（按概率排序）：
1. 标记太少、统计功效不足 — 这是最常见原因
2. 样本量不够
3. 该性状确实无主效 QTL（微效多基因控制）
4. 表型测量误差大
→ 诚实报告，不要过度解读

### 群体结构噪音 (Population Genetics)
- PC1/PC2 解释方差异常高 → 正常的强群体分化
- PC1 解释方差过低 (<3%) → 可能是 LD 未被充分过滤
- Admixture 交叉验证错误不收敛 → 增加 max iterations 或调整 K 范围

### 表型分析

#### 遗传力估计异常低
→ 检查是否有环境/重复效应未建模
→ 多环境数据优先用混合模型估计广义遗传力
→ 单环境单重复数据无法估计遗传力

#### BLUP 值不符合预期
→ 检查是否有缺失值导致估计不稳定
→ 检查模型是否收敛

### 变异检测

#### 比对率异常低 (<70%)
→ 检查参考基因组是否匹配物种
→ 近缘种参考 → 比对率可接受范围 60-80%
→ 检查 FASTQ 质量 (fastqc)
→ 多倍体物种比对率天然较低

#### 大量杂合位点在自交物种中出现
→ 自交物种 (水稻/大豆) 杂合率应 <5%
→ 高杂合率 → 可能混样或污染
→ 异交物种 (玉米) 杂合率高是正常的

### 基因型填充

#### 填充准确率低 (R² < 0.3)
→ 参考面板与目标群体遗传距离太远
→ 标记密度太低 → 先用高密度面板做预填充
→ 某些基因组区域 (着丝粒) 天然难以填充

### QTL 定位

#### 遗传图谱异常 (标记顺序与物理位置不一致)
→ 检查是否有分离畸变标记 (显著偏离 1:1/1:2:1)
→ 检查标记是否有高缺失率
→ 某些区域可能存在倒位/易位

#### LOD 曲线无显著峰
→ 增大 permutation 次数以更准确估计阈值
→ 该性状可能由多个微效 QTL 控制 → 无单一峰
→ 群体太小 (<100) → 只能检测大效应 QTL

### 时间序列表达

#### 聚类结果找不到清晰的模式
→ 检查基因过滤是否太松（太多噪音基因）
→ 增大 fuzzifier m 值（更软的聚类）
→ 检查时间点是否足够（<5 个时间点难以识别趋势）

### 基因调控网络

#### GENIE3 运行极慢
→ 基因太多 → 先按方差过滤到 top 2000-5000
→ 考虑用随机森林参数 ntree=500 加速
#### SCENIC 未检测到显著 regulon
→ 检查物种是否有足够的 TF 注释
→ 非模式植物可能需要从拟南芥同源推断 TF

### 多组学整合

#### MOFA 因子无法解释
→ 检查各数据视图是否经过适当的归一化
→ 尝试调整因子数（先从 5-10 开始）
#### DIABLO 分类性能差
→ 减少特征数（用 loading 筛选 top features）
→ 检查类别是否平衡

### eQTL 分析

#### 大量 cis-eQTL 但没有 trans-eQTL
→ trans-eQTL 需要的统计功效远高于 cis
→ 样本量 <100 时几乎检测不到 trans-eQTL
#### GWAS 和 eQTL 共定位失败
→ 两种分析使用不同的群体 → 遗传结构不同
→ 检查是否使用了相同的 LD 参考

### 环境组学

#### G×E 方差估计为零或极小
→ 检查环境差异是否真的存在
→ 检查基因型×环境交互项是否可估计
#### TPE 定义不清晰
→ 环境数据太少 → 至少需要 5 个环境的主要气候变量
→ 可以考虑用历史气候数据补充

### 小 RNA 分析 (small RNA)

#### miRNA 预测 reads 数过低
→ 检查 sRNA 文库是否经过 size selection (18-30 nt)
→ 测序深度不足 → 建议 small RNA-seq 深度 > 10M reads
→ 降解组数据可辅助验证 miRNA 切割靶点

#### miRNA 靶基因预测假阳性过多
→ 调低 psRNATarget expectation 阈值 (≤ 3.0)
→ 结合降解组数据 (PARE/degradome) 验证切割位点
→ 要求靶位点位于保守区域

### ChIP-seq

#### 峰值数过少或无显著峰
→ 检查抗体质量 — ChIP-grade 抗体是前提
→ 检查 Input control 是否有足够深度
→ 考虑放宽 MACS2 q-value 阈值（但不宜超过 0.1）
→ 宽峰标记 (H3K27me3 等) 需使用 MACS2 --broad 模式

#### FRiP (Fraction of Reads in Peaks) 异常低 (<1%)
→ 信噪比不足 — 可能是抗体效率低或洗涤条件过严
→ 检查是否有大量 PCR duplicate → 使用 Picard MarkDuplicates
→ 考虑更换抗体或优化 IP 条件

### ATAC-seq

#### TSS 富集分数低 (<5)
→ 文库质量差 — 检查插入片段分布 (应有核小体 ladder)
→ 核制备过程中染色质受损 → 优化细胞裂解条件
→ 测序深度不足 → 建议 ≥ 50M reads per sample

#### 峰值集中在非功能区 (intergenic/重复区)
→ 可能是线粒体 DNA 污染 → 比对前过滤 MT 染色体
→ 死细胞比例过高 → 使用新鲜样本或流式分选活细胞

### 甲基化分析 (methylation)

#### 全基因组 CpG 覆盖度不足
→ WGBS 深度不够 → 建议 ≥ 15x 覆盖度
→ RRBS 仅覆盖 CpG 岛区域,不适用于全基因组分析
→ 植物中 CHG/CHH 甲基化也需要关注 (非 CpG 背景)

#### DMR 检测结果过少
→ 检查组间甲基化差异是否确实小 (生物变异)
→ 提高生物学重复数 (建议 ≥ 3 per group)
→ 放宽 DMR 的最小 CpG 位点数要求

### 基因组组装 (genome assembly)

#### BUSCO 完整度 < 80%
→ 杂合度过高导致组装碎片化 → 使用 HiFi 长读长或 trio-binning
→ 多倍体物种天然困难 → 考虑分型组装 (haplotype-resolved)
→ 某些谱系特异性 BUSCO 基因天然缺失是正常的

#### Contig N50 低于技术预期
→ HiFi 组装 N50 应 > 5 Mb → 检查 DNA 提取质量 (HMW DNA)
→ ONT 组装碎片化 → 检查 read N50, 考虑使用 R10.4.1 化学
→ 短读长组装低 N50 → 这是短读长技术的固有局限, 报告 scaffold N50 代替

### 基因组注释 (genome annotation)

#### BUSCO 蛋白完整性 < 80%
→ 检查基因预测模型是否适合该物种 (使用近缘种训练的 Augustus 模型)
→ RNA-seq 证据不足 → 增加组织/发育阶段的 RNA-seq 数据
→ 使用 MAKER 或 BRAKER3 整合多种证据
→ 转座子注释干扰基因预测 → 先做硬遮蔽 (hard masking) 再注释

### 扩增子分析 (amplicon)

#### 稀疏曲线未达到平台期
→ 测序深度不足 → 增加每个样本的 reads 数
→ 群落多样性极高 (如土壤微生物) → 深度 5000 可能不够
→ 检查是否有大量 singleton ASV → 可能是测序错误而非真实多样性

#### ASV 丰度过低 (大量 <0.01% 的 ASV)
→ 考虑过滤低丰度 ASV 以减少噪音
→ 可能是交叉污染 (index hopping) → 使用双端 index
→ 检查阴性对照中是否也出现这些 ASV

### 代谢组学 (metabolomics)

#### 峰检测大量缺失值 (>50%)
→ 信号强度低于检测限 → 优化样品前处理和离子化参数
→ 检查色谱柱性能和流动相新鲜度
→ 考虑使用缺失值填充 (KNN, random forest 等方法)
→ 大量缺失提示代谢物浓度差异过大, 可能具有生物学意义

#### S/N 比低导致峰检测不可靠
→ 检查质谱仪校准状态
→ 样品量不足或降解 → 优化提取方案
→ 提高进样浓度或调整质谱参数 (dwell time 等)

### 泛基因组 (pan-genome)

#### 核心基因组占比异常 (过高或过低)
→ 基因组数 < 5 → 核心基因组占比被高估
→ 物种分化过大 → 同源基因聚类参数需调整 (identity/cov cutoff)
→ 检查是否存在 contamination 导致的错误基因

#### 可变基因组的聚类不收敛
→ 同源性阈值太松 → 提高 protein identity 阈值 (>50%)
→ 部分基因组成员过多 → 可能是转座子或重复序列污染
→ 使用 OrthoFinder 或 PanX 的默认参数作为起点

### 宏基因组 (metagenomics)

#### MAG 完整度低 / 污染度高
→ 组装碎片化 → 增加测序深度或使用 HiFi 长读长
→ 群落复杂度高 → CheckM 可能低估完整度
→ 污染 > 10% → 使用 MAGpurify 或 GUNC 进一步净化
→ 考虑共组装 (co-assembly) vs 单样本组装策略

#### Binning 结果中 MAG 数量过少
→ 增加测序深度 (土壤宏基因组建议 > 10 Gbp)
→ 低丰度物种无法被有效分箱 → 这是方法学局限
→ 尝试多种 binning 工具取并集 (MetaBAT2, CONCOCT, MaxBin2 → DAS Tool)

### 蛋白质组学 (proteomics)

#### 鉴定蛋白数量远低于预期
→ 数据库不完整 → 使用六框翻译基因组或转录组作为搜索空间
→ 酶切效率低 → 检查胰蛋白酶活性和酶切条件
→ FDR 阈值过严 → 蛋白水平 FDR 1% 是标准, 不宜放宽
→ 检查 LC-MS/MS 系统灵敏度

#### 定量结果缺失值过多
→ 使用 DIA/SWATH 而非 DDA 以获得更完整的数据矩阵
→ 考虑 match-between-runs (MaxQuant MBR 功能)
→ 低丰度蛋白的随机缺失是 DDA 的固有局限

### 杂种预测 (hybrid prediction)

#### GCA 预测准确度低 (r < 0.5)
→ 训练群体与预测群体遗传关系远 → 增加训练群体遗传代表性
→ 标记密度不足 → GBLUP 至少 1000 标记, BayesB 至少 5000
→ 训练群体规模太小 (< 100) → 遗传参数估计不稳定
→ 检查表型数据质量 (遗传力是否过低)

#### 预测结果在不同环境中不一致
→ 存在强 G×E 互作 → 按环境分别建模或使用 G×E 模型
→ 检查各环境表型相关性 — 低相关预示预测困难
→ 考虑使用 multi-environment GS 模型 (如 BGGE, EnvRtype)

### CRISPR

#### 预测到的脱靶位点过多
→ 种子区域错配数阈值太松 (建议 ≤ 3)
→ gRNA 序列特异性差 → 重新设计 gRNA 避开重复/低复杂度区域
→ 植物基因组 (尤其是多倍体) 同源序列多 → 这是植物 CRISPR 的天然挑战

#### 编辑效率低或无编辑
→ gRNA 二级结构不良 → 检查 GC 含量 (建议 40-60%)
→ 染色质可及性差 → 靶向开放染色质区域 (结合 ATAC-seq/DNase-seq 数据)
→ 启动子选择影响 Cas 表达时空调控 → 优化启动子 (如使用植物泛素启动子)
→ 多倍体物种所有同源拷贝都需编辑 → sgRNA 设计时考虑同源位点
