# 宏基因组结果解读

## 一句话解释

宏基因组测序：对环境样品中所有微生物的 DNA 进行测序，同时获得分类和功能信息。

## 能回答什么

- 微生物群落中有什么物种？
- 微生物有哪些功能基因（如 N 循环、抗生素抗性）？
- 能否重建微生物基因组 (MAGs)？
- 哪些代谢通路在群落中富集？

## 不能回答什么

- 检测到的基因在采样时正在表达吗？（DNA ≠ RNA，需要宏转录组）
- 所有 MAGs 都是完整的吗？（CheckM 评估：Complete<70% 是部分基因组）
- 低丰度微生物能看到吗？（测序深度的物理限制）

## 典型输出

| 文件 | 含义 |
|------|------|
| `taxonomic_profile.csv` | 各样本中各类群的相对丰度 |
| `functional_annotation.csv` | 功能基因/KEGG 模块丰度 |
| `mags_quality.csv` | MAGs 的完整性和污染度 (CheckM) |
| `metabolic_pathways.csv` | 预测的代谢通路 |

## 常见结果模式

### 植物 DNA >90% → MAGs 数量极少
先去宿主 reads (bowtie2 比对到植物参考基因组，取 unaligned reads)。

### MAGs 质量不高 (Complete<50%)
增加测序深度；用混合组装 (短读长+长读长)；降低 binning 的 contig 长度阈值。
