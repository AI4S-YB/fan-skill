# ATAC-seq 结果解读

## 一句话解释

ATAC-seq 通过 Tn5 转座酶标记开放染色质区域，鉴定全基因组调控元件（启动子、增强子等）的活性状态。

## 能回答什么

- 哪些基因组区域是开放的（可能参与基因调控）？
- 处理前后哪些区域的染色质状态发生变化？
- 哪些转录因子可能在这些区域结合（footprinting）？

## 不能回答什么

- TF 确实结合了吗？（需要 ChIP-seq 验证）
- 开放区域一定在调控基因吗？（相关性 ≠ 功能性）
- 染色质状态变化是原因还是结果？

## 典型输出

| 文件 | 含义 |
|------|------|
| `peaks.narrowPeak` | 开放染色质区域坐标 |
| `peak_annotation.csv` | Peak 在基因组中的位置分类 |
| `differential_accessibility.csv` | 差异可及性区域 |
| `footprint_scores.csv` | TF footprinting 结果 |

## 常见结果模式

### cpDNA reads 比例极高 (>50%)
绿色组织完全正常。比对后去除叶绿体/线粒体 reads 再分析。

### Peak 数量异常少
检查核分离效率、Tn5 酶切效率、MACS2 q 值阈值。
