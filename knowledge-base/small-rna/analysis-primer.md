# 小 RNA 测序结果解读

## 一句话解释

小 RNA 测序：专门捕获和测序 18-30nt 的小 RNA 分子，鉴定 miRNA 和其他小 RNA 种类。

## 能回答什么

- 样品中有哪些已知的 miRNA？
- 是否存在新的 miRNA 种类？
- 哪些 miRNA 在处理条件下差异表达？
- miRNA 可能靶向哪些基因？

## 不能回答什么

- miRNA 真的降解了靶基因吗？（需要降解组测序或 5' RACE 验证）
- 预测的靶基因都是真实靶标吗？（植物中 miRNA 靶标预测精度较高但仍有假阳性）
- miRNA 表达变化一定功能重要吗？（miRNA 的生理效应可能很微妙）

## 典型输出

| 文件 | 含义 |
|------|------|
| `known_mirnas.csv` | 已知 miRNA 的表达量 |
| `novel_mirnas.csv` | 新预测的 miRNA 及其发夹结构 |
| `differential_mirnas.csv` | 差异表达 miRNA |
| `target_predictions.csv` | miRNA 靶基因预测 |

## 常见结果模式

### 非模式物种已知 miRNA 很少
miRBase 中非模式物种收录有限 → 依赖深度测序 + de novo 预测（miRDeep2）。

### miRNA 差异表达但靶基因无变化
植物中 miRNA 引导靶基因切割，但靶基因可能同时受转录调控补偿。
