# 基因调控网络结果解读

## 一句话解释

基因调控网络 (GRN) 分析：从表达数据中推断转录因子和靶基因之间的调控关系。

## 能回答什么

- 哪些转录因子可能调控了目标基因（GENIE3）？
- 哪些 TF 的调控子（regulon）在特定条件下活跃（SCENIC）？
- 表达数据的共表达模块有哪些（WGCNA）？
- 哪些基因是网络中的 hub？

## 不能回答什么

- TF 直接结合靶基因吗？（表达相关 ≠ 直接调控，需 ChIP/DAP-seq 验证）
- 调控方向是怎样的？（表达相关通常是无向的）
- Hub 基因一定重要吗？（高表达变异度的基因倾向于成为 hub）

## 典型输出

| 文件 | 含义 |
|------|------|
| `grn_edges.csv` | TF-靶基因调控关系及权重 |
| `regulons.csv` | SCENIC 调控子（含 motif 支持） |
| `hub_genes.csv` | 网络 hub 基因及其中心性指标 |
| `wgcna_modules.csv` | WGCNA 共表达模块 |

## 常见结果模式

### SCENIC 几乎无调控子通过 motif 过滤
非模式物种 motif 库不完整 → 用 GENIE3 代替（无 motif 过滤），结果标注为"调控假设"。

### Hub TF 的功能已知
如果文献支持该 TF 与所研究性状相关 → 最强的 in silico 验证。

### 多倍体中 cross-homeolog 信号
GRN 边中可能包含跨亚基因组的虚假共表达 → 使用 unique-mapping reads 定量。
