# 基因组组装结果解读

## 一句话解释

基因组组装：将测序 reads 拼接成连续的基因组序列（contig/scaffold/染色体）。

## 能回答什么

- 基因组有多长？N50 是多少？
- 组装包含多少预期的保守基因（BUSCO）？
- 组装的碱基准确度如何（QV 值）？

## 不能回答什么

- 组装是"完全"的吗？（几乎所有基因组都有关键区域缺失）
- 所有 gap 都在正确的位置吗？（Hi-C 挂载可能有 misassembly）
- rRNA/着丝粒区域已完整组装吗？（通常没有）

## 典型输出

| 文件 | 含义 |
|------|------|
| `assembly.fasta` | 最终基因组序列 |
| `assembly_stats.txt` | N50/L50/总长/contig数 |
| `busco_summary.txt` | 完整性评估 |
| `merqury_qv.csv` | k-mer 评估的 QV 和完整性 |

## 常见结果模式

### 组装大小远大于预期
可能存在 haplotype duplication。用 purge_dups 去冗余。或杂合度过高。

### BUSCO Duplicated > 20% (二倍体)
高杂合度导致 redundant haplotypes。使用 purge_dups。

### 组装大小远小于预期
高度重复区域（着丝粒/rDNA）未能组装。植物基因组正常现象。
