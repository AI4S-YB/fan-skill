# 组装质量评估 (QC) -- 工具目录

## 概述

植物基因组组装质量评估需要从多个维度验证：完整性（BUSCO）、碱基准确性（Merqury）、连续性（contig/scaffold N50）、和结构正确性（Hi-C contact map）。推荐使用 BUSCO + Merqury 的双重评估策略。

## 推荐工具

### 1. BUSCO -- 基因组完整性评估

**描述**: BUSCO (Benchmarking Universal Single-Copy Orthologs) 使用一组单拷贝同源基因来评估组装的基因区域完整性。

**安装**:
```bash
conda install -c bioconda busco
```

**基本用法**:
```bash
busco -i assembly.fasta \
  -l embryophyta_odb10 \
  -o busco_output \
  -m genome \
  -c 16
```

**植物 Lineage 选择指南**:
| Lineage | 包含物种 | Busco 基因数 |
|---------|---------|-------------|
| `embryophyta_odb10` | 所有陆地植物 | 1614 |
| `viridiplantae_odb10` | 绿色植物 | 425 |
| `eudicots_odb10` | 双子叶植物 | 2326 |
| `liliopsida_odb10` | 单子叶植物 | 498 |
| `brassicales_odb10` | 十字花科 | 5158 |
| `poales_odb10` | 禾本科 | 4896 |
| `solanales_odb10` | 茄科 | 5950 |
| `fabales_odb10` | 豆科 | 5366 |

**查看可用 lineage**:
```bash
busco --list-datasets
```

**结果解读**:
```
# 输出的 short_summary.txt 示例：
C:90.2%[S:82.1%,D:8.1%],F:4.3%,M:5.5%,n:1614

C = Complete (完全匹配): 90.2%  —— 越高越好
S = Single copy (单拷贝): 82.1%  —— 理想情况应该高
D = Duplicated (重复): 8.1%      —— 对二倍体应 < 10%
F = Fragmented (碎片化): 4.3%     —— 应 < 10%
M = Missing (缺失): 5.5%          —— 应 < 10%
n = Total BUSCOs: 1614
```

**质量基准**:
| 指标 | 优秀 | 良好 | 勉强可接受 |
|------|------|------|-----------|
| Complete (C) | > 95% | > 90% | > 80% |
| Single (S) | > 70% | > 60% | > 50% |
| Duplicated (D) | < 5% | < 10% | < 15% |
| Fragmented (F) | < 5% | < 10% | < 15% |
| Missing (M) | < 5% | < 10% | < 15% |

### 2. Merqury -- k-mer 准确性评估

**描述**: Merqury 使用 k-mer 频谱比较组装置信度。它比较短读长 k-mer 与组装中的 k-mer，输出 QV（碱基准确性）、k-mer 完整性等关键指标。

**安装**:
```bash
conda install -c bioconda merqury
```

**完整流程**:
```bash
# 步骤 1: 使用 meryl 生成 k-mer 数据库
meryl k=21 count output genome.meryl \
  short_R1.fastq.gz short_R2.fastq.gz

# 步骤 2: 运行 Merqury 评估
merqury.sh genome.meryl assembly.fasta output_prefix

# 步骤 3: 查看结果
cat output_prefix.qv
```

**关键输出指标**:
- **QV (Quality Value)**: 组装碱基准确度 (QV = -10 * log10(错误率))
- **Completeness**: 读取组装的 k-mer 比例
- **Spectra-cn plot**: k-mer 频谱图，可视化倍性和重复

```bash
# 提取 QV 和 Completeness
cat output_prefix.completeness.stats
```

**QV 对标**:
| QV 值 | 碱基准确率 | 对应 |
|-------|----------|------|
| QV 20 | 99.0% | 较粗糙 |
| QV 30 | 99.9% | 良好 |
| QV 40 | 99.99% | 优秀 |
| QV 50 | 99.999% | 精准 (HiFi level) |
| QV 60 | 99.9999% | 参考基因组级别 |

### 3. 组装连续性统计

```bash
# 使用 assembly-stats
assembly-stats assembly.fasta

# 或使用 seqkit
seqkit stats assembly.fasta -a
```

**关键连续性指标**:
| 指标 | 说明 |
|------|------|
| N50 | contig 按长度排序，累加到总长 50% 时的 contig 长度 |
| N90 | 累加到 90% 时的 contig 长度 |
| L50 | 达到 N50 的最小 contig 数 |
| Total length | 组装总长度 |
| Longest contig | 最长 contig |
| GC content | 总 GC 含量 |
| # contigs | contig 总数 |

### 4. 组装图评估

```bash
# 统计组装图中的 dead ends 和复杂结构
# 从 Flye 的 GFA 中分析
grep "^S" assembly_graph.gfa | wc -l   # 总 contig 数
grep "^L" assembly_graph.gfa | wc -l   # 总 overlap 边数
```

---

## 植物组装质量评估特殊考虑

### 基因组大小一致性

比较组装总长度与预期基因组大小（来自 flow cytometry 或 k-mer 估计）：

```bash
# 计算组装总长度
TOTAL=$(grep -v "^>" assembly.fasta | wc -c)
EXPECTED=500000000  # 预期 500Mb

echo "Assembly: $TOTAL bp"
echo "Expected: $EXPECTED bp"
echo "Ratio: $(echo "scale=2; $TOTAL/$EXPECTED" | bc)"

# 理想比例: 0.85 - 1.15
# 多倍体可能略高 (haplotype duplication)
```

### 多倍体评估

对多倍体物种，BUSCO D (Duplicated) 比例预计偏高：
- 二倍体: D < 10%
- 同源四倍体: D ~ 40-60%
- 异源四倍体: D ~ 45-55%
- 六倍体: D ~ 65-80%

### 着丝粒和端粒完整性

植物着丝粒和端粒区域极难组装。可通过以下方式评估：
- 搜索端粒重复 (`TTTAGGG` 或 `CCCTAAA`)
- 着丝粒特异性 repeats (如 CentO 在玉米中, CentC/CentA)

```bash
# 搜索端粒重复
grep -c "TTTAGGGTTTAGGGTTTAGGG" assembly.fasta
```

---

## 质量报告清单

生成一份完整的组装质量报告应包括：

1. **GenomeScope 估计**: 基因组大小、杂合度、重复率
2. **BUSCO 结果**: 各类别百分比
3. **Merqury 结果**: QV, k-mer completeness
4. **连续性统计**: N50, L50, 总长
5. **Hi-C contact map** (如有): 染色体结构验证
6. **比对回贴率**: HiFi/ONT reads 回贴到组装的比例
