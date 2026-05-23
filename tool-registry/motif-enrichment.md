# Motif 富集分析 — 工具目录

## 概述

Motif 富集分析用于从 ChIP-seq peaks 中识别富集的 DNA 序列模式，帮助推断转录因子的结合 motif。

## 推荐工具

### 1. MEME-ChIP (推荐)

**描述**: MEME Suite 中专门用于 ChIP-seq 数据 motif 分析的工具，整合了 MEME 和 DREME 两种算法。

**工作流程**:
1. 从 MACS2 输出中提取 peak 峰顶区域 (summit +/- 100bp)
2. 运行 MEME-ChIP 进行 motif 发现
3. 与已知 motif 数据库比对 (JASPAR, CIS-BP)

**关键参数**:
- `-db <motif_database>`: 已知 motif 数据库
- `-meme-mod anr`: MEME 算法模式
- `-nmeme <int>`: 最大 motif 发现数
- `-centrimo`: 中心富集分析

**示例**:
```bash
meme-chip -oc meme_output/ \
  -db JASPAR2022_CORE_plants_non_redundant.meme \
  -meme-mod anr -nmeme 5 \
  peak_summit_sequences.fa
```

**输出**:
- 富集的 motif logos
- Motif 位置分布 (centrimo)
- 已知 motif 匹配结果

### 2. HOMER (findMotifsGenome.pl)

**描述**: 包含从头 motif 发现和已知 motif 富集分析。

**优势**: 与参考基因组集成，可直接使用 BED 文件

**示例**:
```bash
findMotifsGenome.pl peaks.bed tair10 homer_output/ \
  -size given -p 4
```

### 3. STREME

**描述**: MEME Suite 中的快速 motif 发现工具，适用于大规模数据集。

**适用场景**: 大量 peaks 时作为 MEME 的快速替代方案。

---

## 植物 Motif 分析注意事项

1. **植物特异性 motif 数据库**:
   - JASPAR Plants: https://jaspar.genereg.net/
   - PlantTFDB: http://planttfdb.gao-lab.org/
   - CIS-BP: http://cisbp.ccbr.utoronto.ca/

2. **基因组背景模型**: 使用植物基因组序列作为背景

3. **植物 TF 家族特征**:
   - MYB: 保守的 AAC 核心
   - WRKY: W-box (TTGACC/T)
   - NAC: NAC 识别序列
   - bZIP: ACGT 核心 (ABRE, G-box)
   - AP2/ERF: GCC-box, DRE/CRT

4. **中心富集分析 (CentriMo)**: 检测 motif 是否富集在 peak 中心，验证 ChIP 实验质量
