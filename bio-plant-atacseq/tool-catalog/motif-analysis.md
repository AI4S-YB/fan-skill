# Motif 富集分析 -- 工具目录

## 概述

Motif 分析用于从 ATAC-seq peaks 中识别富集的转录因子结合 motif，揭示调控开放染色质的潜在转录因子。

## 推荐工具

### 1. HOMER (Hypergeometric Optimization of Motif EnRichment) -- 首选

**描述**: HOMER 是目前最常用的 de novo 和已知 motif 富集分析工具。特别适合植物 ATAC-seq 数据的启动子/增强子 motif 分析。

**安装**:
```bash
# 通过 conda 安装
conda install -c bioconda homer

# 配置植物基因组数据库
perl /path/to/homer/configureHomer.pl -install athaliana
perl /path/to/homer/configureHomer.pl -install oryza
```

**步骤 1: 提取 peak 序列**:
```bash
bedtools getfasta -fi genome.fa \
  -bed sample_peaks.narrowPeak \
  -fo peaks_sequences.fa
```

**步骤 2: Motif 发现**:
```bash
# De novo motif 发现
findMotifsGenome.pl sample_peaks.narrowPeak \
  genome.fa homer_output/ \
  -size 200 -mask -p 8
```

**参数说明**:
- `-size 200`: 用于 motif 搜索的区域大小（ATAC-seq 推荐 150-200bp）
- `-mask`: 屏蔽重复序列
- `-p 8`: 线程数
- `-len 6,8,10,12,15`: motif 长度范围

### 2. MEME-ChIP (MEME Suite) -- 备选

**描述**: MEME-ChIP 是整合的 ChIP/ATAC-seq motif 分析工具，支持 de novo 和已知 motif 分析。

```bash
meme-chip -oc meme_output/ \
  -db JASPAR2022_CORE_plants_non_redundant.meme \
  -meme-mod anr \
  -nmeme 5 \
  peaks_sequences.fa
```

**关键参数**:
| 参数 | 说明 |
|------|------|
| `-oc` | 输出目录 |
| `-db` | 已知 motif 数据库 |
| `-meme-mod anr` | MEME 的 de novo 发现模式（任意重复次数） |
| `-nmeme 5` | 最多发现 5 个 de novo motifs |
| `-centrimo-score 5.0` | CentriMo 显著性阈值 |

### 3. RSAT plants (备选)

**描述**: Regulatory Sequence Analysis Tools (RSAT) 提供专门针对植物的 motif 分析模块。

在线工具: http://rsat.eead.csic.es/plants/

---

## 植物特有 Motif 数据库

| 数据库 | URL | 描述 |
|--------|-----|------|
| PlantTFDB | http://planttfdb.gao-lab.org/ | 植物转录因子数据库 |
| JASPAR Plants | https://jaspar.genereg.net/ | JASPAR 植物 motif 合集 |
| CIS-BP | http://cisbp.ccbr.utoronto.ca/ | 全面的 TF 结合特异性 |
| PlantPAN | http://plantpan.itps.ncku.edu.tw/ | 植物启动子分析平台 |
| AthaMap | http://www.athamap.de/ | 拟南芥 TF 结合位点 |
| PlantRegMap | http://plantregmap.gao-lab.org/ | 植物转录调控图谱 |
| DAP-seq | https://neomorph.salk.edu/dap_web/pages/index.php | 拟南芥 TF 结合数据（cistrome） |

**提取植物 TF motif 为 HOMER 格式**:
```bash
# 从 PlantTFDB 下载植物 TF motif
# 转换为 HOMER 格式后使用
```

---

## 分析策略

### De novo motif 发现

推荐流程：
1. 取 peak summit 周围 +/- 100bp 序列
2. 使用 HOMER 或 MEME-ChIP 进行 de novo motif 发现
3. 使用 Tomtom (MEME Suite) 将 de novo motif 与已知数据库比对
4. 如果匹配到已知植物 TF，则获得了调控因子的候选

### 已知 motif 富集

1. 取所有 peaks 的序列作为 foreground
2. 取基因组上随机匹配 GC 含量的等长区域作为 background
3. 使用 AME (MEME Suite) 或 HOMER 的已知 motif 扫描
4. 用 Fisher 精确检验评估每种 motif 在 peaks 中的富集程度

**AME 示例**:
```bash
ame --control background_sequences.fa \
  --o ame_output/ \
  peaks_sequences.fa \
  JASPAR2022_CORE_plants_non_redundant.meme
```

---

## 植物 Motif 分析注意事项

1. **植物 motif 短且简并**: 植物 TF 的 DNA 结合域比动物更保守（如 MADS-box、MYB、NAC），但结合位点较短（6-10 bp），假阳性率高
2. **AT-rich bias**: 植物启动子区域 AT 含量高于动物，motif 发现的背景模型需考虑此特性
3. **低复杂度序列**: 植物基因组中简单重复序列（如 (AT)n、(GA)n）较多，必须屏蔽
4. **转座子**: 植物基因组 TE 含量高，TE 中的 TF 结合位点可能是功能性的（如 RdDM 靶向），也可能只是噪音
5. **组织特异性**: 相同 motif 在不同组织中的功能可能完全不同
6. **同源 TF**: 植物 TF 家族通常有数十个成员（如拟南芥 MYB 家族 ~190 个成员），单个 motif 可能对应多个 TF
