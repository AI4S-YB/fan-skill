# TF Footprinting 分析 -- 工具目录

## 概述

转录因子足迹分析（TF Footprinting）利用 ATAC-seq 数据中 Tn5 切割模式的变化来推断特定转录因子在基因组上的结合位置。当 TF 结合 DNA 时，会在结合位点附近产生一个"足迹"——即 Tn5 切割受到保护的信号。

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| --window | 50 | Narrow footprints (ARF, GRF families): reduce to 30; broad binding TFs (MADS-box): increase to 100 | Window size determines footprint detection resolution; plant TF families have different footprint widths depending on DNA binding domain structure |
| --score_method | "bound" | Low coverage data (<30M reads): use "difference"; well-powered data (>50M reads): use "bound" | "bound" uses flanking accessibility to normalize footprint scores; "difference" is simpler but sensitive to coverage noise at low depth |
| motif database | JASPAR Plants (non-redundant) | Non-model crop: supplement with PlantTFDB motifs; well-studied species: include DAP-seq derived motifs | JASPAR Plants covers conserved motifs across species; family-specific and species-specific binding preferences need custom databases for comprehensive analysis |
| --pvalue (BINDetect) | 0.05 | Exploratory hypothesis generation: 0.10; high-confidence regulatory inference: 0.01 | Balance between discovering novel regulators in underexplored plant systems (relax) and minimizing false positives for well-characterized pathways (tighten) |
| --min_score | 0.5 | Divergent species (gymnosperms, ferns): reduce to 0.3; conserved angiosperm crops: keep 0.5 | Motif match stringency; lower scores capture more divergent binding sites in evolutionarily distant plant lineages |

## 推荐工具

### 1. TOBIAS (Transcription factor Occupancy prediction By Investigation of ATAC-seq Signal) -- 金标准

**描述**: TOBIAS 是目前最成熟的 TF footprinting 工具套件，专为 ATAC-seq 设计，包含偏差校正、motif scoring 和差异 footprint 检测功能。

**安装**:
```bash
pip install TOBIAS
```

**完整分析流程**:

#### 步骤 1: ATACorrect (Tn5 切割偏差校正)

```bash
# 对每个样本执行偏差校正
TOBIAS ATACorrect \
  --bam sample1_dedup.bam \
  --genome genome.fa \
  --peaks sample1_peaks.narrowPeak \
  --outdir tobias_atacorrect/ \
  --cores 16
```

**输出**: `<sample>_corrected.bw` -- 偏差校正后的连续信号 track

#### 步骤 2: ScoreMotifs (Motif Footprinting 评分)

```bash
TOBIAS ScoreMotifs \
  --signals condition1_corrected.bw condition2_corrected.bw \
  --motifs JASPAR2022_CORE_plants_non_redundant.meme \
  --genome genome.fa \
  --output tobias_motif_scores/ \
  --cores 16
```

**Motif 数据库选择**:
- 植物: `JASPAR2022_CORE_plants_non_redundant.meme` (JASPAR 植物 motif 集合)
- 广泛: `JASPAR2022_CORE_non_redundant.meme` (全部 motif)
- 自定义: 从 PlantTFDB 或 CIS-BP 提取特异性物种的 motif

#### 步骤 3: BINDetect (差异 Footprint 检测)

```bash
TOBIAS BINDetect \
  --signals condition1_corrected.bw condition2_corrected.bw \
  --condition_names control treatment \
  --motifs JASPAR2022_CORE_plants_non_redundant.meme \
  --genome genome.fa \
  --peaks consensus_peaks.bed \
  --outdir tobias_bindetect/ \
  --cores 16
```

**输出**: 每个 TF motif 的差异结合分数表格

#### 步骤 4: 结果可视化

```r
# R 中读取 TOBIAS 结果
library(TOBIAS)

# 读取 BINDetect 结果
bindetect_results <- read.table(
  "tobias_bindetect/bindetect_results.txt",
  header = TRUE, sep = "\t"
)

# 显著差异结合的 TF
sig_tfs <- bindetect_results[
  bindetect_results$significant == "True",
]
head(sig_tfs[order(sig_tfs$differential_score, decreasing = TRUE), ])
```

### 2. HINT-ATAC (备选)

**描述**: Regulatory Genomics Toolkit (RGT) 中的 HINT-ATAC 工具，可用于 footprint 分析。

```bash
# 使用 RGT-HINT 进行 footprinting
rgt-hint footprinting \
  --atac-seq \
  --organism=arabidopsis_thaliana \
  --output-location=hint_output/ \
  --output-prefix=sample \
  sample_dedup.bam sample_peaks.narrowPeak
```

---

## TOBIAS 参数详解

| 参数 | 说明 | 推荐值 |
|------|------|--------|
| `--cores` | 并行线程数 | 8-32 |
| `--split` | 按染色体分割 (大基因组推荐) | 1 (启用) |
| `--window` | Footprint 滑窗大小 | 50 (默认) |
| `--score_method` | Scoring 算法 | "bound" (默认) |
| `--pvalue` | BINDetect 显著性阈值 | 0.05 |
| `--min_score` | 最小 motif score | 0.5 |

---

## 植物 TF Footprinting 注意事项

1. **使用植物特异性 Motif 数据库**: JASPAR Plants 集合或 PlantTFDB 提取的 motif
2. **大基因组**: 对玉米（2.4Gb）、小麦（17Gb）等大基因组，使用 `--split` 选项
3. **Motif 冗余**: 同一 TF 家族成员（如 MYB、NAC、WRKY）的 motif 可能高度相似
4. **多倍体**: 对异源多倍体，考虑各亚基因组分别分析
5. **覆盖深度**: Footprinting 分析需要较高的测序深度（> 50M reads），太低的深度会导致足迹信号不可靠
6. **拟南芥特异性**: 拟南芥的 TF 绑定位点变异较大，可能比动物 TF motif 匹配度低

---

## 输出解读

### Footprint Score 含义
- **高正分 (> 1)**: 强烈的 TF 结合信号（保护模式明显）
- **接近 0**: 无明显的保护信号
- **负分 (< -1)**: 开放的信号（与 footprint 相反，常为开放的启动子区域）

### 差异 Footprint
- **differential_score > 0**: 处理组中 TF 结合增强
- **differential_score < 0**: 处理组中 TF 结合减弱

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| TOBIAS ATACorrect "signal too weak" or all-zero tracks | Low sequencing coverage or poor Tn5 efficiency | Increase sequencing depth to >50M uniquely mapped reads per sample; verify nuclei isolation quality with positive control tissue |
| Footprint scores near zero for all TFs | Wrong motif database or species mismatch | Use JASPAR Plants non-redundant collection specifically; verify motif format is MEME; confirm genome FASTA matches species |
| No significant differential TFs in BINDetect | Insufficient biological difference between conditions or overly strict threshold | Verify experimental design has true biological contrast; consider lowering pvalue to 0.10 for hypothesis generation |
| TOBIAS ScoreMotifs fails with "motif too long" warning | Plant TF motifs often longer than animal motifs (MADS-box proteins have 10+ bp core motifs) | Pre-process motif file to trim flanking low-information-content bases; set appropriate motif width cutoff |
