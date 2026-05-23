# Peak Calling (峰值检测) — 工具目录

## 概述

Peak calling 是 ChIP-seq 分析的核心步骤，用于识别蛋白质-DNA 结合位点。不同靶标类型需要不同的 peak calling 策略。

## 推荐工具

### 1. MACS2 (Model-based Analysis of ChIP-Seq) — 强烈推荐

**描述**: ChIP-seq peak calling 的金标准工具，支持 narrow peak 和 broad peak 两种模式。

#### Narrow Peak (转录因子 ChIP)

**适用场景**: 转录因子或序列特异性 DNA 结合蛋白的 ChIP-seq。

**关键参数**:
- `-t <treatment>`: 处理组/IP 样本 BAM 文件
- `-c <control>`: 对照组/Input 样本 BAM 文件
- `-f BAM`: 输入文件格式
- `-g <genome_size>`: 有效基因组大小
  - 拟南芥: `-g 1.19e8`
  - 水稻: `-g 3.74e8`
  - 玉米: `-g 2.04e9`
  - 大豆: `-g 9.75e8`
- `-n <name>`: 输出文件前缀
- `--nomodel`: 无模型模式 (植物推荐)
- `--extsize <int>`: 片段延伸长度 (根据超声打断大小)

**示例 (TF ChIP)**:
```bash
macs2 callpeak -t tf_ip_sorted.bam -c input_sorted.bam \
  -f BAM -g 1.19e8 \
  -n TF_sample --nomodel --extsize 200 \
  --outdir macs2_output/
```

#### Broad Peak (组蛋白修饰 ChIP)

**适用场景**: 组蛋白修饰的 ChIP-seq (如 H3K4me3, H3K27ac, H3K27me3, H3K9me2)。

**关键参数**:
- `--broad`: 启用 broad peak 模式
- `--broad-cutoff <float>`: broad peak 阈值 (默认 0.1)

**示例 (组蛋白修饰 ChIP)**:
```bash
macs2 callpeak -t h3k27ac_ip_sorted.bam -c input_sorted.bam \
  -f BAM -g 1.19e8 \
  -n H3K27ac_sample --broad --broad-cutoff 0.1 \
  --outdir macs2_output/
```

### 2. HOMER (备选)

**描述**: 包含 findPeaks 工具，适用于转录因子和组蛋白修饰 ChIP-seq。

**适用场景**: 需要快速初步分析时使用。

**输出文件**:
- `*_peaks.narrowPeak`: narrow peak 结果 (BED 格式)
- `*_peaks.broadPeak`: broad peak 结果 (BED 格式)
- `*_peaks.xls`: 详细的 peak 统计信息
- `*_summits.bed`: peak 峰顶位置 (用于 motif 分析)

---

## 植物 ChIP-seq Peak Calling 注意事项

1. **基因组大小参数**: 必须使用植物的有效基因组大小 (`-g` 参数)，不是总基因组大小
2. **Input 对照**: 强烈推荐使用 Input DNA 作为对照，去除背景噪音
3. **去重复**: peak calling 前需要对 BAM 文件去重复
4. **黑名单过滤**: 去除基因组中的异常高信号区域（目前主要是模式物种有黑名单）
5. **植物特有组蛋白修饰**:
   - H3K9me2: 植物中与转座子沉默相关
   - H3K27me3: 植物中广泛分布于基因体
   - H3K4me3: 活跃启动子标记
