# Hi-C 染色体挂载 (Scaffolding) -- 工具目录

## 概述

Hi-C 染色体挂载利用染色质构象捕获数据（Hi-C）的染色体内接触频率关系，将 contig 级组装提升到染色体级。Hi-C 是基于测序的染色体结构的实验技术，其基本逻辑是：同一染色体上的邻近区域比远处区域的接触频率高，这是挂载的依据。

## 推荐工具

### 1. YAHS (Yet Another Hi-C Scaffolder) -- 首选

**描述**: YAHS 是高性能的 Hi-C scaffolding 工具，使用贪婪算法进行 contig 排序和定向，然后用重叠检测方法合并 contig。

**安装**:
```bash
conda install -c bioconda yahs
```

**标准流程**:
```bash
# 步骤 1: 将 Hi-C reads 比对到组装
bwa mem -t 16 -5SP assembly.fasta \
  hic_R1.fastq.gz hic_R2.fastq.gz \
  > hic_aligned.sam

# 步骤 2: 转换为 BAM 并排序
samtools view -bS hic_aligned.sam | \
  samtools sort -@ 8 -o hic_sorted.bam
samtools index hic_sorted.bam

# 步骤 3: YAHS scaffolding
yahs assembly.fasta hic_sorted.bam \
  -o yahs_output \
  --no-contig-ec
```

**关键参数**:
| 参数 | 说明 |
|------|------|
| `-o` | 输出前缀 |
| `--no-contig-ec` | 禁用 contig 错误纠正 |
| `-l` | 最小 contig 长度 (默认 50000) |
| `--enzyme` | 限制性内切酶 (如 DpnII, HindIII) |

### 2. SALSA2 -- 备选

**描述**: SALSA2 是另一个广泛使用的 Hi-C scaffolding 工具，支持多酶切 Hi-C 数据。

```bash
# SALSA2 scaffolding
python run_pipeline.py \
  -a assembly.fasta \
  -l assembly.fasta.fai \
  -b hic_sorted.bam \
  -e GATC,DpnII \
  -o salsa_output/ \
  -m yes
```

### 3. 3D-DNA / Juicer (备选)

**描述**: 3D-DNA + Juicer 流程是另一个 Hi-C 完整分析流程。

```bash
# Juicer 预处理
juicer.sh -t 16 \
  -g genome_name \
  -s MboI \
  -z assembly.fasta \
  -y restriction_sites.txt \
  -p chrom.sizes \
  -d juicer_work/

# 3D-DNA scaffolding
run-asm-pipeline.sh assembly.fasta \
  juicer_work/aligned/merged_nodups.txt
```

---

## Hi-C 数据预处理

### 限制性内切酶选择

植物 Hi-C 常用酶：
| 酶 | 识别位点 | 特点 |
|-----|---------|------|
| DpnII | GATC | 最常见的 Hi-C 酶（4 碱基切点） |
| MboI | GATC | DpnII 的同切酶 |
| HindIII | AAGCTT | 6 碱基切点，片段更大 |
| Sau3AI | GATC | 通用 |

### 生成酶切位点文件

```bash
# 使用 samtools + seqkit 生成酶切位点
# DpnII = GATC
seqkit locate -p GATC assembly.fasta \
  > dpnii_sites.bed
```

---

## 挂载后的人工校正

Hi-C 自动挂载通常需要手动校正。使用 Juicebox 进行：

### 方法 1: Juicebox 手动校正

1. 在 Hi-C contact map 中识别 misassemblies（表现为对角线的"断点"或"偏移"）
2. 手动划分染色体边界
3. 导出校正后的 assembly

### 方法 2: PretextView 校正

```bash
# 生成 Pretext map
PretextMap -i hic_sorted.bam -o pretext_map.pretext

# 在 PretextView 中手动校正
PretextView -i pretext_map.pretext
```

---

## 挂载后评估

### Hi-C Contact Map 可视化

```bash
# 使用 HiCExplorer 生成 contact map
hicBuildMatrix -s hic_sorted.bam \
  --binSize 50000 \
  -o hic_matrix.h5

hicPlotMatrix -m hic_matrix.h5 \
  -o hic_contact_map.png \
  --log1p --dpi 300 \
  --title "Hi-C Contact Map"
```

### 挂载质量指标

| 指标 | 优秀 | 良好 | 需要检查 |
|------|------|------|---------|
| Scaffold N50 | > 50Mb | > 20Mb | < 10Mb |
| 悬空 contig 百分比 | < 5% | < 10% | > 15% |
| 染色体数匹配 | 匹配核型 | 偏离 < 20% | 严重偏离 |
| Contact map diagonal | 整洁的对角线 | 有模糊但可辨别 | 信号混乱 |
| BUSCO post-scaffolding | 不降 | 轻微降低 < 5% | 大幅降低 |

### 常见 Misassembly 信号

在 Hi-C contact map 中检查以下信号：

1. **对角线偏移**: 可能表示 inversion 或 misorientation
2. **对角线中断**: 表示可能的 translocation
3. **"Christmas tree" 模式**: 着丝粒附近的正常信号（植物中尤其明显）
4. **强 off-diagonal 信号**: 可能表示染色体间的 misassembly 或真实的 inter-chromosomal 互作

---

## 植物 Hi-C Scaffolding 注意事项

1. **植物基因组大小**: 小麦 (17Gb) 或玉米 (2.4Gb) 的大基因组需要更大的 bin size (100-500kb) 以降低噪音
2. **多倍体**: 异源多倍体的亚基因组之间可能存在较强的 Hi-C 信号（homeologous interactions），导致染色体分群错误
3. **端粒-着丝粒信号**: 植物着丝粒通常较大（以 Mb 计），Hi-C 信号在此区域可能有特殊模式
4. **叶绿体和线粒体**: 细胞器不参与 Hi-C 相互作用（没有核小体），在 Hi-C 分析中出现随机信号，在 scaffolding 前需去除
5. **酶切位点密度**: 植物基因组 AT 含量高，DpnII (GATC) 的位点密度通常足够
