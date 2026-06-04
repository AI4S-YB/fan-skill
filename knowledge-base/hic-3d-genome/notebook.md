# 三维基因组分析 -- 分析笔记本

## 分析概览

本笔记本提供植物三维基因组分析的完整流程，涵盖Hi-C数据处理、Loop检测、TAD识别、Compartment分析和差异分析。

---

## 1. Hi-C数据处理流程

### 1.1 原始数据质控

```bash
# FastQC评估原始数据质量
fastqc -t 8 hic_R1.fastq.gz hic_R2.fastq.gz -o qc_output/

# Trimmomatic去除接头
trimmomatic PE -threads 8 \
  hic_R1.fastq.gz hic_R2.fastq.gz \
  hic_R1_trimmed.fq.gz hic_R1_unpaired.fq.gz \
  hic_R2_trimmed.fq.gz hic_R2_unpaired.fq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### 1.2 Juicer流程 (推荐)

```bash
# Juicer完整流程
juicer.sh -t 16 \
  -g genome_name \
  -s MboI \
  -z reference_genome.fasta \
  -y restriction_sites.txt \
  -p chrom.sizes \
  -d juicer_output/

# 参数说明：
# -t: 线程数
# -g: 基因组名称
# -s: 限制性内切酶类型 (MboI/HindIII/DpnII)
# -z: 参考基因组路径
# -y: 酶切位点文件
# -p: 染色体大小文件
```

**Juicer输出文件**：
- `aligned/merged_nodups.txt`: 去重后的比对结果
- `hic_files/inter.hic`: .hic格式contact矩阵
- `aligned/inter_30.hic`: MAPQ≥30的矩阵

### 1.3 使用Cooler构建矩阵

```bash
# Step 1: 比对
bwa mem -t 16 reference_genome.fasta \
  hic_R1.fastq.gz hic_R2.fastq.gz > alignments.sam

# Step 2: 处理比对结果
samtools sort -@ 8 -o alignments.bam alignments.sam
samtools index alignments.bam

# Step 3: 构建像素列表
cooler csort reference_genome.fasta.fai \
  -i 3 -i 7 alignments.bam > pixels.tsv

# Step 4: 构建矩阵
cooler load -f bedpe reference_genome.fasta.fai:10000 \
  pixels.tsv matrix.cool

# Step 5: 标准化
cooler balance matrix.cool
```

---

## 2. Contact Map可视化

### 2.1 Juicebox可视化

```bash
# 命令行生成图片
java -jar juicebox_tools.jar \
  dump observed KR inter.hic BP 1000000 contact_matrix.txt

# 或使用Juicebox GUI交互式查看
# 下载地址: https://github.com/aidenlab/Juicebox
```

### 2.2 HiCExplorer可视化

```bash
# 生成contact map图片
hicPlotMatrix -m matrix.cool \
  -o contact_map.png \
  --log1p \
  --dpi 300 \
  --chromosomeOrder chr1 chr2 chr3

# 添加轨迹
hicPlotMatrix -m matrix.cool \
  -o contact_map_with_tracks.png \
  --track tracks.ini
```

---

## 3. Loop检测

### 3.1 HiCCUPS (高分辨率数据)

```bash
# 需要GPU支持
hiccups -m 5000 -r 5000,10000,25000 \
  -k KR -f .1 -p 4 -i 7 -t 0.02,1.5,1.75,2 \
  inter.hic loop_output/

# 参数说明：
# -m: 最小距离
# -r: 分辨率列表
# -k: 标准化方法
# -f: FDR阈值
```

### 3.2 Mustache (中等分辨率)

```bash
# 安装
pip install mustache-hic

# 运行
mustache -p matrix.cool \
  -o loops.bedpe \
  -res 10000 \
  -dist 2000000
```

### 3.3 Fit-Hi-C (低分辨率/植物)

```bash
# Fit-Hi-C适合植物基因组
FitHiC.py -i interactions.txt \
  -o fithic_output \
  -b 10000 \
  -L 0 -U 10000000
```

---

## 4. TAD识别

### 4.1 Arrowhead (Juicer流程)

```bash
# Arrowhead TAD检测
arrowhead -r 10000 -k KR \
  inter.hic tad_output/

# 输出: tad_output/blocks.txt
```

### 4.2 Insulation Score方法

```bash
# 使用CoolTools计算insulation score
cooltools insulation matrix.cool \
  -o insulation.tsv \
  --window 500000

# 识别TAD边界
python identify_tad_boundaries.py insulation.tsv
```

### 4.3 HiCExplorer TAD检测

```bash
# TAD识别
hicDetectTADs -m matrix.cool \
  -o tad_output/ \
  --minDepth 60000 \
  --maxDepth 1000000 \
  --step 20000

# 可视化
hicPlotTADs -m matrix.cool \
  -t tad_output/tad_score.bedgraph \
  -o tad_plot.png
```

---

## 5. Compartment分析

### 5.1 特征向量分解

```bash
# Juicer方法
eigenvector inter.hic KR BP 1000000 eigenvector.txt

# Cooler方法
cooltools eigdecomp matrix.cool \
  -o eigs.tsv \
  --phasing-track gc_content.bw
```

### 5.2 A/B区室划分

```bash
# 根据特征向量符号划分A/B区室
# 正值 = A区室 (开放染色质)
# 负值 = B区室 (异染色质)

# 与基因密度关联
python correlate_with_gene_density.py eigs.tsv gene_density.bw
```

---

## 6. 差异分析

### 6.1 MultiHiCnorm

```bash
# 多样本标准化和差异分析
multihicnorm -i sample1.cool sample2.cool \
  -o normalized/ \
  -b 10000

# 差异检测
diff_hic -i normalized/ \
  -o differential_results.tsv
```

### 6.2 SELFISH

```bash
# SELFISH差异分析
selfish -m sample1.cool sample2.cool \
  -o selfish_output/ \
  -r 10000 \
  -p 0.05
```

---

## 7. 植物特异性考量

### 7.1 基因组大小影响

| 基因组大小 | 建议测序量 | 推荐分辨率 |
|-----------|-----------|-----------|
| < 500 Mb | 100-200M reads | 5-10 kb |
| 500 Mb - 2 Gb | 200-500M reads | 10-25 kb |
| > 2 Gb | > 500M reads | 25-50 kb |

### 7.2 染色体数量差异

```bash
# 对于高染色体数物种，需要特别处理染色体命名
# 确保chrom.sizes文件正确
awk '{print $1"\t"$2}' reference_genome.fasta.fai > chrom.sizes

# 可视化时指定染色体顺序
hicPlotMatrix -m matrix.cool \
  --chromosomeOrder $(cat chrom_order.txt | tr '\n' ' ')
```

### 7.3 多倍体物种

```bash
# 多倍体物种的同源染色体可能难以区分
# 使用更严格的比对参数
bwa mem -t 16 -a reference_genome.fasta \
  hic_R1.fastq.gz hic_R2.fastq.gz | \
  samtools view -q 30 -b > alignments.bam
```

---

## 8. 质量评估

### 8.1 关键QC指标

| 指标 | 合格标准 | 优秀标准 |
|------|---------|---------|
| 比对率 | > 70% | > 85% |
| 有效配对 | > 50% | > 70% |
| Duplicate率 | < 30% | < 20% |
| 跨片段比例 | > 50% | > 70% |
| 预期距离峰 | 明显 | 锐利 |

### 8.2 QC命令

```bash
# 使用HiC-Pro生成QC报告
hicpro2juicebox.py -i .hicpro_output/ -o juicebox_output/

# 或使用pairtools
pairtools stats alignments.pairs \
  -o stats.txt
```

---

## 常见问题

### Q: Contact map看起来很"模糊"怎么办？
A: 检查测序深度和比对率。提高分辨率或增加测序量。

### Q: TAD边界不明显？
A: 植物基因组的TAD边界可能不如动物清晰。尝试不同的TAD检测参数。

### Q: Loop检测不到？
A: Loop检测需要高分辨率数据（<10kb）。确保测序深度足够。

### Q: Compartment分析与预期不符？
A: 检查GC含量phasing track是否正确。植物基因组可能需要不同的phasing策略。

---

## 参考

- Juicer: https://github.com/aidenlab/juicer
- Cooler: https://github.com/open2c/cooler
- CoolTools: https://github.com/open2c/cooltools
- HiCExplorer: https://github.com/deeptools/HiCExplorer
- Fit-Hi-C: https://github.com/AY-lab/Fit-Hi-C
- Juicebox: https://github.com/aidenlab/Juicebox
