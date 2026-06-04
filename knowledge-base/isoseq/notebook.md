# Iso-Seq全长转录组分析 -- 分析笔记本

## 分析概览

本笔记本提供PacBio Iso-Seq全长转录组分析的完整流程，涵盖CCS生成、Isoform聚类、错误校正、转录本注释和可变剪接分析。

---

## 1. CCS生成

### 1.1 Lima (Sequel II/Revio)

```bash
# 生成CCS
lima movie.subreads.bam barcodes.fasta \
  demux.ccs.bam \
  --peek-guess 0.9 \
  --min-passes 3 \
  --min-rq 0.99 \
  --isoseq

# 提取CCS FASTQ
bam2fastq demux.ccs.bam -o ccs.fastq
```

### 1.2 CCS工具 (Sequel I)

```bash
# 生成CCS
ccs movie.subreads.bam ccs.bam \
  --minPasses 3 \
  --minAccuracy 0.99 \
  --threads 16

# 提取FASTQ
bam2fastq ccs.bam -o ccs.fastq
```

### 1.3 质量检查

```bash
# 检查CCS质量
samtools view ccs.bam | \
  awk '{print length($10)}' | \
  Rscript plot_length_distribution.R

# 统计全长比例
grep -c "full_length" ccs.fastq
```

---

## 2. Isoform聚类

### 2.1 IsoSeq3流程

```bash
# 步骤1: 识别全长和非全长读段
isoseq3 cluster ccs.bam clustered.bam \
  --verbose \
  --threads 16

# 步骤2: 提取高质量转录本
isoseq3 refine clustered.bam \
  hq_transcripts.fasta \
  --min-full-length-passes 3 \
  --min-accuracy 0.99

# 输出：
# hq_transcripts.fasta - 高质量转录本
# lq_transcripts.fasta - 低质量转录本
```

### 2.2 聚类参数说明

| 参数 | 默认值 | 说明 |
|------|--------|------|
| --min-full-length-passes | 3 | 最小全长通过数 |
| --min-accuracy | 0.99 | 最小准确度 |
| --min-aln-coverage | 0.99 | 最小比对覆盖度 |

---

## 3. 错误校正

### 3.1 LoRDEC

```bash
# 构建k-mer数据库
lordec-build -i illumina_R1.fastq -i illumina_R2.fastq \
  -k 21 -s 3 \
  -o kmer_db

# 校正长读长
lordec-correct -i hq_transcripts.fasta \
  -k kmer_db \
  -o corrected_transcripts.fasta \
  -t 16
```

### 3.2 Proovread

```bash
# Proovread校正
proovread --longreads hq_transcripts.fasta \
  --shortreads illumina_R1.fastq,illumina_R2.fastq \
  --threads 16 \
  --output corrected_transcripts.fasta
```

---

## 4. 比对和去冗余

### 4.1 Minimap2比对

```bash
# 比对到参考基因组
minimap2 -ax splice:hq -uf \
  -t 16 \
  reference.fasta \
  hq_transcripts.fasta \
  > alignment.sam

# 转换和排序
samtools view -bS alignment.sam | \
  samtools sort -o alignment_sorted.bam
samtools index alignment_sorted.bam
```

### 4.2 去冗余

```bash
# 使用Cupcake去冗余
# 首先排序和折叠
sort -k 3,3 -k 4,4n alignment_sorted.sam > sorted.sam
collapse_isoforms_by_sam.py \
  -i sorted.sam \
  -o collapsed \
  --min-identity 0.99

# 输出：
# collapsed.filtered.fasta - 去冗余后转录本
# collapsed.filtered.gff - GFF格式
```

---

## 5. 转录本注释 (SQANTI3)

### 5.1 运行SQANTI3

```bash
# SQANTI3完整注释
python sqanti3.py \
  collapsed.filtered.fasta \
  genes.gtf \
  reference.fasta \
  -t 16 \
  -o sqanti3_output \
  --coverage 0.8 \
  --identity 0.9

# 输出文件：
# sqanti3_output_classification.txt - 分类结果
# sqanti3_output_junctions.txt - 剪接位点信息
# sqanti3_output.gff - GFF格式注释
```

### 5.2 分类结果解读

| 类别 | 说明 |
|------|------|
| FSM | Full Splice Match - 完全匹配已知转录本 |
| ISM | Incomplete Splice Match - 部分匹配已知转录本 |
| NIC | Novel In Catalog - 新的组合但使用已知剪接位点 |
| NNC | Novel Not in Catalog - 使用新剪接位点 |
| Antisense | 反义转录本 |
| Intergenic | 基因间区转录本 |
| Genic | 基因内区转录本 |

---

## 6. 可变剪接分析

### 6.1 SUPPA2

```bash
# 生成剪接事件
suppa.py generateEvents \
  -i annotation.gtf \
  -o splicing_events \
  -f ioe -e SE SS MX RI FL

# 定量剪接事件
suppa.py psiPerEvent \
  -i splicing_events.ioe \
  -r transcript_expression.tsv \
  -o psi_values
```

### 6.2 剪接事件类型

- **SE**: Skipped Exon - 外显子跳跃
- **SS**: Alternative Splice Site - 可变剪接位点
- **MX**: Mutually Exclusive Exon - 互斥外显子
- **RI**: Retained Intron - 内含子保留
- **FL**: Alternative First/Last Exon - 可变首/末外显子

---

## 7. 融合基因检测

### 7.1 FusionSeeker

```bash
# 检测融合转录本
python FusionSeeker.py \
  -i alignment_sorted.bam \
  -r reference.fasta \
  -g genes.gtf \
  -o fusion_output
```

### 7.2 过滤假阳性

```bash
# 过滤步骤
# 1. 去除同一基因内的"融合"
# 2. 去除低表达事件
# 3. 检查断点序列

python filter_fusions.py \
  fusion_output.txt \
  --min-reads 3 \
  --min-spanning 2
```

---

## 8. 多聚腺苷酸分析

### 8.1 PolyA尾分析

```bash
# 提取PolyA信息
# 需要保留原始subreads信息
extract_polya.py \
  -b alignment_sorted.bam \
  -o polya_results.tsv
```

### 8.2 APA位点分析

```bash
# 分析可变多聚腺苷酸化
identify_apa_sites.py \
  -i collapsed.filtered.fasta \
  -r reference.fasta \
  -g genes.gtf \
  -o apa_sites.tsv
```

---

## 9. 植物特异性分析

### 9.1 大规模转录本处理

```bash
# 植物转录组通常较大，需要分块处理
# 将转录本分成多个批次
split -l 10000 hq_transcripts.fasta transcripts_chunk_

# 并行处理
for chunk in transcripts_chunk_*; do
  minimap2 -ax splice:hq reference.fasta $chunk > ${chunk}.sam &
done
wait
```

### 9.2 多倍体同源基因区分

```bash
# 多倍体物种的同源基因区分
# 使用位置和序列信息
python distinguish_homologs.py \
  -a alignment_sorted.bam \
  -g annotation.gtf \
  -o homolog_assignments.tsv
```

---

## 10. 质量控制

### 10.1 QC指标

| 指标 | 合格标准 | 优秀标准 |
|------|---------|---------|
| 全长比例 | > 50% | > 70% |
| 平均读长 | > 1 kb | > 2 kb |
| HQ转录本比例 | > 80% | > 90% |
| N50 | > 1.5 kb | > 2.5 kb |

### 10.2 QC命令

```bash
# 计算全长比例
samtools view ccs.bam | \
  awk '/full_length/{full++} END{print full/NR}'

# 计算N50
seqkit stats hq_transcripts.fasta
```

---

## 常见问题

### Q: 全长比例低？
A: 检查RNA质量，增加反转录时间，检查测序参数。

### Q: 转录本太多？
A: 增加聚类严格度，使用去冗余工具。

### Q: FSM比例低？
A: 检查注释质量，可能是物种特异性或发现新转录本。

### Q: 剪接事件过多？
A: 过滤低表达事件，检查比对质量。

---

## 参考

- IsoSeq3: https://github.com/PacificBiosciences/IsoSeq
- Lima: https://github.com/PacificBiosciences/lima
- SQANTI3: https://github.com/ConesaLab/SQANTI3
- SUPPA2: https://github.com/comprna/SUPPA
- LoRDEC: https://github.com/blinlnb/LoRDEC
- Cupcake: https://github.com/Magdoll/cDNA_Cupcake
