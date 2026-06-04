# 基因家族分析 -- 分析笔记本

## 分析概览

本笔记本提供植物基因家族鉴定与进化分析的完整流程，涵盖HMM搜索、BLAST验证、多序列比对、系统发育树构建、Motif分析、基因结构分析、共线性分析、Ka/Ks分析和顺式元件分析。

---

## 1. 基因家族成员鉴定

### 1.1 HMM搜索 (首选方法)

当已知Pfam结构域时，HMM搜索是最敏感的方法：

```bash
# 从Pfam下载HMM profile
wget https://pfam.xfam.org/family/PF00096/hmm -o PF00096.hmm

# 构建HMM数据库
hmmbuild PF00096.hmm <(wget -O- https://pfam.xfam.org/family/PF00096/alignment)

# 搜索基因组蛋白序列
hmmsearch --cpu 8 \
  --tblout hmm_results.tbl \
  -E 1e-5 \
  --domtblout domain_results.tbl \
  PF00096.hmm \
  genome_proteins.fasta

# 提取命中序列
awk '$5 !~ /^#/{print $1}' hmm_results.tbl | sort -u > gene_ids.txt
seqkit grep -f gene_ids.txt genome_proteins.fasta > gene_family_proteins.fasta
```

**HMM搜索参数说明**：
- `-E 1e-5`: E值阈值，控制假阳性
- `--domtblout`: 输出结构域级别结果
- `--cpu`: 线程数

### 1.2 BLAST验证 (补充方法)

```bash
# 使用已知蛋白序列做BLASTP
blastp -query query_protein.fasta \
  -db genome_proteins.fasta \
  -out blast_results.out \
  -evalue 1e-10 \
  -outfmt 6 \
  -max_target_seqs 1000

# 过滤覆盖度
awk '$3 > 50 && $4 > 0.5' blast_results.out > filtered_blast.out
```

### 1.3 迭代搜索策略

对于无Pfam域但有候选序列的情况：

```bash
# Step 1: BLASTP初步搜索
blastp -query seed_sequences.fasta -db genome_proteins.fasta -out initial.hits

# Step 2: 构建HMM profile
mafft --auto initial_hits.fasta > aligned.fasta
hmmbuild custom_family.hmm aligned.fasta

# Step 3: HMM搜索扩大搜索
hmmsearch --tblout expanded.hits custom_family.hmm genome_proteins.fasta
```

---

## 2. 多序列比对

### 2.1 MAFFT (推荐)

```bash
# 小规模 (<50序列): L-INS-i最精确
mafft --localpair --maxiterate 1000 --thread 8 \
  gene_family_proteins.fasta > aligned.fasta

# 大规模 (≥50序列): FFT-NS-2快速
mafft --fft --maxiterate 2 --thread 8 \
  gene_family_proteins.fasta > aligned.fasta

# 超大规模 (>500序列): FFT-NS-1最快
mafft --retree 1 --thread 8 \
  gene_family_proteins.fasta > aligned.fasta
```

### 2.2 比对质量检查

```bash
# 使用trimAl清理比对
trimal -in aligned.fasta -out aligned_trimmed.fasta \
  -automated1

# 检查比对质量
aliscore -s aligned_trimmed.fasta
```

**植物特异性考量**：
- 植物基因家族常有重复结构域，确保比对覆盖完整结构域
- 多倍体物种的旁系同源基因相似度高，比对质量通常较好
- 注意识别和排除片段化序列

---

## 3. 系统发育树构建

### 3.1 IQ-TREE2 (推荐)

```bash
# 自动模型选择 + 超快bootstrap
iqtree2 -s aligned_trimmed.fasta \
  -m MFP \
  -B 1000 --alrt 1000 \
  -T 8 \
  --prefix gene_family_tree

# 输出文件：
# .treefile: 最大似然树
# .iqtree: 完整报告含模型选择
# .ufboot: bootstrap树
```

**IQ-TREE参数说明**：
- `-m MFP`: ModelFinder Plus自动选择最佳模型
- `-B 1000`: 1000次超快bootstrap
- `--alrt 1000`: 1000次SH-aLRT检验
- 支持率≥95%/80%分别表示强/中度支持

### 3.2 SNP数据处理

```bash
# SNP数据必须使用ASC校正
iqtree2 -s snp_alignment.fasta \
  -m MFP+ASC \
  --asc-corr lewis \
  -B 1000
```

### 3.3 树的可视化

```bash
# 使用iTOL在线可视化
# 上传 .treefile 到 https://itol.embl.de/

# 或使用FigTree本地可视化
figtree gene_family_tree.treefile
```

---

## 4. Motif分析

### 4.1 MEME分析

```bash
# 植物基因家族推荐参数
meme gene_family_proteins.fasta \
  -mod anr \
  -nmotifs 10 \
  -minw 6 \
  -maxw 50 \
  -oc meme_output \
  -protein

# 输出文件：
# meme.html: 结果网页
# meme.txt: 文本格式结果
# meme.xml: XML格式结果
```

**MEME参数说明**：
- `-mod anr`: 任意次数重复，适合植物重复结构域
- `-mod zoops`: 每序列最多一次，适合动物
- `-nmotifs`: 搜索的motif数量
- `-minw/-maxw`: motif宽度范围

### 4.2 Motif注释

```bash
# 使用TOMTOM比对已知数据库
tomtom -oc tomtom_output \
  meme_output/meme.txt \
  $MEME_DB/motif_databases/JASPAR_CORE_2022.meme
```

### 4.3 可视化

```bash
# 使用TBtools进行motif可视化
# 或使用R包ggseqlogo
```

---

## 5. 基因结构分析

### 5.1 提取基因结构

```bash
# 使用gffread提取CDS和UTR信息
gffread -g genome.fasta \
  -x gene_cds.fasta \
  -y gene_proteins.fasta \
  -T \
  annotation.gff3

# 提取特定基因结构
gffread annotation.gff3 -T | grep -f gene_ids.txt > gene_structures.gff3
```

### 5.2 GSDS可视化

```bash
# 准备GSDS输入文件
# 格式: GeneID \t ExonStarts \t ExonEnds \t IntronInfo

# 上传到GSDS2.0: http://gsds.gao-lab.org/
```

---

## 6. 共线性分析

### 6.1 MCScanX (推荐)

```bash
# Step 1: 准备输入文件
# 需要BLASTP结果和基因位置文件

# 基因位置文件 (gene_location.txt)
# 格式: Chromosome \t GeneID \t Start \t End
awk '$3=="gene"{print $1"\t"$9"\t"$4"\t"$5}' annotation.gff3 | \
  sed 's/ID=//;s/;//' > gene_location.txt

# Step 2: 运行BLASTP (物种内)
blastp -query proteins.fasta -db proteins.fasta -out blast.out -evalue 1e-5 -outfmt 6

# Step 3: 运行MCScanX
MCScanX gene_location
```

**MCScanX输出文件**：
- `.collinearity`: 共线性区块
- `.tandem`: 串联重复基因
- `.wgd`: 全基因组复制基因对

### 6.2 结果解读

```bash
# 统计复制类型
grep "WGD" gene_location.collinearity | wc -l
grep "Tandem" gene_location.tandem | wc -l
grep "Proximal" gene_location.collinearity | wc -l
```

---

## 7. Ka/Ks分析

### 7.1 准备基因对

```bash
# 从共线性结果提取基因对
awk 'NR>1{print $2"\t"$3}' gene_location.collinearity > gene_pairs.txt
```

### 7.2 计算Ka/Ks

```bash
# 方法1: 使用KaKs_Calculator (批量)
# 准备AXT格式文件
KaKs_Calculator -i gene_pairs.axt -o kaks_results.txt -m MLWL

# 方法2: 使用PAML yn00 (精确但慢)
# 准备codeml配置文件
codeml codeml.ctl
```

### 7.3 结果解读

| Ka/Ks值 | 含义 |
|---------|------|
| < 1 | 负选择(纯化选择)，功能保守 |
| ≈ 1 | 中性进化，功能可能正在变化 |
| > 1 | 正选择，功能快速进化 |

**植物特异性考量**：
- 近期WGD的基因对Ka/Ks通常较低(功能冗余)
- 分化时间久的重复基因可能已功能分化

---

## 8. 顺式元件分析

### 8.1 提取启动子序列

```bash
# 提取起始密码子上游2000bp
bedtools flank -g genome.fasta.fai \
  -b 2000 \
  -i gene_regions.bed > promoters.bed

bedtools getfasta -fi genome.fasta \
  -bed promoters.bed \
  -fo promoter_sequences.fasta
```

### 8.2 PlantCARE分析

```bash
# PlantCARE是在线工具
# 上传FASTA文件到: http://bioinformatics.psb.ugent.be/webtools/plantcare/html/

# 或使用本地PlantCARE脚本
perl plantcare.pl promoter_sequences.fasta output_dir
```

**PlantCARE参数建议**：
- 启动子长度: 2000bp (上游)
- 最小score: 5
- 关注元件类型: ABRE (ABA响应)、G-box (光响应)、MBS (干旱响应)等

### 8.3 结果统计

```bash
# 统计各类型元件数量
grep -c "ABRE" plantcare_results.txt
grep -c "G-box" plantcare_results.txt
```

---

## 9. 表达模式分析

### 9.1 从RNA-seq数据提取表达量

```bash
# 使用featureCounts定量
featureCounts -a annotation.gff3 \
  -o expression_counts.txt \
  -T 8 \
  -p \
  sample1.bam sample2.bam ...

# 提取基因家族成员表达量
grep -f gene_ids.txt expression_counts.txt > family_expression.txt
```

### 9.2 表达热图

```r
library(pheatmap)
expr_matrix <- read.table("family_expression.txt", header=TRUE, row.names=1)
pheatmap(log2(expr_matrix + 1), 
         scale="row",
         clustering_distance_rows="correlation")
```

---

## 常见问题

### Q: HMM搜索结果太多怎么办？
A: 提高E值阈值(如1e-10)，或增加覆盖度要求(>70%)。检查是否有非特异性结构域污染。

### Q: 系统发育树支持率低怎么办？
A: 检查比对质量，可能需要去除高度分歧区域。考虑使用更多的bootstrap重复(2000+)。

### Q: Ka/Ks值异常高怎么办？
A: 检查比对质量，特别是密码子比对。排除比对错误导致的假阳性。

### Q: PlantCARE结果为空？
A: 检查启动子序列是否正确提取。确保序列方向正确(5'到3')。

---

## 植物基因家族分析特有考量

### 多倍体物种

- 同源基因组间需要区分
- 可能需要亚基因组特异性分析
- 复制基因功能冗余是常见现象

### WGD历史

- 植物基因组常有多次WGD历史
- 近期WGD的基因对Ka/Ks低，表达相似
- 古老WGD的基因对可能已功能分化

### 组织特异性

- 植物基因家族常有组织特异性成员
- 表达分析对理解功能很重要
- 考虑取样时间(昼夜节律影响)

---

## 参考

- HMMER: http://hmmer.org/
- MAFFT: https://mafft.cbrc.jp/alignment/software/
- IQ-TREE: http://www.iqtree.org/
- MEME Suite: https://meme-suite.org/
- MCScanX: https://github.com/wyp1125/MCScanX
- PlantCARE: http://bioinformatics.psb.ugent.be/webtools/plantcare/html/
- KaKs_Calculator: https://sourceforge.net/projects/kakscalculator2/
