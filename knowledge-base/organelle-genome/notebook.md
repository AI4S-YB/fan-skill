# 细胞器基因组分析 -- 分析笔记本

## 分析概览

本笔记本提供植物细胞器基因组（叶绿体和线粒体）分析的完整流程，涵盖组装、注释、结构分析和比较分析。

---

## 1. 细胞器基因组组装

### 1.1 GetOrganelle (推荐)

```bash
# 从全基因组数据组装叶绿体和线粒体
get_organelle_from_reads.py \
  -1 reads_R1.fastq.gz \
  -2 reads_R2.fastq.gz \
  -t 16 \
  -o output_dir \
  -F plant_cp  # 叶绿体
  # -F plant_mt  # 线粒体

# 输出：
# output_dir/*.fasta - 组装结果
# output_dir/*.gfa - 组装图
```

### 1.2 NOVOPlasty (叶绿体专用)

```bash
# 准备配置文件
cat > novo_config.txt << EOF
Project:
-----------------------
Project name          = chloroplast_assembly
Type                  = chloro
Genome Range          = 120000-180000
K-mer                 = 39
Max memory            = 16
Extended log          = 0
Save assembled reads  = no
Seed Input            = seed.fasta
Reference sequence    = reference.fasta
Dataset 1:
-----------------------
Read Length           = 150
Insert size           = 350
Platform              = illumina
1                     = reads_R1.fastq.gz
2                     = reads_R2.fastq.gz
EOF

# 运行NOVOPlasty
NOVOPlasty.pl -c novo_config.txt
```

### 1.3 MitoFinder (线粒体专用)

```bash
# 运行MitoFinder
mitofinder -j mito_job \
  -r reference.fasta \
  -1 reads_R1.fastq.gz \
  -2 reads_R2.fastq.gz \
  -t 16 \
  -o mitofinder_output
```

---

## 2. 基因注释

### 2.1 叶绿体注释 (GeSeq)

```bash
# GeSeq是在线工具
# 访问: https://chlorobox.mpimp-golm.mpg.de/geseq.html

# 或使用本地工具CPCAS
python cpcas.py \
  -i chloroplast.fasta \
  -o annotation.gff \
  -r reference.gb
```

### 2.2 线粒体注释

```bash
# 使用MITOS
# 在线: http://mitos.bioinf.uni-leipzig.de/

# 或使用MitoFinder的注释功能
mitofinder -j anno_job \
  -a assembly.fasta \
  -r reference.fasta \
  -o annotation_output
```

### 2.3 功能基因分类

**叶绿体基因组典型基因**：
- 光合作用基因: rbcL, psa, psb, pet, atp
- 遗传系统基因: rRNA, tRNA, rpl, rps
- 其他基因: matK, ccsA, cemA, clpP

**线粒体基因组典型基因**：
- 呼吸链复合物: nad, cob, cox, atp
- 核糖体蛋白: rpl, rps
- RNA基因: rRNA, tRNA

---

## 3. 结构分析

### 3.1 叶绿体IR边界分析

```bash
# 识别IR区域
# 方法1：通过基因排列识别
# rRNA基因通常位于IR区域

# 方法2：使用IRfinder
IRfinder -i chloroplast.fasta -o ir_results.txt

# 可视化IR边界
Rscript plot_ir_boundary.R ir_results.txt
```

### 3.2 SSR分析

```bash
# 使用MISA进行SSR分析
misa.pl chloroplast.fasta

# 输出：
# chloroplast.fasta.statistics - SSR统计
# chloroplast.fasta.ssr - SSR位点列表
```

### 3.3 重复序列分析

```bash
# 使用REPuter分析重复
reputer -seq chloroplast.fasta \
  -f 1 -l 30 -h 3 \
  -o repeats.txt

# 或使用在线工具
# https://www.bioinformatics.nl/reputer/
```

---

## 4. 比较分析

### 4.1 多序列比对

```bash
# 使用Mauve
progressiveMauve --output=alignment.xmfa \
  species1.fasta species2.fasta species3.fasta

# 转换比对结果
java -cp mauve.jar org.gel.mauve.BackboneViewer \
  -output alignment.bb alignment.xmfa
```

### 4.2 基因组可视化比对

```bash
# 使用mVista (在线)
# http://genome.lbl.gov/vista/mvista/submit.shtml

# 或使用本地工具
# Easyfig
easyfig -f genome1.gb -f genome2.gb -o comparison.svg
```

### 4.3 IR边界比较

```bash
# 使用IRscope
# 在线: https://irscope.shinyapps.io/irscope/

# 准备输入文件
# GenBank格式的基因组文件
```

---

## 5. 系统发育分析

### 5.1 单基因建树

```bash
# 提取共享单拷贝基因
python extract_shared_genes.py *.gb

# 每个基因分别建树
for gene in rbcL matK ndhF; do
  mafft --auto ${gene}.fasta > ${gene}_aligned.fasta
  iqtree2 -s ${gene}_aligned.fasta -m MFP -B 1000
done
```

### 5.2 串联基因建树

```bash
# 连接多个基因
cat *_aligned.fasta > concatenated.fasta

# 建树
iqtree2 -s concatenated.fasta -m MFP -B 1000 -T 8
```

### 5.3 全基因组建树

```bash
# 全基因组比对
mafft --auto --thread 8 all_genomes.fasta > aligned.fasta

# 建树
iqtree2 -s aligned.fasta -m MFP -B 1000
```

---

## 6. RNA编辑位点预测

### 6.1 叶绿体RNA编辑预测 (PREP-Cp)

```bash
# PREP-Cp是PREPACT在线工具的一部分
# 访问: http://prep.unl.edu/

# 输入文件准备：
# 1. 叶绿体蛋白编码基因的FASTA序列
# 2. 选择最近的参考物种

# 本地替代方式：使用PREPACT API
wget -O prep_result.txt \
  "http://prep.unl.edu/cp/analysis?cutoff=0.5&species=nearest"

# 解析结果
python parse_prep_result.py prep_result.txt
```

### 6.2 线粒体RNA编辑预测 (PREP-Mt)

```bash
# PREP-Mt用于线粒体基因组RNA编辑位点预测
# 访问: http://prep.unl.edu/

# 准备线粒体蛋白编码基因序列
cat mitochondrial_cds.fasta

# 提交PREP-Mt分析
# 参数：cutoff_value=0.5
# 输出：C-to-U编辑位点列表及位置
```

**RNA编辑位点解读**：
- C-to-U编辑是植物细胞器最常见的RNA编辑类型
- 编辑通常发生在密码子的第1或第2位
- 编辑可改变氨基酸，影响蛋白功能
- 不编辑可能导致蛋白编码基因注释错误（起始/终止密码子偏移）

---

## 7. 选择压力分析

### 7.1 成对Ka/Ks计算 (KaKs_Calculator + ParaAT)

```bash
# Step 1: 准备直系同源基因对
# species1_cds.fa 和 species2_cds.fa

# Step 2: 使用ParaAT批量处理
ParaAT.pl -h species1_cds.fa \
  -h species2_cds.fa \
  -n cds_pairs.fa \
  -a protein_pairs.fa \
  -p proc \
  -m mafft \
  -f axt \
  -k KaKs_Calculator \
  -o output_dir/

# 或手动流程
# 2a. 蛋白比对
mafft --auto protein_pair.fa > protein_aligned.fa

# 2b. 回译为CDS比对
pal2nal.pl protein_aligned.fa cds_pair.fa \
  -output paml > codon_aligned.fa

# 2c. Ka/Ks计算
KaKs_Calculator -i codon_aligned.fa -o kaks_result.txt \
  -m YN  # Yang-Nielsen方法

# 查看结果
cat kaks_result.txt
# 输出列：Sequence, Ka, Ks, Ka/Ks, P-Value
```

### 7.2 基因水平选择压力分析 (PAML codeml)

```bash
# 准备密码子比对文件 (phylip格式)
cat > codeml.ctl << EOF
seqfile = codon_aligned.phy
treefile = species_tree.nwk
outfile = mlc
noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
clock = 0
model = 0
NSsites = 0
icode = 0
fix_omega = 0
omega = 0.4
EOF

# 运行codeml
codeml codeml.ctl

# 提取结果
grep "omega" mlc | head -20

# 批量处理多个基因
for gene in genes/*.phy; do
  sed "s|INPUT|$gene|" codeml_template.ctl > codeml.ctl
  codeml codeml.ctl
  grep "omega" mlc >> all_genes_omega.txt
done
```

**Ka/Ks解读**：
- Ka/Ks < 1：纯化选择（负选择）
- Ka/Ks = 1：中性进化
- Ka/Ks > 1：正选择（适应性进化）
- 叶绿体基因通常Ka/Ks远小于1（高度保守）

---

## 8. 可视化

### 8.1 OGDRAW环形图

```bash
# OGDRAW是在线工具
# https://chlorobox.mpimp-golm.mpg.de/OGDraw.html

# 上传GenBank文件即可生成环形图
```

### 8.2 自定义绘图

```R
# 使用R绘制基因组图
library(ggplot2)
library(gggenes)

# 读取基因注释
genes <- read.table("genes.txt", header=TRUE)

# 绘制基因图
ggplot(genes, aes(xmin=start, xmax=end, y=species, fill=category)) +
  geom_gene_arrow() +
  facet_wrap(~category, scales="free", ncol=1) +
  theme_minimal()
```

---

## 9. 植物特异性分析

### 9.1 叶绿体基因组特点

```
典型结构：
- LSC (Large Single Copy): 80-90 kb
- SSC (Small Single Copy): 15-20 kb
- IR (Inverted Repeat): 20-30 kb × 2

总大小：120-160 kb
基因数：约110-130个
```

### 9.2 线粒体基因组特点

```
特点：
- 大小变异极大：200 kb - 2.7 Mb
- 结构复杂：多染色体、重组
- RNA编辑普遍
- 基因丢失和转移常见

常见植物线粒体大小：
- 拟南芥：367 kb
- 玉米：570 kb
- 南瓜：2.7 Mb
```

### 9.3 母系遗传验证

```bash
# 比较亲本和后代的细胞器单倍型
# 确认母系遗传模式

vcftools --vcf parent_offspring.vcf \
  --indv mother --indv father --indv offspring \
  --hapwindow 1000 \
  --out haplotype_comparison
```

---

## 10. 结果整合

### 10.1 标准输出

```bash
# 组装结果
ls *.fasta

# 注释结果
ls *.gb *.gff

# 比较分析
ls *.xmfa *.svg

# 系统发育
ls *.treefile

# 统计信息
cat genome_statistics.txt
```

---

## 常见问题

### Q: 组装不完整？
A: 检查测序深度和覆盖度。叶绿体需要中等深度，线粒体可能需要更高深度。

### Q: IR边界不清晰？
A: 手动检查rRNA基因位置。使用参考基因组辅助判断。

### Q: 线粒体基因组碎片化？
A: 植物线粒体基因组复杂，可能存在多染色体。检查组装图。

### Q: 注释缺失基因？
A: 检查阈值设置。某些基因可能因变异而未被识别。

---

## 参考

- GetOrganelle: https://github.com/Kinggerm/GetOrganelle
- NOVOPlasty: https://github.com/ndierckx/NOVOPlasty
- MitoFinder: https://github.com/RemiLacroix/Laboratoire_Evol
- GeSeq: https://chlorobox.mpimp-golm.mpg.de/geseq.html
- OGDRAW: https://chlorobox.mpimp-golm.mpg.de/OGDraw.html
- IRscope: https://irscope.shinyapps.io/irscope/
