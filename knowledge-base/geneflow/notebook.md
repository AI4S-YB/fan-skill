# 基因流分析 -- 分析笔记本

## 分析概览

本笔记本提供植物基因流和渐渗分析的完整流程，涵盖群体结构分析、ABBA-BABA检验、TreeMix分析和渐渗方向推断。

---

## 1. 数据准备

### 1.1 VCF文件过滤

```bash
# 使用VCFtools过滤
vcftools --vcf raw.vcf \
  --maf 0.05 \
  --max-missing 0.9 \
  --minDP 3 \
  --maxDP 100 \
  --remove-indels \
  --recode --recode-INFO-all \
  --out filtered

# 使用PLINK转换为bed格式
plink --vcf filtered.recode.vcf \
  --make-bed \
  --out geneflow_data

# 转换为Treemix格式
plink --bfile geneflow_data \
  --freq \
  --within cluster_file.txt \
  --out freq_output
```

### 1.2 数据质控

```bash
# 检查样本缺失率
plink --bfile geneflow_data --missing --out missing_report

# 检查等位基因频率分布
plink --bfile geneflow_data --freq --out freq_report

# 过滤连锁不平衡位点
plink --bfile geneflow_data \
  --indep-pairwise 50 5 0.2 \
  --out pruned
plink --bfile geneflow_data \
  --extract pruned.prune.in \
  --make-bed \
  --out geneflow_pruned
```

---

## 2. 群体结构分析

### 2.1 PCA分析

```bash
# 使用PLINK进行PCA
plink --bfile geneflow_pruned \
  --pca 10 \
  --out pca_results

# R可视化
R
library(ggplot2)
pca <- read.table("pca_results.eigenvec")
ggplot(pca, aes(V3, V4, color=species)) +
  geom_point(size=3) +
  theme_minimal()
```

### 2.2 ADMIXTURE分析

```bash
# 运行ADMIXTURE（K=2到10）
for K in {2..10}; do
  admixture --cv geneflow_pruned.bed $K | tee log${K}.out
done

# 选择最佳K值
grep "CV error" log*.out

# 可视化结果
R
library(pophelper)
qfiles <- list.files(pattern="\\.Q$")
readQ(qfiles)
```

### 2.3 群体结构解读

- PCA第一主成分通常反映最大的群体分化
- ADMIXTURE K值选择：CV error最小的K
- 注意区分群体结构和渐渗信号

---

## 3. 物种树推断

### 3.1 ASTRAL物种树

```bash
# 准备基因树
# 方法1：滑动窗口建树
for i in {1..100}; do
  start=$((i*10000))
  end=$((start+10000))
  vcftools --vcf filtered.vcf --chr chr1 --from-bp $start --to-bp $end --recode --out window_$i
  # 使用IQ-TREE建树
  iqtree2 -s window_$i.recode.vcf -m GTR+ASC -T 1 -pre window_$i
done

# 合并基因树
cat *.treefile > all_gene_trees.tre

# 运行ASTRAL
java -jar ASTRAL.jar -i all_gene_trees.tre -o species_tree.tre
```

### 3.2 SNaQ网络分析

```julia
# 使用Julia运行PhyloNetworks
using PhyloNetworks

# 读取基因树
trees = readMultiTopology("gene_trees.tre")

# 运行SNaQ
net0 = readTopology("species_tree.tre")
net = snaq!(net0, trees, hmax=3, runs=10)
writeTopology(net, "network.nwk")
```

---

## 4. ABBA-BABA检验

### 4.1 Dsuite分析

```bash
# 准备样本分组文件
# 格式: sample_id population_id
cat samples.txt

# 运行Dsuite
Dsuite Dtrios filtered.vcf samples.txt -t species_tree.nwk -o dsuite_output

# 运行f-branch分析
Dsuite Fbranch dsuite_output_tree.txt species_tree.nwk > fbranch.txt

# 可视化f-branch结果
R
library(ggplot2)
fbranch <- read.table("fbranch.txt")
# 绘制热图
```

### 4.2 D统计量解读

| D值 | 解释 |
|-----|------|
| D ≈ 0 | 无基因流 |
| D > 0 | P1和P3间有基因流（ABBA > BABA） |
| D < 0 | P2和P3间有基因流（BABA > ABBA） |
| |Z| > 3 | 统计显著 |

### 4.3 f_d滑动窗口分析

```bash
# 滑动窗口f_d计算
Dsuite Dinvestigate filtered.vcf samples.txt -w 50000,10000 -o fd_windows

# R可视化
R
fd <- read.table("fd_windows.txt")
ggplot(fd, aes(position, f_d)) +
  geom_point(alpha=0.5) +
  geom_smooth() +
  facet_wrap(~chromosome)
```

---

## 5. TreeMix分析

### 5.1 数据准备

```bash
# 使用plink2treemix转换格式
# 需要先计算等位基因频率
plink --bfile geneflow_pruned \
  --freq --within pop_clusters.txt \
  --out treemix_freq

# 转换格式
python plink2treemix.py treemix_freq.frq.strat treemix_input.gz
```

### 5.2 运行TreeMix

```bash
# 运行不同迁移边数
for m in {0..5}; do
  treemix -i treemix_input.gz \
    -m $m \
    -k 1000 \
    -bootstrap \
    -o treemix_m$m
done

# 可视化
R
library(OPTtree)
source("plotting_funcs.R")
plot_tree("treemix_m3")
```

### 5.3 结果解读

- 迁移边数量选择：残差方差最小
- 迁移边方向：箭头指向基因流方向
- 权重：边粗细表示基因流强度

---

## 6. 局部渐渗分析

### 6.1 滑动窗口分析

```bash
# 使用VCFtools计算滑动窗口统计
for chr in chr1 chr2 chr3; do
  vcftools --vcf filtered.vcf \
    --chr $chr \
    --window-pi 50000 \
    --out pi_$chr
done

# 使用自定义脚本计算f_d窗口
python calculate_fd_windows.py \
  --vcf filtered.vcf \
  --samples samples.txt \
  --window 50000 \
  --step 10000 \
  --output fd_windows.txt
```

### 6.2 渐渗区域识别

```bash
# 识别显著渐渗区域
# 方法：f_d > 0.1且Z > 3

# R过滤
R
fd <- read.table("fd_windows.txt")
colnames(fd) <- c("chr", "start", "end", "fd", "zscore")
sig_regions <- fd[fd$fd > 0.1 & abs(fd$zscore) > 3, ]
write.table(sig_regions, "significant_introgression.bed", sep="\t", row.names=F, col.names=F, quote=F)
```

### 6.3 渐渗基因功能分析

```bash
# 提取渐渗区域的基因
bedtools intersect -a genes.gff -b significant_introgression.bed > introgressed_genes.gff

# GO富集分析
# 使用eggNOG-mapper或在线工具
```

---

## 7. 植物特异性分析

### 7.1 自交物种处理

```bash
# 自交物种群体结构可能复杂
# 使用更严格的LD过滤
plink --bfile geneflow_data \
  --indep-pairwise 10 1 0.1 \
  --out strict_pruned
```

### 7.2 多倍体考虑

```bash
# 多倍体SNP calling需要特殊处理
# 使用GATK HaplotypeCaller的ploidy参数
gatk HaplotypeCaller --sample-ploidy 4 -R ref.fa -I sample.bam -o sample.vcf
```

### 7.3 杂交区域识别

```bash
# 对于已知杂交个体，识别杂交区域
# 使用局部祖先推断
rfmix -f query.vcf -r reference.vcf -m mapping.txt -o rfmix_output
```

---

## 8. 结果整合

### 8.1 多方法整合

```R
# 整合D统计量、f_d、TreeMix结果
# 创建综合证据表

results <- data.frame(
  method = c("D-statistic", "f_d", "TreeMix", "f-branch"),
  evidence = c("significant", "significant", "migration_edge", "significant"),
  direction = c("P1->P3", "P1->P3", "P1->P3", "P1->P3")
)

# 综合判断基因流
```

---

## 常见问题

### Q: D统计量不显著？
A: 检查样本量是否足够，外群选择是否合适。尝试其他统计量如f4-ratio。

### Q: 不同方法结果不一致？
A: 检查方法假设是否满足。综合多种证据判断。

### Q: f_d值异常高？
A: 检查群体结构是否被正确建模。可能需要更复杂的模型。

### Q: TreeMix迁移边太多？
A: 选择残差方差下降拐点对应的迁移边数。

---

## 参考

- Dsuite: https://github.com/millanek/Dsuite
- ADMIXTURE: https://dalexander.github.io/admixture/
- TreeMix: https://bitbucket.org/nygcresearch/treemix/wiki/Home
- ASTRAL: https://github.com/smirarab/ASTRAL
- PhyloNetworks: https://github.com/crsl4/PhyloNetworks.jl
- RFMix: https://github.com/slowkoni/rfmix
