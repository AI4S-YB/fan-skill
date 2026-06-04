# 物种分析速查表

分析前先查此表，了解你的物种特性。这些信息会影响方法选择、参数设置和结果解读。

---

## 水稻 (Oryza sativa)

| 属性 | 值 |
|------|-----|
| 倍性 | 二倍体 (2n=24) |
| 繁殖方式 | 自交 |
| 基因组大小 | ~430 Mb |
| 常用标记密度 | 3K array / 44K array / GBS / WGS |
| GWAS 建议 | K 矩阵不可省略；籼粳需分开分析 |
| LD 衰减 | 籼稻 ~100kb, 粳稻 ~200kb |
| 着丝粒注意 | 着丝粒区域 LD 异常高，显著位点需排除着丝粒 |
| 特殊资源 | 水稻3K项目、RiceGFD、RAP-DB、RiceVarMap |
| 参考基因组 | IRGSP-1.0 (Ensembl Plants) / Os-Nipponbare-Reference-IRGSP-1.0 |
| 注释版本 | RAP-DB / MSU7 |

**GWAS 特别注意**：
- 粳稻群体结构通常强于籼稻
- 如果数据混合了籼粳亚群，必须先分亚群分析
- 候选基因优先查 RAP-DB 和 RiceGFD

---

## 玉米 (Zea mays)

| 属性 | 值 |
|------|-----|
| 倍性 | 二倍体 (2n=20) |
| 繁殖方式 | 异交 |
| 基因组大小 | ~2.3 Gb (85% 重复序列) |
| 常用标记密度 | 50K SNP chip / GBS / WGS (>1M SNPs) |
| GWAS 建议 | LD 衰减极快 (<1kb)，高密度标记才能精细定位 |
| LD 衰减 | <1kb (温带), ~2kb (热带) |
| 着丝粒注意 | 玉米着丝粒区域大、重复序列多 |
| 特殊资源 | MaizeGDB、HapMap3、NAM founder lines、GRAMENE |
| 参考基因组 | Zm-B73-REFERENCE-NAM-5.0 (Ensembl Plants) |
| 注释版本 | B73 RefGen_v5 |

**GWAS 特别注意**：
- NAM 群体 (巢式关联群体) 功效远高于自然群体
- 热带/温带种质遗传背景差异大，必须作为协变量
- 玉米泛基因组复杂，PAV (存在/缺失变异) 可能与 SNP 同样重要

---

## 小麦 (Triticum aestivum)

| 属性 | 值 |
|------|-----|
| 倍性 | 异源六倍体 (2n=6x=42, AABBDD) |
| 繁殖方式 | 自交 |
| 基因组大小 | ~17 Gb |
| 常用标记密度 | 90K iSelect / 660K / WGS |
| GWAS 建议 | A/B/D 三个亚基因组分别分析；同源基因可能互补 |
| LD 衰减 | A/B 基因组 ~3-5kb, D 基因组 ~10kb |
| 着丝粒注意 | 三个亚基因组着丝粒位置不同 |
| 特殊资源 | URGI、WheatOmics、T3/Wheat、10+ Wheat Genomes |
| 参考基因组 | IWGSC RefSeq v2.1 (Ensembl Plants) |
| 注释版本 | IWGSC Annotation v2.1 |

**GWAS 特别注意**：
- 六倍体意味着每个位点有 3 个同源拷贝，功能变异可能被其他亚基因组的野生型拷贝掩盖
- 需要区分 A、B、D 亚基因组的 SNP
- 标记命名通常包含染色体信息 (如 IWB12345 在 1A)
- 着丝粒区域在六倍体中占比更大

---

## 大豆 (Glycine max)

| 属性 | 值 |
|------|-----|
| 倍性 | 二倍体 (2n=40)，古多倍体 |
| 繁殖方式 | 自交 |
| 基因组大小 | ~1.1 Gb |
| 常用标记密度 | 50K chip / GBS / WGS |
| GWAS 建议 | K 矩阵 + PC 协变量一起用；古多倍体重复区域需过滤 |
| LD 衰减 | 栽培大豆 ~100-200kb, 野生大豆 ~50kb |
| 特殊资源 | SoyBase、SoyKB、SoyFGB |
| 参考基因组 | Wm82.a4.v1 (Phytozome) / Gmax_ZH13 |
| 注释版本 | Wm82.a4.v1 |

**GWAS 特别注意**：
- 古多倍体历史导致约 75% 的基因是多拷贝的
- 栽培大豆的遗传多样性远低于野生大豆（驯化瓶颈）
- 南方/北方生态型可能构成严重群体分层

---

## 棉花 (Gossypium hirsutum)

| 属性 | 值 |
|------|-----|
| 倍性 | 异源四倍体 (2n=4x=52, AADD) |
| 繁殖方式 | 自交/常异交 |
| 基因组大小 | ~2.5 Gb |
| 常用标记密度 | 63K chip / GBS / WGS |
| GWAS 建议 | At/Dt 亚基因组分别分析 |
| 特殊资源 | CottonGen、CottonFGD |
| 参考基因组 | TM-1 (CottonGen) / ZJU Assembly |

---

## 油菜 (Brassica napus)

| 属性 | 值 |
|------|-----|
| 倍性 | 异源四倍体 (2n=4x=38, AACC) |
| 繁殖方式 | 自交 |
| 基因组大小 | ~1.2 Gb |
| GWAS 建议 | A/C 亚基因组分别分析 |
| 特殊资源 | BnPIR、BRAD |
| 参考基因组 | Darmor-bzh (Ensembl Plants) |

---

## 拟南芥 (Arabidopsis thaliana)

| 属性 | 值 |
|------|-----|
| 倍性 | 二倍体 (2n=10) |
| 繁殖方式 | 自交 |
| 基因组大小 | ~135 Mb |
| 常用标记密度 | 250K chip / WGS |
| GWAS 建议 | 全球群体结构极强，必须严格控制群体分层 |
| 特殊资源 | TAIR、AraGWAS、1001 Genomes |
| 参考基因组 | TAIR10 / Araport11 |

---

## 番茄 (Solanum lycopersicum)

| 属性 | 值 |
|------|-----|
| 倍性 | 二倍体 (2n=24) |
| 繁殖方式 | 自交 |
| 基因组大小 | ~900 Mb |
| GWAS 建议 | 栽培番茄多样性低，考虑包含野生近缘种 |
| 特殊资源 | Sol Genomics Network、Tomato Genome Hub |
| 参考基因组 | SL4.0 / ITAG4.0 |

---

## 马铃薯 (Solanum tuberosum)

| 属性 | 值 |
|------|-----|
| 倍性 | 同源四倍体 (2n=4x=48) |
| 繁殖方式 | 无性繁殖（块茎） |
| 基因组大小 | ~840 Mb |
| GWAS 建议 | 四倍体 GWAS 需专用方法 (GWASpoly)；杂合度高 |
| 特殊资源 | Spud DB、Potato Genome Hub |
| 参考基因组 | DM v6.1 |

**GWAS 特别注意**：
- 同源四倍体的剂量效应：AAAA/AAAa/AAaa/Aaaa/aaaa
- 不能用标准二倍体 GWAS 工具直接分析
- GWASpoly 是目前最好的四倍体 GWAS 工具

---

## 大麦 (Hordeum vulgare)

| 属性 | 值 |
|------|-----|
| 倍性 | 二倍体 (2n=14) |
| 繁殖方式 | 自交 |
| 基因组大小 | ~5.3 Gb |
| 特殊资源 | BAR (Bio-Analytic Resource), eFP Browser |

---

## 高粱 (Sorghum bicolor)

| 属性 | 值 |
|------|-----|
| 倍性 | 二倍体 (2n=20) |
| 繁殖方式 | 自交/常异交 |
| 基因组大小 | ~730 Mb |
| 特殊资源 | SorghumBase |

---

## 甘蔗 (Saccharum spp.)

| 属性 | 值 |
|------|-----|
| 倍性 | 同源多倍体 (2n=100-130) |
| 繁殖方式 | 无性繁殖 |
| 基因组大小 | ~10 Gb |
| GWAS 建议 | 极其复杂的多倍体，GWAS 方法仍在发展中 |
| 特殊资源 | Sugarcane Genome Hub |

---

## 共性规则

### 自交作物 (水稻、大豆、小麦、大麦、番茄)
- GWAS 中 K 矩阵通常不可省略
- 群体结构 (Q 矩阵) 和亲缘关系 (K 矩阵) 一起使用
- LD 衰减较慢 → 显著区域大，精细定位困难

### 异交作物 (玉米、向日葵、黑麦)
- LD 衰减快 → 需要高密度标记
- PCA 单独处理群体结构可能足够
- 有效群体大小远大于自交作物（统计功效不同）

### 多倍体物种
- 必须标注同源 vs 异源
- 异源多倍体：亚基因组必须分别分析
- 同源多倍体：需剂量感知的 GWAS 方法
- 古多倍体 (大豆、白菜)：名义二倍体但大量重复区域

### 无性繁殖作物 (马铃薯、甘蔗、甘薯)
- 杂合度高，基因型复杂
- 品种间亲缘关系可能不明确

### 推断未知物种

如果 `inspect_data.sh` 无法从染色体名推断物种：
1. 检查数据来源说明
2. 染色体数 + 基因组大小可缩小候选范围
3. 让用户确认，不要猜测

---

## 植物特异性参数速查表

以下参数设置基于植物基因组的特殊性，与动物/人类分析有显著差异。

### 差异表达分析 (DESeq2/edgeR)

| 参数 | 通用做法 | 植物特异性调整 | 原因 |
|------|---------|---------------|------|
| `lfcThreshold` | 0 | **1** | 植物小样本(2-rep)策略，减少假阳性 |
| `alpha` | 0.05 | 0.05 | 保持不变 |
| `minReplicatesForReplace` | 7 | 7 | 保持不变 |

```r
# 植物DESeq2推荐参数
dds <- DESeq(dds)
res <- results(dds, lfcThreshold=1, alpha=0.05)
```

### 系统发育树 (IQ-TREE)

| 参数 | 通用做法 | 植物特异性调整 | 原因 |
|------|---------|---------------|------|
| `-m` | GTR+G | **MFP** | 自动选择最佳模型 |
| ASC校正 | 不使用 | **必须使用** (SNP数据) | SNP数据缺失不变位点 |
| `--asc-corr` | - | **lewis** | Lewis校正方法 |

```bash
# SNP数据必须使用ASC校正
iqtree2 -s snp_alignment.fasta -m MFP+ASC --asc-corr lewis -B 1000

# 蛋白质数据标准参数
iqtree2 -s proteins.fasta -m MFP -B 1000 --alrt 1000
```

### lncRNA鉴定 (CNCI/FEELnc)

| 参数 | 通用做法 | 植物特异性调整 | 原因 |
|------|---------|---------------|------|
| CNCI `-m` | ve (脊椎动物) | **pl** (植物) | 植物密码子使用模式不同 |
| FEELnc `--monoex` | 0 | **-1** | 植物单外显子lncRNA常见 |

```bash
# 植物lncRNA分析
CNCI.py -i transcripts.fasta -m pl -o cnci_output
FEELnc_filter.pl -i transcripts.fasta --monoex=-1
```

### Motif分析 (MEME)

| 参数 | 通用做法 | 植物特异性调整 | 原因 |
|------|---------|---------------|------|
| `-mod` | zoops | **anr** | 植物蛋白结构域常重复 |
| `-nmotifs` | 5 | **10** | 植物基因家族更复杂 |

```bash
# 植物基因家族MEME分析
meme proteins.fasta -mod anr -nmotifs 10 -minw 6 -maxw 50 -protein
```

### 甲基化分析 (DMR)

| 参数 | 通用做法 | 植物特异性调整 | 原因 |
|------|---------|---------------|------|
| 窗口大小 | 默认 | **600bp** | 植物甲基化区域特征 |
| 步长 | 默认 | **200bp** | 更精细的分辨率 |
| 最小覆盖度 | 默认 | **4** | 植物样本覆盖度通常较低 |
| \|meth.diff\| | 默认 | **≥20** | 植物甲基化差异阈值 |

```bash
# 植物DMR参数
dmr_caller --window 600 --step 200 --min-cov 4 --meth-diff 20
```

### 启动子分析 (PlantCARE)

| 参数 | 通用做法 | 植物特异性调整 | 原因 |
|------|---------|---------------|------|
| 启动子长度 | 1000bp | **2000bp** | 植物调控元件更分散 |

```bash
# 提取植物启动子序列（2000bp上游）
bedtools flank -g genome.fasta.fai -b 2000 -i genes.bed > promoters.bed
bedtools getfasta -fi genome.fasta -bed promoters.bed -fo promoters.fasta
```

### 共线性分析 (MCScanX)

| 参数 | 通用做法 | 植物特异性调整 | 原因 |
|------|---------|---------------|------|
| GAP_PENALTY | -1 | **-25** | 植物经历多次WGD和基因丢失 |
| MIN_ALIGN | 5 | **5** | 保持不变 |

```bash
# 植物共线性分析参数
# 在MCScanX源码中修改参数后编译
```

### 基因组组装 (hifiasm)

| 参数 | 通用做法 | 植物特异性调整 | 原因 |
|------|---------|---------------|------|
| `--n-hap` | 2 (二倍体) | 根据倍性设置 | 多倍体组装必须指定 |
| BUSCO数据库 | 自动 | **embryophyta_odb10** | 陆地植物专用 |

```bash
# 多倍体植物组装
hifiasm -o assembly --n-hap 4 -t 32 hifi_reads.fastq.gz

# 二倍体植物组装
hifiasm -o assembly -t 32 hifi_reads.fastq.gz
```

### 注意事项汇总

1. **多倍体物种**: 所有涉及基因组/转录组的分析都需要考虑倍性
2. **基因家族扩张**: 植物基因家族通常比动物大，需要调整数量阈值
3. **高重复序列**: 植物基因组重复序列占比高，影响比对和组装
4. **非模式物种**: 注释数据库可能不完整，使用近缘物种数据库
5. **组织特异性**: 植物基因表达组织特异性强，需要考虑采样时间
