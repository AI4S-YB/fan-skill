# 植物数据库接入指南

分析需要的参考数据从哪里获取，以及如何获取。

---

## 参考基因组

| 物种 | 首选来源 | URL 模式 | 备选 |
|------|---------|---------|------|
| 水稻 | Ensembl Plants | `https://plants.ensembl.org/Oryza_sativa/Info/Index` | RAP-DB, RiceGFD |
| 玉米 | Ensembl Plants | `https://plants.ensembl.org/Zea_mays/Info/Index` | MaizeGDB |
| 小麦 | Ensembl Plants | `https://plants.ensembl.org/Triticum_aestivum/Info/Index` | URGI, WheatOmics |
| 大豆 | Phytozome | `https://phytozome-next.jgi.doe.gov/info/Gmax_Wm82_a4_v1` | SoyBase |
| 棉花 | CottonGen | `https://www.cottongen.org/` | Phytozome |
| 油菜 | Ensembl Plants | `https://plants.ensembl.org/Brassica_napus/Info/Index` | BnPIR |
| 拟南芥 | TAIR | `https://www.arabidopsis.org/` | Araport11 |
| 番茄 | Sol Genomics | `https://solgenomics.net/` | Phytozome |
| 马铃薯 | Spud DB | `https://spuddb.uga.edu/` | Phytozome |
| 大麦 | Ensembl Plants | `https://plants.ensembl.org/Hordeum_vulgare/Info/Index` | BAR |
| 高粱 | Phytozome | `https://phytozome-next.jgi.doe.gov/info/Sbicolor_v3_1_1` | SorghumBase |
| 甘蔗 | Sugarcane Hub | `https://sugarcane-genome.cirad.fr/` | — |

### 自动化获取参考基因组

如果用户本地没有参考基因组，建议下载方式：

```bash
# Ensembl Plants (推荐 — 多数作物)
wget ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-XX/fasta/<species>/dna/

# Phytozome (大豆、高粱等)
# 需要 JGI 账号，建议用户提前下载

# NCBI Genome (通用备选)
datasets download genome accession <GCA_ID> --filename genome.zip
```

**Agent 行为**：不要自动下载大型参考基因组（>1GB）。提示用户提供路径，或在 params.yaml 中指定。

---

## 功能注释

### 首选：Plant Reactome

- URL: `https://plantreactome.gramene.org/`
- 覆盖：水稻、玉米、拟南芥等
- 下载：GMT 格式的 pathway 定义文件
- ORA 富集分析优先使用 Plant Reactome，其次 KEGG Plant

### 备选：MapMan / Mercator

- MapMan: `https://mapman.gabipd.org/`
- Mercator: `https://www.plabipd.de/portal/mercator4` — 在线注释，适合非模式植物

### 通用：GO + KEGG

- GO 富集：clusterProfiler (R) / goatools (Python)
- KEGG Plant：`https://www.genome.jp/kegg/pathway.html` → "Plant" 分类
- AgriGO: `http://systemsbiology.cau.edu.cn/agriGOv2/` — 农业物种 GO 富集

---

## 变异与表型

| 资源 | 覆盖 | 用途 |
|------|------|------|
| Gramene Variation | 水稻、玉米、小麦等 | 已知变异注释 |
| European Variation Archive | 所有物种 | 原始变异数据 |
| RiceVarMap | 水稻 | 水稻变异精细注释 |
| SoyBase | 大豆 | 大豆 QTL/GWAS 汇总 |
| WheatOmics | 小麦 | 小麦多组学整合 |

---

## 表达数据

| 资源 | 覆盖 | 用途 |
|------|------|------|
| Expression Atlas (EBI) | 多物种 | 差异表达参考 |
| BAR (Bio-Analytic Resource) | 大麦、拟南芥 | eFP Browser 可视化 |
| Rice eFP | 水稻 | 水稻表达图谱 |
| Maize eFP | 玉米 | 玉米表达图谱 |

### 候选基因优先排序

当 GWAS 显著位点落在基因间区时，按以下优先级查找候选基因：

1. 显著 SNP 上下游 LD 衰减范围 (kb) = 候选区间
2. 候选区间内的基因 → 查表达数据库 (该组织是否表达?)
3. 基因功能注释 (Plant Reactome > GO > KEGG)
4. 同源基因在拟南芥/水稻中的功能 (跨物种推断)

---

## 基因表达数据（RNA-seq 专用）

### 公共表达数据
| 资源 | 覆盖 | 用途 |
|------|------|------|
| Expression Atlas (EBI) | 多物种 | 差异表达参考、基线表达 |
| SRA (NCBI) | 所有物种 | 原始测序数据 |
| GEO (NCBI) | 所有物种 | 公共基因表达数据集 |
| Rice Expression Database | 水稻 | 水稻表达图谱 |
| MaizeGDB Expression | 玉米 | 玉米发育表达图谱 |
| Wheat Expression Browser | 小麦 | 小麦同源基因表达 |
| eFP Browser (BAR) | 拟南芥、大麦等 | 可视化基因表达 |

### 富集分析数据库
| 资源 | 覆盖 | 推荐优先级 |
|------|------|:--------:|
| Plant Reactome | 水稻、玉米、拟南芥 | 1 (首选) |
| MapMan / Mercator | 多物种（在线注释） | 2 (非模式植物) |
| AgriGO v2 | 农业物种 | 2 (GO 富集) |
| KEGG Plant | 多物种 | 3 (备选) |
| clusterProfiler (R) | GO + KEGG | 通用工具 |

---

## 基因组选择资源

### 训练群体数据
| 资源 | 覆盖 | 用途 |
|------|------|------|
| Rice 3K | 水稻 3,024 份种质 | GS 训练群体 |
| Maize NAM | 玉米 25 个 founder × B73 | 多亲本 GS |
| Wheat WatSeq | 小麦多样性面板 | GS 验证 |
| SoyNAM | 大豆 5,600 RILs | 大豆 GS |

### GS 软件与工具
| 工具 | 方法 | 语言 |
|------|------|------|
| rrBLUP | Ridge regression BLUP | R |
| BGLR | Bayesian (A/B/C/Lasso) | R |
| sommer | Mixed models (G×E) | R |
| AlphaPeel | AlphaBayes | CLI |
| DMU | Mixed models (large scale) | CLI |
