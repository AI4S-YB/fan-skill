# 单细胞RNA-seq分析 -- 分析笔记本

## 分析概览

本笔记本提供植物单细胞RNA-seq分析的完整流程，涵盖质量控制、标准化、降维、聚类、细胞类型注释和轨迹分析。

---

## 1. 数据处理流程

### 1.1 Cell Ranger (10x Genomics)

```bash
# 10x Genomics数据处理
cellranger count --id=sample1 \
  --transcriptome=/path/to/reference \
  --fastqs=/path/to/fastq \
  --sample=sample1 \
  --expect-cells=3000 \
  --localcores=16 \
  --localmem=64

# 构建植物参考基因组
cellranger mkref --genome=plant_reference \
  --fasta=genome.fasta \
  --genes=genes.gtf
```

### 1.2 STARsolo (Smart-seq2/其他平台)

```bash
# STARsolo处理单细胞数据
STAR --runThreadN 16 \
  --genomeDir /path/to/genome_index \
  --readFilesIn R1.fastq.gz R2.fastq.gz \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple \
  --soloCBlen 16 \
  --soloUMIlen 12 \
  --soloCBstart 1 \
  --soloCBlen 16 \
  --outFileNamePrefix solo_output/
```

---

## 2. 质量控制

### 2.1 Seurat质控 (R)

```r
library(Seurat)

# 读取数据
sc_data <- Read10X(data.dir = "filtered_feature_bc_matrix/")
seurat_obj <- CreateSeuratObject(counts = sc_data, project = "sample1")

# 计算QC指标
# 植物细胞使用叶绿体基因替代线粒体
chloroplast_genes <- grep("^chloroplast|^cp", rownames(seurat_obj), value = TRUE)
seurat_obj[["percent.chloroplast"]] <- PercentageFeatureSet(seurat_obj, features = chloroplast_genes)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT|^mt")

# 可视化QC指标
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 过滤细胞
seurat_obj <- subset(seurat_obj,
  subset = nFeature_RNA > 100 & nFeature_RNA < 8000 &
           percent.mt < 10 & percent.chloroplast < 5)
```

### 2.2 Scanpy质控 (Python)

```python
import scanpy as sc

# 读取数据
adata = sc.read_10x_mtx('filtered_feature_bc_matrix/')

# 计算QC指标
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['chloroplast'] = adata.var_names.str.startswith('chloroplast-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'chloroplast'], inplace=True)

# 可视化
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# 过滤
sc.pp.filter_cells(adata, min_genes=100)
adata = adata[adata.obs.n_genes_by_counts < 8000, :]
adata = adata[adata.obs.pct_counts_mt < 10, :]
```

### 2.3 植物特异性QC考量

```r
# 植物细胞可能需要调整过滤参数
# 1. 低表达基因：植物组织细胞壁消化可能影响RNA量
# 2. 叶绿体/线粒体：植物细胞需同时考虑两者
# 3. 双细胞：植物细胞大小差异大，双细胞检测更重要

# 使用DoubletFinder检测双细胞
library(DoubletFinder)
seurat_obj <- doubletFinder_v3(seurat_obj,
  PCs = 1:10,
  pN = 0.25,
  pK = 0.005,
  nExp = round(0.05 * ncol(seurat_obj)))
```

---

## 3. 标准化与高变基因

### 3.1 LogNormalize

```r
# 标准对数标准化
seurat_obj <- NormalizeData(seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000)

# 识别高变基因
seurat_obj <- FindVariableFeatures(seurat_obj,
  selection.method = "vst",
  nfeatures = 2000)
```

### 3.2 SCTransform

```r
# SCTransform标准化（推荐用于批次校正）
seurat_obj <- SCTransform(seurat_obj,
  vars.to.regress = "percent.mt",
  verbose = FALSE)
```

---

## 4. 降维与聚类

### 4.1 PCA降维

```r
# PCA降维
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

# 确定PC数量
ElbowPlot(seurat_obj, ndims = 50)

# 植物数据可能需要更多PC
# 检查PC与细胞类型关联
DimHeatmap(seurat_obj, dims = 1:15, cells = 500)
```

### 4.2 UMAP/t-SNE

```r
# 寻找邻居
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)

# 聚类
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)

# UMAP降维
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# 可视化
DimPlot(seurat_obj, reduction = "umap", label = TRUE)
```

---

## 5. 细胞类型注释

### 5.1 Marker基因鉴定

```r
# 寻找cluster marker基因
cluster_markers <- FindAllMarkers(seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25)

# 查看top marker基因
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# 可视化marker基因
FeaturePlot(seurat_obj, features = c("GENE1", "GENE2", "GENE3"))
DotPlot(seurat_obj, features = top_markers$gene) + RotatedAxis()
```

### 5.2 自动注释

```r
# 使用SingleR进行自动注释
library(SingleR)

# 需要参考数据集
ref <- HumanPrimaryCellAtlasData()
predictions <- SingleR(test = seurat_obj@assays$RNA@data,
                       ref = ref,
                       labels = ref$label.main)

seurat_obj$singleR_label <- predictions$labels
```

### 5.3 植物细胞marker基因

```r
# 常见植物细胞类型marker基因
plant_markers <- list(
  epidermis = c("CER1", "CER3", "LTP1"),
  mesophyll = c("RBCS", "CAB1", "PSBA"),
  vascular = c("SUC2", "AAP2", "XYLP"),
  root_tip = c("WOX5", "SCR", "PLT1"),
  guard_cell = c("KAT1", "SLAC1", "OST1"),
  pollen = c("LAT52", "LAT59", "BGP1")
)

# 检查marker表达
DotPlot(seurat_obj, features = unlist(plant_markers))
```

---

## 6. 轨迹分析

### 6.1 Monocle3

```r
library(monocle3)

# 转换为monocle3对象
cds <- as.cell_data_set(seurat_obj)

# 聚类
cds <- cluster_cells(cds, reduction_method = "UMAP")

# 学习轨迹图
cds <- learn_graph(cds)

# 排序细胞
cds <- order_cells(cds)

# 可视化
plot_cells(cds, color_cells_by = "pseudotime")
```

### 6.2 RNA Velocity

```bash
# 使用velocyto计算RNA velocity
velocyto run10x cellranger_output/ output.loom

# 在R中分析
library(velocyto.R)
ldat <- read.loom.matrices("output.loom")
```

---

## 7. 差异表达分析

### 7.1 单样本差异分析

```r
# Cluster间差异分析
markers <- FindMarkers(seurat_obj,
  ident.1 = 1,
  ident.2 = 2,
  min.pct = 0.25,
  logfc.threshold = 0.25)
```

### 7.2 Pseudobulk分析

```r
# 多样本差异分析使用pseudobulk方法
library(DESeq2)

# 聚合为pseudobulk
library muscat
pb <- aggregateData(seurat_obj,
  assay = "RNA",
  fun = "sum")

# 使用DESeq2分析
dds <- DESeqDataSetFromMatrix(
  countData = pb,
  colData = sample_info,
  design = ~ condition)
```

---

## 8. 批次校正与整合

### 8.1 Harmony

```r
library(harmony)

# 运行Harmony
seurat_obj <- RunHarmony(seurat_obj,
  group.by.vars = "batch")

# 使用校正后的embedding
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:30)
```

### 8.2 Seurat整合

```r
# 整合多样本
anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
```

---

## 9. 植物特异性分析

### 9.1 组织特异性分析

```r
# 植物不同组织的细胞类型差异大
# 建议分组织进行分析

# 提取特定组织细胞
tissue_cells <- subset(seurat_obj, subset = tissue == "leaf")
```

### 9.2 细胞周期回归

```r
# 植物细胞周期基因
s_genes <- c("CDC6", "MCM3", "PCNA", "DNA polymerase alpha")
g2m_genes <- c("CDKB1", "CYCB1", "KNOLLE", "Histone H4")

# 计算细胞周期评分
seurat_obj <- CellCycleScoring(seurat_obj,
  s.features = s_genes,
  g2m.features = g2m_genes)

# 回归掉细胞周期效应
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("S.Score", "G2M.Score"))
```

---

## 常见问题

### Q: 细胞数量过少怎么办？
A: 检查质控参数是否过严，考虑降低过滤阈值。检查细胞分离和建库过程。

### Q: 聚类结果不稳定？
A: 尝试不同的resolution参数。检查是否需要批次校正。

### Q: 细胞类型注释不清晰？
A: 收集更多marker基因。考虑使用参考数据集进行标签转移。

### Q: 轨迹分析与预期不符？
A: 检查起点细胞的定义。考虑使用RNA velocity辅助方向判断。

---

## 参考

- Seurat: https://satijalab.org/seurat/
- Scanpy: https://scanpy.readthedocs.io/
- Monocle3: https://cole-trapnell-lab.github.io/monocle3/
- Cell Ranger: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
- Harmony: https://github.com/immunogenomics/harmony
