# 通用排错指南

## 环境问题

### R 包缺失
```
Error in library(GAPIT3): there is no package called 'GAPIT3'
```
→ `devtools::install_github("jiabowang/GAPIT3")` 或使用 singularity 容器

### Conda 环境冲突
```
LibMamba UnsatisfiableError
```
→ 使用 `pixi` 代替 conda，或创建独立环境

### Singularity 镜像不可用
→ 检查 `test_config.yaml` 中 `container_image` 路径
→ 如镜像在 `/share/` 下，确认 `--bind` 参数

## 数据问题

### 染色体名不匹配
```
Warning: Chromosome name 'Chr01' not found in reference
```
→ 植物的染色体命名极不统一。常见模式：
- 水稻: Chr01..Chr12 或 chr01..chr12 或 1..12
- 玉米: chr1..chr10 或 1..10
- 小麦: chr1A..chr7D 或 1A..7D
→ 用 `inspect_data.sh` 的染色体名检测功能

### 编码问题 (中文 Windows 数据)
```
Error: line 1 contains embedded nulls
```
→ `iconv -f GBK -t UTF-8 input.csv > input_utf8.csv`

### 表型文件格式
常见陷阱：
- 表型值用中文标记（"高"/"低"而非数值）
- 缺失值用 "NA"/"."/"-" 不统一
- 品种名含特殊字符（空格、/、括号）
→ `inspect_data.sh` 会自动检测并报告

## 分析问题

### λ 值偏高 (GWAS)
- 加 PCA 协变量 (PC1-PC5)
- 加 K 矩阵 (亲缘关系)
- 检查是否存在极端样本 (>6 SD on any PC)
- 自交作物 λ>1.2 不一定是问题

### 零个显著位点 (GWAS)
可能原因（按概率排序）：
1. 标记太少、统计功效不足 — 这是最常见原因
2. 样本量不够
3. 该性状确实无主效 QTL（微效多基因控制）
4. 表型测量误差大
→ 诚实报告，不要过度解读

### 群体结构噪音 (Population Genetics)
- PC1/PC2 解释方差异常高 → 正常的强群体分化
- PC1 解释方差过低 (<3%) → 可能是 LD 未被充分过滤
- Admixture 交叉验证错误不收敛 → 增加 max iterations 或调整 K 范围
