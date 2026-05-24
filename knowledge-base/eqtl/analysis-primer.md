# eQTL 结果解读

## 一句话解释

eQTL 分析：在全基因组范围内寻找影响基因表达水平的遗传位点。

## 能回答什么

- 哪些 SNP 调控了附近基因的表达（cis-eQTL）？
- 哪些 SNP 远程调控了其他基因的表达（trans-eQTL）？
- GWAS 显著位点是否通过影响特定基因表达来影响表型（coloc）？

## 不能回答什么

- eQTL SNP 是因果调控变异吗？（可能只是与因果变异 LD 关联）
- trans-eQTL 的调控是直接的还是间接的？（需要验证）
- 组织特异的 eQTL 在其他组织中也存在吗？

## 典型输出

| 文件 | 含义 |
|------|------|
| `cis_eqtl_results.csv` | cis-eQTL: SNP-基因对及显著性 |
| `trans_eqtl_results.csv` | trans-eQTL: 远程关联 |
| `egenes_list.csv` | 有显著 eQTL 的基因列表 |
| `coloc_results.csv` | GWAS-eQTL 共定位分析 |

## 常见结果模式

### 零个 eGene
样本量 < 50 时正常。检查协变量是否过多（可能移除了遗传信号）。

### trans-eQTL 热点
一个位点关联数百个基因 → 检查是否多倍体 cross-homeolog 干扰、着丝粒区域 LD。
