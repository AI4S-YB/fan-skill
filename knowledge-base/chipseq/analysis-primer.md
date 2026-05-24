# ChIP-seq 结果解读

## 一句话解释

ChIP-seq 通过抗体富集与目标蛋白结合的 DNA 片段，鉴定该蛋白在基因组上的全部结合位点。

## 能回答什么

- 这个转录因子结合在哪些基因附近？
- 哪些基因组区域有组蛋白修饰（H3K4me3, H3K27ac 等）？
- 处理前后结合模式有什么变化？

## 不能回答什么

- TF 结合一定调控基因表达吗？（需要 RNA-seq 验证）
- 结合位点都是功能性的吗？（有些结合不产生调控效应）
- 非模式植物的抗体从哪里来？（商品化抗体多针对人/鼠）

## 典型输出

| 文件 | 含义 |
|------|------|
| `peaks.narrowPeak` | 显著富集区域 |
| `peak_annotation.csv` | Peak 的基因组注释 |
| `differential_binding.csv` | 差异结合位点 |
| `motif_enrichment.txt` | 富集的 TF 结合 motif |

## 常见结果模式

### Peak 数量太少
降低 MACS2 q 值阈值；检查 IP 效率（ChIP-qPCR 验证）；植物基因组大，peak 数量可以很多。

### Input 信号异常高
检查超声打断效率和 Input DNA 浓度。植物多酚多→提取效率低。
