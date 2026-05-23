# 植物小 RNA 分析 — 特殊注意事项

## 植物特有的小 RNA 特征

### 1. 小 RNA 种类多样性

植物小 RNA 种类比动物更为丰富，包括：
- **miRNA** (microRNA): 21-22nt，由 DCL1 加工
- **siRNA** (small interfering RNA): 24nt，由 DCL3 加工
- **phasiRNA** (phased siRNA): 21nt/24nt，由 miRNA 触发产生
- **tasiRNA** (trans-acting siRNA): 21nt
- **hc-siRNA** (heterochromatic siRNA): 24nt，参与 DNA 甲基化

### 2. miRNA 长度偏好

- 植物 miRNA 以 **21nt** 为主 (动物以 22nt 为主)
- 水稻中 24nt miRNA 也较为常见
- 长度分布图的双峰模式（21nt + 24nt）是植物小 RNA 文库的典型特征

### 3. miRNA-mRNA 互补性

- 植物 miRNA 与靶基因 mRNA 几乎完全互补（通常 < 4 个错配）
- 主要剪切靶位点位于 miRNA 第 10-11 位碱基之间
- 这一特性使得靶基因预测比动物更为准确

### 4. 基因组复杂度

植物基因组通常较大且多倍体化程度高：

| 物种 | 基因组大小 | 倍性 | 注释难度 |
|------|-----------|------|---------|
| 拟南芥 | ~135 Mb | 二倍体 | 低 |
| 水稻 | ~430 Mb | 二倍体 | 低 |
| 玉米 | ~2.4 Gb | 古四倍体 | 中 |
| 大豆 | ~1.1 Gb | 古四倍体 | 中 |
| 小麦 | ~17 Gb | 六倍体 | 高 |

### 5. 重复序列

植物基因组含有大量重复序列，特别是转座子 (TE)。24nt siRNA 来源于 TE 区域，需要注意：

- miRNA 预测时排除 TE 区域
- 比对时允许 multi-mapping reads（`bowtie -k` 参数）

### 6. 物种特异性考虑

#### 模式物种 (推荐方式)
- **拟南芥** (Arabidopsis thaliana): 最完善的 miRBase 注释
- **水稻** (Oryza sativa): 第二大 miRNA 数据库
- **玉米** (Zea mays): 较完善的注释
- **大豆** (Glycine max): 使用 psRNATarget

#### 非模式物种
- 优先使用从头预测方法 (miRDeep2 novel mode)
- 参考近缘模式物种的 miRNA
- 考虑短期 read 比对（允许 0-1 个错配）

### 7. 降解组分析注意事项

植物降解组测序 (PARE/Degradome-seq)：
- 需要 > 10M reads
- 基因间区 reads 可能反映 siRNA 靶向的降解
- 考虑非典型剪切（不在第 10-11 位碱基）

### 8. 常用资源

| 资源 | 链接 | 用途 |
|------|------|------|
| miRBase | https://mirbase.org/ | miRNA 序列数据库 |
| psRNATarget | http://plantgrn.noble.org/psRNATarget/ | 靶基因预测 |
| Plant Non-coding RNA Database | http://structuralbiology.cau.edu.cn/PNRD/ | 植物 ncRNA |
| sRNAanno | http://www.plantsrnas.org/ | 植物小 RNA 注释 |
| PmiREN | https://www.pmiren.com/ | 植物 miRNA 百科全书 |
