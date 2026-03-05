# dist

`pgr dist` 模块提供序列和向量的**距离/相似度计算**功能。它是构建系统发育树（Phylogeny）和聚类分析（Clustering）的核心前置步骤。

## 核心定位

- **定位**：多模式距离计算器（序列、向量）。
- **输入**：FASTA 序列文件、特征向量文件。
- **输出**：Pairwise TSV 格式（`Name1 Name2 Distance ...`），可用于下游分析或矩阵构建。
- **互补**：
  - 上游：`pgr fa`/`pgr fq` (序列处理), `pgr fa count` (生成 k-mer 向量)。
  - 下游：`pgr clust` (聚类/构树), `pgr mat` (矩阵操作)。

## 子命令详解

### 1. `pgr dist seq`: 基于 Minimizer 的序列距离
*利用 Minimizer 采样策略快速计算序列间的 Mash 距离，适合大规模基因组比较。*

- **核心算法**:
  - **Minimizer**: 在窗口 $w$ 内选择哈希值最小的 $k$-mer 作为代表，大幅减少计算量。
  - **Mash Distance**: 基于 Jaccard Index 估算的突变距离（Mutation Distance）。公式：$D \approx -\frac{1}{k} \ln(\frac{2J}{1+J})$。
- **支持指标**:
  - **Mash Distance**: 演化距离估计。
  - **Jaccard Index**: 集合相似度 $J = |A \cap B| / |A \cup B|$。
  - **Containment Index**: 包含度 $C = |A \cap B| / |A|$，适合宏基因组或质粒检测。
- **哈希算法 (`--hasher`)**:
  - `rapid`: RapidHash (默认，速度最快)。
  - `fx`: FxHash。
  - `murmur`: MurmurHash3。
  - `mod`: **Mod-Minimizer**。针对 DNA 序列生成 Canonical k-mers（正反义链统一），避免链的方向影响。
- **主要参数**:
  - `-k`/`--kmer`: k-mer 长度 (默认 21)。
  - `-w`/`--window`: Minimizer 窗口大小 (默认 5)。
  - `--merge`: 将文件内所有序列合并为一个集合计算（例如比较两个基因组）。

### 2. `pgr dist hv`: 基于 Hypervector 的序列距离
*利用超维计算（Hyperdimensional Computing, HDC）技术，将序列映射为固定维度的向量。*

- **核心概念**:
  - 将 k-mer 映射为高维空间（如 4096 维）中的随机向量。
  - 通过向量叠加（Superposition）表示整条序列。
  - 具有全息特性，对噪声鲁棒，且计算速度极快（位运算）。
- **优势**:
  - 维度固定，计算复杂度与序列长度无关。
  - 适合超大规模数据集的快速预筛选。
- **参数**:
  - `--dim`: 向量维度 (默认 4096，需为 32 的倍数)。

### 3. `pgr dist vector`: 通用向量距离
*计算数值向量之间的距离或相似度。*

- **输入格式**: `Name <tab> val1,val2,val3...` (CSV 格式的数值列表)。
- **计算模式 (`--mode`)**:
  - **Euclidean**: 欧氏距离 $L_2$。支持转换为相似度 $S = 1 / e^D$。
  - **Cosine**: 余弦相似度（-1 到 1）。支持转换为距离 $D = 1 - S$。
  - **Jaccard**: 加权 Jaccard 相似度（针对数值向量）。
- **二值化 (`--bin`)**:
  - 将非零值视为 1，零值视为 0。
  - 适合处理 Presence/Absence 矩阵。

## 典型用法

### 场景 A：基因组快速比较 (Mash)
```bash
# 比较两个基因组文件（合并所有 contigs）
pgr dist seq genome1.fa genome2.fa --merge -k 21 -w 10

# 输出: File1 File2 Total1 Total2 Inter Union Mash Jaccard Containment
```

### 场景 B：所有序列两两比较
```bash
# 计算文件中所有序列的两两距离
pgr dist seq genes.fa -k 7 -w 1 > dist.tsv

# 使用 4 线程加速
pgr dist seq genes.fa -p 4 > dist.tsv
```

### 场景 C：向量相似度计算
```bash
# 计算余弦相似度
pgr dist vector features.txt --mode cosine --sim > similarity.tsv

# 计算二值化的 Jaccard 距离 (1 - Jaccard)
pgr dist vector presence.txt --mode jaccard --bin --dis > distance.tsv
```

## 未来规划 (Roadmap)

### Alignment-based Metrics (计划中)
*基于严格比对的距离计算，精度最高但速度较慢。*
- **Kimura 2-Parameter (K2P)**: 区分转换与颠换。
- **Jukes-Cantor (JC69)**: 基础核苷酸替换模型。
- **p-distance**: 简单的 Hamming 距离。

### SciPy 兼容性扩展
*计划借鉴 `scipy.spatial.distance` 支持更多度量：*
- **Bray-Curtis**: 生态学常用。
- **Canberra**: 对小数值敏感。
- **Mahalanobis**: 考虑协方差。
