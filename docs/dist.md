# pgr dist

`pgr dist` 模块提供序列和向量的**距离/相似度计算**功能。它是构建系统发育树（Phylogeny）和聚类分析（Clustering）的核心前置步骤。

## 核心定位

- **定位**：多模式距离计算器（序列 / 索引 / 超向量）。
- **输入**：FASTA/蛋白序列文件、`.pgi` 索引、`.hv` 超向量文件。
- **输出**：Pairwise TSV 格式（`Name1 Name2 Distance ...`），可用于下游分析或矩阵构建。
- **互补**：
  - 上游：`pgr fa`/`pgr fq` (序列处理), `pgr fa six-frame` (蛋白翻译)。
  - 下游：聚类/构树工具。

## 子命令详解

### 1. `pgr dist seq`: 基于 Minimizer / Closed Syncmer 的序列距离
*利用采样策略快速计算序列间的 Mash 距离，适合大规模基因组比较。*

- **核心算法**:
  - **Minimizer**: 在窗口 $w$ 内选择哈希值最小的 $k$-mer 作为代表，大幅减少计算量。
  - **Closed Syncmer** (Edgar 2021): 窗口内最小 s-mer 哈希出现在窗口首端或末端的
    $k$-mer 集合。与 minimizer 不同，syncmer 提供**稀疏但完整**的覆盖（平均深度约 2×），
    且序列与其反向互补得到相同哈希集（Mash/Jaccard 需要）。
  - **Mash Distance**: 基于 Jaccard Index 估算的突变距离（Mutation Distance）。公式：$D \approx -\frac{1}{k} \ln(\frac{2J}{1+J})$。
- **支持指标**:
  - **Mash Distance**: 演化距离估计。
  - **Jaccard Index**: 集合相似度 $J = |A \cap B| / |A \cup B|$。
  - **Containment Index**: 包含度 $C = |A \cap B| / |A|$，适合宏基因组或质粒检测。
- **采样器 (`--sampler`)**:
  - `minimizer` (默认): 见上。
  - `syncmer`: DNA 默认 `-k 8 -w 55`、蛋白默认 `-k 7 -w 5`（syng 默认参数，未显式指定时自动套用）。
- **哈希算法 (`--hasher`)**:
  - `rapid`: RapidHash (默认，速度最快)。
  - `fx`: FxHash。
  - `murmur`: MurmurHash3。
  - `mod`: **Mod-Minimizer**。针对 DNA 序列生成 Canonical k-mers（正反义链统一），避免链的方向影响。
  - 注：`--hasher` 仅对 minimizer 生效；syncmer 使用 2-bit canonical rolling hash (DNA)
    或 RapidHash (蛋白)。
- **主要参数**:
  - `-k`/`--kmer`: k-mer 长度 (默认 7)。
  - `-w`/`--window`: Minimizer 窗口大小 (默认 1)。
  - `--sampler`: `minimizer` (默认) 或 `syncmer`。
  - `--protein`: 声明输入为蛋白序列（影响采样与哈希路径）。
  - `--merge`: 将文件内所有序列合并为一个集合计算（例如比较两个基因组）。
  - `--zero`: 输出 Jaccard 为 0 的结果（默认跳过）。
  - `--sim`: 将 Mash 距离转为相似度输出。
  - `--list-files`: 将输入视为文件列表（每行一个序列文件路径）。
  - `-p`/`--parallel`: 并行线程数。

### 2. `pgr dist hv`: 基于 Hypervector 的序列距离
*利用超维计算（Hyperdimensional Computing, HDC）技术，将序列采样（minimizer 或 closed
syncmer）映射为固定维度的向量。*

- **核心概念**:
  - 将 k-mer 映射为高维空间（如 4096 维）中的随机向量。
  - 通过向量叠加（Superposition）表示整条序列。
  - 具有全息特性，对噪声鲁棒，且计算速度极快（位运算）。
- **优势**:
  - 维度固定，计算复杂度与序列长度无关。
  - 适合超大规模数据集的快速预筛选。
- **`.hv` 文件模式**: 输入为 `.hv` 文件（`pgr pgi to-hv` 产物）时直接比较，
  无需重新采样序列；要求两侧 `k`、维度与稀疏更新数（`--sparse`）一致。
  比较用余弦相似度恢复共享 k-mer 数（`inter = cos × √(n1·n2)`），是
  `pgr dist pgi` 的约 50× 快近似（排序相关性 ρ≈0.97）。
- **参数**（与 `seq` 共享 sampler/hash 参数，含义与默认值相同）:
  - `--dim`: 向量维度 (默认 4096，需为 32 的倍数)。
  - `--sampler`: `minimizer` (默认) 或 `syncmer`（syncmer 默认 DNA `-k 8 -w 55`、
    蛋白 `-k 7 -w 5`）。
  - `--hasher`: 哈希算法（`rapid`/`fx`/`murmur`/`mod`，默认 `rapid`）。
  - `-k`/`--kmer`: k-mer 长度 (默认 7)。
  - `-w`/`--window`: Minimizer 窗口大小 (默认 1)。
  - `--sim`: 将 Mash 距离转为相似度输出。
  - `--list-files`: 将输入视为文件列表。
  - `-p`/`--parallel`: 并行线程数。

### 3. `pgr dist pgi`: 基于 .pgi 索引归并的精确距离
*将两个 `.pgi` 索引的排序 k-mer 流线性归并，计算确定性的 Jaccard/Containment/Mash 距离。*

- **输入**：两个 `.pgi` 索引（`pgr pgi build` 生成）。
- **算法**：两排序流线性归并（O(|K1|+|K2|)），共享 k-mer 集合精确计数。
- **要求**：两侧索引采样参数（k/syncmer/window）必须一致，否则报错。
- **输出格式**：`Name1 Name2 Total1 Total2 Inter Union Mash Jaccard Containment`。
- **注意**：该距离是"采样集合的精确距离"（确定性、零采样方差），但 syncmer
  采样位置随变异漂移，与真实身份率的排序相关性约 0.5；超大规模粗筛建议用
  `dist seq`（k=8，相关性 0.82），`dist hv`（.hv 模式）为 `dist pgi` 的快速近似。

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

### 场景 C：基因组文件间 Hypervector 比较
```bash
# 将每个文件的所有序列合并为一个 hypervector，计算两个基因组之间的距离
pgr dist hv genome1.fa genome2.fa
```

## 未来规划 (Roadmap)

通用向量距离度量（euclid/cosine/jaccard）目前仅在 `libs::linalg` 中提供，
`pgr dist vector` CLI 已于 2026-07 移除，暂无恢复规划。
