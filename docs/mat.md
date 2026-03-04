# mat

`pgr mat` 模块专注于**距离矩阵（Distance Matrix）**的操作与转换。它是 `pgr clust`（聚类与构树）的上游数据准备与预处理工具集。

## 核心定位

- **输入/输出**：主要处理 **PHYLIP** 格式的距离矩阵（稠密）和 **Pairwise TSV**（稀疏/列表）格式。
- **功能**：格式互转、子集提取、矩阵比较、标准化。
- **目标**：为系统发育分析（Phylogenetics）和统计聚类提供标准、高效的数据接口。

## 支持格式与数据结构

### 1. PHYLIP 距离矩阵 (Dense)
`pgr` 内部使用 `NamedMatrix` 结构存储，底层为 `CondensedMatrix`（一维数组存储上三角或下三角），内存占用约为 $O(N^2/2)$。
- **Full (Square)**: 标准 $N \times N$ 矩阵，包含冗余的对称部分。
- **Lower-triangular**: 仅包含下三角部分，文件体积减半。
- **Strict**: 遵循原始 PHYLIP 标准（10字符名称限制），用于兼容老旧软件（如 Phylip 3.695）。
- **Relaxed**: 允许长名称，通常以 Tab 分隔。`pgr` 默认支持读取此格式。

### 2. Pairwise TSV (Sparse-like)
稀疏或列表形式的距离数据，适合存储图结构或仅关注部分配对的情况。
- **格式**：`name1  name2  distance`
- **特点**：
  - 适合稀疏图或作为与其他工具（如 BLAST/MMseqs2）的交换格式。
  - 转换为矩阵时，未列出的配对将被视为缺失值或默认值。

## 子命令详解 (Subcommands)

### 格式转换 (Conversion)

#### `pgr mat to-phylip`
将 Pairwise TSV 转换为 PHYLIP 矩阵。
- **作用**：从比对结果（如 `blast --outfmt 6`）构建距离矩阵，用于后续构树。
- **参数**：
  - `--same <FLOAT>`: 对角线元素（自身到自身）的距离值，默认为 0。
  - `--missing <FLOAT>`: 缺失配对的距离值（如未比上的序列），默认为 1.0（最大距离）。
- **注意**：会自动收集所有出现的 ID 并构建 $N \times N$ 矩阵。

#### `pgr mat to-pair`
将 PHYLIP 矩阵转换为 Pairwise TSV。
- **作用**：将矩阵导出为边列表，用于图聚类（如 `mcl`）或网络可视化（Cytoscape）。
- **输出**：三列格式 `A B 0.123`。通常只输出下三角或上三角部分以避免重复。

#### `pgr mat format`
PHYLIP 格式间的转换与标准化。
- **作用**：清洗矩阵格式，使其符合特定软件的要求。
- **模式 (`--mode`)**：
  - `full`: 输出标准 $N \times N$ Tab 分隔矩阵，保留完整长名称。
  - `lower`: 输出下三角矩阵，节省磁盘空间。
  - `strict`: **截断名称**至10个字符，左对齐填充空格，数值固定宽度。用于兼容原始 Phylip 工具包。

### 操作与分析 (Operations)

#### `pgr mat subset`
基于名称列表提取子矩阵。
- **作用**：从大矩阵中提取特定物种或基因家族的子集进行精细分析。
- **输入**：
  - 矩阵文件。
  - ID 列表文件（每行一个 ID）。
- **行为**：
  - 保持矩阵结构，自动处理行列索引。
  - 若列表中的 ID 在矩阵中不存在，会输出警告并跳过。
  - 输出矩阵的顺序与 ID 列表文件的顺序一致（可用于重排矩阵）。

#### `pgr mat compare`
计算两个矩阵之间的相关性或差异。
- **作用**：评估不同距离计算方法的一致性，或聚类前后的信息损失（Cophenetic Correlation）。
- **前提**：自动取两个矩阵的**公共 ID 交集**进行比较。
- **指标 (`--method`)**：
  - **相关性**：Pearson ($r$), Spearman ($\rho$, 秩相关), Cosine, Jaccard。
  - **距离/误差**：MAE (平均绝对误差), Euclidean (欧氏距离)。
  - `all`: 同时计算以上所有指标。

## 推荐工作流

### 场景 A：从 BLAST 结果构树

```bash
# 1. 解析 BLAST 结果为 Pairwise Distance (假设已计算 distance = 1 - identity)
# 注意：需确保 A-B 和 B-A 均存在，或仅依赖单向
awk '{print $1, $2, 100-$3}' blast.out > pairs.tsv

# 2. 转换为 PHYLIP 矩阵，未比上的设为最大距离 100
pgr mat to-phylip pairs.tsv --missing 100 -o matrix.phy

# 3. 构建 NJ 树
pgr clust nj matrix.phy > tree.nwk
```

### 场景 B：提取子集进行精细分析

```bash
# 1. 准备感兴趣的 ID 列表
cat interesting_ids.txt
# gene_A
# gene_B
# ...

# 2. 从全基因组矩阵中提取子矩阵
pgr mat subset genome_dist.phy interesting_ids.txt -o sub_matrix.phy

# 3. 使用 Ward 聚类分析子集
pgr clust hier sub_matrix.phy --method ward > sub_tree.nwk
```

### 场景 C：评估两种距离计算方法的一致性

```bash
# 比较基于 K-mer (mash) 和基于比对 (ani) 的距离矩阵
pgr mat compare mash_dist.phy ani_dist.phy --method pearson,spearman

# Output:
# Sequences in matrices: 100 and 100
# Common sequences: 100
# Method    Score
# pearson   0.985432
# spearman  0.971234
```

### 场景 D：准备用于 Phylip 软件包的数据

```bash
# 将长名称矩阵转换为严格的 Phylip 格式
pgr mat format modern.phy --mode strict -o input.infile

# 然后运行 neighbor (原始 Phylip 程序)
neighbor < input.infile
```
