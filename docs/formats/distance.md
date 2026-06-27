# Distance and Matrix Formats

pgr uses three matrix structures internally for clustering and distance
analysis, and supports two external file formats for distance matrices:
PHYLIP and Pairwise.

## External formats

### PHYLIP distance matrix

PHYLIP 距离矩阵格式是系统发育分析中的通用格式。`pgr` 在 `src/cmd_pgr/mat`
和 `src/libs/pairmat.rs` 中提供了一系列工具来处理这种格式。

`pgr` 支持标准 (Strict) 和宽松 (Relaxed) 的 PHYLIP 格式。

**Relaxed PHYLIP (默认输入支持)**：
- 第一行（可选）: 序列数量。
- 数据行: 序列名称后跟距离数值。
- 分隔符: 空白字符（空格或制表符）。
- 矩阵形式: 支持全矩阵 (Full Square) 或下三角矩阵 (Lower Triangular)。
- 名称长度: 不受限制。

**Strict PHYLIP (`strict` 模式输出)**：
- 遵循原始 PHYLIP 标准。
- 序列名称: 严格截断为 10 个字符，左对齐并用空格填充。
- 数值格式: 空格分隔，通常保留 6 位小数。

### Pairwise distance

Pairwise 格式是一种简单的三列 TSV 格式，用于表示序列两两之间的距离，
常用于作为中间格式或图数据的输入。

| Column  | Description          |
|---------|----------------------|
| `name1` | Sequence 1 name      |
| `name2` | Sequence 2 name      |
| `distance` | Distance or score |

`pgr` 提供了矩阵与 Pairwise 列表的互转：

- **Matrix to Pair (`pgr mat to-pair`)**：将 PHYLIP 矩阵展平为 Pairwise 列表。
- **Pair to Matrix (`pgr mat to-phylip`)**：将 Pairwise 列表组装回 PHYLIP 矩阵，
  支持 `--missing` 和 `--same` 参数。

## Internal matrix structures

`pgr` 在 `src/libs/pairmat.rs` 中定义了三种核心的矩阵结构。

### ScoringMatrix

- **用途**：稀疏或按需构建的评分/距离矩阵。
- **底层存储**：`HashMap<(usize, usize), T>`。
- **特点**：稀疏存储，仅保存显式设置的值；支持对角线和非对角线的默认值；
  逻辑对称（`get(i,j)` 等价于 `get(j,i)`）。

### CondensedMatrix

- **用途**：高效层次聚类（如 `clust hier`），支持更大规模数据。
- **底层存储**：`Vec<f32>`，仅存上三角（不含对角线），内存占用 $N(N-1)/2$。
- **索引映射**：$(i, j)$ 且 $i < j$ → $k = N \cdot i - i(i+1)/2 + (j - i - 1)$。
- **特点**：强制对称，对角线假定为 0，不存储名称映射，纯数值计算。

### NamedMatrix

- **用途**：稠密且带有行列名称的距离矩阵（如 PHYLIP 内存表示）。
- **底层存储**：`IndexMap`（名称索引）+ `CondensedMatrix` + 可选对角线向量。
- **特点**：组合封装，通过名称索引访问底层 `CondensedMatrix`；
  支持可选对角线存储（`mat transform --normalize` 依赖）；
  $N=10,000$ 时约 200MB。

## pgr commands

| Command                | Description                              |
|------------------------|------------------------------------------|
| `pgr mat format`       | 在 full / lower / strict 格式间转换      |
| `pgr mat subset`       | 按名称列表提取子矩阵                     |
| `pgr mat compare`      | 计算两个矩阵的相关性（Pearson、Spearman）|
| `pgr mat to-pair`      | PHYLIP 矩阵 → Pairwise 列表              |
| `pgr mat to-phylip`    | Pairwise 列表 → PHYLIP 矩阵              |
