# clust boot：层次聚类的多尺度 Bootstrap p-value（pvclust 风格）

本页基于仓库内置的 R 包源码 `pvclust`（见 [pvclust.R](file:///c:/Users/wangq/Scripts/pgr/pvclust/R/pvclust.R)、[pvclust-internal.R](file:///c:/Users/wangq/Scripts/pgr/pvclust/R/pvclust-internal.R)）梳理其算法与数据结构，并给出 `pgr` 侧计划新增命令 `pgr clust boot` 的接口与输出约定。

`pvclust` 的核心价值：给 dendrogram 的每个内部节点（cluster / edge）计算 **BP/AU/SI** 三类支持度（及标准误），用于回答“这个簇是不是稳定/显著”的问题，而不仅是“切成几类”。

---

## 1. 背景与概念

多尺度 bootstrap 为树上的“簇”（内部节点）提供的是一种稳定性度量：在“可重采样的观测单元”扰动下，这个簇是否反复出现。它并不直接提升聚类树的准确性，也不保证“出现频率高”的簇就一定等价于真实分群；它更像是一个诊断工具，用来区分对数据扰动敏感的分支（不稳定、可能受噪声/参数影响）与相对稳健的分支（在当前数据与距离/链接策略下可重复）。

对进化树的适用性取决于重采样单元是否符合系统发育 bootstrap 的假设：若输入矩阵的列是近似独立同分布的位点/特征（例如对齐序列的位点、SNP 位点或可近似为独立的基因家族特征），则对列进行重采样并在每次重采样下重建树，得到的支持度可作为“类 bootstrap”支持；但若列特征强相关、或树的构建依赖显式进化模型（如 ML/贝叶斯），则该方法更适合作为聚类树稳定性评估，而不应直接解读为严格的系统发育分支支持度。

### 1.1 输入数据的方向（非常关键）

`pvclust` 默认假设输入矩阵为 `(n × p)`，即 `Features × Items`：
- **行 (n)**：特征/属性（Attributes/Variables），如基因表达量、SNP 位点。`pvclust` 默认对行进行重采样（Bootstrap）。
- **列 (p)**：样本/条目（Samples/Items），如物种、病人。`pvclust` 计算列与列之间的距离并对列进行聚类。

这与许多生物学数据矩阵（如 `Genome × Domain`，行是物种，列是功能域特征）的方向恰好相反。

为了适配两种场景并避免额外的转置操作，`pgr clust boot` 将提供 `--along` 参数来指定聚类方向：

-   `--along row` (**默认**，推荐用于 `domain.tsv` 等)：
    -   **聚类方向**：行（如物种、基因组）。
    -   **重采样方向**：列（如功能域、特征）。
    -   符合多数生物学特征矩阵的直觉（Samples × Features）。
-   `--along col` (复刻 `pvclust` 原生行为)：
    -   **聚类方向**：列。
    -   **重采样方向**：行。

### 1.3 核心算法流程：多尺度 Bootstrap

`pvclust` 采用 **多尺度 Bootstrap (Multiscale Bootstrap)**，这是其区别于传统 Bootstrap（如 PHYLIP `seqboot`，固定 $r=1.0$）的核心特征。它在每一次 bootstrap 迭代中都会执行完整的流水线：

1. 重采样（默认对行；在 `pgr clust boot --along row` 模式下对应对列重采样）
    - **默认采样比例 ($r$)**：`0.5, 0.6, ..., 1.4` (共10个尺度)。
    - **每个尺度的 bootstrap 次数**：默认 `1000` 次。
    - **总 bootstrap 次数**：`10 (scales) * 1000 (nboot) = 10000` 次。
2. 基于重采样后的数据重新计算 `dist`
    - `pvclust` 会根据 `method.dist`（如 correlation, abscor, uncentered）和 `use.cor`（缺失值处理策略）在重采样后的数据矩阵上重新计算距离。
    - `pgr` 也将复刻此行为，确保每次重采样后的距离计算与原始数据使用相同的度量标准。
3. 基于新的距离矩阵重新做 `hclust`
    - `pvclust` 的实现强依赖于层次聚类（调用 R 的 `hclust`），它通过 `method.hclust` 参数支持 `average/ward/single/complete` 等标准方法，但不支持非层次聚类（如 K-means、DBSCAN）。
    - `pgr clust boot` 也将遵循这一设计，主要支持层次聚类方法（与 `pgr clust hier` 对齐）。
4. 统计簇的出现频率
    - `pvclust` 不直接比较树的拓扑结构，而是将树分解为一组“簇”（Split/Cluster，即内部节点所包含的叶子集合）。
    - 对每棵 Bootstrap 树，将其所有内部节点转换为“成员 Pattern”（例如叶子索引的 0/1 向量或哈希）。
    - 检查原始树中的每个簇（Pattern）是否出现在 Bootstrap 树中，并累加计数。这避免了节点顺序或旋转带来的干扰，只关注“成员集合”的一致性。
5. 拟合得到最终值（BP/AU/SI）
    - 对于原始树中的每个簇，我们得到了一组数据 `(r, BP_r)`，即在不同采样比例 $r$ 下的 Bootstrap 出现频率。
    - 使用加权最小二乘法拟合曲线（模型：$z = -qnorm(BP_r) \approx v\sqrt{r} + c/\sqrt{r}$），源码见 [pvclust-internal.R:L350-L407](file:///c:/Users/wangq/Scripts/pgr/pvclust/R/pvclust-internal.R#L350-L407)。
    - 根据拟合得到的参数 $v$ 和 $c$，计算出最终的三类数值：
        - **AU (Approximately Unbiased)**：**推荐使用**。通过多尺度拟合修正了 BP 的偏差，更接近真实的 p-value。
        - **BP (Bootstrap Probability)**：传统 Bootstrap 值（对应 $r=1.0$），通常有偏差（偏保守）。
        - **SI (Selective Inference)**：选择性推断 p-value。
    - **注意**：如果用户仅指定了一个采样尺度（如 $r=1.0$），则无法进行曲线拟合，此时**无法计算 AU 值**，仅能输出 BP 值（即标准的 Bootstrap 支持度）。

因此 `clust boot` 不适合以“先算一次 PHYLIP 距离矩阵、再反复聚类”的方式实现；为了得到正确的 BP/AU/SI，必须在每次重采样后用同一套距离口径重新计算距离，再构树并统计簇的出现频率。

---

## 2. pvclust 源码结构与关键流程（仓库内 pvclust 目录）

### 2.1 导出 API（NAMESPACE）

`pvclust` R 包对外导出：

- `pvclust()`：主入口（见 [pvclust.R:L1-L63](file:///c:/Users/wangq/Scripts/pgr/pvclust/R/pvclust.R#L1-L63)）
- `msfit()`：多尺度曲线拟合（见 [pvclust-internal.R:L350-L407](file:///c:/Users/wangq/Scripts/pgr/pvclust/R/pvclust-internal.R#L350-L407)）
- `msplot()/seplot()/pvrect()/pvpick()`：诊断与筛选（见 `pvclust/man/*.Rd`）

### 2.2 主入口 pvclust()

`pvclust()` 的职责基本是：

1. 处理并行参数（`parallel` 可为 FALSE / TRUE / 整数 / cluster）
2. 进入 `pvclust.parallel()` 或 `pvclust.nonparallel()`（见 [pvclust.R:L1-L63](file:///c:/Users/wangq/Scripts/pgr/pvclust/R/pvclust.R#L1-L63)）

真正算法实现集中在 `pvclust-internal.R`：

- `pvclust.common.settings()`：计算原始距离与原始 hclust；规范化 `r`（见 [pvclust-internal.R:L3-L37](file:///c:/Users/wangq/Scripts/pgr/pvclust/R/pvclust-internal.R#L3-L37)）
- `boot.hclust()`：对每个 `r`、每次 bootstrap，重采样行、重算距离、重做 hclust，并统计 cluster 出现次数（见 [pvclust-internal.R:L223-L279](file:///c:/Users/wangq/Scripts/pgr/pvclust/R/pvclust-internal.R#L223-L279)）
- `pvclust.merge()`：把多尺度计数合并成 `edges.bp/edges.cnt`，并对每个 edge 调用 `msfit()` 得到 AU/BP/SI 与标准误等（见 [pvclust-internal.R:L281-L332](file:///c:/Users/wangq/Scripts/pgr/pvclust/R/pvclust-internal.R#L281-L332)）

### 2.3 “同一个簇”的判定：hc2split() 的 pattern

`pvclust` 不直接用“内部节点编号”来比较簇，而是把每个内部节点对应的成员集合编码成一个 pattern 字符串：

- `hc2split()` 返回：
  - `member`：每个内部节点的成员索引集合
  - `pattern`：每个内部节点的 0/1 向量拼接成字符串（作为簇 ID）

见 [pvclust-internal.R:L180-L214](file:///c:/Users/wangq/Scripts/pgr/pvclust/R/pvclust-internal.R#L180-L214)。

这意味着：在 `pgr` 侧实现时，最稳健的对齐方式也是“按 leaf-set 比较”，而不是依赖内部节点顺序。

---

## 3. 核心数学：msfit 的拟合模型（实现口径）

`msfit()` 输入：

- `bp`：某个 cluster 在多尺度 `r` 下的 BP（频率）
- `r`：相对样本量 `r = n'/n`
- `nboot`：每个尺度的 bootstrap 次数

关键步骤（见 [pvclust-internal.R:L350-L407](file:///c:/Users/wangq/Scripts/pgr/pvclust/R/pvclust-internal.R#L350-L407)）：

- 变换：`z = -qnorm(bp)`
- 加权最小二乘拟合：
  - 设计矩阵 `X = [sqrt(r), 1/sqrt(r)]`
  - `z ≈ v*sqrt(r) + c/sqrt(r)`（无截距）
- 由 `(v,c)` 计算：
  - `AU = pnorm(-(v - c))`
  - `BP = pnorm(-(v + c))`
  - `SI`：基于 selection probability 的修正版本
- 同时输出 `se.*`、`rss`、`df`、`pchi`

`pgr` 计划复刻该公式与阈值策略（例如 `eps=0.001`、`min.use=3`）以保持与 R 结果可对照。

---

## 4. 计划新增命令：pgr clust boot

### 4.1 定位

`pgr clust boot` 负责“给一棵层次聚类树的每个簇打分（BP/AU/SI）”。

- 与 `pgr clust hier` 的关系：`hier` 只构树；`boot` 在“数据可重采样”的前提下，评估树上簇的统计稳定性
- 与 `pgr clust eval` 的关系：`eval` 比较两份分区/内部指标；`boot` 产出的是**单棵树**的簇置信度（外部一致性并非必要）

### 4.2 使用方式（草案）

```bash
pgr clust boot [OPTIONS] <data.tsv>
```

**输入 `<data.tsv>`（建议默认）**：

- 第一行是列名（对象名，构树叶子名）
- 每一行是一次观测（bootstrap 重采样的单位）
- 可选第一列为观测 ID（通过 `--row-id` 指定）

### 4.3 关键参数（草案）

| 参数 | 取值 | 对应 pvclust | 说明 |
| :--- | :--- | :--- | :--- |
| `--dist` | `correlation/uncentered/abscor/euclidean` | `method.dist` | 距离口径（与 pvclust 同名） |
| `--use-cor` | `pairwise/complete/all` | `use.cor` | 处理缺失值的相关系数策略 |
| `--method` | `average/ward/complete/...` | `method.hclust` | 层次聚类 linkage |
| `--nboot` | int | `nboot` | 每个尺度的 bootstrap 次数 |
| `--r` | `0.5,1.4,0.1` | `r=seq(.5,1.4,.1)` | 多尺度相对样本量参数 |
| `--seed` | int | `iseed` | 随机种子 |
| `--quiet` | flag | `quiet` | 减少进度输出 |
| `-o/--outfile` | path | - | 输出 TSV（edge 指标表） |
| `--out-tree` | path | - | 输出带注释/标签的 Newick（可选） |

说明：

- `--other` 这个命名在 `clust eval` 中用于“另一份分区”；在 `clust boot` 中不需要第二份输入文件，因此不复用该参数名。
- `pvclust` 的 `weight/store/parallel` 先不做强绑定；`pgr` 侧更倾向于用 `--threads`（或继承全局线程设置）控制并行。

### 4.4 输出（草案）

输出建议为 TSV，每行一个 cluster（内部节点），字段对齐 R 包的 `x$edges`/`x$msfit`：

| 列名 | 含义 | 对应 pvclust |
| :--- | :--- | :--- |
| `edge` | 内部节点编号（稳定编号策略见下） | `row.names(x$edges)` |
| `size` | 该簇叶子数 | `length(member[[i]])` |
| `bp` | BP | `x$edges$bp` |
| `au` | AU | `x$edges$au` |
| `si` | SI | `x$edges$si` |
| `se.bp` / `se.au` / `se.si` | 标准误 | `x$edges$se.*` |
| `v` / `c` | 拟合参数 | `x$edges$v/c` |
| `df` / `rss` / `pchi` | 拟合诊断 | `x$msfit[[i]]$df/rss/pchi` |

**内部节点编号策略（建议）**：

- 不依赖构树时的内部节点顺序
- 用 leaf-set 的哈希（或排序后的叶子名拼接）作为稳定 ID
- 同时输出一个“显示用序号”便于用户筛选与可视化

### 4.5 可视化与后处理（建议工作流）

1. 生成 `edge` 指标表：

```bash
pgr clust boot data.tsv --dist correlation --method average --nboot 1000 -o boot.tsv
```

2. 按阈值挑选显著簇（pvclust 的 `pvrect/pvpick` 思路，见 [pvpick.Rd](file:///c:/Users/wangq/Scripts/pgr/pvclust/man/pvpick.Rd)）：

- `au >= 0.95` 常用作强支持阈值（也可用 SI）

3. 将挑选结果映射回树（后续可设计 `pgr nwk label --from boot.tsv --field au` 一类工具）。

---

## 5. 与现有 pgr 命令的组合关系（建议）

- `pgr dist vector → pgr mat to-phylip → pgr clust hier`：现有“向量→距离→树”的通路
- `pgr clust boot`：需要“可重采样的原始观测矩阵”，因此更适合直接吃 `data.tsv`，内部自算距离与树
- `pgr nwk cut`：在 `boot` 给出簇置信度后，再做阈值切割与导出分区
- `pgr clust eval`：对不同切割/不同算法产生的分区做一致性比较（用 `--other`）

---

## 6. 参考与对照

- pvclust R 包版本：2.2-0（见 [DESCRIPTION](file:///c:/Users/wangq/Scripts/pgr/pvclust/DESCRIPTION)）
- 核心实现文件：
  - [pvclust.R](file:///c:/Users/wangq/Scripts/pgr/pvclust/R/pvclust.R)
  - [pvclust-internal.R](file:///c:/Users/wangq/Scripts/pgr/pvclust/R/pvclust-internal.R)
- 关键函数定位：
  - `pvclust()`：[pvclust.R:L1-L63](file:///c:/Users/wangq/Scripts/pgr/pvclust/R/pvclust.R#L1-L63)
  - `boot.hclust()`：[pvclust-internal.R:L223-L279](file:///c:/Users/wangq/Scripts/pgr/pvclust/R/pvclust-internal.R#L223-L279)
  - `pvclust.merge()`：[pvclust-internal.R:L281-L332](file:///c:/Users/wangq/Scripts/pgr/pvclust/R/pvclust-internal.R#L281-L332)
  - `msfit()`：[pvclust-internal.R:L350-L407](file:///c:/Users/wangq/Scripts/pgr/pvclust/R/pvclust-internal.R#L350-L407)
  - `seplot()`：[pvclust-internal.R:L458-L481](file:///c:/Users/wangq/Scripts/pgr/pvclust/R/pvclust-internal.R#L458-L481)
