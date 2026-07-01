# FastK: 高保真 k-mer 计数器

## 简介
FastK 是一个专为处理高质量 DNA 组装数据（如 Illumina 或 PacBio HiFi 数据）而优化的 k-mer 计数工具。它采用了一种新颖的基于 Minimizer 的分发方案，能够利用磁盘存储来处理任意大小的数据集。

其核心优势在于能够直接生成序列的 k-mer 计数档案（Profile），并且在处理低错误率（1% 或更低）数据时，通过两阶段的“先 Super-mer 后加权 k-mer”排序策略，实现了极高的速度。

## 核心概念 (Key Concepts)

### 1. K-mer 与 Canonical K-mer
*   **K-mer**: 长度为 k 的 DNA 序列片段（FastK 支持 k >= 5，默认为 40）。
*   **Canonical K-mer (规范 K-mer)**: 为了处理测序读取方向未知的问题，FastK 将一个 k-mer 及其反向互补序列（Watson-Crick complement）视为同一个 k-mer。在两者中，字典序较小的那个被称为“规范 k-mer”。FastK 的统计表（Table）只记录规范 k-mer。

### 2. Super-mer (超 k-mer)
*   背景: 在传统的 k-mer 计数中，如果一条读长（Read）有 150bp，K=40，那么它包含 111 个 k-mer。如果直接把这 111 个 k-mer 全部拆散单独存下来排序，数据量会膨胀得非常大（111 * 40 bytes），而且丢失了它们原本相邻的信息。
*   定义: Super-mer 是原序列中一段连续的、共享同一个 Minimizer 的子序列。
    *   当 FastK 扫描序列时，只要相邻的 k-mer 的 Minimizer 没有变（或者新出现的 Minimizer 依然是当前这一段里最小的），就把它们连在一起，形成一个 Super-mer。
    *   直到 Minimizer 发生变化，不得不划分到另一个桶时，才切断。
*   例子:
    *   假设序列是 `ABCDE`，K=3。
    *   包含的 k-mer 有: `ABC`, `BCD`, `CDE`。
    *   如果这三个 k-mer 算出来的 Minimizer 都是 `B`。
    *   传统做法: 存 3 个条目 `ABC`, `BCD`, `CDE`。
    *   FastK 做法: 存 1 个条目 `ABCDE` (Super-mer)。
*   优势:
    *   压缩数据: 极大减少了需要写入磁盘和排序的条目数量（通常减少 10-50 倍）。
    *   提升速度: 排序 1 个长条目比排序 10 个短条目要快得多。在最后阶段，程序再从 Super-mer 中还原出具体的 k-mer 进行计数。

### 3. Minimizer (最小标识符)
*   定义: Minimizer 是 k-mer 序列内部的一个特定的、较短的子序列（m-mer，m < k）。通常选择字典序最小的那个 m-mer。
*   分发策略: FastK 根据 k-mer 中包含的 Minimizer 来决定将其放入哪个存储桶（Bucket）。
*   核心作用:
    *   保证同类归并: 如果两个 k-mer 的序列完全相同，它们必然拥有相同的 Minimizer。因此，它们一定会被分发到同一个桶中。这就像把所有姓“李”的书都扔进“L”号箱子，不管这书是从哪里来的，只要是同一本书，它一定在“L”箱子里。
    *   并行独立性: 因为所有相同的 k-mer 都在同一个桶里，我们在统计“L”箱子时，完全不需要去问“Z”箱子有没有漏掉的。这意味着不同的桶可以完全独立地由不同线程并行处理，互不干扰。
    *   内存控制: 无论数据集多大（比如 1TB），我们都可以通过增加桶的数量（比如分成 1000 个桶），让每个桶只有 1GB，从而可以轻松读入内存进行快速排序。

> **注意：与 `pgr` (Sketching) 的区别**
> FastK 中的 Minimizer 与 `pgr/src/libs/hash.rs` 中用于 MinHash Sketch 的 Minimizer 虽然原理相同，但用途截然不同：
> *   FastK (无损路由): 针对每一个 k-mer，找到其内部的 m-mer 标签，用于决定去哪个桶。不丢弃任何数据，目的是全量统计。
> *   PGR Sketch (有损采样): 在长序列的滑动窗口中选出一个 k-mer 代表该窗口。丢弃绝大部分数据，目的是生成稀疏指纹（Signature）用于快速比对。

### 4. 输出文件格式
FastK 生成以下几种核心文件：
*   **直方图 (.hist)**: 记录了每种频次（1次, 2次...）出现的 k-mer 数量。最大计数限制为 32,767。
*   **K-mer 表 (.ktab)**: 一个排序的列表，包含数据集中所有（或满足特定阈值）的规范 k-mer 及其计数。
*   **档案 (.prof)**: (可选) 数据集中每条序列的 k-mer 计数概况。这是一个压缩格式，不仅记录了 k-mer 的出现，还按原序列顺序保留了位置信息。

## 处理流程 (Processing Workflow)

FastK 的内部处理逻辑主要分为四个阶段（Phase）：

### 第一阶段：分片与分发 (Phase 1: Partitioning)
*   输入扫描: 程序首先扫描输入数据集的前 1GB 数据。
*   方案确定: 基于这部分数据，计算 Minimizer 分布，确定如何将 Super-mer 均衡地分发到临时桶中。
*   全量分发: 扫描整个数据集，计算 Super-mer，并根据 Minimizer 方案将它们写入磁盘上的不同临时文件（Buckets）。
*   *代码对应*: `split.c`

### 第二阶段：排序与计数 (Phase 2: Sorting & Counting)
*   并行处理: 对每个临时桶并行执行操作。
*   两级排序:
    1.  首先对桶内的 Super-mer 进行排序。
    2.  然后对 Super-mer 包含的 k-mer 进行加权排序（Weighted Sort）。
*   统计生成: 在排序过程中累积 k-mer 的频次直方图。
*   *代码对应*: `count.c`, `LSDsort.c`, `MSDsort.c`

### 第三阶段：表格合并 (Phase 3: Table Merging)
*   合并: 将第二阶段各线程生成的排序好的 k-mer 片段，按字典序合并成一个单一的、全局有序的 `.ktab` 文件。
*   *代码对应*: `table.c`

### 第四阶段：档案合并 (Phase 4: Profile Merging)
*   合并: (仅当启用 `-p` 选项时) 将分布在各处的 Profile 片段合并成最终的 `.prof` 文件。
*   压缩: Profile 数据采用特定的编码进行压缩（约 4.7 bits/base），以节省空间。
*   *代码对应*: `merge.c`

---
*参考来源: [FastK GitHub Repository](https://github.com/thegenemyers/FASTK)*
