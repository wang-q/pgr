# 方向 D：图构建层 — 设计文档

本文档定义 pgr 从"找到同源"（查询层）到"构建图"（图构建层）的设计方案。 综合三个参考来源：

- **impg-0.4.1**：`graph.rs`（1207 行）— POA → GFA → gfaffix → gfasort
- **minigraph**：`map-algo.c`（502 行）、`ggen.c`（480 行）— k-mer chains → CIGAR → rGFA
- **EKG**：`docs/gfa.md` §6-7 — 变异图哲学 + rGFA anchor 系统

---
## 1. 三种图构建路线

### 1.1 impg 路线：POA → GFA

```
BFS 传递闭包 → 提取序列 → POA (spoa_rs FFI) → GFA → gfaffix → gfasort
```

- **核心算法**：SIMD C++ SPOA + 最长序列优先喂入
- **后处理**：unchop（合并线性节点）+ gfaffix（标准化）+ Ygs 排序
- **输出**：标准 GFA 1.0 with P lines
- **外部依赖**：`spoa_rs`（C++ FFI）、`gfaffix`（外部 binary）、`gfasort`（Rust crate + binary）
- **优势**：全自动，无参数
- **劣势**：依赖重，POA 碎片化严重（per-SNP 分叉）

### 1.2 minigraph 路线：chain → CIGAR → rGFA

```
k-mer seeds → linear chain (DP) → gchain (graph chain) → rGFA
```

- **核心算法**：k-mer 种子收集 → minimizer → 线性链 DP（`mg_lchain_dp`， `gchain1.c:62-134`）→
  图形链合并（`mg_gchain1_dp`）→ path → seq（`mg_path2seq`，`ggen.c:148-268`）

- **minigraph 的 `mg_path2seq`**（`ggen.c:148-268`）是理解 minigraph 如何从 chains 构建线性 MSA
  的核心：

  ```
  while (1) {
  1. 找 rs ≤ r ≤ re 的得分最高 chain
  2. 有 → 写 ref 片段 (g->seg[v].del[voff[0]..voff[1]])
         写 query 序列
         前进 v
  3. 无 → 写剩余 ref 片段，结束
  }
  ```
这个算法本质是"在参考序列上按位置依次插入 query 序列段"—— 不是 all-vs-all MSA，而是
  reference-guided 的线性比对。

- **后处理**：gfa_sort_ref_arc（`--call`）、CIGAR 生成（`mg_gchain_cigar`， `map-algo.c:475-478`）

- **输出**：rGFA 1.0（reference-anchored，`SN:Z:`/`SO:i:`/`SR:i:` tags）

- **外部依赖**：零（纯 C，自包含 `gfa_t` 结构体）

- **优势**：极快（k-mer 免比对）、线性（无分支膨胀）、rGFA 标准化

- **劣势**：依赖 k-mer 参数调优，不适合稀疏比对（如远缘物种间）

### 1.3 pgr 可选路线：PAF → POA → MSA

```
MAF → PAF → PafIndex → BFS 传递闭包 → 提取序列 → POA (纯 Rust) → MSA
```

- **pgr 独特起点**：已有 PAF index + BFS 查询，不需要重新做比对
- **pgr 的 POA**（`libs/poa/poa.rs`）是纯 Rust，无外部依赖
- **输出**：block FASTA/MAF（多序列比对），不是 GFA
- **优势**：零新依赖、代码量少（~100 行）、立即可用
- **劣势**：不输出 GFA（不能直接入 vg/odgi 管道）

---
## 2. 为什么 V1 不做 GFA

三种路线中，minigraph 的 chain → GFA 路线在 pgr 中对应的输入是 **PAF + CIGAR**（PAF 本身就是 chain
的表示）。从 PAF 到 GFA 的直接转换 在理论上是可行的，但需要：

1. **节点定义**：PAF 中每个 CIGAR `M`/`=` block 对应一个 node，还是按 reference 坐标切 node？
   ——minigraph 走前者（CIGAR → segment），impg 走后者（POA 自动决定）
2. **边定义**：PAF 的 strand 信息需要映射到 GFA 的四种边方向
3. **路径定义**：GFA P-line 需要 reference 坐标锚定——pgr 的 PAF target 坐标天然提供 reference 坐标
   （0-based forward-strand），但需要 minigraph 级别的 `SN:Z:`/`SO:i:`/`SR:i:` tag 管理
4. **节点去重**：minigraph 用 `gfa_aux_update_cv` 和 `gfa_sort_ref_arc` 做后处理，pgr
   没有这些基础设施
5. **自环/重复**：minigraph 有 `mg_gchain_set_parent` 和子关系处理 （`map-algo.c:471`），pgr 的 PAF
   来自 multiz MAF，不涉及自环

这些都需要在 Rust 中重建 minigraph 的 `gfa_t` 数据结构（`gfa-priv.h`）， 工作量远超过 100 行。

**结论**：V1 走 PAF → POA → MSA 路线，因为：

- MSA 是 GFA 的前置步骤——先验证 POA 在泛基因组数据上的质量
- MSA 可以用 pgr 已有的 `poa.msa()` 一行调用完成
- 不需要定义 node/edge/path 语义
- 不需要外部依赖

**GFA 推迟到 V3**，届时评估 pgr 是否需要引入 rGFA 标准。

---
## 3. 最小可行实现（V1）：`pgr paf msa`

### 3.1 功能

```bash
# 从 PAF 索引直接生成 MSA
pgr paf msa cohort.paf.idx chr1:1000-5000 -f genomes.fa --transitive -o out.fas
```

内部流程：

```
1. 加载 PAF 索引（.paf.idx 或 .paf）
2. BFS 传递闭包（或单跳）
3. 从 BFS 结果收集所有 query 坐标区间
4. 从 FASTA 提取序列（按坐标切片）
5. POA 多序列比对
6. 输出 block FASTA
```

### 3.2 新增代码

| # | 任务 | 文件 | 行数 | |---|------|------|:--:| | 1 | `collect_intervals()` |
`libs/paf/graph.rs` | ~30 || 2 | `reverse_complement()` | `libs/paf/graph.rs` | ~20 || 3 |
FASTA 序列提取（noodles_fasta） | `cmd_pgr/paf/msa.rs` | ~30 || 4 | POA 调用 + MSA 输出 |
`cmd_pgr/paf/msa.rs` | ~40 || 5 | 集成测试 | `tests/cli_paf.rs` | ~3 || **总计**| | | **~120**|

### 3.3 不做的

- ❌ GFA 输出（推迟到 V3）
- ❌ gfasort/gfaffix（pgr 不做 GFA，不需要）
- ❌ Banded DP post-processing
- ❌ Rayon 并行化
- ❌ minigraph 的 `gfa_t` 在 Rust 中的对应实现

---
## 4. V2/V3 展望

| 阶段 | 内容                                                                 |
|------|----------------------------------------------------------------------|
| V2   | MSA 后处理（裁剪、extract consensus、convert to MAF）                |
| V3   | 评估后决定是否引入 rGFA（参考 minigraph 的 `gfa_t` + `mg_path2seq`） |
| V4   | EKG 的 @tags（`@node_length`、`@align_length`、`IntSpan`）           |

