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
- **输出**：未比对 FASTA / MAF（POA MSA），不是 GFA
- **优势**：零新依赖、代码量少（V1 ~60 行 / V3 ~150 行）、立即可用
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

**GFA 推迟到 V4**，届时评估 pgr 是否需要引入 rGFA 标准。

---
## 3. 最小可行实现（V1）：`pgr paf query -o bed` + 批查

### 3.1 功能

V1 的核心是**坐标输出**，不是 MSA。对照 impg 源码后修订——impg `query` 默认输出 `bed`
（[main.rs#L4873](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L4873)），
README L20-23 明确："It outputs BED / BEDPE / PAF — ready to feed FASTA extraction, multiple
sequence alignment... can also emit GFA directly"。"can also" 表明 GFA/MAF 是**可选附加**，
不是核心。pgr 当前只输出 PAF，缺 impg 的默认 BED，这是错位。

V1 补齐两件事：

```bash
# 1. BED 输出（impg 默认，最 pipe 友好）——3 列：name start end
pgr paf query cohort.paf.idx chr1:1000-5000 --transitive -o bed

# 2. 多 region 批查（impg -b regions.bed）——单 region 限制是 pgr 独有的
pgr paf query cohort.paf.idx -b regions.bed --transitive -o bed
```

PAF 仍作为 `-o paf` 可选输出（含 CIGAR，适合需要完整比对记录的场景）。

### 3.2 为何 V1 不做 fasta/maf

上一轮把 impg 的**可选项**当成了 pgr 的**下一步**，是误判。重新审视：

1. **impg 的 MAF 是可选项**，不是核心——核心是坐标投影 + 传递闭包
2. **pgr 已有 `pgr fas consensus`**——用户要 MSA 时，`pgr paf query -o bed` → `pgr fa range`
   → `pgr fas consensus` 的 pipe 路径已通，不需要在 query 层重复 POA
3. **POA MSA 是图构建层的产物**——按 [[paf-route.md]] §2.4，图遍历和 MSA 是正交步骤，
   不应耦合进 query
4. **V1 的真正用户场景**——"给我 chr1:1000-5000 在 cohort 里的所有同源区段"，输出 BED 即可

`pairwise-selection.md` 变更日志也印证：06-27 曾"BED 成为默认输出"，06-28 又"BED/TSV 删除，
只输出 PAF"。这个回退是错的——把坐标查询和比对记录混为一谈了。BED 三列（`name start end`）
才是 impg 的默认，最 pipe 友好。

### 3.3 新增代码

| # | 任务 | 文件 | 行数 |
|---|------|------|:--:|
| 1 | `-o bed` 输出（3 列，复用现有 results） | `cmd_pgr/paf/query.rs` | ~15 |
| 2 | `-o paf`/`-o bed` 分发逻辑 + `--bed-regions`/`-b` 参数 | `cmd_pgr/paf/query.rs` | ~25 |
| 3 | BED 文件解析（多 region 批查） | `cmd_pgr/paf/query.rs` | ~10 |
| 4 | 集成测试（bed 输出 + 批查） | `tests/cli_paf.rs` | ~10 |
| **总计** | | | **~60** |

### 3.4 不做的

- ❌ `-o fasta`（推迟到 V2，需 `-f` 序列文件）
- ❌ `-o maf`（推迟到 V3，POA 是图构建层产物）
- ❌ GFA/VCF 输出（推迟到 V4，需完整 graph engine）
- ❌ gfasort/gfaffix（pgr 不做 GFA，不需要）
- ❌ minigraph 的 `gfa_t` 在 Rust 中的对应实现

---
## 4. V1/V2/V3/V4/V5 路线

按 impg 各输出格式的依赖链与核心性递进：

| 阶段 | 内容 | 对应 impg | 代码量 |
|------|------|-----------|:---:|
| **V1**（当前缺失） | `-o bed`（默认）+ `-o paf`（完整记录）+ `-b regions.bed` 批查 | impg 默认 `-o bed` + `-b` | ~60 |
| **V2** | `-o fasta`（未比对序列，需 `-f`）| impg `-o fasta` | ~60 |
| **V3** | `-o maf`（POA MSA，需 `-f`）+ `-o fasta-aln` | impg `-o maf`/`-o fasta-aln` | ~150 |
| **V4** | GFA/VCF 物化评估（参考 minigraph `gfa_t` + `mg_path2seq`） | impg `-o gfa`/`-o vcf` | 待评估 |
| **V5** | EKG @tags（`@node_length`、`@align_length`、`IntSpan`） | — | 待评估 |

### 4.1 为何 BED 是 V1 核心

impg 的 11 种输出格式按"是否需要序列文件"分两类：

| 类别 | 格式 | 需 `-f` | 用途 |
|------|------|:---:|------|
| **坐标类**（默认） | `bed`/`bedpe`/`paf` | 否 | "哪些序列的哪些区段同源"——喂给下游工具 |
| **序列类**（可选） | `fasta` | 是 | 提取未比对序列 |
| **MSA 类**（可选） | `maf`/`fasta-aln` | 是 | POA 多序列比对 |
| **图类**（可选） | `gfa`/`vcf`/`gbwt` | 是 | 物化图，需完整 graph engine |

pgr 当前 `paf query` 只输出 PAF，没有 BED。PAF 是**完整比对记录**（含 CIGAR），对"我只想知道
哪些区段同源"的用户是过度输出。BED 三列才是 impg 的默认，最 pipe 友好。

### 4.2 为何 fasta/maf 后移

- **fasta 推到 V2**：impg 也是可选项，需 `-f` 序列文件，依赖 noodles_fasta 索引
- **maf 推到 V3**：POA MSA 是图构建层产物，按 [[paf-route.md]] §2.4 不应耦合进 query 的核心路径。
  `pgr fas consensus` 已提供成熟 POA 后端，`bed → fa range → fas consensus` 的 pipe 路径已通

### 4.3 V4 的能力跃迁

GFA/VCF 需要完整 graph engine（impg 的 `dispatch_gfa_engine` + seqwish + crush + gfaffix +
gfasort），是真正的能力跃迁，不是格式化变种。V4 需要评估是否引入 rGFA 标准（参考
minigraph 的 `gfa_t` + `mg_path2seq`，见 [[minigraph.md]]）。

