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
不是核心。pgr 当前只输出 PAF，缺 impg 的 BED 选项与批查能力，这是错位。

V1 补齐两件事：

```bash
# 1. BED 输出（-o bed，最 pipe 友好）——3 列：name start end
pgr paf query cohort.paf.idx chr1:1000-5000 --transitive -o bed

# 2. 多 region 批查（impg -b regions.bed）——单 region 限制是 pgr 独有的
pgr paf query cohort.paf.idx -b regions.bed --transitive -o bed
```

> **默认输出决策（2026-06-28 修订）**：pgr V1 的**默认输出保持 PAF**，BED 通过 `-o bed` 可选。
> 这与 impg 的"BED 默认"不同，理由：(1) pgr 既有 23 个集成测试已断言 PAF 输出，改默认会破坏；
> (2) PAF 含 CIGAR/gi/bi 标签，对需要完整比对记录的场景更直接，BED 三列是坐标投影的"轻量产物"，
> 用 `-o bed` 显式切换语义更清晰；(3) impg 的"BED 默认"是历史选择，pgr 作为新工具可独立决策。
> 详见 [[pairwise-selection.md]] §5 变更日志。

### 3.2 为何 V1 不做 fasta/maf

上一轮把 impg 的**可选项**当成了 pgr 的**下一步**，是误判。重新审视：

1. **impg 的 MAF 是可选项**，不是核心——核心是坐标投影 + 传递闭包
2. **pgr 已有 `pgr fas consensus`**——用户要 MSA 时，`pgr paf query -o bed` → `pgr fa range`
   → `pgr fas consensus` 的 pipe 路径已通，不需要在 query 层重复 POA
3. **POA MSA 是图构建层的产物**——按 [[paf-route.md]] §2.4，图遍历和 MSA 是正交步骤，
   不应耦合进 query
4. **V1 的真正用户场景**——"给我 chr1:1000-5000 在 cohort 里的所有同源区段"，输出 BED 即可

`pairwise-selection.md` 变更日志也印证：06-27 曾"BED 成为默认输出"，06-28 又"BED/TSV 删除，
只输出 PAF"。本次（06-28 二次修订）的结论是——**PAF 保持默认，BED 通过 `-o bed` 可选**，
两者都保留，由用户按场景选择。BED 三列（`name start end`）是坐标投影的轻量产物，PAF 含完整
CIGAR 适合需要比对记录的场景。

### 3.3 新增代码

| # | 任务 | 文件 | 行数 |
|---|------|------|:--:|
| 1 | `-o bed` 输出（3 列，复用现有 results） | `cmd_pgr/paf/query.rs` | ~15 |
| 2 | `-o paf`/`-o bed` 分发逻辑 + `--bed-regions`/`-b` 参数 | `cmd_pgr/paf/query.rs` | ~25 |
| 3 | BED 文件解析（多 region 批查） | `cmd_pgr/paf/query.rs` | ~10 |
| 4 | 集成测试（bed 输出 + 批查，6 个新测试） | `tests/cli_paf.rs` | ~15 |
| **总计** | | | **~65** |

### 3.4 不做的

- ❌ `-o fasta`（推迟到 V2，需 `-f` 序列文件）
- ❌ `-o maf`（推迟到 V3，POA 是图构建层产物）
- ❌ GFA/VCF 输出（推迟到 V4，需完整 graph engine）
- ❌ gfasort/gfaffix（pgr 不做 GFA，不需要）
- ❌ minigraph 的 `gfa_t` 在 Rust 中的对应实现

---
## 4. V1/V2/V3/V4/V5 路线

按 impg 各输出格式的依赖链与核心性递进：

| 阶段 | 内容 | 对应 | 代码量 |
|------|------|------|:---:|
| **V1**（当前缺失） | `-o paf`（默认）+ `-o bed`（轻量坐标）+ `-b regions.bed` 批查 | impg 默认 `-o bed` + `-b`（pgr 选 PAF 默认，见 §3.1） | ~65 |
| **V2** | `-o fasta`（未比对序列，需 `-f`）| impg `-o fasta` | ~60 |
| **V3** | `-o maf`（POA MSA，需 `-f`）+ `-o fasta-aln` | impg `-o maf`/`-o fasta-aln` | ~150 |
| **V4a** | 粗全局 GFA（`pgr paf graph -o gfa --min-var-len 100`，minigraph 风格）| minigraph `ggen` | 待评估 |
| **V4b** | 区域精细 GFA（`pgr paf query -o gfa -r region`，impg 风格）| impg `query -o gfa` | 待评估 |
| **V5** | 区域 GFA → MAF/VCF（精细分析输出）+ EKG @tags | impg `-o maf`/`-o vcf` | 待评估 |

### 4.1 为何坐标类输出是 V1 核心

impg 的 11 种输出格式按"是否需要序列文件"分两类：

| 类别 | 格式 | 需 `-f` | 用途 |
|------|------|:---:|------|
| **坐标类**（核心） | `bed`/`bedpe`/`paf` | 否 | "哪些序列的哪些区段同源"——喂给下游工具 |
| **序列类**（可选） | `fasta` | 是 | 提取未比对序列 |
| **MSA 类**（可选） | `maf`/`fasta-aln` | 是 | POA 多序列比对 |
| **图类**（可选） | `gfa`/`vcf`/`gbwt` | 是 | 物化图，需完整 graph engine |

pgr V1 同时提供 PAF（默认，含 CIGAR/gi/bi 完整比对记录）与 BED（`-o bed`，3 列轻量坐标）。
PAF 适合需要完整比对记录的场景，BED 三列（`name start end`）是坐标投影的轻量产物，最 pipe
友好——喂给 `pgr fa range` 提取序列。两者由 `-o` 切换，详见 §3.1 默认输出决策。

### 4.2 为何 fasta/maf 后移

- **fasta 推到 V2**：impg 也是可选项，需 `-f` 序列文件，依赖 noodles_fasta 索引
- **maf 推到 V3**：POA MSA 是图构建层产物，按 [[paf-route.md]] §2.4 不应耦合进 query 的核心路径。
  `pgr fas consensus` 已提供成熟 POA 后端，`bed → fa range → fas consensus` 的 pipe 路径已通

### 4.3 V4 的能力跃迁：两段式 GFA

V4 采用**两段式 GFA**——粗全局 + 区域精细，混合 minigraph 和 impg 各自所长：

| 工具 | 全局粗 GFA | 区域精细 GFA | 区域 → MSA/VCF |
|------|:---:|:---:|:---:|
| **minigraph** | ✅（≥100bp SV，rGFA）| ❌（不做）| ❌（小变体用标准工具）|
| **impg** | ❌（不物化全局图）| ✅（`query -o gfa`）| ✅（`query -o maf/vcf`）|
| **pgr V4** | ✅（V4a）| ✅（V4b）| ✅（V5）|

#### 4.3.1 V4a：粗全局 GFA（minigraph 风格）

- **输入**：PAF 索引（已有，无需重新比对——这是 pgr 相对 minigraph 的核心优势）
- **过滤**：`--min-var-len 100`（对齐 minigraph，只保留 ≥100bp SV）
- **输出**：rGFA，含 SN/SO/SR tag（稳定坐标系，见 [[minigraph.md]] §3.2）
- **用途**：可视化（Bandage/odgi）、SV 概览、作为后续 query 的坐标锚
- **数据源**：PAF 索引的显式投影，不引入新的真实源
- **算法骨架**：直接复用 seqwish 的 6 阶段流程（spanning tree → BFS discovery → DSU union-find
  → compact → links → GFA），详见 [[seqwish.md]] §6.2。pgr 输入是 PAF，与 seqwish 一致，
  相对 minigraph（需自跑 minimizer chaining）更天然适配。

#### 4.3.2 V4b：区域精细 GFA（impg 风格）

- **输入**：PAF 索引 + 用户指定 region（或粗 GFA 上定位的区段）
- **流程**：BFS 传递闭包 → 提取序列 → POA/外部 aligner → 局部 GFA
- **输出**：局部 GFA（含 base-level 变异）
- **用途**：特定基因座的精细分析
- **对应**：impg `query -o gfa`
- **完整性保证**：借鉴 seqwish 的 phase 1b orphan recovery 思路——BFS 发现的等价类可能只覆盖部分
  序列，按序列查区间树补漏直到收敛，保证局部 GFA 的等价类完整性。详见 [[seqwish.md]] §3.3、§6.3。

#### 4.3.3 两段衔接

粗 GFA 提供"地图"（哪里有大 SV），用户在粗 GFA 上定位感兴趣的区段，再调用精细 GFA 做碱基级
分析。这比 minigraph（只有粗）和 impg（只有精）更完整。V4a 与 V4b 可独立实现，互不依赖。

#### 4.3.4 局部 GFA 不合并回全局

**设计决策**：V4b 产出的多个局部精细 GFA **独立存在，不合并回全局 GFA**。这是 pgr 与
minigraph 的边界——保持"粗全局 GFA 是不可变投影"的语义，不滑向"可变全局图"。

**为什么不合并**：

1. **哲学一致性** — [[paf-route.md]] §0.1 原则 1 规定"粗 GFA 作为可选投影"。合并局部 GFA
   意味着全局图可变，等于重新实现 minigraph 的 `gfa_t` + augment（[gfa-aug.c](file:///Volumes/ExtHome/Scripts/pgr/minigraph-master/gfa-aug.c)），
   超出 pgr 的边界。
2. **技术难度高** — GFA 合并需解决四个挑战：坐标对齐、边界衔接、节点 ID 冲突、路径完整性。
   无成熟先例，工作量大。
3. **替代方案更优** — 需要"全局精细视图"的场景，用 V4a 粗全局 GFA + 各 region 局部 GFA 配合
   查看即可；需要"合并输出"的场景走 V5 的 VCF/MAF（天然可 concat），不走 GFA 合并。

**按输出格式分治**：

| 输出 | 合并机制 | 说明 |
|------|----------|------|
| **VCF** | `bcftools concat`（按坐标）| VCF 天然按 region 独立 calling，concat 即可 |
| **MAF** | 按坐标拼接 | MAF 本就是分块的，按 region 排序拼接 |
| **GFA** | **不合并** | 需要全局精细 GFA 的用户应直接用 minigraph |

**V4b 的边界**：每个局部 GFA 是独立的、自包含的产物，对应一个用户指定的 region。多个 region
的局部 GFA 之间无依赖、无引用，不形成全局图。

#### 4.3.5 V4 必须引入粗框架过滤（对齐 minigraph）

V4a 物化粗 GFA 时必须加 `--min-var-len`（默认 100）过滤，只把长度差 ≥ 阈值的变异变成图节点。
理由见 [[minigraph.md]] §4.1 引用的论文 L601-609：

1. 不加过滤的图会爆炸——"millions of short segments"
2. minigraph 的 minimizer 索引会失败（pgr 用区间树，但图遍历仍会退化）
3. 小变体用标准方法（VCF/MAF）更易分析
4. 无算法能为数百基因组构建全变体图

**两种正交的过滤维度**（pgr V4a 需同时考虑）：

- **变异长度**（`--min-var-len`，minigraph 风格）— 过滤 < 100bp 的小变体，避免碱基级碎片
- **重复拷贝数**（`--repeat-max` / `--min-repeat-dist`，seqwish 风格）— 限制同一序列在图同一位置
  的拷贝数，避免高拷贝重复把图吹爆。详见 [[seqwish.md]] §3.5 的 `write_graph_chunk` 实现。

两者维度不同：前者按变异长度过滤，后者按重复拷贝数过滤，pgr V4a 可同时启用。

**V4b 不受此约束**——区域精细 GFA 只处理单个区段，序列数通常 < 50，不会爆炸。

#### 4.3.6 边冲突问题

minigraph 用 `mg_ggsimple` 的过滤逻辑（[ggsimple.c:213](file:///Volumes/ExtHome/Scripts/pgr/minigraph-master/ggsimple.c#L213)）
解决"同一区间多个比对"——只保留最优比对，避免图变成一团乱麻。

pgr 从 PAF 索引物化粗 GFA（V4a）时同样面临：
- 同一 query 区间有多个 target 比对（paralog）
- 需要选最优比对（by identity / by length）
- 这就是 impg 的 `--sparsify` 和 minigraph 的"primary chain > 20kb"过滤

**pgr 的解决方案**：复用查询层的 `--min-identity` 等参数，在物化时做同等过滤。不需要新算法，
只是把查询层的过滤逻辑应用到全局物化。

#### 4.3.7 关键区分：查询层 vs 图构建层

粗框架是**图构建层**（V4a）的过滤，不是查询层。V1-V3 查询层全量返回同源区段
（[[paf-route.md]] §2.3），由用户用 `--merge-distance` 等参数控制粗细。V4 物化时才在
graph engine 内部做 `min_var_len` 过滤。这两个层次不能混淆——查询层全量是"让用户决定粗细"，
图构建层粗框架是"避免图爆炸"。

