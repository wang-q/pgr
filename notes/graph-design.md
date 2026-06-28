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
## 3. 最小可行实现（V1）：`pgr paf query` + `pgr paf to-bed` + 批查

### 3.1 功能

V1 的核心是**坐标输出**，不是 MSA。对照 impg 源码后修订——impg `query` 默认输出 `bed`
（[main.rs#L4873](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L4873)），
README L20-23 明确："It outputs BED / BEDPE / PAF — ready to feed FASTA extraction, multiple
sequence alignment... can also emit GFA directly"。"can also" 表明 GFA/MAF 是**可选附加**，
不是核心。pgr 当前只输出 PAF，缺 impg 的 BED 选项与批查能力，这是错位。

V1 补齐两件事：

```bash
# 1. BED 输出（pgr paf to-bed，最 pipe 友好）——3 列：name start end
pgr paf to-bed cohort.paf.idx chr1:1000-5000 --transitive

# 2. 多 region 批查（impg -b regions.bed）——单 region 限制是 pgr 独有的
pgr paf to-bed cohort.paf.idx -b regions.bed --transitive
```

> **默认输出决策（2026-06-28 修订，06-28 二次修订）**：pgr V1 的**`pgr paf query` 默认输出保持
> PAF**，BED 通过独立子命令 `pgr paf to-bed` 提供。这与 impg 的"BED 默认"不同，理由：(1) pgr 既有
> 集成测试已断言 PAF 输出，改默认会破坏；(2) PAF 含 CIGAR/gi/bi 标签，对需要完整比对记录的场景更
> 直接，BED 三列是坐标投影的"轻量产物"，用独立子命令语义更清晰；(3) impg 的"BED 默认"是历史选择，
> pgr 作为新工具可独立决策；(4) 符合 pgr CLI 一贯的 `to-x` 风格（与 `maf to-paf`、`fas to-vcf` 等
> 一致）。详见 [[pairwise-selection.md]] §5 变更日志。

### 3.2 为何 V1 不做 fasta/maf

上一轮把 impg 的**可选项**当成了 pgr 的**下一步**，是误判。重新审视：

1. **impg 的 MAF 是可选项**，不是核心——核心是坐标投影 + 传递闭包
2. **pgr 已有 `pgr fas consensus`**——用户要 MSA 时，`pgr paf to-bed` → `pgr fa range`
   → `pgr fas consensus` 的 pipe 路径已通，不需要在 query 层重复 POA
3. **POA MSA 是图构建层的产物**——按 [[paf-route.md]] §2.4，图遍历和 MSA 是正交步骤，
   不应耦合进 query
4. **V1 的真正用户场景**——"给我 chr1:1000-5000 在 cohort 里的所有同源区段"，输出 BED 即可

`pairwise-selection.md` 变更日志也印证：06-27 曾"BED 成为默认输出"，06-28 又"BED/TSV 删除，
只输出 PAF"。本次（06-28 二次修订）的结论是——**`pgr paf query` 保持 PAF 默认，BED 通过独立
子命令 `pgr paf to-bed` 提供**，两者都保留，由用户按场景选择。BED 三列（`name start end`）是
坐标投影的轻量产物，PAF 含完整 CIGAR 适合需要比对记录的场景。

### 3.3 实际实现（V1 ✅ 已完成）

| # | 任务 | 文件 |
|---|------|------|
| 1 | `pgr paf query` 输出 PAF（默认，含 CIGAR/gi/bi） | `cmd_pgr/paf/query.rs` |
| 2 | `pgr paf to-bed` 输出 BED3（3 列，复用 run_query） | `cmd_pgr/paf/to_bed.rs` |
| 3 | `--bed-regions`/`-b` 参数 + BED 文件解析（多 region 批查） | `cmd_pgr/paf/query.rs` |
| 4 | 共享查询逻辑 `add_query_args` + `run_query`（供 to-bed/to-maf 复用） | `cmd_pgr/paf/query.rs` |
| 5 | 集成测试（PAF + BED + 批查 + 过滤） | `tests/cli_paf_query.rs` |

### 3.4 不做的（V1 时定的边界，部分后续阶段已实现）

V1 立项时定的边界：

- ❌ `to-fasta`（未比对序列提取——`pgr fa range` 已提供独立提取能力）
- ❌ `to-maf` POA MSA（推迟到 V3，pairwise MAF 在 V2 已覆盖大部分需求）
- ❌ GFA/VCF 输出（推迟到 V4/V5，需完整 graph engine）
- ❌ gfasort/gfaffix（pgr 自带 unchop/crush 替代，见 §4.3.2）
- ❌ minigraph 的 `gfa_t` 在 Rust 中的对应实现

V2-V5 已分别实现，下列边界维持：

- ❌ `to-fasta`（用户场景不需要，`pgr fa range` 已提供）
- ❌ indel VCF（V5 仅 substitution，indel VCF 见 §4.3.6 后续）
- ❌ 全局精细 GFA 合并（V4b 产出的局部 GFA 独立，不合并回全局，见 §4.3.4）

---
## 4. V1/V2/V3/V4/V5 路线

按 impg 各输出格式的依赖链与核心性递进：

| 阶段 | 内容 | 对应 | 代码量 |
|------|------|------|:---:|
| **V1** ✅ | `pgr paf query`（默认 PAF）+ `pgr paf to-bed`（轻量坐标）+ `-b regions.bed` 批查 | impg 默认 `-o bed` + `-b`（pgr 选 PAF 默认，见 §3.1） | ~480 |
| **V2** ✅ | `pgr paf to-maf`（pairwise MAF，按 CIGAR 直接还原，需 `-f TSV`）| impg `-o maf` 的 pairwise 子集 | ~560 |
| **V3** ✅ | `pgr paf to-maf --msa`（POA MSA，多序列合并，需 `--transitive` + POA）| impg `-o maf` 的 multi-way | （并入 V2 文件）|
| **V4a** ✅ | 粗全局 GFA（`pgr paf graph -f refs.fa --min-var-len 100`，seqwish DSU 风格）| seqwish `sds`+`links` | ~720（CLI 112 + 引擎 608）|
| **V4b** ✅ | 区域精细 GFA（`pgr paf to-gfa`，impg 风格，含 unchop + LN tag + `--crush`）| impg `query -o gfa` + `crush` | ~500 |
| **V5** ✅ | `pgr paf to-vcf`（POA MSA→VCF，substitution-only）| impg `-o vcf` 的子集 | ~290 |

### 4.1 为何坐标类输出是 V1 核心

impg 的 11 种输出格式按"是否需要序列文件"分两类：

| 类别 | 格式 | 需 `-f` | 用途 |
|------|------|:---:|------|
| **坐标类**（核心） | `bed`/`bedpe`/`paf` | 否 | "哪些序列的哪些区段同源"——喂给下游工具 |
| **序列类**（可选） | `fasta` | 是 | 提取未比对序列 |
| **MSA 类**（可选） | `maf`/`fasta-aln` | 是 | POA 多序列比对 |
| **图类**（可选） | `gfa`/`vcf`/`gbwt` | 是 | 物化图，需完整 graph engine |

pgr V1 同时提供 PAF（`pgr paf query`，含 CIGAR/gi/bi 完整比对记录）与 BED（`pgr paf to-bed`，3 列轻量坐标）。
PAF 适合需要完整比对记录的场景，BED 三列（`name start end`）是坐标投影的轻量产物，最 pipe
友好——喂给 `pgr fa range` 提取序列。两者分别由独立子命令提供，符合 pgr CLI 一贯的 `to-x` 风格
（与 `maf to-paf`、`fas to-vcf` 等一致），详见 §3.1 默认输出决策。

### 4.2 为何 fasta/maf 后移

- **fasta 不做**：用户场景是直接看比对结果，不需要裸序列。`pgr fa range` 已提供独立提取能力。
- **pairwise maf 在 V2**：按 CIGAR 直接还原，不需 POA，是 impg `-o maf` 的子集
- **POA MSA 推到 V3**：多序列合并是图构建层产物，按 [[paf-route.md]] §2.4 不应耦合进 query 的核心路径。
  `pgr fas consensus` 已提供成熟 POA 后端，`-o maf` pairwise → `fas consensus` 的 pipe 路径已通

### 4.2.1 V2 `pgr paf to-maf` 设计（pairwise MAF by CIGAR）

**核心思路**：比对已被 chain/net 等上游程序优化过，**不需要再次 refine**，直接按 CIGAR 还原
pairwise MAF block。每条 query result（PAF record 的坐标投影）输出一个 2 序列 MAF block。

**`-f` 参数语义**（TSV，非单个 FASTA）：

```
# genome_name <tab> bgzf_fasta_path
sample1    /data/cohort/sample1.fa.gz
sample2    /data/cohort/sample2.fa.gz
ref        /data/cohort/ref.fa.gz
```

- 两列 TSV：`genome_name` 与 PAF 中的 query/target name 一一对应
- 每个 FASTA 必须 BGZF 压缩（`.gz`），支持 `pgr fa range` 的随机访问基础设施
- 启动时全量校验：PAF index 中所有 `names` 必须都在 TSV 中，**任一缺失即硬错误退出**
  （用户决策：严格模式，避免静默漏数据）

**MAF block 构建**（按 CIGAR 还原）：

1. 对每条 query result `(query_id, q_iv, t_iv, cigar)`：
   - 从 TSV 查 `query_name` 和 `target_name` 的 fasta 路径
   - 用 [libs/loc.rs](file:///Volumes/ExtHome/Scripts/pgr/src/libs/loc.rs) 的 `fetch_record` 取两条
     record（LRU 缓存最近访问的 genome，复用 `pgr fa range` 的模式）
   - 切片到 `q_iv.first..q_iv.last` 和 `t_iv.first..t_iv.last`
   - 遍历 CIGAR ops 同步推进两条序列，构建带 gap 的对齐字符串：
     - `=`: 两条都取一个碱基
     - `X`: 两条都取一个碱基（不同）
     - `M`: 两条都取一个碱基，逐位比较，相同写 `=` 不同写 `X`（用户决策：需查询两条序列）
     - `I`（insertion in query, gap in target）: query 取碱基，target 写 `-`
     - `D`（deletion in query, gap in query）: target 取碱基，query 写 `-`
   - `-` 链 record：query 序列先反向互补，再应用 CIGAR；MAF block 的 strand 字段标 `-`

2. 输出 MAF block（每条 record 一块，target 在前 query 在后）：
   ```
   a
   s target_name t_start  t_size  +       t_total  aligned_seq_t
   s query_name  q_start  q_size  strand  q_total  aligned_seq_q
   ```
   - `t_start`/`q_start`: 0-based，对齐区段在原序列的起点
   - `t_size`/`q_size`: 不含 gap 的对齐长度（即原序列区段长度）
   - `strand`: `+` 或 `-`（target 恒为 `+`，query 按 PAF record）
   - `t_total`/`q_total`: 原序列总长（从 `loc_of` 或 FASTA record 取）

**与 `pgr fa range` 的关系**：复用 [libs/loc.rs](file:///Volumes/ExtHome/Scripts/pgr/src/libs/loc.rs) 的
`create_loc` + `load_loc` + `fetch_record` 三件套。每个 genome 首次访问时按需创建 `.loc`（若不存在），
后续访问直接走 LRU cache + `fetch_record`。不直接调 `pgr fa range` 子命令（避免 subprocess 开销），
但共享同一套随机访问基础设施。

**不做的事**：
- ❌ 不做 POA 多序列合并（V3 才做）
- ❌ 不做 refine（比对已由上游 chain/net 优化）
- ❌ 不做 `to-fasta`（裸序列提取，用户场景不需要；`pgr fa range` 已提供）

**`-` 链 MAF 处理（V2 已实现）**：
- `PafMetadata` 增加 `strand: char` 字段，`insert_record` 从 `PafRecord` 填充，mirror entry 恒为 `+`
- `QueryResult` 元组第 7 元素为 `strand`，沿 `query` / `query_transitive_bfs` 传出
- `project()` 对 `-` 链把 CIGAR query offset 当作 RC offset，通过 `rc_to_forward()`
  转换回 forward 坐标：RC offset `[rc_lo, rc_hi)` → forward `[query_end - rc_hi, query_end - rc_lo)`。
  全比对时 RC offset = `[0, aligned_q_len)` → forward = `[query_start, query_end)`（与 `+` 链一致）；
  sub-interval 时正确返回 forward 子段。
- `output_maf` 对 `-` 链 record：
  1. 取 forward query[qs:qe]
  2. `reverse_complement` 得到对齐方向序列（CIGAR 列从左到右匹配 RC(query)）
  3. CIGAR 从 offset 0 走（`rec_qs_eff = 0`），`qs_eff = rec_qe - qe`（sub-interval 在 RC
     offset 中的起点；`rec_qe = rec_qs + aligned_q_len`）。这保证 `build_maf_block` 索引
     `q_seq[(cq + skip_t) - qs_eff]` 落在 `[0, qe - qs)` 内，sub-interval 不会越界。
  4. `s` 行 strand 标 `-`，`q_start = q_src_size - qe`（MAF 规范：负链 start 为正向坐标 srcSize - qe）
- 索引版本 bump v3→v4，旧索引需 `pgr paf index` 重建

**已知限制（V2 当前实现）**：
- `M` op 按原样输出两条碱基（MAF 格式本就不区分 `=`/`X`，靠下游逐位比较即可）。

### 4.2.2 V3 `pgr paf to-maf --msa` 设计（multi-way MSA via POA）

**核心思路**：`--transitive` BFS 收集一个 region 的所有同源片段后，用 POA 把它们合并成
**单个多序列 MAF block**（而非 V2 的每条 result 一个 pairwise block）。CIGAR 在此模式下被
忽略——MSA 列由 POA 重新决定。

**入口**：`pgr paf to-maf ... -t --msa -f genomes.tsv`（`--msa` 隐含建议配合 `--transitive`）

**序列准备**（每个 region 一个 block）：
1. **target 序列**：取该 region 第一条 result 的 `t_iv`，作为 block 的第一条 `s` 行（strand `+`）。
2. **query 序列**：遍历所有 result，按 `q_iv` 提取；`-` 链先 `reverse_complement`，使其与 target
   同向进入 POA。MAF `start` 字段：`+` 链为 `qs`，`-` 链为 `srcSize - qe`（与 V2 一致）。
3. **去重**：BFS 经 mirror index 可能返回 target 自身作为 query（自环），用
   `(name, start, strand, src_size)` 元组比对 target key，相同则跳过，避免重复序列进入 POA。

**POA 引擎**：直接复用 [libs/poa/](file:///Volumes/ExtHome/Scripts/pgr/src/libs/poa) 的 `Poa::new` →
`add_sequence` × N → `msa()`。打分参数通过 `--match`/`--mismatch`/`--gap-open`/`--gap-extend` 暴露
（默认 5/-4/-8/-6，与 `fas consensus` 一致）。AlignmentType 固定 Global。

**输出**：一个 `a` block，含 N 条 `s` 行（target 第一，其余按 result 顺序）。每行 `size` =
该行非 gap 字符数。

**与 impg 的差异**：impg 的 multi-way MSA 走 per-bubble POA（先物化 GFA 再平滑），pgr 直接在
BFS 等价类上做 POA，省去图物化步骤。代价是失去 bubble 结构信息——适合"看一眼同源区段的多序列
比对"，不适合全图 genotyping（后者走 V4a→V5 路径）。

**已知限制（V3 当前实现）**：
- POA 不处理循环移位（circular shift）——若两条序列是旋转关系，POA 可能产生末端 gap。
  真实场景中 `-` 链已 RC，循环移位较少见；若确需处理，参考 `fas multiz` 的 banded DP。
- target 序列取自第一条 result 的 `t_iv`——所有 result 共享同一 target_id（run_query 按
  target_id 过滤），但 `t_iv` 可能因 result 不同而区间略异（sub-interval projection）。
  当前取第一条的区间作为 target 序列，与 V2 pairwise 模式下每条 result 各自取 `t_iv` 的行为
  略有不同；多序列场景下统一 target 是必要简化。

### 4.2.3 V5 `pgr paf to-vcf` 设计（substitution-only VCF via POA MSA）

**核心思路**：复用 V3 的 POA MSA 产物，遍历每个 MSA column 调用 substitution 输出 VCF。
**跳过 GFA 物化**——V3 已证明 BFS 等价类 + POA 直接产出 MSA，无需先物化 GFA 再转 VCF。
这是 pgr 与 impg 的关键差异：impg `-o vcf` 走 query→GFA→VCF，pgr 走 query→MSA→VCF，省去图物化。

**入口**：`pgr paf to-vcf ... -t -f genomes.tsv`（推荐配合 `--transitive`）

**复用 V3 基础设施**：
- `build_msa_entries`（target + queries 序列收集，RC，去重）—— V3 提取的 `pub(crate)` 函数
- `FastaStore` / `load_fasta_tsv` / `reverse_complement` —— V3 提取的 `pub(crate)` 项
- POA 引擎与打分参数（`--match`/`--mismatch`/`--gap-open`/`--gap-extend`）

**VCF 生成**（每个 region）：
1. `build_msa_entries` → POA MSA（`Vec<String>`，含 gap，第一条是 target）
2. 遍历 MSA column，**跳过 target 是 gap 的列**（indel column，out of scope）
3. 对 target 非 gap 列：REF = target base，ALT = 其余序列中 distinct 的非 REF ACGT base
4. POS = `target.start + t_aln_pos + 1`（1-based VCF 坐标；`t_aln_pos` 跳过 target gap 列）
5. GT：REF→`0`，ALT[i]→`i+1`，gap/非 ACGT→`.`
6. 仅当存在至少 1 个 ALT 时输出行

**Header**：`##fileformat=VCFv4.2` + `##FORMAT=<ID=GT,...>` + `#CHROM...`（sample 列为 entries 的
name 顺序，target 第一）

**与 `pgr fas to-vcf` 的关系**：逻辑同构（都遍历 MSA column 调 substitution），但输入源不同——
`fas to-vcf` 读 block FA（已 MSA），`paf to-vcf` 自己跑 POA 生成 MSA。未来可考虑提取共享的
"MSA→substitution VCF" 工具函数，但当前两处实现都很短（~30 行核心循环），重复可接受。

**已知限制（V5 当前实现）**：
- **Substitution-only**：indel（target gap 或 query gap 列）被跳过。indel VCF 需把连续 gap 合并成
  indel 事件（左对齐、规范化），工作量显著高于 substitution。后续若需要，可参考 impg `vcf.rs`
  的 indel 处理，或 pipe 到 `fas to-vcf`（若 `fas` 支持 indel）。
- **无 phasing**：GT 是单倍体（0/1/.），不做 phased genotype。多倍体场景需扩展 FORMAT。
- **REF = target**：target 序列作为 REF。若 target 是某个 query 的子集（sub-interval），REF 可能
  不是最合适的参考等位基因，但与 V3 MSA 的 target 选择保持一致。

### 4.3 V4 的能力跃迁：两段式 GFA

V4 采用**两段式 GFA**——粗全局 + 区域精细，混合 minigraph 和 impg 各自所长：

| 工具 | 全局粗 GFA | 区域精细 GFA | 区域 → MSA/VCF |
|------|:---:|:---:|:---:|
| **minigraph** | ✅（≥100bp SV，rGFA）| ❌（不做）| ❌（小变体用标准工具）|
| **impg** | ❌（不物化全局图）| ✅（`query -o gfa`）| ✅（`query -o maf/vcf`）|
| **pgr V4** | ✅（V4a）| ✅（V4b）| ✅（V5）|

#### 4.3.1 V4a：粗全局 GFA（seqwish DSU 风格，✅ 已实现）

- **输入**：PAF 文件 + FASTA 文件（`-f`）
- **过滤**：`--min-var-len 100`（对齐 minigraph，只保留 ≥100bp SV）
- **输出**：GFA v1.0（S/L/P 三类行，节点 1-based）
- **用途**：可视化（Bandage/odgi）、SV 概览、作为后续 query 的坐标锚
- **数据源**：PAF 索引的显式投影，不引入新的真实源
- **实现**：`src/libs/paf/graph.rs`（470 行）+ `src/cmd_pgr/paf/graph.rs`（CLI 包装），
  5 单元测试 + 7 集成测试，覆盖正向/反向/大 indel 切分/小 indel 不切/min_var_len 阈值过滤
- **算法骨架**：seqwish 风格的段级 DSU（CIGAR 切分 → 段对链接 → DSU 传递闭包 → 节点序列
  → 路径构建 + novel 段补全 → 边派生 → GFA 输出）。详见 [[seqwish.md]] §6.2。
- **简化项**（相对 seqwish）：
  - 无 disk-backed interval tree（pgr 内存型，足够 PAF 规模）
  - 无 SparseBitVec（直接用 `Vec<u8>`）
  - 无 lock-free DSU（单线程足够）
  - 路径方向恒为 `+`（反向链段已翻转坐标到正链，rGFA orientation 待 V4b+）
- **rGFA tag 暂缺**：当前输出 GFA v1.0 的 S/L/P，未含 SN/SO/SR（稳定坐标系 tag），
  作为后续可选增强（与 minigraph 兼容性需要时再加）。

#### 4.3.2 V4b：区域精细 GFA（impg 风格，✅ 已实现）

- **输入**：PAF 索引 + 用户指定 region（或粗 GFA 上定位的区段）
- **流程**：BFS 传递闭包 → 提取序列 → POA → 局部 GFA
- **输出**：局部 GFA（含 base-level 变异）
- **用途**：特定基因座的精细分析
- **对应**：impg `query -o gfa`
- **实现**（`src/cmd_pgr/paf/to_gfa.rs`）：
  - 复用 `to-maf` 的 `build_msa_entries` 收集 target/query 序列
  - 调用 `Poa` 生成对齐图，通过 `graph()`/`paths()`/`sequences()` 访问器取得内部结构
  - **Unchop（默认开启）**：合并无分支的连续单碱基节点为多碱基段（节点 v 合并到前驱 p，
    当 p.out_degree=1 且 v.in_degree=1）。节点数减少约一个数量级，bubble 结构保留。
    拓扑序遍历计算段头（segment head），段序列按拓扑序拼接碱基，边/路径映射到段 id 去重。
  - **GFA header + LN tag**：输出 `H\tVN:Z:1.0`；S 行带 `LN:i:<len>` tag
  - **`--crush` 可选**：impg `crush` 风格 bubble 压缩。把 in/out 邻居集合相同的节点组
    合并成一个，保留 weight 最大的 allele。**丢失 base-level ALT 信息**——
    经过 ALT 的路径被改写到保留的 allele。默认关闭，仅用于 SV 概览图
  - SNP/indel 天然形成 bubble（POA 图结构直接表达），无需 MSA→GFA 转换
- **多 region 处理**：每个 region 独立产出 GFA 块（节点 ID 从 1 重启），
  多 region 用 `# region: <name>` 注释行分隔
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

---
## 5. 后续计划（V6-V8 与增量增强）

V1-V5 已完成索引→查询→图构建三层的最小可用闭环。后续按"规模扩展 → 图质量 → 应用层"递进，
对照 impg 与 seqwish 仍可借鉴的能力，结合 pgr 独有优势（UCSC Chain/Net 体系、纯 Rust POA）。

### 5.1 V6：规模扩展（4 万大肠杆菌级）

当前 V4a 在小 cohort（数十基因组）上验证通过。扩到 [[ecoli-cohort.md]] 的 4 万大肠杆菌时，
需引入 seqwish 的工程优化（按收益/复杂度排序）：

| 优化 | 来源 | 何时引入 | 估计收益 |
|------|------|----------|----------|
| **spanning tree 剪枝** | [[seqwish.md]] §3.1 | BFS 边遍历数成为瓶颈时 | N²→N-1 边遍历，BFS 提速 10×+ |
| **lock-free DSU**（`portable_atomic::AtomicU128`） | [[seqwish.md]] §2.4 | V4a rayon 化时 | 并行 union-find，无锁竞争 |
| **`--repeat-max` / `--min-repeat-dist`** | [[seqwish.md]] §3.5 | 高拷贝重复区把图吹爆时 | 限制同序列在同节点的拷贝数 |
| **disk-backed interval tree** | [[seqwish.md]] §2.3 | 全图超 RAM 时 | 兜底可跑性，牺牲速度 |
| **SparseBitVec 序列边界** | [[seqwish.md]] §2.2 | HPRC 规模（数百单倍型、Gb 级） | 边界索引内存 O(m) 而非 O(N) |
| **PosT 单 u64 编码** | [[seqwish.md]] §2.1 | 反链投影成为热点时 | 单棵树同时表达正反链 |

**判断标准**：先跑 4 万大肠杆菌基准测试（无优化版本），按实际瓶颈选 1-2 项引入。不做"预防性优化"。

### 5.2 V7：图质量与归一化（impg graph_report + smoothxg 风格）

V4b 的 `--crush` 是 impg crush 8 阶段流程的最小子集（仅"邻居集合相同的节点合并"）。
完整 crush 有 15 种 ResolutionMethod（[[impg.md]] §6.5），覆盖 median 路由、indel 平滑、
超长 bubble 拆分等。pgr 不打算全抄——按需挑选：

| 能力 | 来源 | pgr 价值 | 优先级 |
|------|------|----------|--------|
| **`pgr paf graph-report`** | impg `graph report`（15+ 维度） | 图质量评估，为参数调优提供依据 | 高（不依赖新算法） |
| **`--smooth`**（smoothxg 风格块归一化） | impg `smooth.rs` + smoothxg | 合并相邻 bubble、填充 gap 列 | 中（需 POA 块分解） |
| **`--gfaffix`**（图归一化） | impg `run_gfaffix` | 重复 bubble 标准化 | 低（外部 binary，pgr 倾向纯 Rust） |
| **完整 crush 8 阶段** | impg `resolution.rs` | 复杂 bubble 处理 | 低（V4b `--crush` 已覆盖主要场景） |

**判断标准**：`graph-report` 优先做（无新算法，纯统计），`--smooth` 视用户需求评估。完整 crush 不做——
V4b 的简化版已覆盖 SNP/简单 indel bubble，复杂 bubble 让用户走 impg/smoothxg 管道。

### 5.3 V8：应用层 — genotyping（impg Genotype/Infer）

impg 能力栈顶端是 `Genotype` / `Infer`（[[impg.md]] §7）。pgr 当前无应用层。基因分型的核心抽象
是 impg 的 **graph feature evidence**：候选 haplotype 编码为特征向量（节点/边/路径 presence），
既可来自物化 GFA，也可来自隐式图查询。

**pgr 的路径**：复用 V4a 粗全局 GFA + V5 VCF，加 `pgr paf genotype`：

```
pgr paf graph cohort.paf -f refs.fa --min-var-len 100 -o graph.gfa  # V4a 物化
pgr paf genotype graph.gfa sample.fa  # 在图上基因分型
```

**为何推迟**：genotyping 需要稳定的图构建（V4-V7）作为前置，且需要真实数据验证。
当前 pgr 的核心价值在"图构建 + 区域查询"，应用层让位于 minigraph-cactus / vg 等成熟工具。

### 5.4 增量增强（不单独成版本）

| 增强 | 当前状态 | 目标 | 优先级 |
|------|----------|------|--------|
| **indel VCF** | V5 仅 substitution | 合并连续 gap 为 indel 事件，左对齐规范化 | 高（用户场景常见） |
| **`--syntenic-filter`** | 无 | 用 UCSC Chain/Net 验证 PAF 边，标低置信度 | 中（pgr 独有优势） |
| **`-m/--max-depth`** | BFS 默认无界 | 控制传递闭包深度（impg 默认 2） | 中（大 cohort 防爆炸） |
| **`--min-transitive-len`** | 无 | 过滤短传递片段（impg 默认 101） | 中（同上） |
| **rGFA tag**（SN/SO/SR） | V4a 缺失 | 给参考路径加 rGFA anchor tag | 低（与 minigraph 工具链互操作时再加） |
| **stage 化管道** | 单命令 | impg `gfa:crush:sort` 风格字符串 | 低（CLAUDE.md "三次相似再抽象"原则，当前单命令够用） |

### 5.5 不做清单（明确边界）

按 CLAUDE.md "不写半成品 / 三次相似再抽象"原则，下列能力**明确不做**：

- ❌ **完整 crush 8 阶段 + 15 ResolutionMethod** — V4b `--crush` 已覆盖主要场景，复杂 bubble 走外部工具
- ❌ **`--sparsify auto`（Mash KNN）** — pgr 走"已有 PAF 覆盖度先验"路线（[[paf-route.md]] §3.2），不引入 mash 依赖
- ❌ **syng 后端（syncmer 免比对）** — pgr 已有 `pgr lav lastz` 全套，不需要免比对后端
- ❌ **1ALN/TPA 格式支持** — pgr 用 PAF/MAF，不需要第三种比对格式
- ❌ **全局精细 GFA 合并** — 见 §4.3.4，保持"粗全局 GFA 是不可变投影"语义
- ❌ **minigraph 的 `gfa_t` 在 Rust 中的对应实现** — pgr 走 seqwish DSU 路线（V4a），不重建 minigraph 数据结构

### 5.6 优先级路线图

```
当前 (V1-V5) → 增量增强 (indel VCF / --syntenic-filter / -m)
            → V6 规模扩展 (4 万大肠杆菌基准 → 选 1-2 项 seqwish 优化)
            → V7 图质量 (graph-report → --smooth 评估)
            → V8 应用层 (genotyping，远期)
```

**近期焦点**：indel VCF（用户场景最直接）+ 4 万大肠杆菌基准测试（验证 V4a/V4b 可扩展性）。
V7/V8 视用户反馈与真实数据验证后再启动，不预先实现。

