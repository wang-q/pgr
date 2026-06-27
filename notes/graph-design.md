# 方向 D：图构建层 — 设计文档

本文档定义 pgr 从"找到同源"（查询层）到"构建图"（图构建层）的设计方案。
参考 impg-0.4.1 的 `graph.rs`（1207 行）、`graph_pipeline.rs`（312 行）和
pgr 已有资产（`libs/poa/`、`libs/fas_multiz.rs`、`docs/gfa.md`）。

---

## 1. impg 怎么做

### 1.1 整体流程

```
BFS 传递闭包 → Vec<Interval<u32>>（同源区段列表）
    ↓
prepare_sequences()          → rayon 并行提取序列（含 reverse-complement）
    ↓
prepare_poa_graph_and_sequences() → SPOA 图（longest-first 喂入）
    ↓
spoa_graph_to_unchoped_gfa() → 原始 GFA（含 strand post-processing + unchop）
    ↓
normalize_and_sort()         → gfaffix（标准化）+ gfasort（拓扑排序）
    ↓
输出 GFA / MAF / FASTA_ALN
```

关键点：
- **per-query POA**：每个 BFS 查询结果跑一次独立的 SPOA，不跨查询优化
- **外部依赖**：`spoa_rs`（SIMD C++ SPOA via FFI）+ `gfaffix` + `gfasort`
- **node 粒度由 SPOA 决定**：SPOA 在变异边界处自然产生分叉节点，
  impg 不做额外的节点拆分

### 1.2 关键数据结构

`SequenceMetadata`（`graph.rs:14-22`）：
```rust
pub struct SequenceMetadata {
    pub name: String,       // 序列名
    pub start: i32,         // MAF-style start（反向链已翻转）
    pub size: i32,          // 片段长度
    pub strand: char,       // '+'/'-'
    pub total_length: usize,// 序列全长
    pub path_name_override: Option<String>,
}
```

这个结构体桥接了"坐标投影结果"和"GFA 路径名"——`path_name()` 方法将
coordinate + strand 转换为 GFA P-line 格式（`name:fwd_start-fwd_end`）。

### 1.3 外部工具链

| 工具 | 作用 | pgr 对应 |
|------|------|---------|
| `spoa_rs` | SIMD C++ POA（FFI） | `libs/poa/`（纯 Rust，无 FFI） |
| `gfaffix` | GFA 标准化（去冗余节点、合并短分支） | **需要**（或跳过，输出原始 GFA） |
| `gfasort` | 拓扑排序 + unchop + SGD 排序 | **需要**（或跳过） |

---

## 2. pgr 有什么

### 2.1 POA 引擎

`libs/poa/poa.rs` — 纯 Rust，无外部依赖。接口：

```rust
let mut poa = Poa::new(params, AlignmentType::Global);
poa.add_sequence(b"ACGT");
poa.add_sequence(b"ACAGT");
let consensus = poa.consensus();  // → b"ACAGT"
let msa = poa.msa();             // → ["AC-GT", "ACAGT"]
let num_nodes = poa.num_nodes(); // 图的节点数
```

**与 impg SPOA 的对比**：

| 维度 | impg `spoa_rs` | pgr `Poa` |
|------|:---:|:---:|
| 实现 | C++ FFI（SIMD） | 纯 Rust（scalar） |
| 外部依赖 | 需要 `spoa_rs` crate | 无 |
| MSA 输出 | ✅ | ✅ |
| GFA 输出 | ✅（`generate_gfa`） | ❌（需新增） |
| 性能 | SIMD 加速 | 单线程 scalar |

**差异**：pgr 的 POA 目前**不能输出 GFA**——只有 `consensus()` 和 `msa()`。
需要新增 `poa.to_gfa()` 方法（约 50 行），将 POA 图转换为 GFA S-line + P-line 格式。

### 2.2 Banded DP 多序列合并

`libs/fas_multiz.rs`（1158 行）— 完整的 banded DP 实现。关键优势：
- 参考序列坐标约束（±radius 对角线）
- 三种 gap model（constant/medium/loose）
- substitution matrix 支持（LASTZ 格式或 preset）
- core/union 两种窗口模式

**对图构建的价值**：如果多个 genome 共享一个 reference，banded DP 可以
利用 reference 坐标约束产生更精确的 MSA，比纯 POA（完全无参考）更可靠。
但代价是——需要指定 reference，而 POA 是 reference-free 的。

### 2.3 MSA 裁剪

`libs/alignment.rs` 的四个裁剪函数（`trim_pure_dash`、`trim_outgroup`、
`trim_complex_indel`、`trim_head_tail`）可以直接在图构建的 POA 阶段之后
清理 MSA 边界。

### 2.4 数据源优势

pgr 已有 UCSC Chain/Net 体系（`pgr chain net`、`pgr chain sort` 等），
可以在图构建前做 syntenic 过滤——只保留 Chain/Net 验证过的同源区域，
减少 POA 阶段的假阳性。

---

## 3. 三个核心设计问题

### 3.1 你要什么样的图？

impg 走 PGGB 路线：SPOA → GFA → gfaffix → gfasort → seqwish/crush。
但 pgr 不需要照搬。

**选项 A：GFA 显式图（跟 impg 一样）**

```
BFS → 提取序列 → POA → GFA → gfaffix → gfasort → seqwish/crush
```

优点：业界标准，下游工具链成熟（vg/odgi/seqwish）。
缺点：引入外部依赖（gfaffix、gfasort），管道长。

**选项 B：隐式图 + 按需 MSA（当前路线延伸）**

```
BFS → 提取序列 → POA → MSA → MAF/FASTA_ALN 输出
```

不物化 GFA，只输出 MSA。查询层继续用 PAF 区间树。
优点：零外部依赖，管道短。
缺点：不能做图统计、图可视化、图比对。

**选项 C：两阶段（先隐式后显式）**

第一步：输出 MSA（选项 B），积累数据和经验。
第二步：评估后决定是否引入 GFA 管道（选项 A）。

**建议**：选项 C。第一步先做 POA → MSA 输出（码量约 100 行），
验证 POA 在真实数据上的质量。若质量过关，第二步引入 GFA。

### 3.2 图的节点粒度怎么定？

**POA 决定节点边界**（impg 方式）：POA 算法自然在变异位点处产生分叉，
节点边界由 POA 的图结构决定。优点是无需额外参数，缺点是节点可能过碎
（POA 对每个 SNP 产生一个分叉）。

**Banded DP 决定节点边界**：先用 banded DP 确定 core 区段（所有 genome 一致），
在 core 边界处切节点。这避免了对高一致性区域的过度分叉。

**混合方式**：先用 banded DP 确定 core 边界 → 在 core 内部不切节点 →
在变异密集区切细节点（但不用 POA 的 per-SNP 粒度）。

**建议**：先用 impg 的 POA 方式（无参数，简单），在真实数据上评估碎片化程度。
如果 per-SNP 碎片化严重，再引入 banded DP 粗细结合。

### 3.3 输入数据够不够？

当前 pipeline 输入是两序列 MAF（ref↔query_i），star topology。
图构建通常需要 denser network（A↔B、B↔C 都有直接比对）。

**够用的场景**：单 locus 查询（HLA/C4），传递闭包已覆盖间接同源。
在同一 locus 上的所有同源片段可以直接跑 POA。

**不够用的场景**：全基因组图构建，需要更密集的 pairwise 网络。
此时需要生成新的 pairwise 比对——但这不属于图构建层的范畴，
属于 `paf-route.md` §3.2 的"第二层挑选问题"，需要单独设计。

**建议**：第一步在 multiz MAF 数据上验证 POA → MSA 输出。
不需要引入新比对，复用已有 PAF 网络上的传递闭包结果。

---

## 4. 最小可行实现（V1）

### 4.1 新增：`PoaGraph::to_gfa()` + `pgr paf msa` 命令

**库层**：给 POA 图新增 GFA 输出能力。不需要 `SequenceMetadata`——
直接用 sequence name 作为 P-line 名，简化 impg 的坐标标注。

```rust
// 新增方法
impl Poa {
    pub fn to_gfa(&self, names: &[String]) -> String;
}
// 用法
let mut poa = Poa::new(params, AlignmentType::Global);
for (name, seq) in &bfs_results {
    poa.add_sequence(seq);
}
println!("{}", poa.to_gfa(&names));
```

**命令层**：`pgr paf msa <region> -i index.paf.idx -f seqs.fa`

```bash
# 单命令：BFS → 提取序列 → POA → MSA
pgr paf msa chr1:1000-5000 -i cohort.paf.idx -f genomes.fa

# 可选输出格式
pgr paf msa ... --gfa -o out.gfa
pgr paf msa ... --msa -o out.fas
pgr paf msa ... --consensus -o consensus.fa
```

这个命令封装了 BFS + 序列提取 + POA 三步，输出 MSA 或 GFA。
与 impg 的 `output_results_gfa` 等价，但更简洁。

### 4.2 不引入外部依赖

pgr 的 POA 是纯 Rust——不需要 `spoa_rs`。gfaffix/gfasort 推迟到第二步。
第一步只输出**原始 POA GFA**——未标准化的图，仅供 pgr 内部使用。

### 4.3 工作清单

| # | 任务 | 文件 | 预估行数 |
|---|------|------|:------:|
| 1 | `Poa::to_gfa()` — GFA 输出 | `libs/poa/poa.rs` | ~50 |
| 2 | `pgr paf msa` — CLI 命令 | `cmd_pgr/paf/msa.rs` | ~80 |
| 3 | 在多序列 BFS 结果上端到端验证 | — | — |
| 4 | 单元测试（POA → GFA 往返） | `libs/poa/poa.rs` | ~20 |
| 5 | 集成测试 | `tests/cli_paf.rs` | ~3 |
