# pgr 泛基因组第一步行动计划

本文档定义 pgr 走向泛基因组的**第一步要做什么**——三个子命令的完整规格与验证标准。
路线决策与理由见 [[paf-route.md]]，PAF 模块的代码级实现参考见 [[paf-implementation.md]]。

核心结论：

> pgr 已有两序列 MAF（可转 PAF），天然避开了 impg 的 all-vs-all 比对。
> 第一步只需实现 PAF 索引 + 区间投影 + 传递闭包 BFS，复用已有 pairwise 基础设施，
> 不做新的比对、不物化 GFA、不引入 wfmash。

---

## 1. 需要实现的三样东西

```
pgr maf to-paf          # MAF → PAF 转换（数据源）
pgr paf index           # PAF → 区间树索引（索引层，不过滤）
pgr paf query <region>  # 区间投影 + 传递闭包 BFS（查询层，在此过滤）
```

对应 impg 的能力栈（[[impg.md]] §1.1.3）：

```
索引层（Index）── pgr paf index    ← 全量装入，不过滤
    ↓
查询层（Query）── pgr paf query    ← 区间投影 / --transitive BFS
    │                                 过滤参数：--merge-distance, --min-identity,
    │                                 --min-output-length, --max-depth
    ↓
（第二步才到图构建层）
```

> **三个命令共用一套 CIGAR 基础设施**（见 §1.1），pgr 当前完全没有 CIGAR 解析能力，
> 需要从零实现。这是跨 `to-paf`/`query` 的共享依赖，独立于任何一个命令。

## 1.1 CIGAR 公共模块 — 跨命令共享基础设施 ✅ 已实现

### 定位

pgr 目前没有 CIGAR 解析能力。CIGAR 是 `maf to-paf`（对齐串 → CIGAR 展开）和
`paf query`（CIGAR 坐标投影）的共同依赖，应实现为 `src/libs/paf/cigar.rs`，
不绑定任何特定命令。

### 需要实现的组件

| 组件 | 类型 | 参考 impg | 说明 |
|------|------|-----------|------|
| `CigarOp` | struct | `impg-0.4.1/src/impg.rs:73-138` | bit-packed：3 位 op code + 29 位 length → `u32` |
| `CigarOp::new(len, op)` | 构造函数 | `impg.rs:79` | `=`/`X`/`I`/`D`/`M` 五种 op |
| `CigarOp::op()` / `len()` | 访问器 | `impg.rs:93-107` | 从 `u32` 解码 |
| `CigarOp::target_delta()` | 方法 | `impg.rs:113-118` | target 坐标推进量：`I`→0, `D`→len, 其他→len |
| `CigarOp::query_delta(strand)` | 方法 | `impg.rs:121-133` | query 坐标推进量，反向链取负 |
| `parse_cigar(s)` → `Vec<CigarOp>` | 函数 | impg 的 `parse_cigar_to_delta` | 字符串解析（如 `"10=5I5D"`） |
| `format_cigar(&[CigarOp])` → `String` | 函数 | 自实现 | 反向格式化 |
| `cigar_from_alignment(ref_aln, qry_aln)` → `Vec<CigarOp>` | **pgr 独有** | 无 | 从两条 MAF `s` 行的 alignment 逐 base 对比生成 CIGAR |
| `gap_compressed_identity(&[CigarOp])` → f64 | 函数 | `main.rs:12042-12061` | `gi = matches / (matches + mismatches + #indel_events)` |
| `block_identity(&[CigarOp])` → f64 | 函数 | 同上 | `bi = matches / (matches + mismatches + indel_bp_total)` |

### bit-packing 设计

直接复用 impg：`u32` 中 bits[31:29] 存 op code，bits[28:0] 存 length（最大 512Mbp）。
`Vec<CigarOp>` 存为 `Vec<u32>`，对齐友好、内存紧凑。
`target_delta`/`query_delta` 是纯位运算，无分支。

### impg CIGAR 源码关键位置

| 函数/类型 | 文件 | 行号 | 说明 |
|-----------|------|------|------|
| `CigarOp` struct | `src/impg.rs` | 73-138 | bit-packed：op()/len()/target_delta()/query_delta() |
| `invert_cigar_ops()` | `src/impg.rs` | 144-160 | I↔D 交换，用作双向索引的反向映射（pgr V1 跳过） |
| `parse_cigar_to_delta()` | `src/impg.rs` | 2884-2898 | 字符串 → `Vec<CigarOp>`，纯字符遍历，**可直接移植** |
| `calculate_gap_compressed_identity()` | `src/impg.rs` | 2901-2922 | `gi` 计算：indel 按事件计数 |
| `read_cigar_data()` | `src/paf.rs` | 68-114 | CIGAR 懒加载 seek+read，pgr V1 不采用 |
| `compact_cigar()` | `src/syng_graph.rs` | 445-464 | byte CIGAR → run-length 字符串，仅当需要时移植 |
| `cigar_stats()` | `src/syng_graph.rs` | 469-483 | 统计 (matches, mismatches, ins, del)，pgr 用 fold 替代 |
| identity 计算（含 bi） | `src/main.rs` | ~12042 | `output_results_paf` 内同时计算 gi 和 bi |

### 外部 CIGAR 解析参考

除 impg 外，以下项目也可作为 CIGAR 设计的参考：

| 项目 | 语言 | 说明 |
|------|------|------|
| **noodles-sam** (`cigar.rs`) | Rust | pgr 已依赖 noodles，其 SAM CIGAR 解析可参考——基于 `u32` 数组的 bit-packing 设计 |
| **rust-htslib** | Rust | htslib 的 Rust 绑定，CIGAR 操作接口成熟（DeepWiki 有专门文档） |
| **rlannescigarparser** | Rust | 专门为 CIGAR 字符串解析发布的轻量 Rust 库（2025 年 LinkedIn 发布） |
| **CIGARStrings.jl** | Julia | CIGAR 解析与操作，设计思路可参考（非 Rust） |

**建议**：pgr 的 CIGAR 以 impg 的 `CigarOp` bit-packing 为蓝本，
`parse_cigar` 参考 impg 的 `parse_cigar_to_delta`。
noodles-sam 作为备选参考（pgr 已有依赖，无需额外引入 crate）。

### pgr 独有：从 MAF 对齐串生成 CIGAR

从 ref 和 query 两条 MAF `s` 行的 `alignment` 向量逐 base 对比：

```
ref=ACG-, qry=ACG  → ref gap → CigarOp(3, '=') + CigarOp(1, 'I')
ref=ACG,  qry=ACG- → qry gap → CigarOp(3, '=') + CigarOp(1, 'D')
ref=ACG,  qry=ACG  → match   → CigarOp(3, '=')
```

**首版折中**：不区分 `=`/`X`，所有非 gap 位统一用 `M`。`M` 的坐标投影行为
与 `=`/`X` 完全一致（都是 target_delta=len, query_delta=len），
不影响查询正确性。

### 测试参考

| 测试 | 来源 | 行号 | 移植方式 |
|------|------|------|---------|
| CIGAR 字符串解析 | `impg.rs` | 3108 | 直接移植 `test_parse_cigar_to_delta_basic` |
| 对齐串 → CIGAR（pgr 独有） | — | — | 新增：ref=ACG-/qry=ACG → `[CigarOp(3,M), CigarOp(1,I)]` |
| identity 计算 | — | — | 新增：已知 CIGAR 的 `gi`/`bi` 与手算一致 |

### impg 测试缺口 & pgr 应补充的测试

impg 的 CIGAR 测试覆盖有限（见下表），pgr 应在此基础上补充边界情况。

**impg 缺失的 CIGAR 解析边界测试**：

| 场景 | 输入 | 预期 | 说明 |
|------|------|------|------|
| 空 CIGAR | `""` | `Ok(vec![])` | impg 未测 |
| 只有数字无 op | `"10"` | 数字被丢弃，`Ok(vec![])` | `parse_cigar_to_delta` 的隐式行为 |
| 0 长度 op | `"0=5I"` | `[CigarOp(0,=), CigarOp(5,I)]` | 0 长度 I 的 target_delta=0, query_delta=0 |
| 超大数字 | `"2147483647="` | i32::MAX 长度 CigarOp | 基因组比对不过此边界，但应验证不 panic |
| 非法 op | `"10Q"` | `panic!`（`CigarOp::new` 中 `panic!`） | impg 有被注释掉的测试，pgr 应改为 `Result` 或 `expect` |

**impg 缺失的 identity 计算测试**：

| 场景 | CIGAR | gi 预期 | bi 预期 |
|------|-------|---------|---------|
| 纯 match | `10=` | 1.0 | 1.0 |
| 含 1 个 I | `10=5I` | 10/(10+0+1) = 0.909 | 10/(10+0+5) = 0.667 |
| 含 1 个 D | `10=5D` | 10/(10+0+1) = 0.909 | 10/(10+0+5) = 0.667 |
| 混合 | `10=2X3I4D` | 10/(10+2+2) = 0.714 | 10/(10+2+7) = 0.526 |
| 空 CIGAR | `[]` | 0.0 | 0.0 |

**pgr 独有的对齐串→CIGAR 测试**：

| 场景 | ref alignment | qry alignment | 预期 CIGAR |
|------|--------------|---------------|-----------|
| 全 match | `ACGT` | `ACGT` | `[CigarOp(4,M)]` |
| ref gap（qry insertion） | `ACG-` | `ACGT` | `[CigarOp(3,M), CigarOp(1,I)]` |
| qry gap（qry deletion） | `ACGT` | `ACG-` | `[CigarOp(3,M), CigarOp(1,D)]` |
| 交错 gap | `AC-TG` | `ACGT-` | `[CigarOp(2,M), CigarOp(1,I), CigarOp(2,M), CigarOp(1,D)]` |
| 首尾 gap | `-ACGT-` | `TACGT` | `[CigarOp(1,I), CigarOp(4,M), CigarOp(1,I)]` |
| 全 gap（全 `-`） | `---` | `---` | `[]`（退化） |

**pgr 应补充的 `compact_cigar` / `format_cigar` 往返测试**：

| 场景 | CigarOp vec | 字符串 | 往返相等？ |
|------|-------------|--------|-----------|
| 单 op | `[CigarOp(10,M)]` | `"10M"` | parse(format(ops)) == ops |
| 多 op | `[CigarOp(3,M), CigarOp(1,I), CigarOp(2,D)]` | `"3M1I2D"` | ✅ |
| 空 | `[]` | `""` | ✅ |

**建议**：将这些测试作为 `src/libs/paf/cigar.rs` 的 `#[cfg(test)] mod tests`，
与 CigarOp 实现放在同一文件中，确保实现和测试紧密绑定。

---

## 2. `pgr maf to-paf` — 数据源 ✅ 已实现

### 现有基础

pgr 已有 MAF 解析器（`src/libs/fmt/fas.rs`），提供以下可直接复用的能力：

- `MafEntry` (line 242)：`s` 行的内存表示——`src`、`start`、`size`、`strand`、`src_size`、
  `alignment: Vec<u8>`（含 gap 的对齐串）
- `MafBlock` (line 292)：一个 MAF block 的所有 `s` 行（`entries: Vec<MafEntry>`）
- `next_maf_block()` (line 298)：从 `BufRead` 流读取下一个 MAF block

`to-paf` 在上述基础上需要补充：
1. 从 `a` 行解析 `score=`——当前 `parse_maf_block()` 忽略了 `a` 行内容
2. 从两条 `s` 行的 `alignment` 对比生成 `Vec<CigarOp>`（`=`/`X`/`I`/`D`）
3. 从 CIGAR 计算 `gi`（gap-compressed identity）和 `bi`（block identity）

> **实现参考**：除 impg 外，[wgatools](https://github.com/wjwei-handsome/wgatools)
> (v1.1.0, Bioinformatics 2025) 也有完整的 MAF→PAF 实现
> (`converter.rs:29-54` + `parser/maf.rs:484-519`)，其 `csv` crate flexible reader
> 和 `AlignRecord` trait 设计值得参考。详见 [[paf-implementation.md]] §13。

### 输入/输出

```bash
pgr maf to-paf ref_vs_query.maf -o ref_vs_query.paf
```

- 输入：两序列 MAF 文件（`a` 行含 ref + query 两条 `s` 记录）
- 输出：标准 PAF 文件（12 列 + 自定义标签）

### MAF block → PAF 行映射规则

```
MAF block:
a score=12345
s ref_name    100 50 + 1000 ACGT...acgt
s query_name   50 50 +  500 ACGT...acgt

→ PAF:
query_name  500  50  50  +  ref_name  1000  100  50  50  50  500  255
  gi:f:0.9500  bi:f:0.9200  cg:Z:10=5X2I3D  ms:i:12345
```

**映射逻辑**：
- PAF col 1-6：从 MAF `s` 行（query 侧）提取 name、length、start、end、strand
- PAF col 7-12：从 MAF `s` 行（target/ref 侧）同样提取
- PAF col 10（matches）：统计对齐串中 match 碱基数
- PAF col 11（block length）：`end - start` 之和（含 gap）
- PAF col 12（mapq）：固定 255（MAF 没有此信息）

**自定义标签**：

| Tag | 类型 | 来源 |
|-----|------|------|
| `gi:f:` | float | gap-compressed identity（从 CIGAR 计算） |
| `bi:f:` | float | block identity（从 CIGAR 计算） |
| `cg:Z:` | string | 从 MAF 对齐串展开的 CIGAR（`=`/`X`/`I`/`D`） |
| `ms:i:` | int | 从 MAF `a` 行 `score=` 提取——pgr 相对原生 PAF 的增值 |

### 需处理的情况

1. **正向 block**（`s` 行 strand = `+`）：坐标直接映射
2. **反向互补 block**（`s` 行 strand = `-`）：PAF 的 query strand 列填 `-`，坐标需翻转
3. **多序列 MAF block**（>2 条 `s` 行）：首版不支持，报错提示"仅支持两序列 MAF"。
   后续可分解成 pairwise 组合（ref↔query1、ref↔query2 各一条 PAF 行）
4. **MAF block 间 gap**：PAF 行之间无 gap 概念，相邻但非连续的 PAF 行由查询层的
   `--merge-distance` 合并（对应 impg 的 `SortedRanges` 合并，[[impg.md]] §3.3）

### 验证标准

- 对已知的 ref_vs_query MAF，生成的 PAF 行数与 MAF block 数一致
- PAF 第 5 列（strand）与 MAF `s` 行一致
- 反向互补 block 的坐标转换正确
- CIGAR 展开后，match base 数 = PAF 第 10 列
- Zero Panic：畸形 MAF 输入返回友好错误

---

## 3. `pgr paf index` — 索引层

### 输入/输出

```bash
# 单文件
pgr paf index sample1_vs_ref.paf -o sample1.paf.idx

# 多文件（cohort）
pgr paf index *.paf -o cohort.paf.idx
```

- 输入：一个或多个 PAF 文件（纯文本或 `.paf.gz`）
- 输出：`.paf.idx` 索引文件——序列化格式包含序列名映射 + per-target 区间树

### 索引结构

名字映射采用 pgr 现有模式——`IndexMap<String, u32>`（与 `libs/loc.rs` 和
`libs/phylo/tree.rs` 的整数索引风格一致），而非 impg 的独立 `SequenceIndex`。

```rust
pub struct PafIndex {
    /// Name → internal ID mapping (reuses pgr's IndexMap pattern).
    pub names: IndexMap<String, u32>,
    /// Per-target interval trees (coitrees).
    pub trees: FxHashMap<u32, Arc<BasicCOITree<PafMetadata, u32>>>,
}
// V1: 纯内存，不序列化；V2: bincode 整体持久化
```

- 每个 target 序列一棵独立的区间树——查询时 O(log n + k) 找到所有重叠 PAF 记录
- 索引时**不过滤**——每个 PAF 行都装入，不论质量（[[paf-route.md]] §2.3）
- 区间树节点 `PafMetadata` 只存 `u32` 坐标和 CIGAR 引用，不存序列名——
  为后续大 cohort / 显式图阶段预留升级空间
- 多文件索引：V1 不做。V2 参考 impg 的 `ImpgIndex` trait + `MultiImpg`
  （[[impg.md]] §3.4），通过 `ForestMap` 实现全局 target_id → 子索引的翻译

### 验证标准

- 索引后，用已知 region 查询，返回的 PAF 记录与 `grep` 原文件一致
- Zero Panic：空 PAF、重复行、超长行均不 panic

---

## 4. `pgr paf query <region>` — 查询层

### 两种模式

```bash
# 模式 1：区间投影（单跳）
pgr paf query chr1:1000-5000 -i cohort.paf.idx -o result.bed

# 模式 2：传递闭包 BFS
pgr paf query chr1:1000-5000 -i cohort.paf.idx --transitive -o result.paf
```

**模式 1（区间投影）**：在 target 序列的区间树上查找所有重叠 PAF 记录，把 target 坐标
lift 到 query 坐标。等价于"找所有直接比对此区间的序列"——类似 `pgr chain lift`
的单链线性投影，但在 all-vs-all 比对网络的并集上做区间树查找。

**模式 2（传递闭包 BFS）**：模式 1 的结果作为 BFS 第 0 层。对每个结果，以其 query 坐标
为新查询，在第 1 层区间树（query 此时变成 target）上继续查找。重复直到深度达到 `--max-depth`
（默认 2）或无新区间。

实现上依赖 `SortedRanges`（[[impg.md]] §3.3）——每轮 BFS 只把"未被现有区间覆盖的新增部分"
加入下一轮，避免重复遍历。参考 impg 实现：
- 索引构建：`Impg::from_multi_alignment_records()` (`impg-0.4.1/src/impg.rs:1549`)
- 单跳查询：`Impg::query()` (`impg-0.4.1/src/impg.rs:1848`)
- 传递闭包 BFS：`Impg::query_transitive_bfs()` (`impg-0.4.1/src/impg.rs:2291`) — 核心参考
- `SortedRanges` 实现：`impg-0.4.1/src/impg.rs:242`

### 过滤参数

所有参数在**查询时**生效，不影响索引：

| 参数 | 默认值 | 含义 |
|------|--------|------|
| `--merge-distance` | 0（不合并） | 同一序列上间距 ≤ D bp 的区间合并为一个（impg `-d`） |
| `--min-identity` | 0.0 | 最低 gap-compressed identity（从 `gi:f:` tag 读取） |
| `--min-output-length` | 0 | 最短输出区间长度（impg `-l`） |
| `--max-depth` | 2 | BFS 最大深度（impg `-m`，仅 `--transitive` 模式） |
| `--subset-sequence-list` | 全部 | 只保留指定序列上的结果 |

Cactus Caf 的过滤维度（Degree / Tree Coverage / Chain Length / Block End Trim）
见 [[paf-route.md]] §4.2 的分析。第一期只支持上述 5 个基础参数，
Caf 维度留作第二期的后处理过滤。

### 输出格式

**BED**（默认，最简输出）：
```
ref_name  1000  5000  query_name:50-450  gi:0.95
```

**PAF**（12 列 + 标签，对齐 impg `output_results_paf`）：
```
query_name  500  50  50  +  ref_name  1000  100  50  50  50  500  255
  gi:f:0.9500  bi:f:0.9200  cg:Z:10=5X2I3D
```

identity 计算（impg `main.rs:12042-12061`）：
```
gap_compressed_identity = matches / (matches + mismatches + #indel_events)
block_identity         = matches / (matches + mismatches + indel_bp_total)
```
- `gi`（gap-compressed）：每个 indel **事件**计 1 个差异——评估"同源性"
- `bi`（block）：每个 indel **碱基**计入差异——评估"序列一致性"

**MAF**（`--output-format maf`）：
- 对传递闭包结果调用 `fas consensus`（`libs/poa/` SPOA）做局部 MSA
- 输出标准 MAF 格式（复用 `libs/fas_multiz.rs`）

### 验证标准

- 区间投影结果与 `grep` 原 PAF 文件一致
- `--transitive` 结果包含所有直接同源 + 间接同源（BFS 深度内）
- `--merge-distance 100` 正确合并间距 ≤100bp 的区间
- `--min-identity 0.9` 正确过滤 `gi:f:` < 0.9 的记录
- 输出 MAF 中的序列与源基因组 FASTA 逐字节一致（path 保真）
- Zero Panic：无效 region 格式、缺失索引等返回友好错误

---

## 5. 已澄清事项与设计决策

以下问题在前序分析中已解决，记录结论供实现时参考。

1. **MAF → PAF 的 CIGAR 展开**（✅ 已实现）

   pgr 现有 MAF 解析器在 `src/libs/fmt/fas.rs`，`MafEntry.alignment: Vec<u8>` (line 246)
   含完整的带 gap 对齐串。从两条 `s` 行的 alignment 逐 base 对比即可生成
   `Vec<CigarOp>`（`=`/`X`/`I`/`D`）。无需引入第三方 MSA 库，纯文本比对即可实现。
   
   参考 impg 的 `CigarOp` bit-packing 设计 (`impg-0.4.1/src/impg.rs:73`)：
   3 位 op code + 29 位 length 压入 `u32`，紧凑且支持投影算术。

2. **多 PAF 文件的统一索引**（⏳ 首版简化，后续扩展）

   cohort 有多个 PAF 文件时，需统一装入同一个索引。impg 的 `ImpgIndex` trait +
   `MultiImpg`（[[impg.md]] §3.4）通过 `ForestMap` 实现全局 target_id → 子索引的翻译。
   
   首版 pgr 只支持单文件索引（类似 impg 的单文件模式），多文件索引留待第二步。
   此时不需要 `ForestMap` 和 trait 抽象。

3. **MAF strand → PAF strand 映射**（✅ 已实现）

   `MafEntry.strand` (fas.rs:257) 直接映射到 PAF 第 5 列（`+`/`-`）。
   反向互补 block 的坐标翻转沿用 `MafEntry::to_range()` (fas.rs:260-288) 的已验证逻辑。

4. **MAF score → `ms:i:` 标签**（✅ 已实现）

   `to-paf` 需从 MAF `a` 行解析 `score=` 字段（当前 `parse_maf_block()` 在 fas.rs:370
   忽略了 `a` 行内容）。转为 PAF 自定义标签 `ms:i:`，保留原始质量注释，
   使查询层可以按 score 过滤。

5. **CIGAR 存储策略**（✅ 已实现：存完整 CIGAR）

   首版 `to-paf` 阶段将 CIGAR 完整展开为 `cg:Z:` 标签写入 PAF，
   `paf index` 阶段将 `Vec<CigarOp>` 直接存入区间树节点，**不采用懒加载**。
   
   pgr 的 PAF 来自 MAF 转换（比对数远小于 all-vs-all），内存压力可控。
   这省去了 impg 的 BGZF file handle 管理 (`thread_local!` 缓存，`impg-0.4.1/src/impg.rs:43`)
   和 CIGAR seek 逻辑 (`paf.rs:68`)，大幅简化代码。

---

## 6. 测试参考：impg 已有测试可移植清单

impg 源码中包含丰富的单元测试，以下按 pgr 子命令分组，标注可直接移植或借鉴的测试。

### 6.1 `pgr maf to-paf` — PAF 行解析与输出验证

参考 impg 测试文件 `impg-0.4.1/src/paf.rs:364-416` 和 `impg-0.4.1/src/impg.rs:3125-3147`：

| 测试 | 来源 | 行号 | 移植说明 |
|------|------|------|---------|
| 有效 PAF 行（无 CIGAR） | `paf.rs` | 368 | 12 列纯坐标，验证坐标映射正确 |
| 有效 PAF 行（含 `cg:Z:`） | `paf.rs` | 394 | 验证 CIGAR 标签的偏移量计算 |
| 无效 PAF（字母污染数字字段） | `paf.rs` | 401 | 验证 zero-panic：返回 `Err` 而非 panic |
| 无效 CIGAR（非法 op `Q`） | `paf.rs` | 409 | 验证 CIGAR 解析的健壮性 |
| 完整 PAF 行解析为 `AlignmentRecord` | `impg.rs` | 3125 | 含 `SequenceIndex` name→id 映射，验证 `data_offset` 和 `data_bytes` |

**pgr 适配**：pgr 的 `to-paf` 输出需经过同样的解析器验证。但 pgr 的 PafRecord 存 12 完整字段（而非 impg 的 8 字段 + 懒加载 offset），验证时要额外检查 col 10-12。

### 6.2 `pgr paf index` — 索引构建验证

参考 impg 的 `impg-0.4.1/src/impg.rs:3125-3147`（同上）。impg 没有独立的索引构建单元测试——索引正确性通过"查询结果与 grep 原文件一致"来间接验证。pgr 可以采用同样的策略：

- 用已知的 PAF 输入构建索引
- 用 `paf query` 单跳模式查询
- 断言返回记录数与 `grep` 一致

### 6.3 `pgr paf query` — 坐标投影与传递闭包

参考 impg 测试文件 `impg-0.4.1/src/impg.rs:2930-3105`，**这是最关键的测试集**。

#### 投影正确性测试（6 组 9 个子测试）

| 测试 | 行号 | 场景 | 断言 |
|------|------|------|------|
| `test_project_target_range_through_alignment_forward` | 2931 | 正向，全 match | query 坐标 = target 坐标平移 |
| `test_project_target_range_through_alignment_reverse` | 2942 | 反向，全 match | query_end < query_start（递减） |
| `test_project_target_range_through_alignment` | 2953 | **混合 CIGAR**（`=`×3 + `I`×2 + `D`×1），**6 个子区间** | 部分重叠、穿越 insertion、跳过 deletion、边界裁剪 |
| `test_forward_projection_simple` | 3008 | 正向 sanity | 100% 覆盖 |
| `test_reverse_projection_simple` | 3022 | 反向 sanity | 坐标反转 |
| `test_forward_projection_with_insertions` | 3036 | `I` 跨越 | query 区段跨越 insertion 正确延长 |
| `test_forward_projection_with_deletions` | 3051 | `D` 跨越 | target 跨越 deletion 时 query 坐标不推进 |
| `test_reverse_projection_with_mixed_operations` | 3066 | 反向 + `D`/`I`/`=` 混合 | 反向递减 + CIGAR 调整 |
| `test_edge_case_projection` | 3087 | CIGAR 含 `=`/`D`/`X`/`I` 混杂，只查前 10bp | CIGAR 正确裁剪 |

**第三个测试是最重要的**——用同一个 CIGAR（`10=5I5D50=50I35=`）和 6 组不同的 target 区间，验证部分重叠/穿越 gap/裁剪等所有情况。pgr 的 `paf query` 单元测试应以此为模板。

#### CIGAR 解析测试

| 测试 | 行号 | 输入 | 验证 |
|------|------|------|------|
| `test_parse_cigar_to_delta_basic` | 3108 | `"10=5I5D"` | `[CigarOp(10,=), CigarOp(5,I), CigarOp(5,D)]` |

**pgr 适配**：pgr 的 CIGAR 来自 MAF 对齐串展开，而非从 PAF `cg:Z:` 标签解析。测试侧重点不同：
- pgr 需新增 `test_maf_alignment_to_cigar` —— 从两条 `s` 行的 alignment 对比生成 `Vec<CigarOp>`
- CIGAR bit-packing 和字符串互转逻辑与 impg 一致，可直接复用 impg 的测试结构

#### CIGAR 反转测试（pgr V1 可跳过，但建议了解）

| 测试 | 行号 | 场景 | 说明 |
|------|------|------|------|
| `test_invert_cigar_forward_strand` | 3150 | I↔D 交换，顺序不变 | 双向索引的 A→B 反向映射 |
| `test_invert_cigar_reverse_strand` | 3172 | I↔D 交换 + 数组反转 | 含反向链的反向映射 |
| `test_invert_cigar_empty` | 3191 | 空数组 | 边界 |
| `test_invert_cigar_matches_only` | 3201 | 纯 match，无 indels | 退化 |

pgr V1 不做双向索引，这些测试不需要移植。但理解其逻辑有助于后续 debug 跨链向的传递闭包边界情况。

### 6.4 pgr 测试策略建议

```
层次 1 — 单元测试（libs/paf/）
  ├── record.rs      : PafRecord 构造 + CigarOp bit-packing（参考 impg.rs:3107+）
  ├── parser.rs      : 纯文本 PAF 解析（参考 paf.rs:364-416）
  ├── cigar.rs       : CIGAR ↔ 字符串互转 + identity 计算
  │   └── + test_maf_alignment_to_cigar（pgr 独有：对齐串 → CigarOp）
  └── index.rs       : 区间树查询 + 坐标投影（参考 impg.rs:2930-3105）

层次 2 — 集成测试（tests/cli_paf.rs）
  ├── maf to-paf     : 输入两序列 MAF → 断言输出 PAF 格式正确
  ├── paf index      : 索引构建 → paf query 对比 grep
  └── paf query      : 单跳 / --transitive 端到端验证
```

---

*last updated: 2026-06-27*
