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

---

## 2. `pgr maf to-paf` — 数据源

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

参考 impg 的 `Impg` struct（[[impg.md]] §3.1）：

```
PafIndex {
    seq_index: SequenceIndex,           // name ↔ u32 ID 双向映射
    trees: FxHashMap<u32, Arc<Coitree<PafRecord>>>,  // target_id → 区间树
    source_files: Vec<String>,          // 源 PAF 文件列表（用于 CIGAR 懒加载）
}
```

- 每个 target 序列一棵独立的区间树——查询时 O(log n + k) 找到所有重叠 PAF 记录
- 索引时**不过滤**——每个 PAF 行都装入，不论质量（[[paf-route.md]] §2.3）
- 多文件索引：参考 impg 的 `ImpgIndex` trait + `MultiImpg`（[[impg.md]] §3.4），
  通过 `ForestMap` 实现全局 target_id → 子索引的翻译

### 验证标准

- 索引后，用已知 region 查询，返回的 PAF 记录与 `grep` 原文件一致
- 多文件索引的 seq_index 中无重复序列名
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
加入下一轮，避免重复遍历。

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

## 5. 待澄清问题

1. **MAF → PAF 的 CIGAR 展开精度**：需要从 MAF `s` 行的对齐串逐 base 判断 `=`/`X`/`I`/`D`。
   pgr 现有 `maf` 模块（`libs/fmt/maf.rs`）已能解析 `s` 行，但不对齐串做逐 base 分析。
   `to_fas` 只提取坐标和 consensus，不做 CIGAR 展开。实现 `maf to-paf` 需要增强解析器，
   或接受用 `M`（不区分 match/mismatch）代替 `=`/`X` 作为首版折中。

2. **多 PAF 文件的统一索引**：cohort 有多个 PAF 文件（每对一个），需统一装入同一个 `PafIndex`。
   impg 的 `ImpgIndex` trait + `MultiImpg`（[[impg.md]] §3.4）通过 `ForestMap` 实现全局
   target_id → 子索引的翻译。pgr 可直接复用此模式。

3. **MAF 的链向信息**：MAF `s` 行 strand `+`/`-` 可直接映射到 PAF 第 5 列。反向互补 block
   的坐标转换需要验证（特别是 query 和 target 同时涉及 `-` strand 时的坐标翻转逻辑）。

4. **MAF score → PAF 自定义标签**：PAF 标准没有 score 列。建议在 `maf to-paf` 转换时把
   MAF `score=` 存入 `ms:i:` 标签。这是 pgr 相对原生 PAF 的增值——保留了原始 MAF 的质量注释，
   使查询层可以按 score 过滤而不仅按 identity。

5. **CIGAR 懒加载的 file handle 管理**：查询层需要跨多个源 PAF 文件读取 CIGAR。
   impg 用 `thread_local!` 缓存 file handle（[[impg.md]] §9.5），pgr 需要类似机制。
   首版可选择"索引时直接存 CIGAR"以跳过此问题（用内存换简单性），后续再优化为懒加载。
