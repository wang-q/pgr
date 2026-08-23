# supermer 质量门控设计（pgr 任务 1.2）

> 目标：让 pgr `kmer::supermer` 支持与 anchr direct 路径一致的质量门控，
> 使 anchr 的 FASTQ 计数可以统一走 supermer（去掉 direct 回退）。
> 对应 `~/Scripts/anchr/notes/design/pgr-tasks.md` §1.2。

## 1. 背景

- anchr 现状：FASTA 默认 supermer，FASTQ 自动回退 direct（`--no-supermer`
  强制 direct），因为 pgr supermer 路径无质量过滤。
- 落地后：FASTQ 走 supermer，重跑 `asm-unitig.md` §12 基准表后决定
  `--supermer` 转正。

## 2. 现状核对

### anchr direct 质量语义（`libs/asm/table.rs::count_read_kmers_packed`）

**窗口级概率硬门控**，非 quorum 混合聚合：

- `prob_error(q)`：q=0 → 0.75，q=1 → 0.7，其余 `10^(-0.1q)`（BBTools 表）。
- 窗口正确概率 `prob = Π (1 - prob_error(q_i))`，滑动维护
  `prob *= prob_correct[q_i]` / `prob *= prob_correct_inv[q_{i-k}]`；
  N 或低质量前导碱基重置 `prob = 1`。
- 发射条件：窗口无 N 且 `prob >= min_prob`（`min_prob <= 0` 或空质量时不过滤）。
- 低质量窗口直接跳过，**没有** `quality.rs` 的 high/low 混合计数。

### pgr supermer 结构（`src/libs/kmer/supermer.rs`）

- stage-1：`pack_sequence_into` 按 N 划分 N-free stretch，`pack_run` 在
  stretch 内按 minimizer 单调 + `max_super` 切分 run，`emit_span` 记录
  (span 碱基, 窗口数)。
- stage-2：展开唯一 span → 加权 k-mer → radix 排序 → 权重累加。
- collapse 收益依赖"run 内窗口全部计入"。

### 与 `kmer/quality.rs` 的区别（重要）

`quality.rs` 是 quorum `hash_with_quality` 移植（`thresh` 硬门控 +
high/low 混合聚合），语义与 anchr direct 路径**不同**，不能直接复用；
supermer 质量门控应对齐 anchr 的 `min_prob` 语义。

## 3. 方案：窗口级布尔门控 + run 切断

质量门控是窗口级布尔判定（每窗口 计入/跳过），与 supermer 的 run 结构
兼容。做法：

1. 在 `pack_sequence_into` 内维护滑动 `prob`（与 anchr 完全相同的计算
   顺序），得到每个窗口起点的有效性。
2. 无效窗口起点当作 stretch 边界：把 N-free stretch 切成
   `[seg_start, i)` / `[i+1, ...)` 的连续有效段，段内所有窗口有效。
3. `pack_run` / `emit_span` / stage-2 完全不变——段只是更短的 stretch。

等价视角：无效窗口起点 = "虚拟 N"；切段只损失 stage-1 collapse 收益，
不改变计数语义。

## 4. API 设计

新增（保持现有无质量 API 不变）：

```rust
/// Super-mer counting with a sliding-window quality gate (anchr `min_prob`
/// semantics); `min_prob <= 0.0` or empty qualities disable the gate.
pub fn build_table_slices_qual(
    seqs: &[&[u8]], quals: &[&[u8]], k: usize, min_prob: f32,
) -> anyhow::Result<KmerTable>

pub fn build_table_slices_qual_with_m(
    seqs: &[&[u8]], quals: &[&[u8]], k: usize, m: usize, min_prob: f32,
) -> anyhow::Result<KmerTable>
```

- `prob_correct` / `prob_correct_inv` 表（128 项）放 `kmer` 模块内部私有
  （anchr 已有自己的表，不需要导出）。
- 外层 `quals.len() != seqs.len()`：`anyhow::bail!`（Zero Panic；anchr 的
  `count_read_kmers_packed` 会越界 panic，pgr 不复制该行为）。
- 单条 read 质量为空：该 read 不过滤（与 anchr per-read 的
  `!quals.is_empty()` 规则一致）；非空质量必须与序列等长，否则 bail。

## 5. 逐字节一致性要点

- `prob_error` 表与 BBTools 完全一致（0/1 特例 + `10^(-0.1q)`）。
- f32 运算顺序与 anchr 一致：`1.0 - p`、`1.0 / c`、滑动乘除、`prob >= min_prob`
  比较——f32 非结合，任何重排都会改变窗口集合。
- 输出表（keys + counts）应与"质量发射 + `count_keys`"逐字节一致。

## 6. 验证

- pgr 单测：随机 (seq, qual) 输入，`build_table_slices_qual` 与
  "质量发射 + `count_keys`"（reference 实现）逐字节对比；`min_prob = 0`
  与 `build_table_slices` 一致；空质量输入一致；qual 过短 bail。
- anchr 侧（落地后）：bump rev → `verify-migrate.sh` → `asm-gate.sh`
  smoke/single；FASTQ 去 direct 回退，重跑 `asm-unitig.md` §12 基准表。

## 7. 性能权衡

- 有效窗口数不变 → stage-2 展开量不变；run 变短 → stage-1 去重粒度变细、
  记录数可能增加。低质量数据上 supermer 相对 direct 的优势缩小，但不会
  劣于 direct（折叠退化到逐窗口）。
- 质量计算成本：每碱基 2 次 f32 乘除 + 查表，叠加在打包热路径上；anchr
  direct 路径已有同样成本，两者可比。

## 8. 待决策点

1. 语义用 anchr `min_prob`（推荐，替换目标就是它）还是 quality.rs
   `thresh`。
2. API 形态：新增 `_qual` 变体（推荐）vs 并入现有函数破坏签名。
3. 是否同时提供 Vec 版 `build_table_qual`——anchr 只用 slices 版，暂不做。
