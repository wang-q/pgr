# pgr pgi align vs FastGA 端到端基准

> 目的：对比 `pgr pgi align` 全流程（索引构建 ×2 + 扩展比对）与 FastGA 单命令
> 的端到端耗时。FastGA 是 C 的极致优化参照（GIX 归并 + wave aligner）。

## 环境与输入

- 输入：MG1655 vs Sakai（E. coli K-12 vs O157:H7，各 ~4.6/5.5 Mb）
- 依赖：`hyperfine`、FastGA（PATH 中）、pgr release 二进制
- 机器：32 核

## 执行（可直接复制）

```bash
cd "$(git rev-parse --show-toplevel)"
cargo build --release
PGR="$PWD/target/release/pgr"
B="$(mktemp -d /tmp/pgialign.XXXX)"
REF="$PWD/tests/genome/mg1655.fa.gz"
QUERY="$PWD/tests/genome/sakai.fa.gz"
hyperfine --warmup 1 --runs 3 --export-markdown "$B/bench.md" \
  -n "pgr full (build 2x + align ext)" \
    "$PGR pgi build $REF -o $B/r.pgi 2>/dev/null && $PGR pgi build $QUERY -o $B/q.pgi 2>/dev/null && $PGR pgi align $B/r.pgi $B/q.pgi --ref-seq $REF --query-seq $QUERY -o $B/o.psl 2>/dev/null" \
  -n "FastGA -psl (one-shot)" \
    "FastGA -psl -T8 -P$B $REF $QUERY > $B/f.psl 2>/dev/null"
cat "$B/bench.md"
```

## 结果（2026-08-02，MG1655 vs Sakai，3 次）

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|---|---:|---:|---:|---:|
| `pgr full (build 2x + align ext)` | 1.322 ± 0.006 | 1.317 | 1.328 | 1.08 ± 0.01 |
| `FastGA -psl (one-shot)` | 1.220 ± 0.010 | 1.212 | 1.232 | 1.00 |

优化后与 FastGA **基本持平**（1.08×）。优化前为 2.63s（2.24× 慢）。

## 优化历程（align 扩展阶段）

| 版本 | 自比对 | 跨株 | 说明 |
|---|---:|---:|---|
| v2（单窗口、全列扫描） | 37.8s | 2.0s | 自比对 332 窗口的主链在单线程串行 |
| v3（分窗） | 37.8s | 2.0s | 主链窗口仍在一线程内串行 + 每窗 DP 扫全列 |
| v3.1（带限列循环） | ~12s* | ~2s | DP 只扫 `|j-i-diag0|≤band` 的 65 列而非 16000 列 |
| v3.2（窗口摊平并行） | **0.84s** | **0.66s** | 所有链的所有窗口进入同一 rayon 流，负载均衡 |

\* 估算：带限列循环修复后未单独计时（与并行修复同批落地）。

**两个关键发现**：

1. **banded DP 内层必须按带限列迭代**：原实现对每行扫描全部 m 列再做带外
   过滤，band=32 时 16000 列里只有 65 列有效，浪费 246×。
2. **窗口级负载均衡**：rayon 按"链"分片时，自比对的主链（332 窗口）成为
   单线程长尾任务（37s），其余 744 条小链瞬间完成；把 (链, 窗口) 摊平成
   单一并行流后全部窗口均匀分布。

两者修复后自比对 45×、跨株 3× 加速，输出与修复前逐字节一致。

## 对照说明

- 身份率：pgr 98.42% vs FastGA 97.83%（pgr banded 局部取精确核心，略高）；
- 覆盖（**真实并集**，2026-08-02 复核）：MG1655 vs Sakai pgr 75.8% /
  FastGA 78.2%；vs Nissle 两者均为 77.3%——**基本打平**。早前的
  "95.7% vs 99.7%" 是块区间 span 求和（重叠重复计数）的假象；未覆盖的
  ~22% 是株系特异序列（O157/Nissle 特有岛），任何比对器都无法映射到
  MG1655。**lcp/adaptamer 变长种子因此不是当前优先级**；剩余小差距
  （sakai +2.4%）来自分歧区的 wave 式补齐（banded 窗口对低分区间跳过）。
- 时间可比：两流程均从 FASTA 开始（pgr 建索引 ×2 + 比对；FastGA 内部
  FAtoGDB + GIXmake + 比对）。
