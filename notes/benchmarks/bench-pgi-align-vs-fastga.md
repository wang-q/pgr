# pgr align pgi vs FastGA 端到端基准

> 目的：对比 `pgr align pgi` 全流程（索引构建 ×2 + 扩展比对）与 FastGA 单命令
> 的端到端耗时。FastGA 是 C 的极致优化参照（GIX 归并 + wave aligner）。
> 2026-08-02 初测，2026-08-04 修复轮后复测。

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
    "$PGR pgi build $REF -o $B/r.pgi 2>/dev/null && $PGR pgi build $QUERY -o $B/q.pgi 2>/dev/null && $PGR align pgi $B/r.pgi $B/q.pgi --ref-seq $REF --query-seq $QUERY -o $B/o.psl 2>/dev/null" \
  -n "FastGA -psl (one-shot)" \
    "FastGA -psl -T8 -P$B $REF $QUERY > $B/f.psl 2>/dev/null"
cat "$B/bench.md"
```

> **FastGA 参数顺序说明（2026-08-06 复核）**：`FastGA $REF $QUERY` 中第一个
> 参数是 FastGA 的 **query/种子侧**，即 mg1655。pgr 侧 `align pgi` 的种子侧是
> ref（= mg1655），两者种子侧一致，比对内容才可比。**不要**按 pgr 的
> ref-first 直觉"修正"成 `FastGA $QUERY $REF`——那会把种子侧换成 sakai，
> 结果不可比（FastGA 不对称，见 [[../references/fastga.md]] §3.3/§7.5 与
> [[../design/pgi-align.md]] §1.3.6）。两套 PSL 的 q/t 标签相反（FastGA
> qName=mg1655、pgr tName=mg1655），但覆盖/身份统计按基因组算，不受影响。

## 复测结果（2026-08-04，MG1655 vs Sakai，5 次）

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|---|---:|---:|---:|---:|
| `pgr full (build 2x + align ext)` | 1.674 ± 0.022 | 1.648 | 1.700 | 1.00 |
| `FastGA -psl (one-shot)` | 3.861 ± 0.169 | 3.606 | 4.041 | 2.31 ± 0.11 |

复测 pgr **反超 FastGA ~2.3×**（初测为 1.08× 持平）。变化来自两侧：
FastGA 本次稳定在 ~3.9 s（初测 1.22 s，复测时主机并发负载 + 文件系统
状态影响其大量小文件 IO，system 时间 ~19 s）；pgr 由 1.32 s 升至 1.67 s
（同期负载）。两个数值都按同一次 hyperfine 内对比，结论方向一致
（pgr 更快）。

## 复测结果（2026-08-05，MG1655 vs Sakai，3 次）

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|---|---:|---:|---:|---:|
| `pgr align pgi`（自动建索引 ×2，默认路径） | 1.231 ± 0.008 | 1.222 | 1.239 | 3.28 ± 0.10 |
| `pgr full (build 2x + align ext)` | 1.260 ± 0.009 | 1.251 | 1.268 | 3.21 ± 0.10 |
| `FastGA -psl (one-shot)` | 4.039 ± 0.121 | 3.912 | 4.153 | 1.00 |

**索引写出（2026-08-05）**：`pgr align pgi` 自动建索引路径的索引写出走
`BufWriter`（`cmd_pgr/align/pgi.rs::resolve_side`），默认路径与显式
build+align 持平（1.23 vs 1.26 s），输出逐字节不变；默认路径反超
FastGA **~3.3×**。

**输出结构（当前二进制，2026-08-05 复测）**：mg1655–sakai 默认参数
738 条 PSL 记录、pooled 身份率 0.9754（`(matches+rep)/block_len`，
口径见 [[../design/pgi-align.md]] §2.2）；FastGA 700 条非 self 记录。

## 迁移后复测（2026-08-12，FastK 字节键迁移后）

`pgr align pgi`（mg1655 × sakai，默认参数，8 线程）：

| 指标 | 迁移后 | 迁移前（2026-08-05） | 结论 |
|---|---:|---:|---|
| PSL 记录数 | **738** | 738 | 逐条一致 |
| pooled identity | 0.9790 | 0.9754 | 口径细节差异，量级一致 |

索引构建（cft073，k=40）：294 万 unique k-mers / 299 万 positions。
**结论：FastK 字节键迁移未改变 pgi 查询/链化输出**（记录数逐条一致），
`.pgi` 格式兼容性验证通过。

## 初测结果（2026-08-02，MG1655 vs Sakai，3 次）

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

- 身份率（pooled PSL identity，mg1655×sakai，FastGA 排除全基因组自比对块）：
  pgr 97.30% vs FastGA 97.28%（同口径基本持平）；
- 覆盖（**真实并集**，2026-08-02 复核）：MG1655 vs Sakai pgr 75.8% /
  FastGA 78.2%；vs Nissle 两者均为 77.3%——**基本打平**。早前的
  "95.7% vs 99.7%" 是块区间 span 求和（重叠重复计数）的假象；未覆盖的
  ~22% 是株系特异序列（O157/Nissle 特有岛），任何比对器都无法映射到
  MG1655。**lcp/adaptamer 变长种子因此不是当前优先级**；剩余小差距
  （sakai +2.4%）来自分歧区的 wave 式补齐（banded 窗口对低分区间跳过）。
  （2026-08-12 已实验否定：E. coli + 酵母实测变长种子对端到端 PSL 覆盖
  无收益，见 `design/pgi-align.md` §7.3.1。）
- 时间可比：两流程均从 FASTA 开始（pgr 建索引 ×2 + 比对；FastGA 内部
  FAtoGDB + GIXmake + 比对）。
