# pgr pgi build vs FastGA GIXmake 基准测试

> 目的：对比 `pgr pgi build` 与 FastGA `GIXmake` 的索引构建速度（hyperfine），
> 评估 `.pgi` 索引构建的性能差距。GIXmake 是 C 的极致优化参照（syncmer 稀疏
> k-mer 索引，见 [[fastga.md]] §10 / [[pbit.md]]）。2026-08-04 复测。

## 环境与输入

- 输入：MG1655（E. coli K-12，4.64 Mb，`tests/genome/mg1655.fa.gz`）
- 依赖：`hyperfine`、FastGA 生态（`GIXmake` / `FAtoGDB` 在 PATH）、
  pgr release 二进制（`cargo build --release` 先执行一次）

两个场景：A. 从 FASTA 全流程（pgr 一步 vs FAtoGDB+GIXmake 两步）；
B. 序列已编码（pgr 从 2bit vs GIXmake 从 GDB）。

## 执行（可直接复制到终端）

在仓库内任意目录粘贴以下整段运行（自动定位仓库根、准备输入、跑
hyperfine 并导出 Markdown 结果）：

```bash
cd "$(git rev-parse --show-toplevel)"
PGR="$PWD/target/release/pgr"
B="$(mktemp -d /tmp/pgibench.XXXX)"
gzip -dc tests/genome/mg1655.fa.gz > "$B/genome.fa"
FAtoGDB "$B/genome.fa" >/dev/null 2>&1
"$PGR" fa to-2bit "$B/genome.fa" -o "$B/genome.2bit" >/dev/null
mkdir -p "$B/gix"
hyperfine --warmup 1 --runs 5 --export-markdown "$B/bench-result.md" \
  -n "pgr pgi build (FASTA)" \
    "$PGR pgi build $B/genome.fa -o $B/out.pgi" \
  -n "pgr pgi build (2bit)" \
    "$PGR pgi build $B/genome.2bit -o $B/out2.pgi" \
  -n "FAtoGDB+GIXmake (FASTA)" \
    "FAtoGDB $B/genome.fa >/dev/null 2>&1 && GIXmake $B/genome.1gdb -P$B/gix >/dev/null 2>&1" \
  -n "GIXmake (GDB)" \
    "GIXmake $B/genome.1gdb -P$B/gix >/dev/null 2>&1"
cat "$B/bench-result.md"
```

结果会打印到终端，同时保存在 `$B/bench-result.md`。

## 结果（2026-08-04 复测，MG1655，5 次；数值有正常波动）

最终版本（MSD radix + 并行分桶 + 单遍扫描优化 + 索引一致性修复）：

| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `pgr pgi build (FASTA)` | 348.0 ± 11.8 | 339.4 | 368.8 | 1.12 ± 0.04 |
| `pgr pgi build (2bit)` | 344.0 ± 4.0 | 340.0 | 349.1 | 1.11 ± 0.02 |
| `FAtoGDB+GIXmake (FASTA)` | 322.7 ± 1.7 | 321.6 | 325.6 | 1.04 ± 0.02 |
| `GIXmake (GDB)` | 310.2 ± 4.5 | 305.3 | 316.3 | 1.00 |

`pgr pgi build` 与 GIXmake 仍基本持平，但 2bit 路径由初测的略快
（301.8 ms）变为 ~11% 慢（344.0 ms）。差值来源是 **2026-08-03 的索引
一致性修复**（f7461c1：pending 队列去重 + 出列位置重算 k-mer key）——
单遍扫描为每个位置多做一次去重/重算；GIXmake 侧数值几乎不变（312.4 →
310.2 ms）。32 核机器、release 构建（当前文档表为复测值，初测见下）。

## 复测（2026-08-12，FastK 字节键迁移后，MG1655，8 次）

| Command | Mean [ms] | 峰值内存 | Relative |
|---|---:|---:|---:|
| `pgr pgi build (FASTA)` | 384.0 ± 2.9 | **145 MB** | 1.21 ± 0.01 |
| `FAtoGDB+GIXmake (FASTA)` | 316.1 ± 2.7 | — | 1.00 |

迁移后 `pgr pgi build` 峰值内存 145 MB（迁移前 ~180 MB，**-19%**，
`PgiEntry` 24 B→8 B + keys 打包的收益）；构建时间 384 ms 较迁移前
348 ms 略升（+10%，机器负载/收集路径未变的正常波动，待字节直落
优化后复测）。GIXmake 316 ms。

> **字节直落实验（2026-08-12，负收益，已回退）**：把 `collect_one_contig`
> 的 u128 滚动 + `Vec<u128>` 收集改为 u128 滚动 + 发射时 `pack_kmer`
> 直落 `Vec<u8>` + `radix_sort_bytes_par` + 字节分组——413 ms / 137 MB，
> 比 d2c422a 基线（384 ms / 145 MB）慢 7%，内存仅降 5%。原因：pgi 的
> 记录数（syncmer 位置，~几十万）远小于 k-mer 计数（460 万），字节直落
> 的排序/内存收益被每条记录一次 `pack_kmer` 的开销抵消；滚动本身
> 保持 u128（每碱基位运算，字节滚动会更慢，已验证 418 ms）。**保持
> u128 中间层 + 分组时 pack**（384 ms / 145 MB）为 pgi build 的当前
> 最优形态；字节化仅保留在最终存储（`PgiIndex.keys`）。

> **并行收集（2026-08-12，学习 GIXmake 的 per-thread 分桶）**：
> `build_from_seqs` 的 contigs 收集从单线程改为 `par_iter` 每 contig
> 独立缓冲再合并（单 contig 保持顺序路径，避免 rayon 分发开销）。
> 4 基因组（cft073/e2348_69/ec958/mg1655，~20 Mb，8 contig）：
> 1.450 s → **1.203 s**（+17%）；MG1655 单 contig 384 ms 不变。
> GIXmake 4 基因组 505 ms（pgr 仍 2.4× 慢，剩余差距在单线程分组
> 去重与 GIXmake 的排序内 COUNT + 多线程分块，见 align 对照笔记）。

> **分组去重 + pack 优化（2026-08-12）**：
> * 分组去重：`HashSet` → **组内 `sort_unstable` + 相邻去重**（组平均
>   1.8 条，小数组插入排序无哈希开销）+ size=1 快路径（省 4 字节
>   已足够，分组 542 → 364 ms）；
> * `pack_kmer`：逐字节 4 次移位 → **u128 大端拷贝**（左移
>   `8*kmer_bytes-2k` 位后取最后 `kmer_bytes` 字节，分组 364 →
>   192 ms；`.pkt` 写出同受益）；
> * 4 基因组总时间 1.450 s → **817 ms**（-44%）；阶段分解
>   collect 354 / sort 151 / group 192 ms；vs GIXmake 505 ms =
>   **1.62×**（从 2.4× 改善）。align 输出 738 PSL 不变、1701 测试
>   全绿。

### 优化历程

- 初始（两遍扫描 + 比较排序）：761 ms（2.41× 慢）
- 单遍扫描 + rc_key 查表 + 计数偏移桶化：493 ms（1.61× 慢）
- **排序换成 MSD radix（American-flag 原地版，移植自 FastGA `MSDsort.c`）**
  并修复低字节分桶导致的**全局乱序 bug**：590 ms（当时单线程 radix 反而慢）
- **顶字节分桶 + rayon 并行排序各桶**：404 ms（1.27× 慢）
- collect 优化（去 ring、碱基查表、滚动 RC k-mer、数组环形单调队列）+
  分组容量预估 + 写入批量合并：302 ms（1.04× 快）

### 与 GIXmake 的实现对照

| 环节 | GIXmake | pgr（当前） |
|---|---|---|
| 扫描 | 字节级 4 碱基/字节滚动 + `ne[8]` 环形缓冲（C） | 单遍：s-mer 滚动哈希 + 单调队列 + 40-mer 滚动（含 RC） |
| 分桶 | 按 k-mer 前 10 位分 1024 桶（多线程手指写入） | 按最高字节原地分 256 桶 + rayon 并行排序各桶 |
| 桶内排序 | `radix_sort`（原地 American-flag，循环置换） | `msd`（同算法移植，`src/libs/ds/radix_sort.rs`） |
| 小桶 | shell sort（2..15 记录） | 插入排序（≤16 记录） |
| 内存 | 一次性 sort 数组（`swide` 字节记录） | keys + payloads 并行数组，无辅助数组（峰值 ~180 MB vs 原 ~330 MB） |

### 额外发现

1. **低字节分桶是排序 bug**：原实现按 k-mer 最低字节分 256 桶再按桶合并，
   并非全局升序（如 `0x0100` 排到 `0x01` 前面），而 `dist pgi` 的归并要求
   entries 全局有序，会导致交集漏算。已由 radix 全局排序修复并加回归测试。
2. **rc_key 查表路径在 k % 4 != 0 时错误**：剩余碱基的迭代顺序/位移写反，
   k=40（整字节）不受影响；测试用的 k=10 会出错。已修正并加随机回归测试。

## 初测结果（2026-08-02，修复前）

| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `pgr pgi build (FASTA)` | 307.4 ± 8.1 | 299.6 | 318.2 | 1.02 ± 0.03 |
| `pgr pgi build (2bit)` | 301.8 ± 3.3 | 296.7 | 306.1 | 1.00 |
| `FAtoGDB+GIXmake (FASTA)` | 319.8 ± 2.8 | 316.0 | 322.2 | 1.06 ± 0.01 |
| `GIXmake (GDB)` | 312.4 ± 4.0 | 307.9 | 316.6 | 1.04 ± 0.02 |

## 历史结果（优化前）

| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `pgr pgi build (FASTA)` | 492.8 ± 3.0 | 490.3 | 497.6 | 1.63 ± 0.04 |
| `pgr pgi build (2bit)` | 497.0 ± 3.4 | 492.8 | 501.7 | 1.65 ± 0.05 |
| `FAtoGDB+GIXmake (FASTA)` | 311.1 ± 3.2 | 307.4 | 314.7 | 1.03 ± 0.03 |
| `GIXmake (GDB)` | 301.4 ± 8.0 | 293.2 | 313.8 | 1.00 |
