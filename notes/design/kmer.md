# kmer 模块设计：原生 k-mer 计数 / profile / run 提取（rept s-kmer / e-kmer）

> 状态：**已实现（2026-08-09）**。本笔记记录 `libs/kmer/` 的设计与行为契约：
> 让 `pgr rept s-kmer` / `pgr rept e-kmer` 不再依赖外部 FastK / Profex。
> 参考源码：仓库内 `FASTK-master/`（2025-09-13 下载的快照，
> **等于 FASTK-1.2**，README 标注 Current: April 18, 2021，**不是上游当前
> master**，见 §2.3 版本核对）；行为契约以本机安装的 FastK / Profex
> （CBP 安装，2025-03 构建，上游 commit ddea6cf，**无源码补丁**，见 §2.3）
> 实测为准。
> **旧管线视为实验性参考**，迁移不要求输出与旧管线字节级一致，不保留已知
> 缺陷（尾 run quirk，见 §2.2）。

> **2026-08-09 后续（`pgr kmer` 命令组 + 三种格式）**：在 rept 用途之上新增
> 顶级命令组 `pgr kmer`（table/profile/hist）与三种持久化格式——
> `.pkt`（表，原 `.pgrk` 改名）、`.pkp`（profile，pgr 自有单文件）、
> `.hist`（直方图，**FASTK 字节兼容**）。兼容性决策：直方图单文件、代价
> 小（~50 行），做兼容（外部可比较：Histex/KatGC/GenomeScope 直接可读，
> 已实测与 FastK 输出逐行一致）；profile 因 FASTK `.prof` 是 stub + 分片
> 多文件布局（`.pidx.N`/`.prof.N` + RLE），**不做兼容**，用 `.pkp` 单文件
> 自有格式（header + raw u16）；`.ktab` 维持不做（§3.0）。详见 §10。

## 1. 现状与目标

现管线（`src/libs/pl/repeat.rs::run_repeat_pipeline`）：

```text
s-kmer: FastK -p -k17 -Ngenome → genome.prof → Profex -z genome <sn> → .rg → runlist
e-kmer: FastK -t（repeat 库建表）→ FastK -p:repeat -k17 → genome.prof
        → Profex -z genome <sn> → .rg → runlist
```

目标：用 `src/libs/kmer/` 原生实现替换 FastK（k-mer 计数 + profile 生成）和
Profex（profile → run 提取），其余（runlist cover/fill/excise/fill、点号染色体名
映射、`--keep-index` 缓存机制、tempdir 机制）保持现状（缓存文件格式另见
§4.3）。磁盘分桶与 `.ktab/.prof` 外部格式兼容不做（那是 FASTK 为 TB 级数据
设计的，pgr 场景用不上）；**super-mer/minimizer 两段计数 2026-08-14 起实现**
（`unitig-bucket.md` 阶段 B，pgr 侧计数层，见 §3.6）。

## 2. 行为契约（已实测）

> 契约以本机安装的 FastK / Profex（CBP 安装，pgr 现有管线实际依赖的版本，
> 对应上游 commit ddea6cf）实测为准。源码核对：
> 仓库内 `FASTK-master`（= FASTK-1.2 快照）与 CBP 源码包
> `/tmp/cbp_fastk_check.*/fastk/`（= ddea6cf）逐文件比对，结论见 §2.3：
> 两版**计数核心语义一致**，差异仅限 Profex 输出细节与周边小修。

### 2.1 FastK 计数语义

* k-mer 编码为 2-bit，**canonical**（正向与反向互补取字典序小者），大小写不敏感。
* 窗口含 N 则该 k-mer 无效（FastK 按 gap 拆分），profile 对应位置值为 0
  （作为 run 分隔符，见 §2.2）。
* `-p`（s-kmer）：profile 第 i 个值 = 该 k-mer 在整个数据集（全部染色体）
  中的频次（≥1，跨 read 累计，实测确认）。
* `-p:<table>`（e-kmer）：profile 第 i 个值 = 该 k-mer 在 repeat 表中的 count，
  不在表中为 0（已用双拷贝探针验证：值取表内 count 而非基因组 count）。
* profile 长度 = n - k + 1。

### 2.2 Profex 输出与 pgr 解析

实测输出（`Profex -z genome 1`）：

```text
Read 1:
     0 -   316 (1)
   300 -   420 (2)
   404 -   536 (1)
   520 -   640 (2)
   624
```

* run 行：`<start> - <end> (<depth>)`；`start` 为 0-based k-mer 位置，
  `end` = start + run 长度（k-mer 数）+ k - 1，即 **1-based 基因组闭合端**。
* run 边界 = profile 值**恒定**的连续段（值 > 0 才输出；0 是分隔符）。
* **尾 run 不闭合**（最后一行只有裸 start，无 end/depth）。这是该 Profex
  版本的输出缺陷：pgr 现有处理（`run_profex_per_chr`）只能猜——e-kmer 用
  染色体长度闭合，s-kmer 因 depth 未知而丢弃，导致染色体末端的重复区间
  被漏掉。
* pgr 写 `.rg`：`chr:start+1-end`（start 转 1-based，end 原样）。

**原生实现不保留该 quirk**：profile 向量在内存中是完整的，尾 run 的 depth
已知，直接按普通 run 输出即可，顺带修复旧管线的漏区间问题。

### 2.3 源码核对记录（版本关系 + 差异清单）

**版本关系**（2026-08 核对）：

* 仓库内 `FASTK-master/` 是 2025-09-13 下载的快照，与 FASTK-1.2 完全一致
  （README `Current: April 18, 2021`；Profex.c/count.c/libfastk.c 等
  md5 与行数逐一相同）。**它不是上游当前 master**——上游提交历史显示
  2021 之后仍有 2022-12 / 2023-06 / 2024-10-23 的提交。
* 本机 `/home/wangq/.cbp/bin/` 的 FastK / Profex 由 CBP 安装（2025-03-10）：
  配方 `~/Scripts/cbp/packages/fastk.json` 指向上游 commit
  `ddea6cf254f378db51d22c6eb21af775fa9e1f77`（提交标题 "Logex space
  consumption issue fixed"，GitHub 已查证该提交同时改了 Profex.c）。
  CBP 构建脚本 `scripts/fastk.sh` 只改 Makefile 链接（libdeflate / libhts /
  -lz）并用 zig cc 编译，**没有源码补丁**；CBP 源码包解压目录
  `/tmp/cbp_fastk_check.*/fastk/` 与 `/tmp/fastk_src/fastk/` 的
  Profex.c / FastK.c / count.c / libfastk.c md5 全部一致。
* 因此下面"本地补丁"的说法修正为：**ddea6cf 相对 FASTK-1.2 的上游改动**，
  不是本地修改。

| 语义 | 源码依据 | 结论 |
| :--- | :--- | :--- |
| canonical = 2-bit 字典序较小者 | `count.c::kmer_list_thread`（`kb<hb` 取正向）+ `Comp` 表（逐 2-bit 互补，索引倒序实现序列反转） | 与 §2.1 一致 |
| count 上限 32767 | `count.c`（`ct>=0x8000` 时 cap `0x7fff`） | 对应 §3.3 u16 cap |
| N 拆分、profile 对应位置为 0 | `split.c` 按 gap 拆分；实测 N 段在 `-z` 输出中无内容 | 与 §2.1 一致 |
| profile 值 = 跨 read 的数据集级频次 | 实测：read1 重复段 count=3（read1×2 + read3×1） | 与 §2.1 一致 |
| 相对 profile 值 = 表内 count，缺省 0 | `count.c::cmer_merge_thread`（命中取表 count，否则 0） | 与 §2.1 一致 |
| `-t` 无参数 = cutoff 1（全量入表） | `FastK.c`（`flags['t']` → `DO_TABLE=1`） | e-kmer 建表全量 |
| `-p:<table>` 要求 k 一致 | `FastK.c`（`PRO_TABLE->kmer != KMER` 报错） | 对应 §4.3 header 校验 |

**ddea6cf 相对 FASTK-1.2 的差异**（`/tmp/fastk_src/fastk/` vs 仓库内
`FASTK-master/` 逐文件 diff）：

* `Profex.c`：`-z` 语义两版相同（非 ASCII 分支输出 run 形式）；差异是
  run 闭合的 end 由 1.2 的 `i-1`（0-based k-mer 位置）改为
  `i + kmer - 1`（1-based 基因组闭合端）、**不再闭合尾 run**，并移除
  `-A`（ASCII 输出）。
* `libfastk.c`：删除 histogram 读写；`Fetch_Profile` 的 `plen`/返回类型
  int64 → int（profile 长度上限 2^31，pgr 场景无影响）。
* `count.c`：`RUN_BYTES` → `PLEN_BYTES` 修正与 `Runer_Reload` 简化
  （profile 编码小修，不涉及 canonical/count 语义）。
* `merge.c` / `table.c`：仅文件权限位；`FastK.c`：`-P` 默认 /tmp；
  `split.c`：编译清理。

**计数核心语义（canonical、count 上限 32767、N 处理、相对 profile = 表内
count）在两版中一致**，上表核对结论对 1.2 与 ddea6cf 均成立。Profex 输出
契约（§2.2）以本机 ddea6cf 版实测为准——`run_profex_per_chr` 解析的正是
该行为；原生实现直接生成 §2.2 语义，不受上游 `Profex.c` 影响。

## 3. 模块设计：`src/libs/kmer/`

新目录，注册到 `src/libs/mod.rs`。三层职责：

```text
src/libs/kmer/
  mod.rs     公共类型 KmerTable、模块文档、re-export
  count.rs   计数表构建（canonical key 收集 → radix sort → 分组计数）+ 持久化
  profile.rs profile 向量生成（自计数 / 相对表）+ RLE 编码（备用）
  extract.rs profile → run 提取 → 写 .rg（替代 run_profex_per_chr 的核心）
```

### 3.0 格式决策：不扩充 PGI

`KmerTable` 是**独立格式**，不复用、不扩充 `.pgi`：

* `.pgi` 的语义绑定 syncmer 采样（header 含 `smer/window`，entries 只含
  syncmer k-mer，positions 存出现位置），`align pgi` / `dist seq` / `sd`
  都依赖这套语义；塞入"全量计数模式"会让一个格式承载两种语义。
* 需求只需 `canonical k-mer → count` 的查表，不需要位置信息；全量 k-mer
  存位置会爆炸且无用。
* 实现层面复用：2-bit 滚动编码（`kx/kxr` + `nt::rc_key`）、
  `ds/radix_sort`、bincode 持久化模式（magic/version header，风格对齐
  `PgiIndex::write/read`），但文件格式互相独立。
* **持久化紧凑编码**：内存用 u128 key（与 pgi 一致），落盘时把 key 打包成
  `ceil(2k/8)` 字节 + count，避免裸 bincode 序列化 u128 的 3 倍浪费；
  bincode 只当容器，不直接序列化 `Vec<u128>`（详见 §4.3）。

### 3.1 核心类型

```rust
/// Sorted canonical k-mer table with parallel counts.
pub struct KmerTable {
    pub k: usize,
    pub keys: Vec<u128>,   // 升序、去重、canonical
    pub counts: Vec<u32>,  // 与 keys 平行
}
```

### 3.2 count.rs：构建计数表

```rust
pub fn build_table(seqs: &[Vec<u8>], k: usize) -> anyhow::Result<KmerTable>;
pub fn save(table: &KmerTable, path: &Path) -> anyhow::Result<()>; // 紧凑编码
pub fn load(path: &Path, k: usize) -> anyhow::Result<KmerTable>;   // header 校验
```

1. 逐序列滚动 2-bit key（复用 `pgi/build.rs` 的 kx/kxr 滚动与 `nt::rc_key`，
   但只保留 canonical = min(正, 反)，N 清零重滚，含 N 的窗口无 key）；
2. rayon 按序列并行收集 `Vec<u128>`；
3. `ds::radix_sort::radix_sort_u128_par` 全局排序（与 pgi 一致）；
4. 一趟分组得 `(keys, counts)`。

直接路径仍是默认：内存中 `KmerTable` 用打包字节 key（§12），
~5 B/唯一 k-mer（k=17）：5 Mb 细菌 ~5 M key ≈ 25 MB。FastK 式
super-mer 两段计数实现在 `supermer.rs`（§3.6），以 `pgr kmer table
--supermer` **显式选项**接入（不做自动判断，选择权交给用户；默认仍是
直接路径）。超大输入的后备（分块计数）不在本期范围，接口上
`build_table` 与 `profile` 分离即可，将来可换实现。持久化用紧凑编码，
与内存表示无关（§4.3）。

### 3.3 profile.rs：生成 profile

```rust
pub fn self_profiles(seqs: &[Vec<u8>], k: usize, table: &KmerTable) -> Vec<Vec<u16>>;
pub fn relative_profiles(seqs: &[Vec<u8>], k: usize, table: &KmerTable) -> Vec<Vec<u16>>;
```

* 逐 k-mer 在 `keys` 上 `partition_point` 二分查表（无额外内存）。
* self：查得 count（≥1）；relative：查得表内 count，缺省 0。
* 含 N 的 k-mer 位置没有 key，profile 值为 0（run 分隔符，见 §2.1）。
* 用 `u16` 对齐 FastK 的 32767 上限（真实场景不触发，超限 cap）。

### 3.4 extract.rs：run 提取（Profex 等价）

```rust
/// 把每条染色体的 profile 写成 prof.<sn>.rg 文件（1-based 闭合区间）。
pub fn write_rg(
    profiles: &[Vec<u16>],
    chrs: &[String],
    k: usize,
    min_depth: Option<u16>,
    rg_files: &mut Vec<String>,
) -> anyhow::Result<()>;
```

逻辑 = Profex `-z` + `run_profex_per_chr` 语义：

* 扫描 profile，切分**恒定值 > 0** 的 run；
* 每 run（含尾 run）：`start = 0-based k-mer 起点 + 1`，
  `end = start0 + len + k - 1`；
* `min_depth` 过滤（s-kmer = 2）。

profile 完整时尾 run 的 end 自然正确（最后一个 k-mer 覆盖到序列末尾），
不需要染色体长度，故 `write_rg` 无 `lens` 参数。

只写 `.rg`，不复刻 Profex 文本；染色体名映射（点号 → `cN`）仍在管线层做。

### 3.5 FastK 查表结构对照（2026-08-09，profiling 驱动）

§3.3 的"全局排序数组 + `partition_point`"是迁移时的简化取舍（语义一致、
零额外内存），但 **FastK 原版查表不是全局二分**，而是前缀索引 + 分桶：

* `.ktab` 参考表是"前缀压缩 + 索引"结构（`libfastk.c` `_Kmer_Stream`）：
  `ibyte` 字节前缀 → `index[1 << (8*ibyte)]` 偏移表 + `inverse_index`
  （反向前缀索引）。
* `Split_Table`（`split.c`，相对 profile 时）把参考表按 1-byte 前缀
  （`Pro_File::index[257]`）分发到 NPARTS×NTHREADS 个块表；序列侧 k-mer
  同前缀路由到块，排序合并时**块内**取 count——即 O(1) 前缀定位 +
  桶内小范围查找（几十条，可入 cache），且多线程独立。

pgr 迁移丢失了这层结构。profiling 实测（[[../benchmarks/bench-profile-hotspots.md]]）：
`rept s-kmer` 中 `table_profiles` 占 78.5%，每窗口 ~23 次 u128 比较全是对
~73 MB 全局表的 DRAM 随机访问（cache-miss 41%，每窗口 ~28 次）；
mg1655 为单 contig，"按序列 par_iter" 无效。

**已实现（2026-08-09，排序合并替代逐窗口查表）**：前缀索引方案先实测证伪
（隔离基准：全局二分 1.125 s vs 前缀桶 1.195 s，见
[[../benchmarks/bench-profile-hotspots.md]]——73 MB 表随机访问延迟主导，
比较次数 23→7 无收益）。最终采用 FastK 的**排序合并**路线：收集全部窗口
key（并行）→ `radix_sort_u128_par` 排序 → 与 `table.keys` 线性归并一次
写回。基准：self_profiles 1.43 s → 250 ms（5.2×）、relative_profiles
1.41 s → 270 ms（5.4×），`rept s-kmer` 整命令 1.67 s → 0.50 s（3.4×）。
接口 `self_profiles`/`relative_profiles` 不变；语义由
`sort_merge_matches_binary_search` 对照测试保证。

### 3.6 supermer.rs：FastK 式两段计数（2026-08-14，阶段 B 原型）

> 分工（`~/Scripts/anchr/notes/design/unitig-bucket.md` §3.1）：**pgr 实现计数
> 层**（本模块），anchr 待接口更新后接入（分桶表 k-way 合并 + `TadpoleTable`
> 构建入口）。本条目记录 pgr 侧实现与"先原型验证收益"的实测结论。

**算法**（对照 `FASTK-master/split.c`/`count.c`，见 `notes/references/fastk.md`
§3）：

1. **m-mer 值**：每条序列按 N-free run 滚动 canonical m-mer（fwd/rc 两个
   u32 滚动，取小者），m = min(12, k-1)（FastK 自适应 `PAD_LEN` 典型 10–13，
   固定值简化实现）。
2. **run 划分**：FastK 式——run 持续到出现严格更小的 m-mer 或距定义 m-mer
   位置 ≥ MAX_SUPER = k-m+1 窗口（force cut）。pgr 用**无重叠**版本（切点
   窗口归新 run；FastK 把边界窗口计入前后两个 span，靠第二段加权合并兜底，
   计数语义等价但记录有 ~1 窗口/切的冗余）。
3. **记录**：span 按定义 m-mer 的 canonical 方向打包（`flip` 使正反链同一
   区域产出字节相同的 span，可在第一段合并），固定尺寸 = `ceil((2k-m+1)/4)`
   字节 + u16 窗口数。
4. **第一段**：整条记录做 `radix_sort_bytes_par` 排序，折叠相同 span
   （多重度 ct）。
5. **第二段**：每个唯一 span 展开 canonical k-mer（沿用 `canonical_keys`
   的滚动 + 半长比较），每条以权重 ct 入数组（u32，超 u16 不封顶），
   `radix_sort_bytes_par`（key=packed k-mer，payload=weight）后按 key 累加。
   输出与 `count::count_keys` **逐字节一致**（`matches_direct_*` 系列测试）。

**正确性测试**：随机数据（k=5/8/17/31/64/100）、含 N、大小写、重复读、
70,000 重（>u16 权重）、k 全扫 3..=40、m 边界（k=3/m=2、m=k-1）、正反链合并、
空/短输入、非法参数。

**基准结论**（`benches/supermer_benchmark.rs` 内部对比 +
`notes/benchmarks/bench-supermer-vs-fastk.md` 端到端 FastK 对照，
mg1655 + 合成 reads，release，多次均值）：

| 数据 | k | direct | supermer | 结论 |
| :--- | ---: | ---: | ---: | :--- |
| genome（单拷贝，无冗余） | 17 | 148 ms | 404 ms | **慢 2.7×** |
| genome | 31 | 191 ms | 339 ms | 慢 1.8× |
| genome | 100 | 288 ms | 367 ms | 慢 1.3× |
| 150 bp reads ×20 覆盖（唯一起点，无重复读） | 17 | 303 ms | 211 ms | 快 1.4× |
| 同上 | 31 | 321 ms | 209 ms | 快 1.5× |
| 同上 | 100 | 254 ms | 483 ms | **慢 1.9×** |
| 同 reads ×10 重复（极端高冗余） | 17/31/100 | — | — | 快 2.1–3.8× |

阶段拆分（genome k=17）：gen ~154 ms（≈2–3× 直接生成）、stage1 排序 ~16 ms、
**expand 单线程 ~145 ms（最大头）**、stage2 排序 ~48 ms。

**与 FastK 端到端对照**（99.5 M bp / 663k reads，32 核）：
FastK k=31 **0.74 s / 411 MB**、k=100 **0.95 s / 907 MB**；pgr 直接路径
端到端 1.74 s / 1.46 GB（k=31）、1.60 s / 1.95 GB（k=100）；pgr supermer
lib 0.87 s（k=31）/ 1.85 s（k=100）。同一输入下 pgr supermer 与 FastK 的
span 实例数与 stage-2 加权 k-mer 数**几乎逐一对上**（k=31 均 ~4.5× 折叠、
k=100 均 ~1.1×）——算法已同构。FastK k=100 依然快是因为窗口总量少
（150 bp 读只有 51 窗口）+ **工程效率**（C 位打包、907 MB vs pgr 1.95 GB、
32 线程），而非 super-mer 折叠（FastK 自己告警
`Too much of the data is in reads on the order of the k-mer size`）。

**结论**：super-mer 的收益**只来自 span 级冗余**（同一基因组区域被多条读
覆盖时，内部 span 字节相同可在第一段折叠）；几何条件是 **span 长度 << 读长**
（span ≈ 2k-m，k ≤ 读长/3 左右时折叠明显）。**k 接近读长时 span ≈ 整条读，
无折叠，两段式变成纯开销**——这正是 `unitig-bucket.md` §3.1 想解决的长 k
场景（k=100、150 bp 读），原型验证未达预期（慢 1.9×），**FastK 本体在该
场景同样无折叠收益**（savings 1.1×），其领先来自工程效率。无冗余长序列
（单拷贝基因组）同样无收益。可能的出路（未实现）：① 按
自适应阈值切换路径；② bcalm 式"窗口 minimizer 相同才连段"的短 span
变体（k=100 时 span ~100 bp，跨读折叠仍有戏，但 FastK 的 k=100 折叠实测
只有 1.1×，预期有限）；③ 只用于中低 k（≤31）的读数据（k=31 时 pgr
supermer 0.87 s 已接近 FastK 0.74 s）。**已定（2026-08-14，用户决定）：
不做自动判断，接入为 `pgr kmer table --supermer` 显式选项（默认直接
路径，输出逐字节一致，CLI 测试锁定）。**

## 4. 集成改动

### 4.1 `src/libs/pl/repeat.rs`

* 数据流：`pgi::build::read_fasta` 一次性读入 `(names, seqs)`；
  `has_sequences` 预检改为检查内存序列（友好报错语义不变）；
  `chr.sizes` / `pgr fa size` 调用删除（`write_rg` 不需要染色体长度，
  名字直接从内存取）。
* `run_repeat_pipeline`：`FastK -p / -t / -p:<prefix>` 三个 `run_cmd!` 分支
  替换为：
  * s-kmer：`kmer::count::build_table(seqs)` →
    `kmer::profile::self_profiles` → `kmer::extract::write_rg`；
  * e-kmer：`build_table(库)`（缓存命中则 `load`）→
    `kmer::profile::relative_profiles` → `write_rg`。
* `RepeatOpts` 删除 `re_prof` 字段，其余字段不变；命令层
  （s_kmer.rs / e_kmer.rs）同步删除 regex 构造与传参。
* 删除 `run_profex_per_chr`。
* `-P` 排序目录逻辑删除；tempdir（`PipelineCtx::enter`）保留（`.rg` 中间文件）。
* 日志去掉 FastK 字样，沿用 `==>` 风格：`==> Counting k-mers`、
  `==> Building k-mer table`、`==> Extracting repeats`。

### 4.2 命令层与文档

* `src/cmd_pgr/rept/s_kmer.rs` / `e_kmer.rs`：CLI 参数不变；`after_help` 删除
  "External dependencies: FastK / Profex"。
* `README.md`、`docs/rept.md`、`docs/usage_examples.md`：依赖说明改为无外部
  依赖；`docs/rept.md` 中 FastK 并行 SIGSEGV、`-P` 目录等注意事项删除/改写。
* `CHANGELOG.md` 记录迁移与旧 `.ktab` 缓存作废。

### 4.3 `--keep-index` 缓存

* 新格式：`<库>.pkt` 单文件（`lib.fa` → `lib.pkt`，`lib.fa.gz` →
  `lib.fa.pkt`），**紧凑编码**：
  header（magic/version/k/条目数/key 字节数）+ 每条目
  `packed key（ceil(2k/8) 字节，复用 pgi 的 pack_kmer）+ u32 count`；
  k=17 时约 9 B/条目，5 Mb 库 ~45 MB（对比裸 bincode u128 的 ~100 MB）。
  bincode 只当容器；不再需要 FastK 的隐藏分片，也去掉 `.complete` 标记——
  完整性由 header 校验兜底（损坏即重建），写入用原子 rename（临时文件 +
  rename）。`cache_is_fresh` 保留 mtime 检查，判断对象从 `.ktab/.complete`
  变为单个 `.pkt` 文件。
  命名遵循项目 sidecar 惯例（替换扩展名，同 `.pgi`，实现参考
  `align/pgi.rs::sibling_pgi_path` 的 `.gz` 分支）：文件名不带 k 或
  场景限定，k 存在 header 里，读取时校验——k 与命令行不一致或 mtime 旧
  则重建（对齐 `align pgi` 的 sibling index 检查；缓存是纯加速，重建比
  报错友好）。`KmerTable` 是通用格式，当前只有 e-kmer 用它。
* **旧 FastK `.ktab` 缓存不兼容**：升级后首次运行自动重建（README 注明
  一次重建成本）。

## 5. 验证计划

1. **单元测试（kmer 模块）**：
   * canonical 编码与 `nt::rc_key` 一致性、N/gap 拆分；
   * 小序列手工核对计数；重复序列频次；
   * relative profile：表内 count / 0；
   * run 提取边界：恒定值切分、min_depth 过滤、尾 run 正常输出
     （旧管线漏区间修复项）；
   * `KmerTable` save/load roundtrip、截断文件判脏、header k 与命令行
     不一致判 stale（沿用 `cache_is_fresh` 测试）。
2. **集成测试（tests/cli_rept.rs）**：现有 e2e 用例去掉
   `FastK/Profex in $PATH` 跳过条件；新增"无外部工具可跑通"断言。
3. **合理性复核（一次性脚本，不进 CI）**：MG1655 上新管线结果与
   FastK+Profex 粗略对照（旧管线实验性，仅作参考），人工复核重复区间
   覆盖合理、染色体末端尾 run 不再漏。
4. **边界输入**：全 N 序列、单染色体、超短序列（< k）、空库（沿用预检报错）。
5. **基准**：`benches/` 下计数 + profile 生成（MG1655 级，参照现有
   `pgi build` bench）。

## 6. 工作量估算（净增行数）

| 部分 | 估算 |
| :--- | :--- |
| `libs/kmer/count.rs`（含持久化） | 350–450 |
| `libs/kmer/profile.rs` | 250–350 |
| `libs/kmer/extract.rs` | 150–250 |
| `libs/kmer/mod.rs` + lib.rs 注册 | ~50 |
| `pl/repeat.rs` 集成（净） | 100–150 |
| 命令层 / 文档（净） | -30 |
| 测试 | 400–600 |
| **合计** | **1,300–1,800** |

实现顺序：count → profile → extract → 管线集成 → 合理性复核 → 文档清理。

## 7. 依赖与风格

* 复用：`fmt/fa`（读 FASTA/gz）、`pgi::build::read_fasta`（返回
  `(name, seq)` 列表）、`nt::rc_key`、`ds/radix_sort`、rayon、
  `bincode + serde`（已有依赖，不新增）、`pgi` 的
  `pack_kmer`/`unpack_kmer`（pub(crate)，同 crate 直接可用）；
  持久化写读模式参考 `pgi/mod.rs` 的 magic/version header 实现。
* 新代码全部在 `libs/`；`cmd_pgr` 保持薄壳；公共 API 写一行英文 doc comment。
* 不引入新依赖。

## 8. 风险与决策点

1. **尾 run quirk**：已决定不保留，原生输出完整 run（修复旧管线漏区间）。
   与旧管线的差异仅限染色体末端区间。
2. **内存**：单次全量计数无磁盘分桶。若未来目标变成 Gb 级基因组，在
   `build_table` 内部加分块，接口不变。
3. **旧缓存失效**：`--keep-index` 的 FastK 格式缓存作废，重建一次。
4. **profile u16 上限**：对齐 FastK 32767，超限 cap（不 panic）；
   `KmerTable` 内部 count 用 u32，不受此限。
5. **e-kmer profile 值语义**：e-kmer 不读 depth 值本身，但 profile 值参与
   恒定值 run 的切分（§2.2），因此表内 count 的具体值影响 run 边界——原生
   实现必须生成真实表内 count，不能退化成 0/1 存在性标记（与 FastK/Profex
   语义一致，也解释了 §2.1 的相对 profile 定义）。

## 9. FASTK 功能对照与缺口（2026-08 补充）

> 配套算法机制分析（Super-mer/Minimizer/分桶）见 [fastk.md](../references/fastk.md)。
> 下表是"FASTK 能提供的能力" vs "pgr `libs/kmer` 现状"的功能对照。

| FASTK 能力 | pgr 现状 | 说明 |
|---|---|---|
| `-p` profile | ✅ `profile.rs self_profiles` | rept s-kmer 用 |
| `-t` k-mer 表 | ✅ `count.rs KmerTable`（.pkt 缓存） | 内存版；不做 .ktab 磁盘分桶 |
| `-p:<table>` 相对 profile | ✅ `relative_profiles` | rept e-kmer 用 |
| Profex `-z` run 提取 | ✅ `extract.rs write_rg` | 已实现（修复尾 run quirk） |
| **直方图 `.hist` + Histex（含 `-G` 基因组大小格式）** | ✅ `hist.rs`（.hist 兼容写） | **已实现（2026-08-09）**：`from_table` 聚合 + FASTK 二进制布局写（固定 low=1/high=32767，含 ilowcnt/max_inst）；实测 Histex 读 pgr 输出与 FastK 自产逐行一致；峰值/GenomeScope 模型拟合仍是第二步（R 侧，未做） |
| **profile 落盘（`.prof` 等价）** | ✅ `.pkp`（自有格式） | **已实现（2026-08-09）**：`save_profiles`/`load_profiles`，header + raw u16；**不做 FASTK `.prof` 兼容**（多文件分片 + RLE，代价中等且 pgr 无外部消费者） |
| `-c` homopolymer 压缩 | ❌ 无 | PacBio/HiFi 专用（homopolymer 错误率高），低成本 |
| Logex（表逻辑运算 + 计数阈值过滤） | ❌ 无 | 两个 k-mer 库的 AND/OR/NOT 集合运算 |
| KmerMap（k-mer → .bed 区域） | ❌ 无 | k-mer 在目标序列上的覆盖区域 |
| Tabex / Symmex | ❌ 无 | 表查看/导出、canonical → 对称表 |
| Fastmerge/Fastcat/Fastrm 等 | ❌ 不做 | TB 级磁盘分桶/分布式设计，pgr 场景用不上（§8.2 已声明） |

**结论**：`libs/kmer` 已覆盖 FASTK 的 **rept 用途**（profile/相对 profile/
run）与 **直方图**（.hist 兼容）。剩余缺口优先级：
1. 峰值/基因组大小估计（GenomeScope 模型拟合，R 侧，未立项）；
2. homopolymer 压缩（PacBio/HiFi 场景）；
3. Logex / KmerMap（按需求再议）。

### 9.1 anchr 2_fastk 的实际用法（2026-08 补充）

anchr `templates/2_fastk.tera.sh`（对 R/S/T 三个样本，各自单端或双端 fq.gz）：

```bash
FastK -v -T<threads> -t1 -k<21|51|81> <S>1.fq.gz [<S>2.fq.gz] -NTable-<k>
Histex -G Table-<k> | Rscript ../../0_script/genescopefk.R -k <k> -p 1 ...
KatGC -T<threads> -x1.9 -s Table-<k> <P>-Merqury-KatGC-<k>
Fastrm Table-<k>
```

- **只用 FASTK 的 k-mer 表（`-t1`，cutoff=1 全量）+ 默认直方图（`.hist`）+
  `Histex -G`**（GenomeScope 2.0 ASCII 格式）；**没用 profile/相对 profile**
  （那是 rept 的活，pgr 已实现）。
- 下游 `genescopefk.R`（R 脚本，anchr 模板 `templates/genescopefk.R.gz`）做
  GenomeScope 模型拟合，输出 `summary.txt`/`model.txt`（2_fastk 解析
  `model.txt` 的 `kmercov` 字段并汇总成 statFastK.tsv/md）。KatGC 是
  Merqury 家族的外部工具（k-mer 覆盖度 vs GC），不在 FASTK 范围内。
- **结论**：直方图 + `Histex -G` 格式正是 §9 标注的最大缺口，且是
  anchr 2_fastk 的直接消费者。
- **范围提示**：直方图聚合（`.hist` 文本输出）成本低；GenomeScope 的
  基因组大小/杂合度估计是**模型拟合**（由 R 脚本承担），pgr 若完整替代
  需原生实现拟合（k-mer 谱 → 泊松/负二项模型 → kmercov/基因组大小），属
  第二步、成本更高。可先做直方图 + ASCII 输出（对齐 `Histex -A/-G` 格式），
  拟合后续再议。

## 10. `pgr kmer` 命令组与三种格式（2026-08-09 定稿）

### 10.1 命令归属

新增顶级命令组 `pgr kmer`（table/profile/hist/gc/qhist/qcheck/gsize），**`rept s-kmer` /
`e-kmer` 保留不动**：它们是完整"重复提取管线"（计数只是第一步），语义归
rept；`pgr kmer` 是通用 k-mer 分析（表/谱/直方图），消费者是
GenomeScope/MerquryFK 场景。`libs/kmer` 共享，命令层互不依赖。

### 10.2 三种格式

| 格式 | 内容 | 策略 | 实现 |
|---|---|---|---|
| `.pkt` | canonical k-mer 计数表（原 `.pgrk` 改名，magic `PKTT`） | pgr 自有单文件 | `count.rs`（紧凑编码，同 §4.3） |
| `.pkp` | 逐序列 profile（magic `PKPP`） | pgr 自有单文件 | `profile.rs save/load_profiles` |
| `.hist` | 频次直方图 | **FASTK 字节兼容** | `hist.rs`（新增） |
| `.kgc` | GC×覆盖度矩阵 | **KatGC 兼容**（实测逐行一致） | `gc.rs`（新增） |
| qhist 输出 | 质量偏置直方图 | **quorum `histo_mer_database` 格式兼容** | `quality.rs`（新增） |

`.hist` 布局（与 FastK `count.c` 写侧一致，实测 Histex 读 pgr 输出 =
FastK 自产，diff 为空）：

```text
int32 k | int32 low=1 | int32 high=32767 | int64 ilowcnt(=bin1 数)
int64 max_inst | int64 hist[1..=32767]     # 28B 头 + 32767×8B = 262164B 固定
```

`.pkp` 布局：`magic PKPP | u32 version | u32 k | u64 n_seqs`，每条
`u64 length | u16[length]`（raw，未压缩；RLE 属内部编码，将来优化不动
header）。

### 10.3 输入形态

* 序列输入：FASTA/FASTQ（复用 `fmt/seq` FAFQ reader），支持 gz/stdin。
* `table`：多输入合并计数（FastK `-t1` 语义，含 singleton），`-o .pkt`。
* `profile`：序列必填；无 `-t` = self（内部建表不落盘），有 `-t` =
  relative（表内 count，缺省 0）。`k` 解析规则：有表用表 k（header 读取，
  `count::k_of`），命令行给 k 则校验一致；无表时 `-k` 必填。
* `hist`：序列直算或 `-t` 表二选一（`required_unless_present`）。
* `gc`：同 hist 输入形态；`-X`（绝对 x 上限，兼作 count cap）/`-x`（倍数，
  默认 2.1）对齐 KatGC 峰值语义；`--tex` 渲染 heat 图（复用
  `plot hh` 的 heatmap.tex + 自适应轴，pgr 无 R 依赖）。
* `qhist`：FASTQ 必需（质量阈值判定）；阈值默认 = 自动检测 Phred 偏移
  （复用 `fq trim` 的 BBDuk 检测）+ 5（quorum 默认），`-q` 可显式覆盖。
* `qcheck`：FASTQ 必需；建质量表（同 qhist 阈值/bits）+ 逐 read 判定
  （anchor + 双向 extend），`-o` 保留 / `--discard-file` 丢弃。
* `gsize`：同 hist 输入形态；输出 peak_coverage/total_distinct/total_kmers/
  genome_size（= total_kmers / peak_cov，简单单倍体估计）；`--model` 跑
  GenomeScope 完整拟合（见 §10.7）。

### 10.7 GenomeScope 完整迁移（2026-08-09）

`gsize --model` = **genescopefk.R（GenomeScope 2.0）原生移植**
（`libs/kmer/genomescope.rs`）。范围：

* **模型**：p=1（unique + repeat 两负二项分量）、p=2（AA/AB/BB 四类混合，
  `predict2_1`）；公式逐行对齐 R（alpha 系数 + `dnbinom(x, size=kmercov*i/
  bias, mu=kmercov*i)`），`predict` 返回 mixture，公式层乘 `x^1 * length`
  （R 的 `y_transform ~ x*length*predict`）；`dnbinom` 对 size=Inf 退化为
  泊松（R 语义，bias→0 时 LM 不会产生 NaN）。
* **拟合**：**minpack `lmdif` 完整移植**（genescopefk.R 的 nlsLM 无解析
  雅可比 → lmdif）：`fdjac2`（前向差分，eps=sqrt(machine eps)）+
  `qrfac`（Householder QR 带列置换）`qrsolv`（Givens 消去阻尼对角）+
  `lmpar`（trust-region 阻尼参数二分）+ `lmdif` 主循环（外循环雅可比/QR/
  梯度判据 + 内循环 trust-region 步长 + ftol/ptol/gtol 收敛判据 +
  info 1-8）；参数边界用 R 的投影法（fcn 调用前 clamp）。拟合后从最终
  R 算 hessian（P^T R^T R P）求逆得真实 SE（对齐 R summary）。
* **输出**：`summary.txt`（property/min/max 表）+ `model.txt`（参数表含
  Estimate + Std. Error，`kmercov` 行可被 anchr `2_fastk` 的
  `grep '^kmercov' | cut -f 2` 解析）。
* **不做**：p>2 多拓扑、错误分量、端部修正、R 绘图（明确记录）。
* **验证（本机 R 4.4.2 + minpack.lm 端到端对照，2026-08-10）**：真实
  60× 1 kb reads 同一直方图喂 R `genescopefk.R` 与 pgr `gsize --model`：
  kmercov 55.7 vs 55.73、bias 0 vs 0、d 0 vs 0、length 1018 vs 988、
  kmercov SE 0.745 vs 0.787（全部高度一致）；Model Fit 62.1% vs 64.4%
  （2% 差，score_model 的 first_zero 边缘细节，记录为已知偏差）。
  无噪声合成谱精确还原参数（单测）。
* **排错记录**：`predict` 误含 x*length 因子（R 的 predict 只返回
  mixture）→ LM 无法收敛；稀疏/稠密直方图语义（R 读文件只含非零行，
  稠密数组把 count=0 的 x 也算残差）→ m/score 错误；`qrsolv` 的 sdiag
  只清 j..n（清全部会破坏前几列 S 对角）→ lmpar 零步长；`dnbinom`
  size=Inf 需退化泊松 → bias 无法到 0 边界；est_length 初值稀疏后已是
  正确量级（x 因子，非 x³），用 est 附近多初值。

### 10.5 QUORUM 质量偏置计数（2026-08-09 补充）

quorum `hash_with_quality` 的 nval 编码 `(count<<1)|quality` 更新语义经
推导**与出现顺序无关**：key 只要出现过一次高质量（窗口 k 碱基全 ≥ 阈值），
最终 count = 高质量出现次数、quality = 1（低质量出现不参与计数）；从未
高质量则 count = 低质量出现次数。因此 pgr 用"收集 (key, is_high) → 排序 →
分组聚合"实现，无需保序哈希表。质量判定：非 ACGT 拆分两条 stretch，
碱基 < 阈值只断 high stretch，窗口 high 当且仅当 high_len ≥ k（对齐
`quality_mer_counter`）。端到端对照受本机无 Jellyfish 2.0 限制未做，
语义经源码逐行核对 + 单测（顺序无关性/N 拆分/阈值/bits 封顶），并新增
**独立逐事件仿真对照**（随机 reads k=4/8/17 × 20 组：按 quorum `add()`
的 nval 更新规则保序仿真 vs 顺序无关聚合，直方图完全一致）。

### 10.6 QUORUM read 判定器（2026-08-09 补充）

`pgr kmer qcheck` 只判定不修正：重建 quorum 的**错误信号**（quorum.md
§6.1 确认的 anchr 场景 = 检测有错 read 直接丢弃）。关键语义：

* **anchor**：`get_val` 对低质量 k-mer 返回 0（`v.second ? v.first : 0`），
  故只有高质量 k-mer（quality=1）能作 anchor；连续 `good` 个 count ≥
  anchor-count。
* **extend 判定**（双向，backward 用反转序列镜像为 forward）：
  `get_best_alternatives` = 把当前（最新）碱基替换为 4 种，取最高质量
  等级候选的 count；count==0 → truncation，count==1 且唯一候选 ≠ 当前
  碱基 → substitution，count>1 且当前碱基不满足保留条件（count ≤
  min-count 或 Poisson 检验不通过）→ 事件。任何事件 → read 丢弃。
* **Poisson 检验**：`p = Σcounts × (apriori_error_rate/3)`，
  `poisson_term < threshold(1e-6)` 保留（Stirling 近似 i≥11）。
* 未移植：修正输出、err_log 窗口限制（判定器只看有无事件）、homo_trim、
  污染库。参数默认对齐 quorum（skip=0/good=1/anchor-count=1/min-count=1/
  cutoff=4/error-rate=0.01）。

### 10.4 兼容性结论（决策记录）

* `.hist` 兼容：单文件固定布局，实现 ~50 行，**做**——外部可比较
  （Histex/KatGC/GenomeScope 直接读）。
* `.prof` 兼容：FASTK `.prof` = stub + `.<root>.pidx.N` + `.<root>.prof.N`
  多文件分片 + RLE 编码（`merge.c`/`libfastk.c`），且 profile 按 read id
  索引；实现需 300–500 行 + 分片语义。**不做**（用户拍板：不喜欢分片，
  与 `.pgrk` 选自有格式同一理由）。
* `.ktab` 兼容：stub + 分片 + 前缀索引，维持不做（§3.0）。
* profile 落盘需求来源：MerquryFK 的 QV/completeness/trio 分型是
  self + 相对 profile 逐位置比较（`MerquryFK.c scan_asm`、`hap_plotter.c`），
  pgr 若做组装质量评估需要 `.pkp` 数据流；当前直方图不经过 profile。

### 10.8 真实数据验证（2026-08-10，Lambda SRR5042715）

用 BBTools 迁移引入的 Lambda 真实双端 reads（`tests/bbtools/Lambda/golden/
filter.fq.gz`，36384 reads，trim/filter 后）对 kmer 命令做真实数据验证，
并修复了 gsize 的 peak 估计：

* **table/hist（k=31）**：87941 unique 31-mers 与 BBTools `#unique_kmers`
  一致；khist-text/peaks 与 golden 逐字节一致（已有测试，M7）。
* **gsize 修复**：原 `estimate()` 取全局频次众数，真实数据下 count=1 的
  错误 kmer 主导 → peak=1、genome_size 267 万 bp（真实 48.5 kb 的 55 倍）。
  改为复用 `khist::call_peaks`（BBTools 完整移植）取体积最大峰 center →
  peak=56（= R.peaks `#main_peak`）、genome_size=47786 bp（误差 1.5%）。
  合成 30× 回归不变。新增真实数据回归测试
  `command_kmer_gsize_real_lambda_matches_bbtools_peak`（引用现有 golden，
  不新增测试料）。
* **gsize --model**：kmercov=55.3、genome_size=46873 bp（误差 3.4%）、
  bias=0.713、converged=true；与 R `genescopefk.R` 同直方图对照：
  summary 除路径行外逐字节一致、参数千分位级吻合
  （`design/genescopefk.md` §5.1）。
* **gc**：peak 53（31 GC × 111 count 矩阵；KatGC 峰值语义，与 CallPeaks
  的 56 因算法不同略有差异，合理）。
* **qhist**：auto 检测 Phred+33 → threshold=38 正确；输出 quorum
  `histo_mer_database` 格式；~80× 覆盖下所有 kmer 至少一次全高质量窗口
  → low=0（quorum 高质量主导语义，合理）。
* **qcheck**：kept 35221 / flagged 1163（3.2%，真实测序错误量级合理）。
* **profile**：self/relative 各 36384 profiles 正常。
* **测试覆盖改进**：原 `--model` 合成测试（60× 1kb reads）拟合落在
  病态 bias=0/d=0 泊松退化边界（length 高估 4.8×，断言区间 (500,10000)
  仅为容忍该病态），只证明"能跑通"。新增真实数据测试
  `command_kmer_gsize_model_real_lambda`（kmercov≈55/bias≈0.7/
  length≈46789，覆盖正常拟合区域）及 `command_kmer_gc_real_lambda`（peak
  ∈ 50..60）、`command_kmer_qhist_real_lambda`（threshold 38 + depth-1
  计数 = golden 38961）、`command_kmer_qcheck_real_lambda`（flagged
  2–5%）、`command_kmer_profile_real_lambda`（36384 profiles，self/
  relative）；全部引用现有 golden，不新增测试料，合成行为测试保留。

## 11. k 的范围与表示（2026-08-11 记录）

> 背景：排查"项目里哪些命令能用 k=81"时发现，k 的上限不是全局统一的，
> 而是由**三种 k-mer 表示**各自决定。此前文档只零散提及（`pgi.md`
> "at most 64"、`fq-assemble.md` §7 的 k>64 radix 搁置），没有一处
> 完整记录，故补本节。

### 11.1 三条表示路线

| 表示 | 上限 | 使用方 | 校验点 |
|---|---|---|---|
| u128 单字（2 bit/碱基） | **k ≤ 64** | `libs/kmer`（`KmerTable.keys: Vec<u128>`）、`libs/pgi`（`PgiEntry` u128）、`libs/map`（`MapIndex.keys`） | `count::build_table`、`pgi::build_from_seqs`、`map::build_index` 均 `ensure!(1..=64)`；`.pkt`/`.pgi` 读时另有 header 校验 |
| tadpole 多字 `Kmer`（`libs/asm/tadpole.rs`，`Vec<u64>` 共 2k 位，镜像 BBTools long array） | **无上限（k≥1）** | `asm contig`/`asm unitig`、`fq extend`/`fq ec-kmer`、`fq merge` extend2（硬编码 k=81） | CLI 仅 `range(1..)`；排序用 `cmp_bases` 比较排序，与 k 无关 |
| FastK 参考实现（`FASTK-master/`） | **无硬上限** | 参考实现本身 | `ARG_POSITIVE` 只要求 k>0 |

要点：

* **u128 = 128 bit ÷ 2 bit/碱基 = 64**，这是 map/pgi/kmer 三族上限
  的共同来源；CLI 层大多只收 `usize`，越界在 lib 层报友好错误。
* **tadpole 多字 Kmer 是项目里唯一为 k>64 设计的表示**，`asm contig`/
  `unitig` 用 `cmp_bases` 比较排序（与 k 无关），k=81 可用；
  `fq-assemble.md` §7 的"k>64 多 word radix 泛化搁置"只是排序优化的
  取舍，不是能力上限。
* **FastK 的 40 是默认值不是上限**：`FastK.c` `KMER=40` + `ARG_POSITIVE`
  （仅 >0）；k-mer 字节打包 `KMER_BYTES=(2k+7)>>3`，k 多大都能存，只是
  每条记录内存线性增长。**pgr 的 64 上限是移植选型**（u128 排序/radix/
  查表方便），不是 FastK 的约束。

### 11.2 实际场景对照（为什么重要）

* anchr `unitigs.tera.sh`：tadpole `k ∈ opt.kmer`（如 "31 81"）→
  anchr `asm contig`/`unitig` 已覆盖（多字 Kmer，k=81 OK）。
* anchr `2_fastk.tera.sh`：`FastK -t1 -k<21|51|81>` → **k=81 超出
  `pgr kmer table` 当前能力（u128 上限 64）**，是已知缺口（k=21/51
  没问题）。若将来替代 2_fastk 需要 k=81，得给 `libs/kmer` 扩表示
  （参考 FastK 字节打包，或 u128 双字），当前未做。
* anchr `asm map`/`pgr pgi`：默认 31/40，均 ≤64，anchors/GIX 场景无缺口。

## 12. 长 k-mer 落地：统一到 FastK 表示（2026-08-12 修改版，**已实施**）

> 承接 §11.2 的 k=81 缺口。目标：`pgr kmer table`（及同族
> profile/hist/gc/extract）支持 k=81，替代 anchr `2_fastk.tera.sh`
> 的 `FastK -k<21|51|81>`。**2026-08-12 用户裁定：项目只保留一套
> k-mer 实现，且以 FASTK-master 的表示为准**——不用 tadpole 的
> `Vec<u64>` 方案，也不新增独立 KmerKey；u128 键族与 tadpole Kmer
> 全部迁移到 FastK 风格字节键，最终项目里只有一个 k-mer 键类型。

### 12.1 关键事实（源码核实，2026-08-12）

* **FastK 编码**（`FASTK-master/`）：
  * 2-bit/碱基字节打包，`KMER_BYTES=(k+3)>>2`（`count.c`
    `kmer_list_thread`），**k 无硬上限**（`ARG_POSITIVE` 仅 >0）；
  * **字节序**：字节间 5'→3'（`Fetch_Kmer`/`Current_Kmer` 先输出
    最高字节），**字节内 5' 端碱基在高 2 位**（`setup_fmer_table`：
    `fmer[i]` 解码 = `dna[l3] dna[l2] dna[l1] dna[l0]`，l3 是 bit6-7；
    `kclip = {0xff,0xc0,0xf0,0xfc}` 裁剪低位、保留高位）——即"字节内
    大端"，**与 pgr `pgi::pack_kmer` 一致**（2026-08-12 实测修正：
    k=5 表条目 `acgta` = 字节 `1b 00`，Tabex 对照；上一版"字节内小端"
    记录有误，以实测为准）；
  * **canonical**：正链 f 与反链 r 逐字节比较取小（`count.c:495-530`，
    `Comp` 互补表），k%4 余位用 `KCLIP` 裁剪；**比较起点（最高/最低
    字节）不影响结果**——正/反链字节互为"反转+互补"（自逆映射），从
    两端做字典序选出的代表元相同，与 u128 `min(fwd,rc)` 大端比较一致；
  * 表：`Kmer_Table` = 前缀字节索引（`ibyte`/`inver`/`index`，
    `ixlen=1<<(8*ibyte)`）+ 后缀 + u16 count（32767 cap）；
  * 排序：Myers 风格 MSD/LSD radix（`MSDsort.c`/`LSDsort.c`）——
    与 pgr `ds/radix_sort.rs`（FastGA MSDsort 移植）同源。
* **pgr `.pkt` 磁盘格式已经是"打包字节 + count"**（`libs/kmer/count.rs`）：
  header（`PKTT`/v1/k/n_entries/key_bytes）+ 每条目 `ceil(2k/8)` 字节 key
  + u32 count。**key_bytes 由 k 决定，格式与 k 无关**——k=81 → 21 字节/条，
  结构天然支持长 k。**字节序已实测与 FastK 一致（字节内大端）**——迁移
  不需要翻转字节语义，`.pkt` 格式不变、**不 bump**，旧缓存直接兼容
  （上一版"方向相反、条目翻转、PKT_VERSION bump"记录有误，2026-08-12
  实测修正）。
* **内存表示是唯一瓶颈**：`KmerTable.keys: Vec<u128>`（`count.rs:35`
  `ensure!(1..=64)`）；同样 u128 的还有 `libs/pgi`（PgiEntry）、
  `libs/map`（MapIndex）、`libs/kmer/quality.rs`。
* **`kmer::n` 是唯一窗口发射函数**（`libs/kmer/mod.rs`），被
  count/profile/norm/map 共用——统一后全部发射 FastK 字节键。
* **tadpole Kmer**（`libs/asm/tadpole.rs:227`，`Vec<u64>` 小端窗口）：
  **不采用**（用户裁定，不信赖该实现），组装侧一并迁移到 FastK 字节键。
* **FASTGA（pgi 参考项目）的 k-mer = FastK 同一套**（2026-08-12 确认）：
  `FASTGA-main/` 自带 `libfastk.c`/`ONElib.c`/`gene_core.c`（与
  FASTK-master 同作者共享库）；`GIXmake.c` 建 GIX 索引用 2-bit 压缩
  （`KBYTES=KMER/4`，KMER 默认 40）+ `Comp` 互补表 + canonical 取小
  （`TMap[y] < TMap[z]`）；`FastKS.c`/`FastGA.c` 的 seed 匹配直接读
  `libfastk.h` 的 `Kmer_Stream`（merge 两个基因组的排序 k-mer 表 →
  adaptamer 变长种子）。pgr 的 `align pgi` 是 FastGA 移植/对照方向——
  **统一到 FastK 表示与 pgi 方向一致**；pgr `radix_sort_u128` 本就移植自
  FastGA `MSDsort.c`。

### 12.2 方案（唯一 k-mer 键 = FastK 字节编码 + 打包存储）

**唯一 k-mer 表示 = FastK 字节编码**（2-bit/碱基，字节间 5'→3'、
**字节内 5' 端碱基在高 2 位**，canonical = 正/反链逐字节取小
（`Comp` 互补表），k%4 余位 KCLIP 裁剪，`key_bytes=(k+3)>>2`），
放 `libs/kmer/key.rs`：

* 操作（一套实现）：`push_right`/`push_left`（字节滚动，组装侧延长）、
  `rc`、`canonical`、`base_at(i)`、`byte_at(i)`（radix 用）、
  `to_bytes`/`from_bytes`（`.pkt`/`.pgi` 直接承载）；
* `Ord` 按字节序（= FastK 表序 = radix 升序 = `partition_point` 二分
  一致）；canonical 语义与现有 u128 版一致（字典序取小），k≤64 结果
  逐 k-mer 相同。

**存储形态 = FastK 式连续打包**（关键：参考项目实测，**无 per-key
对象头**）：

* `KmerTable`：`keys: Vec<u8>`（连续打包区，每条 `key_bytes` 字节）+
  `counts: Vec<u32>`——与 FastK 表条目（`tbyte=kbyte+2`，
  `libfastk.c:419-420`）同一形态；
* 组装侧 HashMap（tadpole 迁移）：`Vec<u8>` 值键（kbyte 字节），同一套
  编码与操作；
* pgi/map 索引：同为排序表，同样打包存储（默认 k=31/40，打包后比现状
  u128 更小）。

**内存核算**（每条键，2026-08-12 对照参考项目实测）：

| k | 现状 u128（16+4 B） | ~~定长 [u8;32]（32+4 B）~~ | FastK 打包（key_bytes+4 B） | FastK 原版（tbyte=kbyte+2） |
|---|---|---|---|---|
| 21 | 20 | 36 | 10 | 8 |
| 31 | 20 | 36 | 12 | 10 |
| 51 | 20 | 36 | 17 | 15 |
| 81 | 20 | 36 | 25 | 23 |

> 定长对象键（32 B）会让内存翻倍——**用户顾虑，已否决**；打包存储下
> k≤64 比现状**更小**，k=81 仅 1.25 倍（能力扩展的必要成本）。FastK/
> FASTGA 均为变长打包（无 per-key 对象头），本方案与其一致。

**radix 泛化**：`radix_sort_u128` → 按打包条目排序（每条 `key_bytes`
字节，counts 并行交换），k≤64 排 8-16 字节，k=81 排 21 字节。

**窗口发射统一**：`kmer::n` 直接发射 Kmer 字节键（u128 快路径与 n_long
分叉都不要），count/profile/norm/map 全部随之适配。

**消费方全部迁移**（消灭 u128 键族 + tadpole Kmer）：

| 消费方 | 现状 | 迁移 |
|---|---|---|
| `libs/kmer`（KmerTable 系列） | u128 | 换 Kmer（本计划主体） |
| `libs/kmer/quality.rs` | u128 | 换 Kmer |
| `libs/pgi`（PgiEntry / build） | u128 | 换 Kmer；`.pgi` 字节布局已一致，格式不变（不 bump） |
| `libs/map`（MapIndex） | u128 | 换 Kmer |
| `libs/fq/norm.rs` | `kmer::n`（u128） | 跟随 n 自动适配 |
| `libs/asm/tadpole.rs` | `Vec<u64>` Kmer | **替换为 Kmer**（用户裁定） |
| `libs/nt.rs`（rc_key 等） | u128 辅助 | 随迁移收敛/删除 |

### 12.3 里程碑与验证

| 里程碑 | 内容 | 验证 |
|---|---|---|
| M1 | `key.rs`（Kmer 字节编码 + 打包存储，FastK 字节序）+ radix 泛化 | 现有测试绿；**与 FastK 二进制对照的 canonical 字节序列单测**（k=21/51/81，本地可编译 FASTK-master） |
| M2 | kmer 族换 Kmer（count/profile/hist/gc/extract/quality） | k≤64 与 u128 语义逐 k-mer 一致（golden 回归）；`.pkt` 格式不变、旧缓存直接兼容（字节已实测一致） |
| M3 | `kmer::n` 统一发射 Kmer + norm/map 适配 | map/norm 现有测试绿；k=81 表可建可读 |
| M4 | pgi 换 Kmer；`.pgi` 字节布局已一致，版本按实际变更定稿 | pgi 测试绿；旧 `.pgi` 兼容或明确重建策略 |
| M5 | tadpole 组装侧迁移 + 收尾 | asm contig/unitig k=31/81 与迁移前对照；FastK `-p` 端到端对照（-k 21/51/81）；radix/哈希基准；全项目 grep 无 u128 键/`Vec<u64>` Kmer 残留；更新本文档 + todo |

> **实施完成（2026-08-12，测试 1701 全绿，fmt/clippy 干净）**：
> M1–M5 全部落地。Kmer 字节键（`key.rs`）+ radix 字节泛化 +
> FastK golden 对照（k=21/51/81，`tests/kmer/fastk_k*.golden.gz` +
> `tests/kmer/m1.fa.gz`，逐条一致）；KmerTable/quality/norm/map/pgi/
> tadpole 全部打包字节化，qcheck 查询侧经 `key_to_kmer` 转 Kmer；
> `.pkt`/`.pgi` 格式不变、旧缓存直接兼容（字节实测与 FastK 一致）；
> k=81 端到端 2939 条与 FastK 一致；tadpole 的 `Vec<u64>` 小端 Kmer
> 由 `key::Kmer` 薄包装替代（`base_at` 索引反转保持旧语义）。
> **FastK `-p` 端到端对照（2026-08-12 实测）**：`FastK -p` + `Profex -z`
> 对 `tests/kmer/m1.fa.gz` 的 profile run（start/end/depth）与 pgr
> `self_profiles` 在 k=21/51/81 逐 run 一致；唯一差异是 Profex 尾 run
> 不闭合（旧版缺陷，§2.2），pgr 完整闭合尾 run（`2101-3145` 等），
> 与 §8 裁定一致。
> 基准（MG1655 4.6 Mb）：count_mg1655 219 ms；canonical 发射
> 25.6 ms（双窗口滚动，k=17 每窗口 ~20 ops，字节表示的必要成本，
> 相对 u128 键基线：count 整体 +29%（u128 ≈170 ms，由初版字节
> criterion change +75% 反推），canonical 纯发射 +780%
> （25.6 vs 2.9 ms）。

> **性能修复（2026-08-12，向 FastK 学习后的优化）**：初版字节实现
> count 297 ms（+75%）→ 双窗口滚动 219 ms（+29%）→ **158 ms**。
> 拆解（MG1655 k=17）：发射 ~27 ms、radix 排序 45 ms（独立测，
> 含 clone）、分组 ~47 ms；158 ms 已快于 u128 基线（~170 ms）与
> FastK 8 线程实测（188 ms）。关键优化：
> * **收集直落字节**：`build_table`/`table_profiles` 不再收集
>   `Vec<Kmer>`（40 B/条）再转字节，emit 直接 `extend` 打包字节
>   （省 460 万 × 40 B 拷贝）；
> * **emit 传引用**：`canonical_keys` 闭包参数 `&Kmer`，消除每窗口
>   40 B 值拷贝；
> * **`append` 移动缓冲**（-27%）：`build_table` 把 per-seq 字节
>   缓冲 `append` 移动而非 `extend` 拷贝（单 contig 23 MB 零拷贝）；
> * **canonical 半长比较**（学 FastK `KMd2`）：正/反链镜像对称，
>   比较前 `ceil(kb/2)` 字节即可——实测非瓶颈，保留语义；
> * 尝试过 u64 块滚动（k=17 时块转换开销 > 字节循环），回退。
> 结论：字节表示的真实开销集中在发射（~27 ms）与 radix（45 ms），
> 通过消除自身拷贝浪费（`Vec<Kmer>` 中间层、23 MB 复制）即可
> 反超 u128 与 FastK；k=81 能力与打包内存由此无性能代价获得。

> **线程对比（2026-08-12 实测，MG1655 k=17 建表）**：
>
> | 配置 | pgr | FastK | 差距 |
> |---|---|---|---|
> | 单线程 | 347 ms | 481 ms（`-T1`） | pgr 快 28% |
> | 8 线程 | 158 ms | 188 ms（`-T8`） | pgr 快 16% |
>
> pgr 单线程 347 ms 分解：发射 27 ms（8%）、radix 排序 240 ms（69%）、
> 分组/收集 ~80 ms。可比性：pgr 基准从内存序列开始，FastK 含文件读 +
> 写表（481 ms 中一部分），真实差距略小于 28% 但方向不变。多线程
> 加速均受内存带宽限制：pgr 2.2×、FastK 2.56×。**下一步优化方向 =
> radix 排序本身（单线程 240 ms）**，对照 FastK `MSDsort.c`/`LSDsort.c`
> 的分块 + MSD/LSD 混合 + shell 阈值策略。

> **radix 排序对照 FastK 的结论（2026-08-12）**：
> * 精确拆解（MG1655 k=17，单线程）：发射 27 ms、radix 排序
>   240 ms（par 单线程 242 与顺序版 243 相当，无 par 开销）、
>   分组扫描仅 12.6 ms、其余 ~67 ms 为收集/分配。
> * FastK `MSDsort.c` 结构 = American-flag MSD（`radix_sort` 递归
>   分桶 + cycle 置换）+ shell 小段（`S_thr0/1/2`、`GAP1/2`）+
>   **排序叶子段直接 `COUNT`**（`hist_kmers`/`invert_kmers`，省二次
>   分组扫描）+ 按第一字节分块多线程（`PARTS`/`sort_thread`）；
>   `LSDsort.c` 用于 profile（迭代 LSD）。pgr `radix_sort_bytes`
>   与之同源同构（同一 Myers MSD 家族），单线程 240 ms 与 FastK
>   排序同量级。
> * **尝试过且无收益、已回退**：u64 块滚动（k=17 块转换开销 > 字节
>   循环）、置换路径 `copy_within` 替代逐字节拷贝（编译器已合成为
>   memcpy，copy_within 反而 +1.6%）。半长比较（学 FastK `KMd2`）
>   非瓶颈但语义正确，保留。
> * 结论：单线程 radix 240 ms 是 American-flag MSD 的随机访问
>   本质成本；FastK 的进一步手段（super-mer 加权、LSD）依赖其
>   数据模型（2026-08-14 起 pgr 已实现 super-mer 两段计数原型，
>   收益有条件，见 §3.6；LSD 方向不迁移）。当前 pgr
>   单线程/多线程均快于 FastK，此方向收尾。

### 12.4 风险与决策点

* **tadpole 迁移是最大工程量**：组装热路径（asm contig/unitig、fq
  extend/ec-kmer/merge）换键类型，`push_right`/`push_left` 从 u64 位移
  改字节位移——M5 基准保性能（字节滚动预期不劣于 u64 版，需实测）。
* **字节序必须与 FastK 精确一致**：已实测一致（字节内大端，与现有
  `pack_kmer` 相同）——M1 用 FastK 二进制对照锁定，M2 无需处理旧
  `.pkt` 缓存（格式不变、直接兼容）。
* **存储形态统一**：排序表打包（FastK 式，`Vec<u8>` 连续区）与组装侧
  HashMap 值键（`Vec<u8>`）共用同一套字节编码/操作；**不做定长对象键**
  （32 B 内存翻倍，用户否决）。k 无上限由 `key_bytes` 参数决定。
* **k≤64 性能与内存**：radix 有效字节数不变；打包存储下 pgi/map
  （k=31/40）内存**下降**（10-12 B vs 16 B）；M5 基准确认（预期 ≤10%）。
* **canonical 语义**：FastK 逐字节取小与现有 u128 `min(fwd,rc)` 都是
  字典序取小，k≤64 集合不变——M2 用 golden 回归锁定。

### 12.5 参考实现对照记录（2026-08-12 补）

> 实施时 `push_right`/`push_left` 字节滚动与 radix 字节版为推导实现、
> 以 golden 锁定；此处对照 FASTK-master 源码确认等价性，作为依据追溯。

* **表条目布局**（`count.c` `kmer_list_thread`）：按 k-mer 最高字节值分桶
  （`fours[kf]`），桶内条目 = 占位 `0` + 后续字节（最高字节由桶身份
  提供，最终落在前缀索引）+ u16 count；`fill[-1] &= KCLIP` 裁余位。
  pgr `KmerTable` 直接存完整字节（`keys` 打包区）+ u32 count，与 FastK
  表字节逐条一致（golden 实测：`acgta` = `1b 00`，k=5）。
* **canonical**（`count.c`）：正/反链逐字节比较取小，k%4 余位 KCLIP
  裁剪；pgr `Kmer::canonical` 字节比较等价（与比较起点无关，§12.1）。
* **排序**（`MSDsort.c`）：American-flag MSD radix，直接操作字节数组，
  从最高字节（digit=0）向低位推进；pgr `radix_sort_bytes` 同构（分离
  key/payload 并行交换、无限 cycle 栈、insertion 小段）。排序结果
  均为字节升序，等价。
* **窗口滚动**（`count.c`）：FastK 在 super-mer 连续位流中滚动（
  `fptr`/`rptr` 位偏移 `fs`/`rs`），k-mer 字节按偏移提取——super-mer
  位流编码本身 pgr 不复刻（§3.6 用固定尺寸记录 + `Kmer::from_bytes` +
  滚动展开，语义等价）；pgr 的 Kmer 值类型直接维护
  打包字节（`push_right`/`push_left` 整体移位），语义与位流滚动等价
  （golden 逐条一致），差异是内存模型配套的选择。
