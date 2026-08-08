# FastK / Profex 原生迁移设计（方案 B）

> 状态：**已实现（2026-08-09）**。目标：让 `pgr rept s-kmer` / `pgr rept e-kmer` 不再依赖外部
> FastK / Profex。参考源码：仓库内 `FASTK-master/`（2025-09-13 下载的快照，
> **等于 FASTK-1.2**，README 标注 Current: April 18, 2021，**不是上游当前
> master**，见 §2.3 版本核对）；行为契约以本机安装的 FastK / Profex
> （CBP 安装，2025-03 构建，上游 commit ddea6cf，**无源码补丁**，见 §2.3）
> 实测为准。
> **旧管线视为实验性参考**，迁移不要求输出与旧管线字节级一致，不保留已知
> 缺陷（尾 run quirk，见 §2.2）。

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
§4.3）。不做 super-mer、磁盘分桶、`.ktab/.prof` 外部格式兼容 —— 那是 FASTK
为 TB 级数据设计的，pgr 场景用不上。

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

不做 super-mer：内存中 `KmerTable` 用 u128 key（与 pgi 一致），
u128 + u32 ≈ 20 B/唯一 k-mer：5 Mb 细菌 ~5 M key ≈ 100 MB；
50 Mb 真菌 ~50 M key ≈ 1 GB，可接受（用完即释放）。超大输入的
后备（分块计数）不在本期范围，接口上 `build_table` 与 `profile`
分离即可，将来可换实现。持久化用紧凑编码，与内存表示无关（§4.3）。

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

* 新格式：`<库>.pgrk` 单文件（`lib.fa` → `lib.pgrk`，`lib.fa.gz` →
  `lib.fa.pgrk`），**紧凑编码**：
  header（magic/version/k/条目数/key 字节数）+ 每条目
  `packed key（ceil(2k/8) 字节，复用 pgi 的 pack_kmer）+ u32 count`；
  k=17 时约 9 B/条目，5 Mb 库 ~45 MB（对比裸 bincode u128 的 ~100 MB）。
  bincode 只当容器；不再需要 FastK 的隐藏分片，也去掉 `.complete` 标记——
  完整性由 header 校验兜底（损坏即重建），写入用原子 rename（临时文件 +
  rename）。`cache_is_fresh` 保留 mtime 检查，判断对象从 `.ktab/.complete`
  变为单个 `.pgrk` 文件。
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
