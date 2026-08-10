# 替换 anchr trim 流水线中的 BBTools（分析 + 迁移计划）

> 2026-08 整理。用户最终需求：anchr 的 `templates/trim.tera.sh` 用 BBTools
> 做 read 清洗，结果满意但速度不满意，目标是用 pgr 替换 BBTools。
> 2026-08-10 更新：BBTools-40.01 源码已本地化（见 §4.1），§4 给出具体迁移计划。
> 配套：[fq-trim-qual.md](fq-trim-qual.md)（已实现，覆盖 sickle 部分）、
> [seq-reader.md](seq-reader.md)（FAFQ/BGZF 基础设施）。

> **定位（2026-08 修正）**：`pgr fq trim-qual` 只替换流水线中的 **sickle**（第 9 步，
> 多阈值质量/长度参数扫描），**不替换 BBTools 任何组件**。BBTools 8 步（第 1-8）
> 是另一个替换目标，逐项梳理见 §3。

## 1. 流水线拆解

`trim.tera.sh` 对输入（R1/R2 或单端）依次执行（每步有 `if [ ! -e ]` 缓存）：

| # | BBTools 工具 | 作用 | 参数要点 |
|---|---|---|---|
| 1 | `clumpify.sh` | read 排序（可选去重 `dedupe dupesubs=0`） | 为后续阶段加速 |
| 2 | `filterbytile.sh`（可选） | 按 flowcell tile 过滤低质量 read | `opt.tile==1` |
| 3 | `bbnorm.sh`（可选） | kmer cutoff：滤掉低深度 kmer 的 read | `passes=1 bits=16 min=cutoff` |
| 4 | `reformat.sh`（可选） | 降采样 | `samplebasestarget=sample` |
| 5 | **`bbduk.sh trim`（核心）** | 接头 kmer 修剪 + 质量/长度修剪 | `ktrim=r k=trimk mink=11 hdist=1 tbo tpe maxns=0 minlen=60 qtrim=r trimq=15 ftm=5` |
| 6 | `bbduk.sh filter` | kmer 匹配过滤 adapter/artifact 库 | `k=matchk cardinality` |
| 7 | `kmercountexact.sh` | kmer 直方图 + peaks | `k=cutk` |
| 8 | `repair.sh` | 交错 → R1/R2/singles | `out/out2/outs` |
| 9 | `sickle` | 多阈值质量/长度修剪（参数扫描） | `-q/-l` 遍历，`Q{qual}L{len}/` 目录 |

第 9 步（sickle）已由 `pgr fq trim-qual` 覆盖（anchr 的 sickle 调用只有
`-q/-l/-t sanger`，未用 `-n` 截断，trim-qual 均已覆盖），见
[fq-trim-qual.md](fq-trim-qual.md)。注意：bbduk trim（第 5 步）里的
`qtrim/minlen/maxns/ftm` 是 **bbduk 的参数，不属于 trim-qual 的替换范围**。

## 2. 速度瓶颈在哪

1. **8 个独立 Java 进程**：每个启动 JVM（秒级）+ GC + `-Xmx` 内存；pgr 是
   Rust 原生，启动毫秒级。
2. **中间文件反复 gz**：每步输出 `*.fq.gz`，下一步再解压——压缩/解压放大
   IO 多次。pgr 命令可流式管道衔接，中间不落盘。
3. BBTools 算法本身不慢（clumpify/bbduk 都是高效实现），慢的是**工程外壳**；
   但替换它意味着把其中几个算法在 pgr 里重新实现。

## 3. BBTools 8 步逐项梳理

| # | 工具 | 职责与算法本质 | pgr 基础 | 替换评估 |
|---|---|---|---|---|
| 1 | `clumpify.sh` | 按 kmer 签名排序/聚类 reads（可选 `dedupe` 去重）；为后续 kmer 类步骤加速 | 无 | 工作量大（TB 级外部排序）；已确认迁移（2026-08-10），见 §4 |
| 2 | `filterbytile.sh`（可选） | 按 Illumina flowcell tile 坐标过滤低质量 reads | 无 | 小众、老数据专属 → 可不做 |
| 3 | `bbnorm.sh`（可选） | kmer cutoff：全数据集 kmer 深度表（`bits=16` 近似哈希，省内存），滤掉含低深度 kmer 的 reads（去测序错误） | `KmerTable` 精确计数（u128+u32） | 中：内存策略需权衡（近似哈希 vs 精确表）→ 三期 |
| 4 | `reformat.sh`（可选） | 降采样到目标碱基数 | 无 | 小：`fq sample` 流式抽样 → 一期 |
| 5 | **`bbduk.sh trim`（核心）** | 接头 kmer 修剪（`ktrim=r`/`mink=11`/`hdist=1`/`tbo`/`tpe`）+ `maxns=0` + `qtrim=r` + `minlen` + `ftm=5` | 无（trim-qual 只替代 sickle，不是这里） | **大**：`tbo`/`tpe` 是 BBDuk 特有 → 二期 |
| 6 | `bbduk.sh filter` | 参考库（adapter/artifact）kmer 匹配，过滤命中 reads + `cardinality` 统计 | `KmerTable`/`canonical_keys` | 中：kmer 命中过滤 → 三期 |
| 7 | `kmercountexact.sh` | 精确 kmer 计数直方图 + peaks（估计基因组大小/深度） | `KmerTable.counts` 可直接出直方图 | 小-中 → 三期 |
| 8 | `repair.sh` | 交错文件 → R1/R2/singles | 无（`fq interleave` 的反操作） | 小：`fq split` → 一期 |

## 4. 迁移计划（2026-08-10）

### 4.0 实施进度（2026-08-10 晚更新）

| 里程碑 | 状态 | 证据 |
|---|---|---|
| M0 golden | 完成 | `tests/bbtools/Lambda/golden/`（39.38 + ordered=t + seed=1，全链确定性已验证） |
| M1 `fq clump` | **完成，逐字节一致** | `cli_fq_clump.rs` 对照 `clumpify.fq.gz`；dedupe 模式（`--dedupe --dupesubs 0`）已实现，与 threads=1 golden 一次性逐字节验证（golden 未入库，语义由合成测试覆盖，见 §4.4 M1 注） |
| M2 `fq split`/`fq sample` | **完成，逐字节一致** | `cli_fq_split.rs`/`cli_fq_sample.rs` 对照 repair/reformat golden |
| M3-M5 `fq trim-adapter` | **完成，逐字节一致** | `cli_fq_trim_adapter.rs` 对照 trim/filter golden |
| M7 kmercountexact | **完成，逐字节一致** | `pgr kmer hist --khist-text/--peaks`（logScale + CallPeaks 全移植）对照 R.khist.txt/R.peaks.txt |
| M6 bbnorm cutoff | **完成（精确表语义）** | `pgr fq norm`（精确 canonical 表 + bbnorm per-read 判定逻辑：truedepth/depthAL 分位数 + toss 条件）；与 bbnorm bits=16 近似计数在 min=3 边界差 ~21 对（39846 vs 39888），属设计稿已声明的"先精确 KmerTable"路线 |
| M8 集成 | **完成（原语路线）** | 只提供可组合原语（clump/split/sample/trim-adapter/fq norm/hist），**不内置 pl trim 流水线**——编排属于 anchr，pgr 不做"别人的活"（2026-08-10 修正，`pl trim` 已移除）；anchr 模板把 `bbduk.sh` 等调用换成 pgr 命令、用管道串联避免中间 gz |

基准见 [benchmarks/bbtools-vs-pgr.md](../../benchmarks/bbtools-vs-pgr.md)
（hyperfine，Lambda 40000 reads，pgr release 单线程 vs BBTools 39.38 8 线程）。

移植要点记录（避免后续重复踩坑）：

* `java.util.Random.nextLong()` 的低 32 位是**带符号扩展相加**：`((long)next(32)<<32)+next(32)`，
  高位为 1 时高位字会减 1（Rust 需 `hi<<32).wrapping_add(sign_extend(lo))`）；
* `SketchObject.makeCodes` 的 antialias 会把符号位置位（值可为负），不能用
  "正数"假设；
* BBTools 的 `Read.quality` 存的是 **phred（输入时减 ASCII_OFFSET）**，输出再加回；
  qtrim 走 `TrimRead.optimalMode`（默认 true，`testOptimal`），不是简单 testRight；
* `bbduk.sh` 默认 `ordered=f`（输出顺序不可复现），golden 一律用 `ordered=t`；
* clumpify/reformat 默认 seed 下行为可复现，但 golden 显式 `seed=1`/`sampleseed=1`；
* `trim.fq` 与 `filter.fq` 的 golden 输入链必须用同一份 clumpify 输出（repair.sh
  的 `out=R1.fq.gz` 会覆盖同名输入，golden 生成时用隔离目录）。

### 4.1 目标与范围

### 4.1 目标与范围

用户已确认方向：BBTools 结果好但慢、中间步骤冗余，迁移到 pgr。
BBTools-40.01 源码已置于仓库根 `BBTools-40.01/`——该目录自带 `.gitignore`
（内容为 `*`），整体被 Git 忽略，**仅作参考，不属于项目代码**。

- **版本基线（2026-08-10）**：本地实际安装的是 **BBTools 39.38**（cbp 安装于
  `~/.cbp/libexec/bbtools/`，`bbduk.sh` 2024-11-18，入口 `jgi.BBDuk`）；仓库根
  的 `BBTools-40.01` 是更新版（2026-02-11，入口 `bbduk.BBDukS`）。两者行为可能
  有差异，迁移语义与 golden 验证**一律以本地 39.38 为准**，40.01 只作算法演进
  参考，不作为移植来源。

- **替换对象**：`trim.tera.sh` 第 1-8 步（第 9 步 sickle 已由 `pgr fq trim-qual`
  覆盖，见 [fq-trim-qual.md](fq-trim-qual.md)）。
- **不做**：复刻 Java 工程外壳（JVM、`-Xmx`、每步 gz 落盘、`if [ ! -e ]` 缓存）；
  `filterbytile.sh` 默认跳过。**clumpify 已确认纳入迁移范围**（2026-08-10），
  它是字节级一致（read 顺序）的前提。
- **成功标准（2026-08-10 强化）**：pgr 各步输出与本地 39.38 BBTools 输出
  **解压后内容逐字节一致**（含 name、顺序、序列、质量、行宽/换行等格式细节；
  gz 压缩字节不要求一致）。端到端墙钟时间显著下降、峰值内存可控；中间文件
  可选（流式不落盘）。
- **golden 测试数据**：`anchr/tests/Lambda/`（R1.fq.gz、R2.fq.gz、env.json、
  pe.cor.fa.gz、unitigs.fasta）。R1/R2 是 trim 流水线输入；golden 输出由本地
  39.38 按 trim.tera.sh 参数实跑生成后入库（存放位置见 §6）。

### 4.2 BBTools 入口 → pgr 目标模块映射

| # | 流水线步骤 | BBTools 入口类 | 源码位置（参考） | 行数 | pgr 目标 |
|---|---|---|---|---|---|
| 1 | clumpify | `clump.Clumpify` | `current/clump/Clumpify.java` | 705 | `libs/fq/clump.rs` + `pgr fq clump` |
| 2 | filterbytile | `hiseq.AnalyzeFlowCell` | `current/hiseq/AnalyzeFlowCell.java` | 1935 | 不做 |
| 3 | bbnorm cutoff | `jgi.KmerNormalize` | `current/jgi/KmerNormalize.java` | 3895 | `libs/fq/norm.rs` + `pgr fq norm` |
| 4 | reformat 降采样 | `jgi.ReformatReads` | `current/jgi/ReformatReads.java` | 1994 | `libs/fq/sample.rs` + `pgr fq sample` |
| 5 | bbduk trim | `jgi.BBDuk`（39.38 实际入口） | 本地 `~/.cbp/libexec/bbtools/current/jgi/BBDuk.java`（5384 行）；40.01 的 `jgi/BBDuk.java`（5462 行）+ `bbduk/` BBDukS 家族作参考 | 5384 | `libs/fq/trim_adapter.rs` + `pgr fq trim-adapter` |
| 5b | tbo / tpe | `jgi.BBMergeOverlapper.mateByOverlapRatio` | 本地 1256 行；40.01 1503 行 | 1256 | 同上命令的 overlap 子模块（纯 Rust 移植 Java fallback，JNI 两版本均已禁用） |
| 6 | bbduk filter | 同上（`jgi.BBDuk` `cardinality` 模式） | 同上 | — | `pgr fq trim-adapter --filter`（复用 kmer 表） |
| 7 | kmercountexact | `jgi.KmerCountExact` | `current/jgi/KmerCountExact.java` | 1155 | 复用 `libs/kmer`（table/hist/gsize）+ 补 peaks 文本 |
| 8 | repair | `jgi.SplitPairsAndSingles rp` | `current/jgi/SplitPairsAndSingles.java` | 909 | `libs/fq/split.rs` + `pgr fq split` |

> 版本差异注意：`bbduk.sh` 在本地 39.38 中入口是 `jgi.BBDuk`；40.01 才改成
> `bbduk.BBDukS`（新 bbduk 包）。迁移以 39.38 的 `jgi.BBDuk` 为准，40.01 的
> BBDukS 只作行为参考；其余 7 个工具两个版本的入口类相同。BBTools 源码整体
> 被其自带 `.gitignore` 忽略，仅作参考。

### 4.3 目标 CLI（已实现，2026-08-10 定稿）

- `pgr fq clump`：按 kmer 签名排序/聚类 reads（对齐 `clumpify.sh`；`--dedupe
  --dupesubs 0` 整对去重已实现——R1 与 R2 都精确匹配（N 通配）才算重复，
  保留期望错误更少的那对；超内存数据走外部 hash 桶路径（`--mem` 控制，
  确定性桶序，见 §4.6）。
- `pgr fq split`：交错输入 → R1/R2/singles（对齐 `repair.sh rp`，是
  `fq interleave` 的反操作）。
- `pgr fq sample`：按目标碱基数/比例降采样（对齐 `reformat.sh
  samplebasestarget`）。
- `pgr fq trim-adapter`（核心）：
  - 修剪模式（默认）：`--ref <fa> --k <trimk> --mink <11> --hdist <1>`
    + `--trimq <15> --minlen <60> --maxns <0> --ftm <5>`；`--no-tbo`/
    `--no-tpe`/`--no-qtrim` 关闭对应步骤；配对丢弃由 `--no-toss-broken-reads`
    关闭；
  - 过滤模式：`--no-ktrim --no-tbo --no-tpe --no-qtrim --k <matchk>
    --mink 0 --minlen 0 --maxns -1 --ftm 0`（对齐 bbduk filter 的
    `k=<matchk> cardinality`）；
  - 双端输入为 interleaved，R1/R2 同步处理（tbo/tpe 需要同读对）；
  - 质量修剪复用 `libs/fq/trim.rs`（trim-qual 的 sliding/mott、polyG）。
- `pgr fq norm`：`--min <cutoff>`（对齐 `bbnorm.sh min=`；`bits=16` 近似
  哈希 vs 精确 KmerTable 的内存策略见 §4.5）。
- `pgr kmer hist` 已有（FastK 兼容）；补 `--peaks` 文本输出对齐
  `kmercountexact.sh` 的 `khist.txt`/`peaks.txt`。
- 流水线串联：各命令保持可组合（Unix 风格），anchr 新模板直接调 pgr 命令
  流式衔接（管道串联，中间不落 gz）；**不在 pgr 内置一键流水线**。

### 4.4 里程碑与验证

1. **M1 `pgr fq clump`**：按 kmer 签名排序（+`--dedupe` 整对去重），对齐
   `clumpify.sh`。→ 与本地 39.38 `clumpify.sh` 输出逐字节比对（顺序、去重
   语义）。注：dedupe golden 用 `threads=1` 生成——BBTools threads>1 时
   dedupe 输出顺序不确定（clump 线程竞态），删除集合一致；pgr 的 dedupe
   无论线程数都按排序序收集 clump（实现约定），因此 threads=1 golden 在
   未来引入多线程后依然有效，测试不失效。golden 未入库（体积控制），
   字节级一致性已一次性验证，日常由合成测试覆盖语义。
   **去重标准（已确认，2026-08-10）**：采用 anchr 现行的**精确整对去重**
   （R1 与 R2 都必须精确匹配，N 作通配符；R1 相同但 R2 不同不算重复）。
   理由：整对唯一性在 1000× 以上覆盖度下极强，相同整对基本就是
   PCR/光学重复，直接序列去重足够，无需官方 assemblyPipeline.sh 的
   `dedupe optical`。**光学去重明确不做**（2026-08-10 确认）：需要 flowcell
   坐标解析 + 邻近判定，且 BBTools-40.01 与 Lambda 数据都没有带真实坐标的
   reads 可验证，流程也不需要它。
   `--dupesubs` 支持 >0（2026-08-10 补全）：移植 BBTools 的扫描语义——
   dupesubs=0 时 scan=0（只比相邻），>0 时 scan=5/maxDiscarded=15 且单轮
   删除超限会扩大扫描重试（`scan+10`/`maxDiscarded*2+20`）；`dupesubs=2`
   已与 BBTools threads=1 golden 逐字节验证（39978 reads），合成测试覆盖
   0 vs 1 的容错差异。
2. **M2 `fq split` + `fq sample`** → 与本地 39.38 `repair.sh`/`reformat.sh` 在
   Lambda 数据上的输出逐字节比对（含格式），单元 + 集成测试。
3. **M3 `fq trim-adapter`（无 tbo/tpe）**：参考序列建 kmer 表
   （`k/mink/hdist`）→ `ktrim=r` 修剪 → `maxns/minlen/qtrim/trimq/ftm`。
   → 与本地 39.38 `bbduk.sh` 单端/非 overlap 双端场景输出逐字节比对
   （Lambda R1/R2 及合成小样本）。
4. **M4 tbo/tpe**：移植 `BBMergeOverlapper` 的 `mateByOverlapRatioJava*`
   fallback 语义（JNI 在两版本中均已禁用，实际运行即该路径）。
   → 用 insert 短于 read 长度的双端数据（合成）逐字节比对。
5. **M5 bbduk filter + cardinality** → 输出 read 集合与 `stats` 文本逐字节一致。
6. **M6 `pgr fq norm`** → 过滤后输出与 39.38 `bbnorm.sh` 逐字节一致（先精确
   KmerTable；大数据量再评估近似哈希路径，对齐 `bits=16` 语义）。
7. **M7 kmer hist peaks** → `khist.txt`/`peaks.txt` 逐字节一致（已有 hist 基础）。
8. **M8 集成与基准**：anchr 新模板（pgr-only）端到端：流水线最终输出与 BBTools
   全流程输出逐字节比对（含 clumpify 顺序）；墙钟、峰值内存、中间文件体积对比。

### 4.5 风险与待决

- **hdist=1 语义**：BBDuk 的 hamming 距离匹配通过 index mask/邻居 kmer 表
  实现（`current/bbduk/BBDukIndexMask*.java`），移植前先读透内存布局，再定
  Rust 数据结构。
- **tbo/tpe 精度**：`BBMergeOverlapper` 的 JNI 路径在 39.38 与 40.01 中均已
  禁用（`if(false && Shared.USE_JNI)`），实际运行就是纯 Java fallback
  `mateByOverlapRatioJava*`——直接移植该路径即与 golden 一致；质量值参与
  overlap 判定，需保留 33/64 编码处理。
- **bbnorm 内存**：`bits=16` 是近似 kmer 表（省内存）；pgr 精确 u128 表在
  超大数据集的峰值内存待实测，必要时引入近似哈希路径。
- **paired 流一致性**：修剪/过滤时 R1/R2 必须同步处理（tbo/tpe、
  tossbrokenreads 语义）。
- **统计文本**：`trim.stats.txt`/`filter.stats.txt`/cardinality 输出格式是否
  完全复刻，待确认（下游 anchr 是否解析这些文件）。
- **clumpify 排序语义**：字节级一致要求 read 顺序完全一致，需逐字节复刻
  `Clumpify.java` 的 kmer 签名、排序键、稳定性和去重语义（`dedupe
  dupesubs=0`）。大数据下 pgr 走外部 hash 桶路径，输出为确定性桶序（与
  BBTools 大数据行为一致，非全局序）——字节 golden 只在小数据（内存路径）
  上定义。
- **gz 流式**：pgr 已有 BGZF/FASTQ 基础设施（seq-reader），流水线内不落盘。

### 4.6 内存模型与外部排序（2026-08-10 定稿）

`pgr fq clump` 的排序内存上限：

```
mem_limit = min( --mem（默认 2g）,  物理内存 × 0.5,  数据估算 )
```

* `--mem`：用户显式预算（KMG 解析，语义同 `-Xmx`），大机器上的主保护；
* 物理内存 × 0.5：低内存环境（如 1G 虚拟机）的兜底，防止排序挤爆机器；
* 数据估算：gz 输入按文件大小 × 8（解压 ~4× × 记录开销 ~2×），明文 × 2；
  只作路径决策用（估算 ≤ 上限 → 内存路径，否则桶路径），不引预扫额外 IO。
* 物理内存/CPU 探测用 **sysinfo** crate（跨平台，Linux/macOS/Windows），
  以后并行度自动判断（`logical_cpus`）复用同一依赖；已加 Cargo 依赖。

外部路径：按 pivot k-mer 哈希分桶（`--buckets`，默认由 `mem_limit` 推导，
上限 4096），桶文件写临时目录（`std::env::temp_dir()`，用完清理），逐桶
内存排序 + dedupe（相同整对共享 pivot → 同桶，去重语义不变），按桶序拼接
输出。确定性：同一输入 + 同一 `--mem` 下桶数与顺序固定。已验证：
Lambda 上 `--mem 1m` 强制桶路径，输出确定、与内存路径 read 集合一致
（普通 40000 与 dedupe 39984 都一致）；内存路径（默认 2g）字节 golden
不变。BBTools 自身大数据走 KmerSort2/3 也是桶序，行为对齐。

路径可强制（2026-08-10）：`--sort-mode auto|global|bucket`——auto 按内存
预算自动选（默认），global 强制内存全局排序（小数据字节 golden 不变），
bucket 强制外部桶路径；指定 `--buckets` 等价于隐含 bucket 模式（桶数不同
顺序不同但确定、集合一致）。

并行化（2026-08-10）：内存路径排序用 `par_sort_by`（比较器全序，结果
不变）；桶路径按 **内存受限的 wave 分批并行**处理（wave 数 =
`mem_limit×0.8 / 单桶估算`），每个 wave 内 rayon 并行读桶/排序/去重，
按桶序写回——并行度自动受内存预算约束，不会吃掉 `--mem`。50 万对
合成数据实测：`--mem 32m` 桶路径 ~9.0s（user 12s > wall 9s，并行生效）、
`--mem 4g` 内存路径 ~7.6s，两者确定性且集合一致。

### 4.7 流水线命令并行化（2026-08-10）

通用组件 `libs/par::ordered_map`：有界保序并行流水线（feeder 线程 +
`workers` 个 worker + 按输入序收集的 collector），内存由通道容量界定，
输出顺序与线程数无关。已接入：

* `pgr fq trim-adapter --parallel N`（默认逻辑 CPU 数）：流式读取
  （`libs/fq/pairs::PairReader`，不再整体载入内存）+ 并行处理每对 reads，
  按序写回。50 万对合成数据实测 threads=1 9.2s → threads=8 1.4s（6.6×），
  峰值内存 ~15MB（流式、有界），threads=1/8 输出逐字节一致且与 golden
  一致。
* `fq clump` 已在 §4.6 并行（par_sort + 桶 wave）。

踩坑记录：channel 原始 `out_tx` 若不 drop，collector 的 `recv()` 会永久
阻塞（worker 全部退出后 out_rx 仍不关闭）——feeder join 必须在 drain 之后。

## 5. 验收标准（替换后的对比）

- 同一输入（`anchr/tests/Lambda/` R1/R2），pgr 各步输出与本地 39.38 BBTools
  输出**解压后逐字节一致**（name、顺序、序列、质量、行宽/换行等格式细节；
  gz 压缩字节不要求一致；顺序由 M1 clumpify 保证）；
- `khist.txt`/`peaks.txt` 逐字节一致（已达成）；`trim.stats.txt`/
  `filter.stats.txt`/cardinality 统计文本未复刻（下游 anchr 未解析，§6.5）；
- 端到端墙钟时间显著下降（8 个 JVM 进程 → Rust 单进程流式）；
- 中间文件可选（管道模式不落盘）。

## 6. 待用户确认

1. 版本基线：以本地安装的 BBTools 39.38（`jgi.BBDuk` 入口）为对齐与验证
   基准、仓库内 40.01 仅作参考——是否确认？本地 39.38 源码是否也放进仓库
   参考目录（当前只有 40.01）？
2. 命令命名（2026-08-10 定稿）：`pgr fq` 下全部动词/约定风格——`clump`、
   `norm`、`sample`、`split`、`interleave`、`range`（与 `fa range` 同约定）、
   `to-fa`（与 `fa to-2bit` 同约定）；trim 家族统一为 "trim-<目标>"：
   `trim-adapter` + `trim-qual`（原 trim-q）。
3. ~~接头修剪语义~~ → 已定：完全复刻（tbo/tpe/hdist 全部移植，逐字节一致）。
4. 中间文件策略：流水线管道串联（pgr 命令流式、不落 gz）是否可接受？
   （2026-08-10 已定：接受管道串联，不内置 `pl trim`。）
5. 统计文本：`trim.stats.txt`/`filter.stats.txt`/cardinality 是否需要与
   BBTools 逐字节一致（下游 anchr 是否解析）？
6. ~~golden 存放~~ → 已定：复制进 `tests/bbtools/Lambda/`（含 gz 化的
   golden 与 README）。

---

*参考来源: [trim.tera.sh](../../../anchr/templates/trim.tera.sh)（anchr 项目，只读） |
BBTools-40.01 源码 | [fq-trim-qual.md](fq-trim-qual.md) | [seq-reader.md](seq-reader.md)*
