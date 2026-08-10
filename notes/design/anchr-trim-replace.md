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
| M2 `fq split`/`fq sample` | **完成，逐字节一致** | `cli_fq_split.rs`/`cli_fq_sample.rs` 对照 repair/reformat golden；2026-08-10 晚复核 sample：6 组 target/seed（含超总量边界）+ 单端输入，与 39.38 `reformat.sh samplebasestarget` 逐字节一致（FastRandomXoshiro/allowUpsample=false/per-pair 决策全对上） |
| M3-M5 `fq trim-adapter` | **完成，逐字节一致** | `cli_fq_trim_adapter.rs` 对照 trim/filter golden；`--stats` 输出 bbduk `stats=` 3 列格式，与 39.38 逐字节一致（#File 路径行除外，见 §6.5）；2026-08-10 晚复核：19 组 trim 变体（k/mink/hdist/minlen/trimq/ftm/maxns/tbo/tpe/qtrim/组合）+ 3 组 filter k 变体 + 质量边界，与 39.38 `bbduk.sh ordered=t` 逐字节一致；修复 changequality 与 qtrim 空 read 边界（见下） |
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
* bbduk 默认 `changequality=t`（读入即做）：ACGT 碱基质量 <2 的提到 2、
  N 碱基强制 0——发生在 minbasequality/maxNs 等判定**之前**；pgr 原实现
  只加回 ASCII 偏移，默认参数下被 maxns=0 掩盖，`maxns=-1` 时暴露
  （修复：`make_read_buf` 读入时钳制）。
* qtrim 的边界语义（`TrimRead.trimByAmount(r, 0, right, 1)`）：
  `right >= len` 时剪到 `max(1, len-1)`——len=1 的 read 被剪成空 read
  （minlen=0 时保留输出），len≥2 保留 1bp；pgr 原实现 `saturating_sub`
  在 len=1 时错误保留 1bp、len≥2 全剪时错误剪空（修复：复刻钳制公式）。
* maxNs 检查在 qtrim **之后**（bbduk 顺序 qtrim → minlen → maxNs）——
  qtrim 剪掉 N 后剩余 N 数才参与 maxNs 判定；交叉验证用 N×100 read 确认。

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
- **bbnorm 内存**：精确 vs 近似的完整分析记录见 §4.8（2026-08-10，
  **尚未定论**）。
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

### 4.8 bbnorm norm：精确 vs 近似分析记录（2026-08-10，未定）

**原始用法（anchr `templates/trim.tera.sh`）**：

```bash
bbnorm.sh in=temp.fq.gz out=highpass.fq.gz \
    passes=1 bits=16 min={{ opt.cutoff }} target=9999999 \
    threads=... -Xmx...
```

注释 "Remove reads without high depth kmer"——**纯 highpass filter**，这是
bbnorm 在 anchr 里的唯一用途（khist/peaks 由后续 `kmercountexact.sh` 单独
完成，与 bbnorm 无关）。

**判定语义（BBTools-40.01 `current/jgi/KmerNormalize.java`）**：

* `target=9999999`：`coin>target` 分支永不触发（`coin=1..depthAL`，
  depthAL ≤ 65535），降采样/归一化完全禁用；`passes=1` 且无 ECC
  （`USE_ECC1/ECCF` 默认 false）→ 整个调用退化为 filter。
* toss 条件只有两条：`depthAL<0`（read 中 ≥15 个 k-mer 的计数
  ≥ `max(min, high/125)` 不成立，`MIN_KMERS_OVER_MIN_DEPTH=15`、
  `ERROR_DETECT_RATIO=125`、`HIGH_PERCENTILE=0.90`）或
  `maxTrueDepth<min`（truedepth = 46 分位 k-mer 计数，
  `DEPTH_PERCENTILE=0.54`，取一对中较大者）。`MIN_LENGTH=1` 不生效。
* pair 级决策：任一 mate 满足上述条件即整对保留
  （`USE_LOWER_DEPTH=true` 取 minAL 有值者、`REQUIRE_BOTH_BAD=false`）。
* 表由**全量 reads** 构建（`tablereads=-1`），runPass 内部两遍
  （`makeKca` 全量建表 → `count` 全量过滤），确定性、与顺序无关——
  **不是对 reads 采样**（khmer 式单遍在线会顺序相关，已排除）。

**精确 vs 近似 分析（结论未定，2026-08-10 保留）**：

1. 判定只需要阈值区分（≥min 与 ≥max(min, high/125)），不需要精确大计数。
2. CMS（min-of-tables）只会高估不会低估 → 误差单向：只可能把 <min 的
   read 误判为保留，不可能误杀；而本 filter 的用途是滤掉含低深度（错误）
   k-mer 的 reads，单向高估与用途相反。
3. bbnorm 自身是近似表（bits=16 + 哈希 + minprob），结果依赖 `-Xmx`
   （装载率变 → 碰撞率变 → 边界判定变），"与 bbnorm 字节级一致"是移动靶；
   Lambda min=3 的 ~21 对差异（精确 39846 vs bbnorm 39888）是定义差异，
   不是实现缺陷。近似路径也无法保证追上（除非克隆 KCountArray 内部并锁死
   内存配置）。
4. 精确外部桶路径已实现（§4.6 同款 mem_cap 约束，输出与内存路径一致）；
   CMS "固定内存"≠"小内存"（1B unique k-mer 低装载需 ~10GB 量级，且随
   unique 数线性涨；1G 虚拟机下 CMS 判定失真，精确外部路径只是慢但正确）。
5. 速度上精确外部桶 = 多遍 I/O + 排序，与 bbnorm 自身两遍全扫描同量级。
6. 工程上精确路径已实现/已并行/已测试；CMS 是新代码 + 任意参数
   （bits/表数/哈希数/饱和），结果随配置漂移。

**倾向：精确**（外部桶路径即 1TB 答案）。转向近似的唯一场景：单机、
极小内存、无大磁盘、接受判定噪声——非 pgr/anchr 语境。待用户定稿。

**2026-08-10 复核（多参数交叉验证）补充**：

* **changequality 缺失**（与 bbduk 同源）：bbnorm 读入时 N 质量强制 0、
  ACGT 质量提到 2，输出质量随之变化；pgr 原样输出导致 Lambda 上数百处
  `#` vs `!` 质量差异。已修复（读入时钳制 + 输出一致）。合成数据上
  验证 N@2 → `!`（0）与 bbnorm 逐字节一致。
* **minq=6 建表过滤缺失**：bbnorm 的 KmerCount（bits=16 → KmerCount4）
  跳过含质量 <6 碱基的 kmer（`quals[i]<minQuality` 即重置）。已实现；
  Lambda 输入无 ACGT<6（仅 N@2，N 的 kmer 本就被跳过），故无行为变化，
  但真实低质量数据上正确。
* **minprob=0.5 名义存在但不生效**：KmerCount4.addRead 只有 minQuality
  检查、无概率乘积逻辑（minProb 在 bloom 包的 KmerCountAbstract 定义但
  KmerCount4 不用；KmerTableSet 的 addKmersToTableAA 才有 prob 逻辑，非
  bbnorm 路径）。最初误实现 minProb 导致差异从 21 对膨胀到 35 对，已回退。
* **剩余差异纯判定**（精确表 vs bits=16 近似表的边界）：质量行 0 处
  差异（无 c 型 diff）；整 read 差异 min3=21 对、min5=27、min10=31、
  min20=37 对——近似表碰撞影响随 min 阈值增大而增多。外部桶路径与
  内存路径输出一致（--mem 1k 实测）。

**可选的收尾项（待确认）**：① 把 21 对差异正式定义为"精确语义 vs bbnorm
近似语义"；② `.pkt` count 字段按 bits 截断（对齐 bbnorm bits=16 → u16，
更激进可 u8）缩小落盘体积，判定字节不变（阈值在低端，截断到 65535
不影响任何判定）。

khmer（Count-Min Sketch + 在线 diginorm）源码分析见
`notes/references/khmer.md`。fairy（FracMinHash 稀疏采样 + 宏基因组
coverage）源码分析见 `notes/references/fairy.md`。

### 4.9 bbnorm kmer 深度分箱（2026-08-10 讨论记录，未定、暂不实现）

**结论**：bbnorm 原生支持按 kmer 深度把 reads 分箱（"Depth binning
parameters"，39.38 与 40.01 一致），但这是**新功能**，不在 anchr 替换
范围内（anchr 只用 highpass filter 一项），用户暂不打算做，先记录。

**参数与语义**（`bbnorm.sh` usage + `KmerNormalize.java` 主循环）：

* `lowbindepth=10`（lbd）/ `highbindepth=80`（hbd）+
  `outlow=<file>` / `outmid=<file>` / `outhigh=<file>`，一次运行 3 个 bin，
  更多层需多次运行不同阈值取交集。
* 判据是 `depthAL`——kmer 深度数组的稳健分位（约 46-54 百分位，只统计
  深度 ≥ `max(min, high/125)` 的 kmer），与 §4.8 的 toss 判定同一套量；
  help 里 "median" 是口语化说法。阈值控制的是 kmer 深度而非 read 深度。
* low：两条 mate 的 `depthAL` 都 < lbd；high：两条都落在 [lbd, hbd) 之外
  （源码 read1 用 `>HBD`、read2 用 `>=HBD`，小不对称，疑似笔误）；
  其余进 mid。双端按对分类。
* 分箱与 keep/toss 正交：任何 read 都会恰好进入一个 bin，可只分箱不
  归一化；两遍式（全量建表 → 分类），确定性，与已实现的 norm 语义一致。

**若将来要做**：`pgr fq norm` 的扩展——建表与 depthAL 计算全为现成，
加 `lbd/hbd` + 三个输出通道即可，成本不大。当前无动作。

### 4.10 reformat.sh 功能全景盘查（2026-08-10）

**结论**：anchr 实际用 reformat 做两件事——① `trim.era.sh` 的
`samplebasestarget`（已核对 ✅，§M2）；② **`2_insert_size.era.sh` 的
`ihist`**（SAM → 插入片段直方图，pgr 目前无对应，**这是真缺口**）。

**ihist 细节（39.38 实测）**：

* 调用：`reformat.sh in=<sam.gz> ihist=<file>`，输出格式：
  `#Mean/#Median/#Mode/#STDev/#PercentOfPairs` 五行 + `#InsertSize Count`
  分布表（Lambda bbmap sam 实测：Mean 421.7 / Median 407 / STDev 98.4）。
* 实现：`tracker/ReadStats.java` 的 `addToInsertHistogram(SamLine)`（取
  `|TLEN|`，`pairedOnSameChrom && x>0` 计入 paired，`MAXINSERTLEN` 截断）
  与 `insertSizeMapped(r1,r2)`（proper pair + 同染色体，cigar 长度换算，
  重叠对近似）；输出汇总在 `BBMerge.writeHistogram`。39.38/40.01 的
  `.class` 均含 `MAKE_IHIST`，两版本都支持。
* **重要认知修正**：BBTools 仓库里的 `.java` 与发布的 `.class` **不同步**
  ——40.01 的 `ReformatReads.java` 全文只有 1 处 "hist"（无任何 histogram
  实现），但 `ReformatReads.class` 里有 `MAKE_IHIST`/`addToInsertHistogram`
  （实测也真的输出）。后续源码分析以 `.class` strings + 黑盒实测为准，
  `.java` 只作提示性参考。
* 移植成本：pgr 没有 SAM reader，需先能读 sam（`|TLEN|` 或双线
  insertSizeMapped）+ 直方图输出；属于新命令（如 `pgr fq ihist` 或挂在
  `paf`/`fa` 侧）。**是否纳入迁移待用户定**（ihist 属于 2_insert_size
  模板，不在 trim 流水线 8 步内）。

**reformat 其余功能与 pgr 现状对照**（anchr 未用到）：

| reformat 功能 | pgr 现状 |
|---|---|
| 格式转换 fq↔fa、sam→fq、qual/scarf/oneline | `fq to-fa` ✅；sam→fq 无（pgr 无 SAM 命令） |
| 命名/序列操作（addslash/underscore/tuc/rcomp/uniquenames/remap/fixjunk/utot/pad） | `fa rc` ✅；fq 侧无专门命令 |
| 采样家族（samplerate/samplereadstarget/upsample/prioritizelength/reads/skipreads） | `fq sample` 只做 samplebasestarget（anchr 其他降采样用 hnsm split，不用 reformat） |
| 过滤家族（qtrim/trimq/minlen/maq/maxns/barcode/GC/forcetrim） | `fq trim-qual`/`fq trim-adapter` 覆盖（与 bbduk 交叉，已移植） |
| 直方图（bhist/qhist/lhist/gchist、sam 系列 ehist/idhist...） | `kmer hist`/`kmer qhist`/`fa n50`；fq lhist/gchist 无（anchr 未用） |
| sam/bam 过滤（mappedonly/mapq/flag/cigar） | pgr 无 SAM 命令（anchr 未用） |
| k/cardinality（loglog 唯一 kmer 数） | `kmer table` 精确计数 ✅；loglog 无（bbduk filter 的 cardinality 已随 trim-adapter 移植） |

**交叉关系**（BBTools 家族特点）：reformat 的 qtrim/minlen/maxns ↔ bbduk
（已移植）；cardinality ↔ bbduk filter（已移植）；ihist ↔ BBMerge/Tadpole
（共享参数）；sample 家族在 reformat 是主入口。**pgr 真正的缺口只有
ihist 一项**（anchr 用到且 pgr 无）。

**为什么 reformat 里有 ihist（设计动机，2026-08-10 反推）**：reformat 的
定位是"通用流式 read 处理器"而非纯格式转换器（官方 ReformatGuide 第一句：
"designed for generic streaming read-processing tasks... such as format
conversion, subsampling, and various filtering operations"）。ihist 在其中的
理由：① 反正要逐条扫 reads，直方图是循环内 O(1) 增量，零边际成本且保持
低内存；② 输入为 paired sam 时 TLEN/proper pair 信息现成，sam→fq 转换是
reformat 常见用途，顺路提取；③ 实现不在 ReformatReads 而在共享
`tracker.ReadStats`，谁消费 paired reads 谁挂上（BBMerge/Tadpole/BBWrap/
RQCFilter 同款 ihist）——"命令=输入流类型"的组织模式；④ 基础 histogram
与 bbduk 共享（官方明说，bbduk 更快更耗资源，reformat 低资源），高级变体
才单独出工具（如 readlength.sh）。这也是"工具交叉"现象的成因。

**insert size 的两种来源（2026-08-10 修正）**：reformat 的 ihist **必须
输入 sam（比对结果）**——取 `|TLEN|`/坐标差，未比对则无从得知，这是
"输入为比对流"的顺路统计，不是 reads 自带信息。未比对时只有另一条路：
BBMerge 在 overlap merge 时从 read 对重叠/缺口反推 `bestInsert`
（`hist[bestInsert]++`），可输出同格式 ihist，但仅对能 merge 的短插入
文库（insert < 2×read length）有效且为估算。同名 `ihist` 在两个工具里是
两条不同实现——"工具交叉"的另一面：**同名参数在不同工具里可能不是一回事**。

### 4.11 bbduk 参数全景盘查（2026-08-10）

**anchr 用 bbduk 做三件事**：

1. trim.era.sh trim（ktrim=r/k/mink/hdist/tbo/tpe/minlen/qtrim/trimq/ftm/
   maxns/stats/tossbrokenreads）——已移植并核对 ✅；
2. trim.era.sh filter（k/cardinality/stats/tossbrokenreads）——已移植 ✅；
3. **merge.era.sh 纯 qtrim**（`bbduk.sh qtrim=r trimq={{opt.qual}}
   minlen={{opt.len}}`，无 ref）——**缺口**：`pgr fq trim-adapter` 的
   `--ref` 目前必填；`pgr fq trim-qual` 是 sickle 语义（sliding/Mott），
   与 bbduk 的 optimalMode（testOptimal）输出不一致，不能顶替。

**结论：bbduk 全部参数里，对 anchr 有用的只剩 merge.era.sh 的纯 qtrim
一项**。已实现（2026-08-10）：`trim-adapter --ref` 改为可选——无 ref 时
跳过 kmer 匹配/ktrim/tbo/tpe，只做 qtrim/minlen/maxns/ftm/toss（qtrim
语义已与 39.38 逐字节核对）。实测 5 组无 ref 参数（含 merge 的 q25l60、
minlen=0、ftm=5 组合）与 `bbduk.sh qtrim=r trimq=... minlen=...` 逐字节
一致；新增回归测试 `no_ref_quality_trim_only`。merge.era.sh 替换时用
`--no-ktrim --no-tbo --no-tpe --maxns=-1 --ftm 0 --trimq <qual> --minlen <len>`。
注意：trim-adapter 输入语义为交错双端（1 文件）或 R1/R2（2 文件），
不支持单端多条（单端文件会按交错解析，第二条被当作 mate）。

**其余参数分类与用途标记**（anchr 均未用，pgr 不主动加）：

| 类别 | 参数 | 备注 |
|---|---|---|
| kmer 容错/判定 | maskmiddle、edist、qhdist、minkmerhits/fraction、mincovfraction、findbestmatch、forbidn | pgr 已实现 k/mink/hdist/rcomp |
| trim 模式 | ktrim=l、qtrim=rl/l/w、kmask、ksplit、ktrimtips、tp、ftl/ftr/ftr2、trimclip | anchr 只用 r 模式 |
| 过滤家族 | minavgquality、minbasequality、mcb、mingc/maxgc、tossjunk、swift、entropy、chastity/barcode/tag | 通用 QC，anchr 未用 |
| 质量处理 | quantize、recalibrate、mincalledquality/maxcalledquality | changequality 已按默认行为实现 |
| 输出流 | outm/outs/out2、refstats、rpkm、dump、rename、statscolumns=5 | anchr 只用 stats 3 列 |
| 直方图 | bhist/qhist/lhist/gchist/enthist/ihist 等 | anchr 的 khist 走 kmercountexact |
| sam/bam | sam、trimclip、varfile/vcf、ehist/idhist | pgr 无 SAM 命令 |
| 其他 | ecco（BBMerge 纠错）、amino、literal、samplerate、copyundefined | — |

## 5. 验收标准（替换后的对比）

- 同一输入（`anchr/tests/Lambda/` R1/R2），pgr 各步输出与本地 39.38 BBTools
  输出**解压后逐字节一致**（name、顺序、序列、质量、行宽/换行等格式细节；
  gz 压缩字节不要求一致；顺序由 M1 clumpify 保证）；
- `khist.txt`/`peaks.txt` 逐字节一致（已达成）；`trim.stats.txt`/
  `filter.stats.txt` 统计文本已复刻（`pgr fq trim-adapter --stats`，3 列
  格式与 39.38 逐字节一致，见 §6.5）；
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
5. 统计文本（2026-08-10 已实现）：`pgr fq trim-adapter --stats <file>` 输出
   bbduk `stats=` 3 列格式（`#File`/`#Total`/`#Matched`/`#Name` + 每参考序列
   行），排序 = StringCount（bases 降序、reads 降序、name 升序），每 read
   记首个命中 kmer 的 scaffold（ktrim 的 `id0` / countSetKmers 的
   maxBadKmers 命中），bases = 命中时 read 全长；Lambda trim/filter 两模式
   与 39.38 逐字节一致（`#File` 行是输入路径，天然路径相关）。**注意**：
   anchr `2_trim.tera.sh` 确实解析 `trim.stats.txt`/`filter.stats.txt`
   （保留 `#Matched`/`#Name` 行 + 第 3 列 >0.1/0.01 的数据行），早期
   "下游未解析"的判断有误——这也是本项必须做字节级一致的原因。
6. ~~golden 存放~~ → 已定：复制进 `tests/bbtools/Lambda/`（含 gz 化的
   golden 与 README）。

---

*参考来源: [trim.tera.sh](../../../anchr/templates/trim.tera.sh)（anchr 项目，只读） |
BBTools-40.01 源码 | [fq-trim-qual.md](fq-trim-qual.md) | [seq-reader.md](seq-reader.md)*
