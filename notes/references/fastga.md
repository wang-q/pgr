# FastGA 源码与论文分析

> 整理于 2026-08，源自对 `FASTGA-main/` 目录源码（约 4.6 万行 C）及 README 的通读。
> 目的：理解 FastGA 的快速全基因组比对算法（adaptive seeds + wave aligner + trace points），
> 为 pgr 的 pangenome 上游比对（verify-pangenome.sh 已用 `FastGA -psl/-pafx`）与对齐算法
> 提供参考。

## 1. FastGA 概览

- **工具定位**: 在两个高质量基因组之间（或一个基因组自比对）寻找全部局部 DNA 比对，
  默认输出 PAF，也可输出 PSL 或 ONEcode `.1aln` 格式。
- **作者/版本**: Gene Myers（daligner 作者）与 Chenxi Zhou；2023-05 首次发布，
  V1.5（2025-12-30）为当前版本。
- **核心假设**: 输入为近完整组装（至多几千 contig），序列质量 Q40+。
- **性能**: 2 Gbp 蝙蝠基因组 vs（8 核）约 5 分钟找到几乎所有 >100 bp、≥70% 相似区域；
  63.5 万个比对压缩到 44.5 MB `.1aln`。
- **算法来源**: adaptive seed（Martin Frith 的 adaptamer 思想）+ 首个 wave-based
  local aligner（源自 daligner 2012）；数据编码用 Gene Myers 的 ONEcode 框架。

**与 pgr 的关系**：pgr 的 pangenome 路线把 FastGA 当作上游比对器
（`FastGA -v -psl/-pafx A B` → `pgr pl chainnet --syn`），FastGA 负责产出 pairwise
alignment，pgr 负责 chain/net 精修与下游 PAF 图。FastGA 自己的 chaining 在 pgr 路线中
不使用（统一由 UCSC chain/net 承担，见 [[biser.md]] §6.8）。

## 2. 整体架构与数据流

```
FASTA/ONEcode
    │  FAtoGDB
    ▼
GDB (.1gdb 元数据 + 隐藏 .bps 2-bit 序列)     ← 随机访问、4 倍省 IO
    │  GIXmake
    ▼
GIX (.gix 稀疏 k-mer 索引, k=40 + (12,8) syncmer) ← 两个索引直接互查
    │  FastGA 主流程（种子扫描 → 链 → wave 对齐）
    ▼
.1aln（ONEcode trace point 编码，按 contig1→contig2→start 排序）
    │  ALNtoPAF / ALNtoPSL（线性时间）
    ▼
PAF / PSL
```

- 所有步骤可由 `FastGA` 一次触发；`FAtoGDB` / `GIXmake` / `ALNtoPAF` 等子进程允许
  分步控制，GDB/GIX 可持久化复用（-k 保留，多基因组重复比对时显著省时）。
- GIX 很大（每 Gbp 约 14 GB），建议批量比对前构建、之后用 `GIXrm` 清理，保留 GDB。
- 自比对模式（`FastGA A`）可检测基因组内部重复/单倍型间同源。
- **方向不对称**：adaptamer 依赖 source1 的种子，`FastGA A B` ≠ `FastGA B A`；
  `-S` 用两个基因组的 adaptamer 做对称（更慢，重复结构分析用；synteny 场景不建议）。
- **Soft mask**（V1.3+）：FASTA 小写=掩码，存入 GDB 的 `.1ano` 文件；默认忽略，
  `-M` 或 `#mask.1ano` 参数启用。

## 3. 核心算法

### 3.1 GDB：genome database（GDB.c）

- **两级结构**：scaffold → contig。`GDB_SCAFFOLD` 记录 scaffold 长度、首/末 contig、
  header 偏移；`GDB_CONTIG` 记录 contig 长度、scaffold 内起始、`.bps` 文件字节偏移。
- **序列存储**：2-bit 压缩（每个碱基 2 bit），存于隐藏文件 `.foo.bps`；元数据与序列分离，
  不需要序列的应用只读轻量 `.1gdb`。
- **N 处理**：FASTA 中 N 默认为 contig 间 gap；`-n` 指定阈值，短于阈值的 N 视为未知碱基
  （按 'a' 处理）。
- 派生自 daligner 的 GDB 代码。

### 3.2 GIX：syncmer 稀疏 k-mer 索引（GIXmake.c）

- 对每个 GDB 构建 k-mer 索引（`-k` 默认 40），但**不是全后缀数组**：只索引"以
  (12,8) syncmer 起始"的 40-mer（GIXmake.c: `TMER=12, SMER=8, SOFF=4`；
  `is_syncmer` 对 12-mer 的 8 个 s-mer 窗口取 canonical（正反链）最小值）。
- 每个索引条目 = 40-mer + 位置 + 掩码前缀信息；排序表存为 `-T` 个隐藏
  `.ktab.<int>` 分片（`.gix` 只是代理文件）。
- **体量**：README 实测约 14 GB / Gbp（人类 ~3 Gbp ≈ 42 GB）；但 **FastGA 默认在
  退出时自动删除自己创建的 GIX/GDB**（`Clean_Exit` 调 `GIXrm`，`-k` 才保留），
  运行时索引落在 `TMPDIR` 或 `-P` 指定目录。
- 关键设计：FastGA **直接比较两个 GIX**（两个排序的 k-mer 位置流线性归并找相同
  40-mer），而不是把一方的序列在另一方的索引中逐条查询——这是速度来源之一。
- 索引用 2-bit 编码 + canonical 方向，正反链统一。

### 3.3 Adaptive seeds（adaptamer，libfastk.c / FastGA.c）

- **定义**：位置 p 的 adaptive seed = 从 p 开始、在另一基因组中也出现的最长字符串。
- **频率过滤**：若该字符串在 source2 出现次数 > `-f`（默认 10），视为重复、不作为种子。
- **最小化**：`is_minimal`（libfastk.c:590）把种子与其反向互补做字典序比较，保留更小者
  （canonical 方向），正反链统一——与 pgr 的 canonical minimizer/syncmer 思路一致。
- **种子命中**：adaptamer 在 source2 的每个出现位置 (p, q) 都是一个 seed hit。
- `-S` 对称模式取两个基因组的 adaptamer 并集。

### 3.4 种子链（chaining，FastGA.c align_contigs）

种子 hit 按 **anti-diagonal（反对角线）空间** 排序扫描（`print_seeds` / `align_contigs`
中维护按 `(ipost, apost)` 归并的种子流）。合法链需满足：

1. 所有种子落在宽度 128 的对角线带内；
2. 相邻种子间距 < `CHAIN_BREAK`（源码 2000 = 2×`-s` 1000，anti-diagonal 空间 2 倍）；
3. 链在两侧覆盖 ≥ `CHAIN_MIN`（源码 170 = 2×`-c` 85）个 anti-diagonal。

满足条件的链在"tube"（`alow..ahgh` × `dgmin..dgmax`）内触发 wave aligner
（`Local_Alignment`）。self 比对时跳过完全相同的对角线段。

### 3.5 Wave-based local alignment（align.c）

源自 daligner 的 wave-front 对齐：

- **forward_wave / reverse_wave**（align.c:336 起）：沿对角线扩展 wave，维护
  - `V[k]`：对角线 k 的最远到达点（furthest reaching point）；
  - `M`：最近 TRIM_LEN 列的匹配数（位向量 1-bit 计数）；
  - `T`：隐含对齐最后列的位向量（用于轨迹）；
  - `Pebble` cells：wave cell 记录（ptr/diag/diff/mark）。
- **Local_Alignment**（align.c:1423）：在给定对角线带与 anti-diagonal 区间内做局部对齐，
  长度 ≥ `-l`（默认 100）、相似度 ≥ `-i`（默认 70%）。
- **Compute_Alignment**（align.c:5426）：divide-and-conquer trace（`dandc_nd` /
  `trace_nd` / `middle_np` / `iter_np`），用 sparse DP 在 wave 之间回溯完整比对路径，
  按 `tspace`（trace spacing）压缩轨迹。
- **Gap_Improver**（align.c:6714）：对 gap 区域做二次精修。

### 3.6 Trace points 与 .1aln 编码（alncode.c / ONEaln.c）

- 每条比对记录为**轨迹点（trace points）**：按 tspace 采样的比对路径编码，配合
  diff/长度信息，在 ONEcode 二进制中极紧凑（63.5 万比对 → 44.5 MB）。
- `.1aln` 头：`1 3 aln 2 1`、`!` 记录 FastGA 版本与参数、`<` 引用两个 GDB。
- **排序**：按 source1 contig # → source2 contig # → source1 start 排序，便于线性扫描。
- ALNtoPAF / ALNtoPSL（多线程）在**线性时间**把轨迹展开为 PAF/PSL（含 CIGAR：
  `-pafx` = `=`/`X`，`-pafm` = `M`；`-pafs/S` = CS 字符串）。
- ONEaln.c 提供 C 库读取 .1aln（依赖 GDB/ONElib/alncode/align/gene_core 一起编译）。

## 4. 源码模块结构

| 模块 | 职责 |
|------|------|
| `FastGA.c`（~4k 行）| 主流程、参数解析、种子扫描、anti-diagonal 链、调用 aligner |
| `align.c`（~6.7k 行）| wave aligner（forward/reverse_wave、Local_Alignment）、Compute_Alignment、trace、Gap_Improver |
| `libfastk.c` / `FastKS.c` | FastK 生态 k-mer 计数库：读写 GIX 的 `.ktab` 表（Histogram / Kmer_Table / Kmer_Stream / Profile_Index）|
| `GDB.c` / `GDB.h` | genome database：scaffold/contig 两级结构、2-bit 序列随机访问 |
| `GIXmake.c` | syncmer 稀疏 k-mer 索引构建（k=40 + (12,8) syncmer，含 mask 支持）|
| `ONElib.c` / `ONEaln.c` | ONEcode 数据编码框架、.1aln 读取 C 库 |
| `alncode.c` | trace point 编解码 |
| `ALNtoPAF.c` / `ALNtoPSL.c` | 轨迹 → PAF/PSL（多线程、含 CIGAR 生成）|
| `ALNchain.c`（Chenxi Zhou）| 按局部链过滤 .1aln 比对（-c/-s 阈值）|
| `ALNreset.c` | 重设 .1aln 对 GDB 的内部引用 |
| `select.c` | 基因组选择表达式解析（只比对选定的 contig/区间）|
| `PAFtoALN.c` / `PAFtoPSL.c` | 反向转换（PAF 带 X-CIGAR → .1aln/.psl）|
| `FAtoGDB.c` / `GDBtoFA.c` | FASTA/ONEcode ↔ GDB 互转 |
| `GDBshow`/`GDBstat`/`ALNshow`/`ALNplot`/`ANOshow`/`ANOstat` 等 | 查看/统计/绘图工具 |

> **libfastk 移植评估**：libfastk 是 FastK/GIX 私有格式（`.ktab.<int>` 分片）的访问库，
> 含 Histogram（k-mer 频率直方图）、Kmer_Table（加载排序表 + Fetch/Find）、
> Kmer_Stream（流式遍历 + GoTo 定位）、Profile_Index（raw reads profile）。pgr
> **不需要移植任何实现**：格式绑定 FastK 生态、pgr 不做 k-mer 计数/raw reads；
> 其中 `is_minimal`/`compress_norm`/`compress_comp`（2-bit 编码 + canonical 最小化）
> 与 pgr `nt.rs`/`syncmer.rs` 等价。仅 Kmer_Stream 的"排序流迭代 + 定位"接口形态
> 值得未来原生 `sd search --mode kmer` 的两流归并借鉴（Rust 自研实现）。

## 5. 关键参数（FastGA main 默认值，源码与 README 的对应）

| 参数 | 默认 | 含义 | 源码位置 |
|------|------|------|----------|
| `-f` | 10 | 最大种子频率（超过视为重复，不作为 adaptamer）| `FREQ = 10` |
| `-c` | 85 | 最小链覆盖 bp（源码 `CHAIN_MIN` 存 2×=170，anti-diagonal 空间）| `CHAIN_MIN = 170; <<= 1` |
| `-s` | 1000 | 相邻种子最大间距（源码 `CHAIN_BREAK` 存 2×=2000）| `CHAIN_BREAK = 2000` |
| `-l` | 100 | 最小局部比对长度 | `ALIGN_MIN = 100` |
| `-i` | 0.7 | 最小比对相似度（源码 `ALIGN_RATE = 1.-sim`，默认 .3；合法 [0.55,1)）| `ALIGN_RATE = .3` |
| `-k` | 40 | GIX k-mer 大小（GIXmake）| — |
| `-T` | 8 | 线程数 | `NTHREADS = 8` |
| `-S` | off | 对称 adaptamer（两个基因组种子）| flags |
| `-M` | off | 使用 GIX 中的 soft mask | flags |
| `-v` / `-L` | — | 详细模式 / 日志文件 | flags |

## 6. 输出格式

- **PAF**（默认）：12 列标准 PAF；`-pafx` 追加 `cg:Z:`（`=`/`X`/`I`/`D`），
  `-pafm` 用 `M`，`-pafs/S` 追加 CS 字符串。
- **PSL**（`-psl`）：UCSC PSL 格式，可直接喂给 `pgr pl chainnet` / `pgr psl chain`。
- **.1aln**（`-1:path`）：ONEcode 二进制（须指定输出文件），可用 ALNtoPAF/ALNtoPSL
  按需转换。

## 7. 对 pgr 的启示

1. **Adaptive seeds vs pgr 的 syncmer/minimizer**：FastGA 的 adaptamer 是"最长共享字符串"
  （长度自适应），pgr 的 closed syncmer（`src/libs/syncmer.rs`）是固定 k 的有界间隔采样。
  FastGA 的 `is_minimal` canonical 判断与 pgr 的 canonical rolling hash 同思路；
  若 pgr 未来做原生 k-mer search（如 [[biser.md]] 的 `pgr sd search --mode kmer`），
  adaptive seed 的"频率自适应 + 最长共享"是比固定 k-mer 更灵敏的候选。
2. **Wave aligner vs pgr 的 ScalarAlignmentEngine**：pgr 的 POA 对齐是标量 O(nm) 矩阵 DP；
  FastGA 的 wave-front（V/M/T 位向量 + Pebble cells）在线性空间内扩展，与 WFA 同族。
  若 pgr 需要更快的小片段 pairwise 对齐（如 pbit 的 CIGAR 精修、SD refine），
  wave-front 是明确的优化方向（远优于当前 O(nm)）。
3. **Trace point 编码 vs pgr 的 PAF/MAF 存储**：FastGA 用轨迹点 + ONEcode 压缩比对集合
  （63.5 万比对 → 44.5 MB），支持线性时间重放为任意格式。pgr 的 paf index 已做
  CIGAR 懒加载（BGZF vpos），但完整比对的紧凑存储可借鉴 trace point 思路。
4. **GDB 2-bit 序列库 vs pgr 的 twobit**：FastGA 把 2-bit 序列存隐藏文件 + 元数据分离，
  与 pgr `TwoBitFile` 类似；pgr 的 `.loc`/BGZF 随机访问已覆盖同等需求。
5. **对称性语义**：`FastGA A B` ≠ `FastGA B A`（adaptamer 不对称）对 pgr 有直接影响——
  pangenome 管线里 FastGA 的 query/target 顺序会影响找到的比对集合；
  verify-pangenome.sh 固定 `FastGA(b,a)` 方向后 chainnet 统一精修，顺序影响被下游
  chain/net 部分吸收（但重复区域仍可能不对称）。

## 8. 版本与许可

- 当前 FASTGA-main 对应 V1.5（2025-12-30），含 ONEcode ANO 文件支持。
- LICENSE：MIT（ALNchain 单独标注 Chenxi Zhou，MIT）。
- 参考：https://github.com/thegenemyers/FASTGA ；ONEcode:
  https://github.com/thegenemyers/ONEcode ；daligner:
  https://github.com/thegenemyers/DALIGNER

## 9. GDB 与 pgr 存储格式对比

### 9.1 对比对象

| 格式 | pgr 侧实现 | 定位 |
|------|-----------|------|
| FastGA GDB（`.1gdb` + `.bps`）| — | 组装基因组数据库：元数据与 2-bit 序列分离 |
| pgr 2bit（`pgr fa to-2bit`）| `src/libs/fmt/twobit.rs` | 标准 UCSC 2bit：单级 contig + 内嵌 mask/N block |
| pgr loc 索引 FASTA | `src/libs/loc.rs` | 原样 FASTA/BGZF + `.loc` 偏移索引（`fa range`、paf `FastaStore`）|
| pbit 参考层 | `src/libs/pbit/` | 复用 2bit 记录格式（`read_2bit_record` / `write_2bit_record`）|

### 9.2 序列编码与空间效率

- **GDB**：2-bit 压缩（`COMPRESSED_LEN = ceil(len/4)`），编码表 **A=0, C=1, G=2, T=3**
  （libfastk.c 的 `code[128]`），N 在 FASTA 阶段即按 gap 拆分（`-n` 阈值）。
- **pgr 2bit**：同为 `ceil(len/4)` 字节，但编码表是 **UCSC 标准 T=00, C=01, A=10, G=11**
  （twobit.rs:118）。N 保留为 n_blocks（长度不变）。
- **空间结论**：两者序列密度等价（0.25 B/bp）。pgr loc+FASTA 是 1 B/bp（4 倍），
  但保留原文、无转换成本。
- **编码表差异是硬约束**：GDB 的 packed bytes 与 pgr 2bit 的 packed bytes 不能互读，
  且负链互补映射不同（GDB `comp` 表 vs pgr 的 T↔A、C↔G）。若要互操作必须走文本/ASCII
  再重新编码，不能直接搬 packed 数据。

### 9.3 结构与元数据模型

| 维度 | GDB | pgr 2bit | pgr loc |
|------|-----|----------|---------|
| 层级 | scaffold → contig 两级（N 即 gap）| 单级 contig 平铺 | 单级 FASTA record |
| 元数据 | 与序列**分离**（轻量 `.1gdb` 不需 `.bps`，只读骨架免载序列）| 一体（记录头含 dna_size + block 表）| 一体（原文）|
| N/gap | N 拆分 contig（组装语义，scaffold 保留 N 的 gap 长度）| N 保留为 n_block（序列语义，长度不变）| 原样保留 |
| mask | **外置** `.1ano` 区间文件（可多个 mask union，改 mask 需重建 GIX）| 内嵌 mask_blocks（保留 soft-mask 语义）| 原文大小写 |

- **GDB 的 scaffold 语义**对组装输入（contig + gap 长度估计）更友好；pgr 2bit 面向
  "序列就是序列"的通用场景，N 当未知碱基。
- **GDB 元数据/序列分离**是实际优势：统计 scaffold 数、长度、名称时不用碰 `.bps`。
  pgr 2bit 单文件一体化，读元数据要扫整条记录。
- **mask 模型**：GDB 外置 `.1ano` 可 union 多个 mask 且不改序列本体（但改 mask 必须
  重建 GIX）；pgr 2bit 内嵌（一次打包、不可变）。pgr 的 `pgr fa mask`（runlist 硬/软
  mask，长度保留）接近"外置改字符"，但语义不同。

### 9.4 随机访问

- **GDB**：`Get_Contig_Piece` 按 `boff + beg/4` 做 `fseeko`，读 `ceil((end-beg)/4)` 字节后
  `Uncompress_Read`；也支持整库 mmap（`seqstate != EXTERNAL` 时 `Get_Contig` 直接
  memcpy）。区间读取是 **O(区间长)**。
- **pgr 2bit**：`read_2bit_record` 同样 seek 到 `packed_dna_start + first_byte_idx`，
  只读区间字节后解压，**O(区间长)**——与 GDB 等价。
- **pgr loc**：`fetch_range_seq` 先 `fetch_record` **读整条 record** 再切片
  （除非 BGZF 虚拟位置做细分）。对长 contig 的短区间访问，loc 明显低效；
  这是 pgr 内部"区间提取优先用 2bit"（chain ScoreContext、pbit）的原因。
- **构建成本**：三者都是单遍扫描 FASTA；loc 最便宜（只记偏移），GDB/2bit 需编码。

### 9.5 输出形态与工具生态

- **GDB**：`Get_Contig` 支持 COMPRESSED / NUMERIC(0-4) / LOWER_CASE / UPPER_CASE 四种
  形态按需转换；ONEcode 生态（.1seq 输入、.1aln 输出）；`GDBtoFA` 可逆回 FASTA。
- **pgr 2bit**：`read_sequence(no_mask)` 控制大小写；标准 UCSC 2bit 可被 kent-tools /
  其他工具直接读取（pgr 与 UCSC 字节级一致，见 [[ucsc.md]]）；pbit 参考层复用同一记录格式。
- **pgr loc**：底层是标准 FASTA/BGZF，任何工具可读原文。
- **标准性**：pgr 2bit 互操作最强（UCSC 标准）；GDB 是 FastGA 生态私有格式，外部工具
  无法直接消费，属于生态锁定。

### 9.6 结论与建议

- **核心等价**：GDB 与 pgr 2bit 在"2-bit 压缩 + 字节偏移随机访问"上是同构设计，
  空间与区间读取性能相当。真正的差异在语义层（scaffold 两级 vs 平铺、mask 外置 vs
  内嵌、元数据分离 vs 一体）和互操作（私有 vs UCSC 标准）。
- **GDB 值得借鉴的两点**：
  1. **元数据/序列分离**：pgr 若需要"只看骨架不读序列"的场景（如大量基因组的名/长
     统计、`fa size` 的快速版），可参考 `.1gdb` 轻量件设计；
  2. **mask 外置（.1ano）**：pgr 的 mask 内嵌 2bit 后不可变，若要支持"同一参考不同
     mask 重复比对"，外置区间文件 + 读取时过滤（类似 GDB 的 mask 参数）更灵活。
- **不建议**：为互操作而模仿 GDB 的编码表/ONEcode——pgr 2bit 已是 UCSC 标准，
  与 kent-tools 字节级兼容是现有资产（ucsc pipeline 验证依赖它），不应为 FastGA
  私有格式放弃。
- **pgr loc 的定位**：适合"保留原文 + 便宜索引"（fa range、FastaStore），区间密集
  随机访问场景应继续用 2bit。

## 10. GIX 分析：好处与 pgr 借鉴评估

### 10.1 GIX 是什么

GIX（GIXmake.c）是每个基因组的**syncmer 稀疏 k-mer 索引**：只取"以 (12,8) canonical
syncmer 起始"的 k=40 的 k-mer（2-bit 压缩 10 字节），按字典序桶排序（首字节 1024 桶 +
`Ksplit` 均衡分割 + 多线程），排序后每个条目附位置信息（contig #、contig 内偏移、
方向位）与 lcp。README 实测体量约 **14 GB / Gbp**（`.gix` 代理 + `-T` 个 `.ktab.<int>` 隐藏
分片），但 **FastGA 默认退出时自动删除**（`Clean_Exit` → `GIXrm`，`-k` 才保留）；
构建/运行的临时分片在 `TMPDIR` 或 `-P` 目录。

> **与 pgr 的直接关联**：GIX 的 (12,8) syncmer 就是 [[syng.md]]/pgr 已实现的 closed
> syncmer 同族采样——FastGA 用 syncmer 稀疏化 40-mer 索引（密度约 2/(w+1)，大幅低于
> 全后缀数组），这正是 §10.4 建议 pgr"用 syncmer 稀疏化借鉴"的依据：**FastGA 自己
> 就是这么做的**，pgr 不必重新发明。

> **为何用户"排人类基因组没感觉生成多大的数据"**：FastGA 默认（不加 `-k`）在
> 结束时 `GIXrm` 清理自己创建的 GDB/GIX，只留下输出 PAF/PSL；GIX 分片只存在于
> 运行期（`TMPDIR`/`-P`，人类 ~42 GB），结束后即被删除。只有显式 `-k` 或预建
> GIX（`GIXmake`）才会留下这几十 GB。

### 10.2 GIX 的核心好处

1. **两个索引线性归并找同源（最重要）**：FastGA 不把 A 的序列在 B 的索引里逐条查询，
   而是把两个**已排序的 k-mer 位置流做一次归并**（FastGA.c 的 PAIR 文件流），相同
   40-mer 的两侧位置对 (i,j) 在一次 O(|A|+|B|) 线性扫描中全部发现。逐查询是
   O(|A|·log|B|)，归并消除了常数级和随机 IO 开销——这是 2 Gbp vs 5 分钟的主要来源。
2. **lcp 连续传播 → 种子长度自适应**：排序流中相邻相同 k-mer 的 lcp（最长公共前缀）
   直接给出共享长度，40-mer 命中自然扩展为任意长度的最长共享字符串（adaptamer），
   无需对不同 k 反复查询，且天然支持"频率过滤后最长种子"的语义。
3. **anti-diagonal 坐标（在种子流中计算，不在 GIX 里）**：`.ktab` 条目只存
   contig + 位置；两索引归并产种子命中 (i,j) 时**即时**算 `diag = i−j`、
   `anti = i+j` 写入种子流条目，链扫描（`align_contigs` 的 tube 逻辑）直接使用。
   "预编码"是归并阶段的 O(1) 实现选择，不是 GIX 的数据内容。
4. **流式 + 定宽条目**：`Post_List` 按块流式读入（POST_BLOCK），固定宽度条目
   （swide）顺序扫描，内存驻留可控。
5. **桶排序近线性构建**：MSD radix 风格（首字节 1024 桶 + Ksplit 负载均衡），
   多线程并行，构建接近线性；GIX 持久化后多轮比对复用（-k）。

### 10.3 代价与限制

- **空间巨大**：14 GB/Gbp（40-mer 表 + 位置 + 桶），细菌规模（5 Mb）约 70 MB 可接受，
  但大规模集合不可行。
- k=40 固定（可调但需 ≥12 且被 4 整除）；依赖完全匹配种子，靠 lcp 扩展与频率过滤
  （-f 10）控制重复区域。
- 私有格式，构建/读取绑定 FastGA 生态。

### 10.4 pgr 是否需要借鉴

**当前不需要**，理由：

1. **路线不符**：pgr 的 pangenome 路线明确"复用 pairwise 资产、不重新做比对"
  （[[paf-pangenome.md]] §1），de novo 同源检测由 FastGA/lastz 等外部比对器承担；
  pgr 的索引（paf index）消费的是**已比对的 PAF**，不是序列索引。
2. **规模不匹配**：GIX 14 GB/Gbp 的代价对 4 万大肠杆菌（总 ~200 Gbp）完全不可行；
  pgr 现有的 Mash/syncmer sketch 才是这一规模下的"近似同源过滤"正确工具。

**未来值得借鉴的三个子项**（若 pgr 需要原生序列同源检测）：

1. **两排序 k-mer 流归并找命中**——如果未来实现 [[biser.md]] 的 `pgr sd search
   --mode kmer` 或自比对找重复（SD/重复家族），"两个排序流归并"比逐查询更优；
   FastGA 自己就用 (12,8) syncmer 把索引稀疏化到 14 GB/Gbp；pgr 可沿用这一
   模式，用 **closed syncmer**（已落地于 `src/libs/syncmer.rs`）或更大的 w 进一步
   降密度，代价是灵敏度（syncmer 只保采样点）。
2. **lcp 连续传播 → 种子长度自适应**——固定 k-mer 的"最短种子"思路（如 minimap2）
   会漏掉短于 k 的同源；lcp 扩展天然得到"该位置的最长共享字符串"，灵敏度更高。
3. **anti-diagonal 坐标变换 + 桶排序 + 流式**——种子命中 (i,j) 即时算 diag/anti 进
   种子流、MSD 桶排序、定宽流式扫描，都是实现高效种子检测的成熟工程模式，可整体
   迁移到 Rust（`libs/sd/kmer_index.rs` + `plane_sweep.rs` 蓝图）。

**一句话结论**：GIX 的"归并 + lcp + anti-diagonal 坐标"是高效的序列索引设计，但 pgr 的
架构（外部比对 + PAF 图）和规模（4 万细菌）决定了它**现在不引入**；等有原生
同源检测需求时，以 syncmer 稀疏化 + 两流归并的形式借鉴其思想，而不是照搬
14 GB/Gbp 的截断后缀数组。

## 11. 从 GIX 到 Wave align 的完整算法管线

> 本节把 FastGA 从索引到比对的完整数据流串起来，标注源码位置与关键参数，
> 便于整体理解或移植其中的算法模式。

```
FASTA
  │ 1. GDB 构建（FAtoGDB / GDB.c）
  ▼
GDB（.1gdb 元数据 + .bps 2-bit 序列，scaffold→contig 两级）
  │ 2. GIX 构建（GIXmake.c）
  ▼
GIX（.gix 代理 + N×.ktab.<int> 分片 = (12,8) syncmer 起始的 40-mer 排序表）
  │ 3. 归并找种子（FastGA.c new_merge_thread）
  ▼
PAIR 流（种子位置对：ipost/icont/jpost/jcont/lcp，按前缀面板归并）
  │ 4. 链扫描（FastGA.c align_contigs）
  ▼
"tube" 命中（对角线带 × anti-diagonal 区间，含链覆盖判定）
  │ 5. Wave 局部对齐（align.c Local_Alignment / forward_wave / reverse_wave）
  ▼
比对路径（start/end + diff + trace cells）
  │ 6. Trace 回溯（align.c Compute_Alignment / dandc_nd / trace_nd）
  ▼
trace points → .1aln（ONEcode，按 contig1→contig2→start 排序）
  │ 7. 格式输出（ALNtoPAF.c / ALNtoPSL.c，线性展开）
  ▼
PAF / PSL（含 CIGAR）
```

### 11.1 步骤 1-2：GDB + GIX（离线，可持久化复用）

- **GDB**（GDB.c）：FASTA → 2-bit 压缩（`COMPRESSED_LEN = ceil(len/4)`），
  scaffold/contig 两级，N 按 `-n` 阈值拆 gap；元数据（.1gdb）与序列（.bps）分离，
  可 mmap 整库。
- **GIX**（GIXmake.c）：只索引"以 (12,8) canonical syncmer 起始"的 k=40 k-mer
  （`TMER=12, SMER=8, SOFF=4`；`is_syncmer` 对 12-mer 的 8 个 s-mer 窗口取
  正反链最小值）。排序：首字节 1024 桶 → `Ksplit` 均衡分片 → 多线程桶排序，
  输出 `-T` 个 `.ktab.<int>` 分片（Kmer_Stream 流式读取，libfastk.c）。
- 参数：`-k 40`（k 大小）、`-f 10`（种子频率阈值，GIX 构建时同时算频率）。

### 11.2 步骤 3：归并找种子（两个排序 k-mer 流的 join）

`new_merge_thread`（FastGA.c:610，逐前缀面板并行）：

1. 两个 GIX 各自以 `Kmer_Stream` 流式迭代（字典序）；线程各处理 256 个前缀
   （`[pbeg<<8, pend<<8)`），按 `Kmer_Stream.index` 定位起点。
2. 对 T1 的每个前缀 `cpre`：T2 跳过前缀更小的条目，把前缀 == cpre 的 T2 条目
   载入小缓存（`cache`），然后做**前缀面板内的归并**。
3. 面板内相同的 40-mer（T1 条目 vs T2 缓存）产出种子位置对，写 PAIR 流；
   **频率过滤**：某 40-mer 出现次数 > `kfreq = FREQ × kbyte`（-f 10）则跳过。
4. 种子条目 = (ipost, icont, jpost, jcont, lcp)：两侧位置 + contig + 与下一个
   相同 k-mer 的**最长公共前缀**（lcp 连续传播 → adaptamer 可扩展到任意长度）。

复杂度：O(|A|+|B|) 线性归并（每个 k-mer 恰好处理一次），而非逐查询。

### 11.3 步骤 4：链扫描（anti-diagonal 空间）

`align_contigs`（FastGA.c:2973）把 PAIR 流转成 anti-diagonal 坐标
（`diag = i-j`、`anti = i+j`），按对角线桶（`diag >> BUCK_SHIFT`，桶宽 64）
组织后扫描：

1. 按对角线分三段（b/m/e）归并相邻对角线，保证链不跨对角线带。
2. 维护链的 `alow..ahgh`（anti 区间）与 `dgmin..dgmax`（对角线区间）"tube"；
   种子间距超过 `CHAIN_BREAK`（2000，=2×-s 1000）时结束当前链。
3. 链的 anti 覆盖 ≥ `CHAIN_MIN`（170，=2×-c 85）则触发"tube"处理
   （否则丢弃）。self 比对时跳过完全相同的对角线。
4. tube 处理：加载两个 contig 序列，按 `amid = alow + BUCK_ANTI` 分块调用
   `Local_Alignment`，每次对齐一个 anti-diagonal 子区间。

### 11.4 步骤 5：Wave 局部对齐（Myers wavefront）

`Local_Alignment`（align.c:1423）：

1. 分配 wave 数组（V/M/HA/NA/T，5×vlen）与 trace cells 空间。
2. `forward_wave` 从 mid-line 正向扩展；`reverse_wave` 从低端反向扩展；
   自比对时用 `minp=1/maxp=-1` 防止与自身完全重合。
3. `fshort`/`rshort`：若正向或反向扩展太短（< `DUB_TRIM`），调整边界后
   只重跑短的一侧。

`forward_wave`（align.c:336）核心：

- **0-wave 初始化**：对每个对角线 k，从 `x=(mida+k)/2` 做 snake（同向延伸
  匹配），`V[k]`=最远点、`T[k]`=轨迹位、`M[k]`=匹配数；每 `tspace` 个匹配
  生成一个 `Pebble` cell（ptr/diag/diff/mark）。
- **wave 推进**（`while more && lasta >= besta - TRIM_MLAG`）：
  - 每轮 `dif += 1`，对角线带 ±1 扩展；
  - 每个 k 按 Myers 三分支更新波前（`ac`/`am`/`ap` 三条候选取 max，
    mismatch +1、双 gap +2），再 snake 延伸并更新 M/T；
  - `TRIM_MLAG` 提前终止：最优波前推进超过一定滞后即停止。
- 阈值：长度 ≥ `ALIGN_MIN`（-l 100）、相似度 ≥ `1-ALIGN_RATE`（-i 70%）。

### 11.5 步骤 6：Trace 回溯与编码

`Compute_Alignment`（align.c:5426）根据任务类型组装：

- **DIFF_ONLY**：`split_nd` 只算差异数（用于种子阶段评估）。
- **PLUS_ALIGN**：`dandc_nd`（Hirschberg 风格分治）——`split_nd` 找中间点，
  递归左右，D==1 时输出单个 I/D/S 操作，得到完整路径。
- **DIFF_TRACE / PLUS_TRACE**：`trace_nd` 在中间点按 `tspace` 采样，把比对
  压缩为 trace points（每点记录到下一 trace point 的 diff 与坐标增量）。
- `Gap_Improver`（align.c:6714）对 gap 区域二次精修。

结果按 contig1 → contig2 → start 排序写入 `.1aln`（ONEcode 编码，
alncode.c）；ALNtoPAF/ALNtoPSL 多线程线性展开 trace → CIGAR（`-pafx` 的
`=`/`X` 或 `-pafm` 的 `M`）。

### 11.6 关键设计点总结

| 阶段 | 核心技巧 | 复杂度 |
|------|----------|--------|
| GIX | (12,8) syncmer 稀疏 + 首字节桶排序 | 近线性构建 |
| 归并 | 两排序流前缀面板 join + lcp 传播 | O(|A|+|B|) |
| 链 | anti-diagonal 坐标 + tube 扫描 | 线性（链稀疏）|
| Wave | Myers wavefront（V/M/T + Pebble cells）| 与差异数成正比，优于 O(nm) |
| Trace | Hirschberg 分治 + tspace 采样 | 线性于路径长 |

## 12. 简化移植方案与代码量估算

> 场景：pgr 原生实现"简化 FastGA"——输入基因组序列 → 输出 PSL 块（不做
> chaining），接入 `pgr pl chainnet` 做 UCSC 链化。以下估算基于 Rust 实现 +
> pgr 现有资产复用（nt.rs 2-bit 编码、syncmer.rs canonical hash、fmt/psl.rs
> PSL 写入、ScalarAlignmentEngine 局部比对）。

### 12.1 边界

只移植 §11 的"种子发现 + 局部扩展"段，输出成块 PSL；FastGA 的种子链（tube）
简化为"独立扩展 + 重叠块合并"（或直接去掉，由 chainnet 的块链化兜底）。

### 12.2 两档估算（行数，含测试）

| 组件 | 极简版（细菌可用）| 完整版（接近 FastGA）|
|------|------------------|---------------------|
| k-mer 索引（2-bit + canonical）| ~200（HashMap，5 Mb 约 70 MB）| ~400（(12,8) syncmer 稀疏 + 桶排序）|
| 种子发现（频率过滤 + lcp 扩展）| ~100 | ~250（两流归并 / adaptamer）|
| 种子链 / 归组 | 可去掉（扩展 + 重叠合并）| ~300（anti-diagonal tube）|
| 局部扩展 | ~250（复用 `ScalarAlignmentEngine::Local`）| ~600（Myers wavefront：V/M/T + snake + 三分支 + trace cells）|
| PSL 输出 | ~50（复用 `Psl::from_align`）| ~200（trace → PSL）|
| 命令层 + 帮助 + 测试 | ~350 | ~500 |
| **合计** | **~950** | **~2250** |

复用 pgr 资产可省约 **400-600 行**（2-bit / canonical / PSL 写入 / 局部 DP）。

### 12.3 三个关键决策

1. **扩展器**：wavefront 是最难写的 ~600 行（V/M/T 位向量 + snake + 三分支更新 +
   trace cells，调试成本高）。务实路径：先用 `ScalarAlignmentEngine::Local`
   （O(nm)）验证正确性，再按 §11.4 把 wavefront 作为性能优化加入。
2. **索引**：细菌级 HashMap 足够（~5 M 条目）；人类级才需要 syncmer 稀疏 + 流式
   归并。第一版直接 HashMap，"稀疏化"留作未来开关。
3. **chaining 交给 pgr**：简化版对种子直接扩展 + 合并重叠块，PSL 交给
   `pgr pl chainnet`（非 --syn 或 --syn）做 UCSC 链化——这正是 pgr 已字节级验证
   的主场（[[ucsc.md]]）。

### 12.4 建议路线

1. 极简版（~1000 行）：k-mer 种子 + ScalarAlignmentEngine 扩展 + PSL 输出，
   E. coli 端到端验证（chainnet 接手链化）。
2. 性能版（+~1200 行）：syncmer 稀疏索引 + wavefront 扩展器，替换步骤 1 的两个
   热点，按需启用。

### 12.5 索引选型结论

- **第一版不需要 GIX 式索引**：pgr 是单基因组逐个比对，E. coli 5 Mb 的全 40-mer
  HashMap 索引约 100-200 MB，逐查询种子发现 ~0.3 s；GIX 的两流归并/流式分片是
  人类规模（14 GB/Gbp）才需要的工程，细菌规模收益可忽略，却要付几百行复杂度。
- **若要优化，做"syncmer 稀疏 + 排序数组"**：closed syncmer 密度 ~2/(w+1) 把条目
  降到 ~1/4（5 Mb 约 40 MB），`Vec<(u64,u32)>` 排序后二分替代 HashMap（省内存、
  cache 友好），排序流天然支持相邻条目 lcp（种子扩展/归并所需）。
- **完整 GIX（桶排序 + .ktab 分片 + Kmer_Stream 流式）不移植**：那是 14 GB/Gbp
  规模的工程复杂度，见 §10.4。
