# pgr 迁移设计：`.1aln`（ONEcode 轨迹点）读写

> 定位：把 FastGA 的 `.1aln`（ONEcode 二进制轨迹点比对存储）读写迁入 pgr 的设计与
> 实施参考。动机：`.1aln` 文件比展开后的 PAF/PSL 小很多（63.5 万比对 → 44.5 MB）。
> 引入它是为了**后续 40k E. coli 计划**：大规模全基因组比对会产生巨量 PAF，用紧凑的
> `.1aln` 落盘能**有效减少总文件大小**。
>
> 状态：**2026-08-05 起设计**，本文档为实施参考。用户规划确定为**"读 + 写"**
> （§3 第三档，~2700–3800 行），见 §6。
>
> 结构：§0 结论 → §1 格式是什么 → §2 迁移拆解与工作量 → §3 三个范围档位 →
> §4 关键注意点 → §5 与既有 pgr 资产/决策的关系 → §6 决策与写侧要点 →
> §7 实施细节（精读参考代码后补全）。

## 0. 结论

按 **读 + 写** 规划（§3 第三档），迁移 `.1aln` 到 pgr 约需 **~2700–3800 行 Rust**
（读侧 1900–2600 + 写侧 800–1200）。
最难的比对回溯部分（`dandc_nd`/`split_nd`/wave）pgr 的 `[[pgi-align.md]]`
wave.rs 已移植，主要新工作是 ONEcode 二进制容器 + 整数 codec + schema/GDB
骨架解析，以及写侧的轨迹点重采样（`cigar2tp`，见 §6.1.1）。

核心权衡：`.1aln` 是 FastGA **私有格式**。读侧展开成 PSL/PAF 必需源序列；写侧可由
带路径的比对（PAF-with-CIGAR / MAF / wave `EditOp`）经 `cigar2tp` 生成，无需源序列。
紧凑存储在 pgr 内部已有 `[[pbit.md]]`（CIGAR delta）可替代——所以 `.1aln` 更偏向
"与 FastGA 生态互操作"（读其输出、写标准 `.1aln` 供其消费），而非"pgr 自己的紧凑
格式"。

## 1. 格式是什么

`.1aln` 是 FastGA 的二进制比对存储格式，套在 **ONEcode 容器**里。每条比对记录是
"轨迹点（trace points）"：按 `tspace`（默认 100）采样的比对路径 + 每段差异数，
而不是逐碱基的 CIGAR。这正是它比 PAF 小很多的原因——细菌级别 PAF 展开后很大，
但轨迹点 + 差分整数编码能压到非常紧凑。

读取它解压成 PAF/PSL 需要 6 个独立组件（依赖 `ONEaln.c` + `ONElib.c`/`alncode.c`/
`GDB.c`/`align.c`）。

## 2. 迁移拆解与工作量

| # | 组件 | 来源 | 说明 | 估算 Rust 行数 |
|---|---|---|---|---|
| A | **ONEcode 二进制容器读取** | `ONElib.c`（~4900 行） | 头部 schema/provenance/reference 解析、二进制行码表、标量与列表字段读取、footer/对象索引。只需读路径 | ~550–700 |
| B | **整数 codec（vc）+ DNA codec** | `ONElib.c` `vcCreateCodec/vcEncode/vcDecode` | 轨迹点 T/X 行是 `INT_LIST`，用"字节差分 + Huffman"压缩；DNA 字段用 2bit。这是文件小的核心 | ~270–350 |
| C | **`.1aln` schema + 记录读取** | `alncode.c` + `ONEaln.c` 前半 | schema 文本、`A/D/R/T/X/U` 行→`AlnRecord`、GDB skeleton（scaffold/contig） | ~300–400 |
| D | **序列访问** | `GDB.c Get_GDB` | FastGA 从 `.1gdb/.bps` 读源序列；pgr 直接复用现有 2bit/loc 或 `--ref-seq/--query-seq` 模式（同 `align pgi`） | ~100–200 |
| E | **轨迹展开 → CIGAR/PSL/PAF** | `align.c Compute_Trace_PTS`/`iter_np`/`cigar_core` | **pgr 已移植大半**：`wave.rs` 已有 `split_nd`/`dandc_nd`/`EditOp`/`local_alignment`，且**轨迹点间盒内 DP 已由 `dandc_nd`/`banded_edit_ops` 覆盖**（§7.10）。真正要新写的是：轨迹点解码 → 还原比对路径的编排 + CIGAR/indel 生成 + PSL/PAF 输出 | ~400–500（另有数百行直接复用） |
| F | **CLI + 集成测试** | — | 新命令 `pgr 1aln`（§6.4），输出 PSL/PAF/stats | ~250–400 |

**净新增合计（读 + 展开成 PSL/PAF）：约 1900–2600 行 Rust**，其中 E 复用现有
`wave.rs`，是最大的一块已经省掉。

## 3. 三个范围档位

| 档位 | 内容 | 估算行数 |
|---|---|---|
| 最小（只读 + 统计） | 只解析头部、记录、spans/diffs，不碰序列 | ~800–1200 |
| 完整（读 + 展开成 PSL/PAF） | 上表全部 | ~1900–2600 |
| **读 + 写 `.1aln`（紧凑再输出）【已定】** | 读侧（§2 全表，1900–2600）+ 写侧（§6.2，800–1200） | 总计 ~2700–3800 |

> **2026-08-05 用户规划定为"读 + 写"**——既能消费 FastGA 的 `.1aln`，也能把
> pgr 的比对结果写成紧凑 `.1aln`（pgr align pgi 输出端可改用 `.1aln` 落盘）。

## 4. 关键注意点

1. **展开必须依赖源序列**：`.1aln` 头部只存源基因组文件名引用，要把轨迹展开成
   PSL/PAF 必须有原始序列（就像 `pgr align pgi` 的 `--ref-seq/--query-seq`）。
   只做统计则不需要。
2. **Gap_Improver 不用移植**：与 [[pgi-align.md]] §3.1 结论一致——它只对 tspace
   采样轨迹补全，而 pgr 的 wave 是精确回溯。展开端**轨迹点间的盒内 DP 已由
   `wave.rs` 覆盖**（`dandc_nd`/`banded_edit_ops`，§7.10），CIGAR 生成也有
   `cigar_from_alignment`；真正要新写的是把轨迹点解码、还原路径的编排 + PSL/PAF 输出。
3. **格式是 FastGA 私有**：ONEcode + GDB 编码表与 pgr 的 UCSC 2bit 不互通
   （[[fastga.md]] §9.2 已记录）——pgr 读不了 FastGA 的 `.bps`，FastGA 也读不了
   pgr 的 2bit。但若 pgr 写出**标准 `.1aln`**，FastGA 生态工具（`ALNchain`/
   `ALNtoPSL`/`ALNshow`/`ALNplot`）能直接消费（见 §6.3 互操作价值）。

## 5. 与既有 pgr 资产/决策的关系

- **wave.rs 已覆盖对齐回溯**：`split_nd`/`dandc_nd`/`EditOp`/`local_alignment`
  是展开端最难的算法部分，已移植，避免了重复实现。
- **设计稿 §0.4 曾把 `.1aln` 定为"不做"**，理由是"人类规模才需要"。本次因
  **存储体积**重新评估，动机已变。
- **pbit 是 pgr 自己的紧凑格式**：若目标只是"紧凑存储"，`[[pbit.md]]`（CIGAR
  delta 压缩）已存在，更贴合"pgr 内部紧凑格式"，而非互操作 FastGA 私有格式。
- **`.paf.idx` 是 pgr 已有的二进制 PAF 先例**（`pgr paf index`，`PGRI`+bincode，
  [persist.rs](../../src/libs/paf/persist.rs)）：它是**查询索引**（target 区间树 +
  双向 BFS），不是紧凑归档。CIGAR 逐 op 位打包（4B/op），比文本 PAF 小，但每条
  op 4 字节、达不到 `.1aln` 的轨迹点级压缩，且带区间树开销。它与 `.1aln` **正交**：
  前者服务 `paf query` 复用，后者服务 40k E. coli 的落盘体积。其 `Lazy` 模式
  （vpos 指向 BGZF 源文件）代表"引用源文件、不展开 CIGAR"的思路，与 `.1aln` 的
  "重采样"是两种不同的紧凑策略。

## 6. 决策与写侧要点

### 6.1 决策（2026-08-05）

**目标范围 = 读 + 写 `.1aln`**（§3 第三档，~2700–3800 行）。读侧即 §2 全表
（1900–2600）；写侧（§6.2，800–1200）在读侧基础上新增 ONEcode 容器写（A'）、
codec 编码（B'）、轨迹点重采样（G）、写 CLI 集成（H）。

### 6.1.1 写侧接入点修正（PAFtoALN 确认）

`PAFtoALN.c` 证明：**从带碱基级路径的比对格式 → `.1aln` 是可行的**，且无需像早先
假设的那样"挂在 wave 回溯（`EditOp`）里"。PAFtoALN 从带 `cg:Z`
X-CIGAR 的 PAF 出发，用 `cigar2tp` 把 CIGAR 重采样成轨迹点 + 差异，再
`Write_Aln_Overlap`/`Write_Aln_Trace` 落盘。

**关键判据：源格式是否保留碱基级路径。**

| 源格式 | 带碱基级路径 | 能否 → `.1aln` |
|---|---|---|
| PAF 带 `cg:Z` X-CIGAR | ✅ | ✅ 即 PAFtoALN |
| 两序列 MAF（`s` 行带显式 gap） | ✅ | ✅ 逐位推 X-CIGAR 后同路 |
| pgr `align pgi` 内部 wave 的 `EditOp` | ✅ | ✅ |
| PSL（只存坐标+计数） | ❌ | ❌ |
| chain / net（更压缩、丢非 syntenic） | ❌ | ❌ |

**因此写 `.1aln` 的正确接入点是"对齐输出端"**（保留 PAF-with-CIGAR / MAF /
`EditOp` 任一路径），而不是 chainnet 之后（chain/net 已丢路径，无法反推）。
推荐做成**独立转换器**：`align pgi ──► PAF(with cg:Z) ──► pgr aln-write
──► .1aln`，解耦自对齐器；`cigar2tp`/`cigarCheck`/`Write_Aln_Overlap/Trace`
即写侧组件 B'/A'/G 的完整参考。

**MAF 写 `.1aln` 的骨架说明（contig 级 MAF）**：

我们的 MAF 是 **contig 级**（无 scaffold 层）——每个 `src` 就是顶层 source
sequence（contig），`s` 行的 `srcSize` 即该 contig 全长。因此：

- **骨架信息 MAF 自带**：`srcSize`（contig 全长）+ 每个块的 `start/size/strand`
  足以直接构建 `.1aln` 头部 GDB skeleton（`Write_Skeleton`），每个 contig 一个
  顶层序列 + 全长，**无需额外源基因组**（`.1gdb` 或 fasta）。
- **无需 scaffold→contig 分组**：不要这层，故也不必从 `i` 行 contiguity 推断
  上级结构。
- **路径层无损**：MAF 的 gap 已在 `s` 行显式，逐位推出 X-CIGAR（`=`/`X`/`I`/`D`）
  后走 `cigar2tp`，与 PAFtoALN 同路。
- **唯一注意**：MAF 常是过滤后的子集（syntenic/分层），转出的 `.1aln` 通常不是
  FastGA 原始全量比对集——数据选择问题，非格式可转换性问题。

即：**contig 级原始比对 MAF → `.1aln` 完全可行，无骨架限制**。

### 6.2 写侧新增组件

| # | 组件 | 来源 | 说明 | 估算 Rust 行数 |
|---|---|---|---|---|
| B' | **vc 整数编码 + DNA 压缩** | `ONElib.c vcEncode/Compress_DNA` | 读侧 vcDecode 的逆操作；字节差分直方图累加 + Huffman 码表构建 | ~+150–200 |
| A' | **ONEcode 二进制写** | `ONElib.c` 写路径 | 头部 prolog、二进制行写出、对象索引（`&` 行）与 footer（计数/压缩器）写 | ~+350–500 |
| G | **轨迹点重采样** | `PAFtoALN.c cigar2tp`；`align.c` 逆 `Compute_Trace_PTS` | 把带路径比对（PAF-with-CIGAR / MAF / wave `EditOp`）按 `tspace` 重采样为 `tpoints/tdiffs`（含每段差异数）；`cigar2tp` 是现成参考 | ~+200–300 |
| H | **写 CLI 集成** | `open_Aln_Write`/`Write_Aln_Overlap`/`Write_Aln_Trace` | 源基因组名引用、`t` 行 tspace、`A/D/R/T/X/U` 行写出 | ~+150–200 |

写侧合计约 **+800–1200 行**。

### 6.3 写侧关键点

1. **轨迹点采样是新增算法**：把带路径的比对（PAF-with-CIGAR / MAF / wave `EditOp`）
   按 `tspace` 重采样为轨迹点 + 每段差异数。`PAFtoALN.c` 的 `cigar2tp` 是现成参考
   （CIGAR → 轨迹点），pgr 需移植它或等价的 `Compute_Trace_PTS` 逆操作——这是写侧
   最需要新写的部分。
2. **互操作价值**：写成标准 `.1aln` 后，FastGA 生态工具（`ALNchain`/`ALNtoPSL`/
   `ALNshow`/`ALNplot`）能直接消费，实现跨工具互操作。
3. **源基因组引用**：写侧需在头部 `reference` 里记录源基因组文件名（count 1/2/3），
   FastGA 读侧靠它定位 `.bps`/源序列——与 pgr 的 `--ref-seq/--query-seq` 语义对应。
4. **字节级一致性**：若要 pgr 写出与 FastGA 逐字节一致（类似 chain/net 的 UCSC
   验证标准），需对照 `vcEncode`/`Write_Aln_Trace` 的格式细节（含 endian 标志、
   footer 计数、对象索引布局）。

### 6.4 子命令命名方案（定稿）

命名对齐项目惯例：统计类用 `stat`（同 `paf stat`），转换类一律 `to-<输出>`（同
`maf to-paf`、`psl to-paf`、`axt to-psl`）。

**读侧（`.1aln` → 其他）**——新命令 `pgr 1aln`（`src/cmd_pgr/1aln/`）：

| 子命令 | 作用 | 需源序列 |
|---|---|---|
| `pgr 1aln stat` | 头部/统计（tspace、来源、计数） | ❌ |
| `pgr 1aln to-paf` | `.1aln` → PAF | ✅ |
| `pgr 1aln to-psl` | `.1aln` → PSL | ✅ |

**写侧（其他 → `.1aln`）**——挂在源格式命令里，与现有 `to-*` 并列：

| 子命令 | 作用 | 模块 |
|---|---|---|
| `pgr paf to-1aln` | PAF(cg:Z) → `.1aln`（PAFtoALN 等价） | `src/cmd_pgr/paf/to_1aln.rs` |
| `pgr maf to-1aln` | contig 级 MAF → `.1aln` | `src/cmd_pgr/maf/to_1aln.rs` |

说明：
- `pgr 1aln to-paf` 与 `pgr paf to-1aln` 互为镜像，命名规则统一。
- 各 `mod.rs` 只需各加一个转换子模块（或 `1aln` 命令的三个子模块），符合目录结构。
- 写侧 `to-1aln` 的骨架：`paf to-1aln` 用 PAF `srcSize`/源序列，`maf to-1aln` 用
  contig 级 MAF 自带 `srcSize`（§6.1.1）。

## 7. 实施细节（精读参考代码后补全）

> 本节基于 `ONElib.c`/`alncode.c`/`ONEaln.c`/`ALNtoPAF.c`/`PAFtoALN.c`/`GDB.c`
> 精读整理，作为迁移实现参考。所有函数名/行类型均以 `FASTGA-main/` 为准。

### 7.0 建议实现顺序（依赖排序）

算法部分已确认由 `wave.rs` + `paf/cigar.rs` 覆盖，故实现主体是**格式管道**。按依赖
自底向上构建，每阶段可独立验证：

| 阶段 | 内容 | 依赖 FastGA 参考 | pgr 侧重 | 验证 |
|---|---|---|---|---|
| P1 | 整数 codec（LTF） | `intPut`/`intGet` §7.3 | 新写 | `intPut`→`intGet` 往返单测 |
| P2 | vc Huffman + DNA codec | `vc*`/`Compress_DNA` §7.4 | 新写 | `vcEncode`→`vcDecode` 往返单测 |
| P3 | ONEcode 容器（头部/行 I/O/对象索引/footer） | `writeHeader`/`one[Read\|Write]Line`/`oneWriteFooter` §7.2 | 新写 | 往返字段比对 |
| P4 | `.1aln` 记录层（schema + A/D/R/T/X/U + 骨架） | `alncode.c` §7.5/7.6/7.7 | 新写 | 读 gold `.1aln` 字段对齐 |
| P5 | 读侧展开编排 | `ALNtoPAF.c gen_paf` §7.8 | **改可见性** + 薄编排 | `pgr 1aln to-paf` vs `ALNtoPAF` |
| P6 | 写侧 `cigar2tp` + 编排 | `cigar2tp` + `open_Aln_Write` §7.9 | **`cigar2tp` 新写** + 编排 | `ALNshow`/`ALNtoPAF` 回读 |
| P7 | CLI 命令 | §6.4 | 新写 | 集成测试 |

> 关键：P1–P4 是全新代码（容器 + codec + schema），P5–P6 的**算法已被现有资产覆盖**，
> pgr 侧只补可见性（`pub`）+ 编排 + `cigar2tp`。P5/P6 是"组装"，不是"新算法"。

### 7.1 参考文件与依赖

> 行号随文标注（`文件:行`），供实现时直接定位。`ONElib.c` 有两处声明+定义，
> 行号取**定义处**。

| 文件 | 作用 | pgr 对应 |
|---|---|---|
| `ONElib.c` | ONEcode 容器（header/schema/行编码/footer）+ vc Huffman + LTF 整数 + DNA 2bit | 需全移植（读+写） |
| `alncode.c` | `.1aln` schema 常量 + `open/read/write` 记录骨架 | 直接复刻常量 |
| `ONEaln.c` | 对齐器，产出 `.1aln`（tspace 采样 + 写轨迹） | 复用 wave.rs，只取写侧 |
| `ALNtoPAF.c` | 读侧展开 `.1aln`→PAF（CIGAR/CS） | 读侧主参考 |
| `PAFtoALN.c` | 写侧 PAF(cg:Z)→`.1aln`（`cigar2tp`） | 写侧主参考 |
| `align.c` | `Compute_Trace_PTS`/`Gap_Improver`/`Decompress_TraceTo16` | 大部已由 wave.rs 覆盖 |
| `GDB.c` | `Read_Skeleton`/`Write_Skeleton`/`Get_GDB` | 骨架 + 序列访问 |

**关键函数行号速查**：

| 符号 | 位置 |
|---|---|
| `.1aln` schema 常量 `alnSchemaText` | `alncode.c:19` |
| `open_Aln_Read` / `open_Aln_Write` | `alncode.c:59` / `alncode.c:239` |
| `Read_Aln_Overlap` / `Read_Aln_Trace` / `Skip_Aln_Trace` | `alncode.c:136` / `168` / `210` |
| `Write_Aln_Overlap` / `Write_Aln_Trace` / `Copy_Aln_Trace` | `alncode.c:272` / `288` / `307` |
| `oneFileOpenRead` / `oneFileOpenWriteNew` | `ONElib.c:1350` / `1817` |
| `oneReadLine` / `oneWriteLine` | `ONElib.c:1077` / `2355` |
| `writeHeader` / `writeInfoSpec` / `writeCounts` | `ONElib.c:2211` / `455` / `2186` |
| `oneWriteFooter` / `oneFinalizeCounts` / `oneFileClose` | `ONElib.c:2617` / `2678` / `2803` |
| `intPut` / `intGet` / `ltfWrite` / `ltfRead` | `ONElib.c:3778` / `3737` / `3835` / `3804` |
| `vcAddToTable` / `vcCreateCodec` / `vcSerialize` / `vcDeserialize` | `ONElib.c:2971` / `3012` / `3298` / `3338` |
| `vcEncode` / `vcDecode` / `vcMaxSerialSize` | `ONElib.c:3479` / `3621` / `3293` |
| `Compress_DNA` / `Uncompress_DNA` | `ONElib.c:3443` / `3577` |
| `Read_Skeleton` / `Skip_Skeleton` / `Write_Skeleton` / `Get_GDB` | `GDB.c:1952` / `2062` / `2070` / `1655` |
| `Compute_Trace_PTS` / `Gap_Improver` / `Decompress_TraceTo16` | `align.c:6171` / `6714` / `3912` |
| `gen_paf` / `main`（读侧展开） | `ALNtoPAF.c:102` / `638` |
| `cigar2tp` / `main`（写侧） | `PAFtoALN.c:215` / `745` |

### 7.2 ONEcode 容器二进制布局

**头部（ASCII，`writeHeader`）**：

```
1 <len> <fileType> <major> <minor>          # 主类型名 + 版本（major 必须相等，minor ≤ 当前）
[2 <len> <subType>]                          # 可选子类型
[. <headerText>]                             # 头部正文
[! 4 <plen> <prog> <vlen> <vers> <clen> <cmd> <dlen> <date>]   # provenance
.                                            # 空行分隔
[< <len> <filename> <count>]                 # reference（count 1/2/3 = db1/db2/cpath）
[> <len> <filename>]                         # deferred
.
[~ O|D <type> <nField> [<strlen> <TYPE> <strlen> <TYPE> ...]]  # schema，每行类型一条
[~ G <type> 0]                               # group 行
$ <isBig>                                    # endian 标志（二进制）
```

**二进制行编码（`oneWriteLine`）**：每条数据记录一行，由一个字节引导，接着是字段区
和可选的列表区：

- 引导字节 `x = binaryTypePack`（每行类型一个唯一编码，高位 0x80 置位以区分 ASCII）；
  若该行列表用 codec 压缩，则 `x |= 0x01`。
- 字段区：`INT`/`CHAR` 用 LTF（§7.3）写；`STRING`/`DNA` 先写长度再写字节。
- 列表区：
  - `INT_LIST`：先用 LTF 写首元素；若 `listLen==1` 结束；否则再接一个字节
    `intListBytes`（后续每元素字节数），然后要么按 `intListBytes` 定宽写剩余元素，
    要么（压缩时）写 `ltfWrite(nBits)` + `⌈nBits/8⌉` 字节的 Huffman 数据（§7.4）。
  - `STRING_LIST`/`REAL_LIST`：按 ASCII 逐元素写。
- 行与行之间用 `\n` 分隔；数据区结束前有一个空行（`\n`）作为 end-of-data 标记。

**对象索引 + footer（`oneWriteFooter`）**：文件末尾对每个计数 >0 的行类型：

- 写计数行 `# <type> <count>`、`@ <type> <max>`、`+ <type> <total>`（ASCII）。
- 若该类型是对象（object），写 `&` 行：`& <count+1> [<byteOffset>...]`——即每个对象
  的字节偏移索引，供 `oneGoto` 随机访问。
- 若该类型用了 list codec，写 `;` 行：`; <n> [serialized codec]`（§7.4）。
- 结束标记 `^\n`，随后 8 字节写入 footer 起始的 `off_t` 偏移。

### 7.3 整数 codec（LTF，`intGet`/`intPut`）

值域分级、首字节高位编码：

| 首字节 | 含义 | 字节数 |
|---|---|---|
| `0x40–0x7f` | 单字节正数（6 bit） | 1 |
| `0xc0–0xff` | 单字节负数（符号扩展 8 bit） | 1 |
| `0x20–0x3f` | 双字节正数（13 bit） | 2 |
| `0x00–0x07` | 正数，低 3 位 = 后续字节数（1–8） | 3–9 |
| `0x80–0x87` | 负数，低 3 位 = 后续字节数（1–8），符号扩展 | 3–9 |

对应 `intPut`/`intGet` 完全对称。这是 `.1aln` 里所有标量字段（坐标、长度、diffs）的
基础编码，也是 T/X 列表首元素/`nBits` 的编码。

### 7.4 vc Huffman codec（`vc*`）

**建码**：`vcAddToTable` 累加字节直方图 → `vcCreateCodec(partial)` 用
Larmore-Hirschberg（`JACM 73,3 (1990)`）长度受限 Huffman 建码。`partial=1` 时若存在
零计数字节则保留一个 escape 码（代码长度 `esp_len`，码值 `esc_code`），用于流中途出现
未见字节时的转义。

**编码 `vcEncode`**：`ilen` 输入字节 → 压缩到 `obytes`，返回**比特数**。逐字节查
`codelens`/`codebits`，未见的字节走 escape（写 `esc_code` 后再写原始 8 位）。若中途
压缩收益为负（`tbits > ilen*8`）则放弃，首字节写 `0xff` 后跟原始字节（明文 fallback）。

**解码 `vcDecode`**：`ilen` 比特 → 解出字节。首字节 `0xff` 即明文 fallback。用 16 位
前缀查找表 `lookup[0x10000]` 快速解码；`GET` 宏维护 64 位位缓冲。endian 安全：压缩流
首字节 `0x40` 位记录 endian，读写不一致时按 64 位翻转。

**序列化 `vcSerialize`**：写 `isbig`(1B) + `esc_code`(int) + `esc_len`(int) + 256 个
`{codelen(1B), [code(2B) 若 len>0 或为 esc]}`。最大 `vcMaxSerialSize()` 字节。这是
footer 里 `;` 行存的内容。

**DNA codec（`Compress_DNA`/`Uncompress_DNA`）**：2bit/碱基，**小端**打包（`Number`
表 A/C/G/T=0/1/2/3，非 ACGT 记为 0），4 碱基/字节，余数按尾部处理。`vc==DNAcodec`
时 `vcEncode`/`vcDecode` 直接走此路径。

### 7.5 `.1aln` schema（`alnSchemaText`，头部内嵌）

```
P 3 seq SEQUENCE
  O s 2 3 INT 6 STRING       scaffold: length + id
  G S                        scaffold(s) 组 sequence 对象 S
  D n 2 4 CHAR 3 INT         scaffold 内非 acgt 部分
  O S 1 3 DNA                sequence
  D I 1 6 STRING             sequence 标识符
P 3 aln ALIGNMENTS
  D t 1 3 INT                trace point spacing（tspace，全局）
  O g 0                      GDB skeleton（可缺失）
  G S                        scaffold 集合
  O S 1 6 STRING             scaffold id
  D G 1 3 INT                gap 长度
  D C 1 3 INT                contig 长度
  O a 0                      A 的共线链
  G A                        chain(a) 组 alignment(A)
  D p 2 3 INT 3 INT          a/b 相邻比对间距
  O A 6 3 INT *6             aread abpos aepos bread bbpos bepos
  D L 2 3 INT 3 INT          a、b 序列长度
  D R 0                      反向互补标记（b）
  D D 1 3 INT                diffs = 置换 + 插删数
  D T 1 8 INT_LIST           轨迹点（b 坐标，差分）
  D X 1 8 INT_LIST           每个轨迹区间内差异数
  D Q 1 3 INT                质量（未用）
  D E 1 3 INT                match 数（未用）
  D Z 1 6 STRING             CIGAR（未用）
  D U 1 3 INT                TR 比对单元长度
```

> 关键：schema 全文内嵌于头部，**读侧无需硬编码**，直接解析头部即可。pgr 只需用一个
> 常量复刻 `alnSchemaText`（写侧用）或从头部反解析（读侧用）。

### 7.6 记录布局（`Read_Aln_Overlap`/`Read_Aln_Trace`）

每条比对按序为：

- `A` 行：6 个 INT（`aread abpos aepos bread bbpos bepos`），object 起始。
- 可选 `R` 行：置 `COMP_FLAG`（b 反向互补）。
- `D` 行：1 个 INT = `diffs`。
- `T` 行：INT_LIST，长度 `k`，内容 → 轨迹点。
- `X` 行：INT_LIST，长度 `k`，内容 → 每区间差异数。
- 可选 `U` 行：1 个 INT = TR 单元长度。

**轨迹内存布局（`Read_Aln_Trace` 后）**：`tlen = 2k`，交错存放
`trace[0..tlen)`：`trace[2i] = diff_i`（来自 X），`trace[2i+1] = tpoint_i`（来自 T）。
其中 `tpoint_0` 是相对 `bbpos` 的坐标，`tpoint_i (i>0)` 是相对前一点的差分。

### 7.7 骨架（`Read_Skeleton`/`Write_Skeleton`）

对象 `g` 内按序：

- `S` 行（STRING）：scaffold id（名字）。
- `C` 行（INT）：contig 长度。
- `G` 行（INT）：两个 contig 之间的 gap 长度。

骨架语义：每个 scaffold 含一组连续 contig（`fctg..ectg`），contig 之间允许 gap；每个
contig 有 `sbeg`（scaffold 内偏移）、`clen`（长度）、`boff`（压缩序列偏移）。**无
scaffold 层时**（我们的 contig 级场景）即每个 scaffold 恰好一个 contig、无 gap。

### 7.8 读侧数据流（`.1aln` → PAF/PSL，参考 `ALNtoPAF.c gen_paf`）

1. `open_Aln_Read`：解析头部，读 `t` 行得 tspace，`<` 行得源基因组名（count 1/2/3）。
2. 骨架：若 `g` 行在比对前，`Read_Skeleton` 建 GDB；否则 `Get_GDB` 从源 FA/GDB 读。
3. 每条比对：`Read_Aln_Overlap` + `Read_Aln_Trace` 得 `aread/abpos/aepos/bread/bbpos/
   bepos/diffs/trace`。
4. `Decompress_TraceTo16`：把 uint8 轨迹扩宽为 uint16（`tspace > TRACE_XOVR` 时需要）。
5. 取两侧 contig 序列（`Get_Contig`/`Get_Contig_Piece`，COMP 时补全互补）。
6. **盒内 DP（复用 `wave.rs`）**：FastGA `Compute_Trace_PTS` 在相邻轨迹点之间做盒内
   DP 重建碱基层路径。pgr 用 `banded_edit_ops`（`wave.rs:235`，需 `pub`）实现：相邻两
   轨迹点作锚点，取对角带 `[k_lo,k_hi]` 调用，得到 `Vec<EditOp>`。
7. **对齐列（复用 `wave.rs`）**：`ops_to_columns(..., ops)`（`wave.rs:332`，需 `pub`）
   把 `EditOp` 展开为 `(q_aln, t_aln, matches)`。`Gap_Improver` 跳过（见 §4 注意点 2）。
8. **CIGAR（复用 `paf/cigar.rs`）**：`cigar_from_alignment(q_aln, t_aln)`（`cigar.rs:332`）
   生成 X-CIGAR；`-m` 合并 `=`/`X` 为 `M` 时自行合并，或直接用 `=`/`X`/`I`/`D`。
9. **输出 PAF/PSL**：坐标换算（COMP 时逆转 b），`blocksum`、`iid=blocksum-diffs`、
   mapq 255、`dv:f:` 恒等分数、`df:i:` diffs、可选 `cg:Z:`/`cs:Z:`。

> 步骤 6–8 的**算法全部复用**，仅需给 `wave.rs` 三个私有函数加 `pub`；pgr 侧新写的是
> 步骤 1–5 的轨迹点解码/骨架编排，以及步骤 9 的 PSL/PAF 落盘。

### 7.9 写侧数据流（→ `.1aln`，参考 `PAFtoALN.c` + `alncode.c`）

写侧**输入是带碱基层路径的比对**，来源三选一：`WaveAlign`/`LocalAlign`（`wave.rs`，
PAF 带 `cg:Z`，或 MAF `s` 行）。pgr 侧链条：

1. **CIGAR（复用 `paf/cigar.rs`）**：`cigar_from_alignment(q_aln, t_aln)` 从对齐列生成
   X-CIGAR（`=`, `X`, `I`, `D`）。若输入已是 PAF 带 `cg:Z`，直接用其 `CigarOp`。
2. **轨迹点采样（`cigar2tp`，唯一新写算法）**：把 X-CIGAR 按 tspace 重采样为
   `tpoints`（差分）+ `tdiffs`（每区间差异数）。参考 `PAFtoALN.c:215`。
3. **记录落盘（复用/新写的容器层）**：
   - `open_Aln_Write`：写头部（类型 `aln`、provenance、reference、`t` 行 tspace）。
   - `Write_Skeleton`：写 `g` 对象（若源提供骨架）。
   - `Write_Aln_Overlap`：写 `A` 行 + 可选 `R` + `D` 行。
   - `Write_Aln_Trace`：交错写 `T`（tpoints）+ `X`（tdiffs）+ 可选 `U`。
4. `oneFileClose` → `oneFinalize` → `oneWriteFooter`：合并线程索引、写 `&`/`;`/`^`。

> 步骤 1 复用 `cigar_from_alignment`，步骤 2 是**唯一新写的写侧算法**，步骤 3–4 是
> P3/P4 阶段已建好的容器/记录层。接入点：`align pgi` 输出端（`WaveAlign`→CIGAR）或
> `pgr paf to-1aln`/`pgr maf to-1aln`（读入 CIGAR）。

### 7.10 与 pgr 现有资产对应

| 组件 | FastGA | pgr 现状 |
|---|---|---|
| 轨迹点间 DP | `Compute_Trace_PTS`/`iter_np` 盒内 DP | `wave.rs:banded_edit_ops` 已覆盖 |
| 精确回溯 | `dandc_nd`/`split_nd` | `wave.rs:184/42` 已覆盖 |
| 对齐列表示（写侧源） | `align->tpoints`/`tdiffs` | `wave.rs:WaveAlign{q_aln,t_aln}` 已覆盖 |
| CIGAR 生成（写侧） | `cigar2tp` 的输入 | `paf/cigar.rs:cigar_from_alignment`/`cs_from_alignment` 已覆盖 |
| 轨迹点采样（写） | `cigar2tp`（PAFtoALN） | **唯一需新写的写侧算法**（输入 CIGAR 已现成） |
| 整数 codec | `intGet`/`intPut`（LTF） | 需新写 ~120 行 |
| Huffman codec | `vcEncode`/`vcDecode`/`vcSerialize` | 需新写 ~400 行 |
| DNA 2bit | `Compress_DNA`/`Uncompress_DNA` | 复用/对照 2bit 模块 |
| 容器行 I/O | `oneReadLine`/`oneWriteLine` | 需新写 ~500 行 |
| 骨架 | `Read/Write_Skeleton` | 需新写 ~150 行 |
| PAF/PSL 输出 | `ALNtoPAF`/`ALNtoPSL` | `fmt/psl.rs`+`align.rs` 已覆盖 |

**已核实的接口签名（2026-08-05 精读 wave.rs / paf/cigar.rs）**：

读侧（`.1aln`→PAF/PSL，还原轨迹点间路径）：
- `banded_edit_ops(q,t,q_abs,t_abs,k_lo,k_hi,ops)->usize` `wave.rs:235` ——对角带
  `[k_lo,k_hi]` 内、锚点对角线上重建路径，产 `Vec<EditOp>`；**即轨迹点间盒内 DP**。
- `dandc_nd(q,t,q_abs,t_abs,ops)->usize` `wave.rs:184` ——精确 Myers D&C（diff 大时）。
- `ops_to_columns(...,ops)->(q_aln,t_aln,matches)` `wave.rs:332` ——`EditOp`→对齐列。
- `cigar_from_alignment(r#ref,qry)->Result<Vec<CigarOp>>` `cigar.rs:332`；`cs_from_alignment` `cigar.rs:363` ——对齐列→X-CIGAR / CS。

写侧（→`.1aln`）：
- `WaveAlign{q_aln,t_aln,q_start,t_start,matches}` `wave.rs:12`；`LocalAlign{...,diffs}` `wave.rs:722` ——对齐列（碱基层路径）。
- `cigar_from_alignment(q_aln,t_aln)` → `Vec<CigarOp>` → `cigar2tp`（新移植）。
- `CigarOp` 位打包 `cigar.rs:9`：`bits[31:29]`=0/1/2/3/4 → `=`/`X`/`I`/`D`/`M`，
  与 FastGA `cigar2tp` 消费的字符一一对应，解包转换平凡。

**可见性缺口**：`banded_edit_ops`/`dandc_nd`/`ops_to_columns` 在 `wave.rs` 中是**私有
`fn`**，读侧需先改 `pub`（各一行）。读侧还需一段**薄编排**：给定两个轨迹点，决定对角带
`[k_lo,k_hi]` 再调 `banded_edit_ops`（对应 FastGA `Compute_Trace_PTS` 的 box 选择）。

> 结论：迁移主体是 **ONEcode 容器 + codec + schema 的格式管道**；比对逻辑（读侧盒内
> DP）与写侧 CIGAR 生成均已被 `wave.rs` + `paf/cigar.rs` 覆盖，真正新写的写侧算法仅
> `cigar2tp`（轨迹点采样）。

### 7.11 验证/测试策略

测试料：FASTGA 目录内仅 `EXAMPLE/H1vH2.1aln`（22MB，人类、多 contig 骨架），不适合做
常规回归。黄金 `.1aln` 由 FastGA 工具链在 E. coli MG1655/Sakai（`tests/genome/`）上
**自产**，脚本：`scripts/gen-1aln-golden.sh`（约 117KB，单 contig，贴合 40k E. coli
用途）。对应命令：

```bash
export PATH="$PWD/FASTGA-main:$PATH"          # FastGA 经 system() 调 FAtoGDB，需在 PATH
./FASTGA-main/FastGA -k -T8 -1:mg1655-sakai.1aln \
    tests/genome/mg1655.fa.gz tests/genome/sakai.fa.gz
./FASTGA-main/ALNtoPAF mg1655-sakai.1aln | head   # 回读验证
```

- **读侧**：FastGA `ALNtoPAF`/`ALNtoPSL` 输出作 golden，pgr 输出与其逐字段一致
  （坐标、strand、CIGAR、CS）。
- **写侧**：pgr 写出的 `.1aln` 用 FastGA `ALNshow`/`ALNtoPAF` **回读**，比对还原结果；
  若需字节级一致，逐字节对照 `vcEncode`/`Write_Aln_Trace`/`oneWriteFooter`。
- **纯容器**：`ONElib` 自带的 `oneFile` 往返测试思想可移植为 Rust 单元测试（写→读→
  比对字段）。
- **codec 性质**：`vcEncode`→`vcDecode` 往返恒等；`intPut`→`intGet` 往返恒等；跨端序
  读（小端/大端标志）测试。
- **冒烟**：`EXAMPLE/H1vH2.1aln`（22MB）仅作容器/codec 冒烟 + 多 contig 骨架覆盖，
  不入常规回归。