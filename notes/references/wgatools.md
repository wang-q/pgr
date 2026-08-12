# wgatools：Rust PAF/MAF 工具库（源码分析）

> 2026-08-12 整理，纯源码分析（`wgatools-1.1.0/`，Rust，BSD-2）。wgatools
> 是**全基因组比对格式（PAF/MAF/Chain）工具库**：三种格式互转 + 20 个已注册
> 子命令（`maf2paf`/`paf2maf`/`paf2chain`/`chain2maf`/`maf-index`/`maf-ext`/
> `chunk`/`call`/`tview`/`stat`/`dotplot`/`filter`/`rename`/`pafcov`/
> `pafpseudo`/`gen-completion`/`validate` 等，见 `src/cli.rs` 的 `Commands`
> 枚举，另有 DEV 态 `maf2sam`）。**与 pgr 的关系**：pgr `paf` 模块是自研 PAF
> 隐式图（index/query/to-gfa/to-vcf/stat/validate），wgatools 是最接近的
> 外部参照，可补 pgr 缺失的 PAF 通用操作（filter / end-trim / 解析器对照）。

> **注意（1.1.0 版本事实核查）**：`trimovp` 与 `pileup` **并非已注册子命令**——
> 二者在 `src/cli.rs:285-291`、`src/cli.rs:315-324` 中均被注释掉（`main.rs:194-196`
> 的 dispatch 与 `utils.rs:643-653` 的 wrap 同样被注释），`tools/trimovp.rs` 全文件
> 为注释。下文 §4.4 记录了其**注释掉的意图算法**，供 pgr `--end-trim` 参考。

## 1. 工具对照（pgr paf vs wgatools）

| wgatools | 作用 | pgr 现状 |
|---|---|---|
| `stat` | PAF/MAF 统计（identity/similarity/indel，按 ref 聚合） | ✅ `paf stat` |
| `validate` | PAF 校验，`-f` 可选修复末端坐标 | ✅ `paf validate`（2026-08-06，仅报告不修复） |
| `maf-index` | **MAF** 索引（JSON：name→区间+字节偏移），非 PAF 索引 | ✅ `paf index`（区间树 + `.paf.idx`，设计更进阶） |
| `pafcov` | PAF → 逐碱基覆盖度 BED | ✅ `pgr paf coverage`（2026-08-12 新增，cg:Z 扫描线） |
| `filter` | 按 block/query 大小、pair 累计比对大小过滤 | ❌ pgr 无 PAF 过滤子命令 |
| `trimovp` | overlap 末端修剪 | ⚠️ wgatools 1.1.0 中已注释停用，仅余意图算法 |
| `call` | MAF/PAF → VCF 变异（SNP/INS/DEL/INV） | 部分与 `paf to-vcf` 重叠 |
| `maf-index`+`maf-ext` / `chunk` / `tview` / `rename` / `dotplot` / `pafpseudo` | MAF 索引提取 / 分块 / 终端查看 / 重命名 / dotplot / 伪 MAF | 部分与 `paf to-fas/to-maf`、`plot dot` 重叠 |

## 2. 解析器与 CIGAR（`src/parser/`）

### 2.1 PAF 解析 —— 用 `csv` crate 解析制表符格式

`PAFReader`（`src/parser/paf.rs:13`）**直接复用 `csv` crate 做制表符行解析**：

```rust
ReaderBuilder::new()
    .flexible(true)      // 允许可变列数（PAF 12 列 + 任意 tag 数）
    .delimiter(b'\t')
    .has_headers(false)
    .comment(Some(b'#')) // 自动跳过 `#` 注释行
```

`PafRecord`（`paf.rs:50`）12 个必填字段 + `#[serde(default)] tags: Vec<String>`。
`AlignRecord` trait（`common.rs:142`）统一 PAF/MAF/Chain 三类的存取接口，让
stat/filter/call/dotplot 等工具泛型复用。

### 2.2 CIGAR —— `nom` 分词器（零分配、无正则）

CIGAR 用 `nom` 手写分词而非正则（`cigar.rs:268-274` 注释明确说"正则比 nom 慢 3 倍"）：
- `parse_cigar_str_tuple`（`cigar.rs:59`）：`take_while(dec_digit)` 取长度 +
  `take_till(dec_digit)` 取 op，返回 `(op, len)` 元组；空输入返回 Eof 以终止 `fold_many1` 循环。
- `cst2cu`（`cigar.rs:43`）：校验 op 恰好为单字符（`chars.next()` 再取一次必须为 None）。
- `fold_many1(parse_cigar_str_tuple, null, ...)` 迭代累加，错误沿 `IResult` 传播（`cigar.rs:277`）。

核心 CIGAR 解释函数（均接受 `&dyn AlignRecord`，自动从 `cg:Z:`/`cs:Z:` 取串）：
- `parse_paf_to_cigar`（`cigar.rs:629`）→ `Cigar` 结构（match/mismatch/ins/del 的事件数+碱基数，
  负链时 ins/del 计入 `inv_*`）；`RecStat::from(Cigar)`（`common.rs:116`）再算
  `aligned_size = matched+mismatched+del+inv_del`。
- `cs_to_cigar`（`paf.rs:159`）：把 `cs:Z:` 紧凑格式（`:N`/`-seq`/`+seq`/`*x`）转成标准
  CIGAR（`=`/`D`/`I`/`X`），正则 `(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)` 逐段归并。
- `parse_cigar_to_chain`（`cigar.rs:251`）：CIGAR → Chain `ChainDataLine`，`cigar_unit_chain`
  （`cigar.rs:460`）：`M/=X` 累加 size，`I` 累加 target_diff，`D` 累加 query_diff，非零 diff
  时写出 dataline。
- `update_cov_vec`（`cigar.rs:710`）：覆盖度扫描线——仅 `M`/`=` 递增 `cov_vec[pos]`，
  `I`/`S` 不动（query-only），其余（`D` 等）只推进 pos。即覆盖度只计比对/target 消耗碱基。
- `parse_maf_seq_to_cigar`（`cigar.rs:344`）：双序列逐位分类（`cigar_cat_ext`：相等`=`、
  目标`-`→`I`、查询`-`→`D`、否则`X`），用 `itertools::group_by` 聚成 run；`with_h=true`
  时在头尾加 `H`（用于 SAM 输出，pgr 不需要）。

### 2.3 MAF 解析

`MAFReader`（`maf.rs:15`）按行读：跳过非 `s` 开头的空行，连续 `s`-line 聚成一个 `MAFRecord`
（`maf.rs:374-421`）。`MAFSLine::get_col_coord`（`maf.rs:81`）做**区域坐标→列索引**映射
（跳过 `-`），供 `slice_block`（`maf.rs:223`）按区间截取序列，其余 s-line 同步截列并重算
`align_size = 列长 - gap 数`。负链 query_start 用 `size - start - align_size` 换算
（`maf.rs:438`）。MAF 索引 `build_index`（`tools/index.rs:14`）记录每条 s-line 的
`start/end/strand/byte-offset` 到 JSON，并做**重复名检测**与"同一序列不能既 ref 又 query"检查。

## 3. 子命令语义（`src/tools/`）

- **`stat`**（`stat.rs`）：`rayon` 并行，按 `Pair`（ref/query 名+size）分组。`identity =
  matched/aligned_size`，`similarity = (matched+mismatched)/aligned_size`；合并模式下
  `unaligned_size = ref_size - aligned_size`（`stat.rs:216`）。输出 TSV（`csv::WriterBuilder`
  tab 分隔），用 `natord`（自然序）排序 ref_name（`stat.rs:116`）。`-e/--each` 逐 block 不聚合。
- **`validate`**（`validate.rs`）：`exp_query_end = query_start + matched+mismatched+ins(+inv_ins)`
  （`validate.rs:80`）、`exp_ref_end = target_start + matched+mismatched+del(+inv_del)`
  （`validate.rs:98`），不一致即标记；`-f FILE` 时把末端坐标改为期望值输出修复版 PAF。
  ⚠️ 它 `rec.get_stat().unwrap()`（`validate.rs:77`）——缺 CIGAR 会 **panic**，
  违背 pgr 的 zero-panic 准则（对比见 §4.3）。
- **`filter`**（`filter.rs`）：`filter_alignrec`（`filter.rs:91`）用
  `block_length = target_align_size() = target_end - target_start`、`query_length` 双阈值
  （注意用**按位或 `|`** 非短路）。`filter_paf_align_pair`（`filter.rs:108`，即 `-a`）为
  all-to-all 设计：并行按 `(q_name,t_name)` 累加 `align_size` 之和，≥ 阈值才保留该 pair 的
  全部记录；`-a` 生效时忽略 `-b/-q`（`utils.rs:559` warn）。
- **`pafcov`**（`pafcov.rs`）：`HashMap<target, Vec<usize>>` 全分辨率覆盖向量，rayon
  `try_fold/try_reduce` 并行累加，输出**逐碱基 BED4**（`pafcov.rs:57`）。注释里留了一段
  "collapsed 恒定段"备选实现（`pafcov.rs:62-81`）——pgr `coverage` 正是用后者（省内存更高效）。
- **`call`**（`caller.rs`）：MAF/PAF → VCF，`-s` SNP、`-i` INV、`-l` SVLEN 阈值。
- **`chunk`**（`chunk.rs`）：按长度切 MAF，`recount_align_size`（`common.rs:179`）数非 `-` 碱基
  以保持各 s-line 对齐。

## 4. pgr 可借鉴的工程点

### 4.1 PAF 解析用 `csv` crate 做制表符行解析（`paf.rs:22-30`）

`flexible(true)` + `delimiter(b'\t')` + `comment(b'#')` 一个 ReaderBuilder 即搞定
"可变列数 + 注释行"两个 PAF 痛点。pgr `paf::parser` 目前手写逐列 split，可考虑改用
`csv`/`noodles` 的制表符 Reader 以省去边界处理（不引新依赖，pgr 已有 `noodles`）。

### 4.2 `cs:Z:` → CIGAR（`paf.rs:159 cs_to_cigar`）

pgr `maf_import`（`maf to-paf`）已产出 FastGA 风格 `cs:Z:`。若未来 pgr 的 PAF 工具需消费
第三方 `cs:Z:` 输入，可对照 `cs_to_cigar` 的逐段归并逻辑（`=`/`D`/`I`/`X`）实现反向解析。

### 4.3 validate 的"期望末端"公式（`validate.rs:80,98`）与 pgr 对齐

pgr `libs/paf/validate.rs`（`ValidationReport::validate`，行 38-74）用**同一套公式**
（query = +ins_bp，target = +del_bp），且对 missing/malformed CIGAR **计数跳过、绝不 panic**，
比 wgatools 的 `.unwrap()` 更稳健（符合 pgr zero-panic 准则）。wgatools 独有的 `-f` **修复模式**
（改写末端坐标为 CIGAR 期望值）是 pgr validate 可考虑补的能力（目前仅报告）。

### 4.4 `--end-trim` 的现成语义与意图算法（`cigar.rs:155,202` + `trimovp.rs` 注释）

- `parse_cigar_to_trim`（`cigar.rs:202`）/ `parse_maf_seq_to_trim`（`cigar.rs:155`）：
  沿 CIGAR/双序列扫描，统计**头尾**的 INS/DEL 碱基数（遇 `M/=X` 置 `head_indel=false`、
  清空 tail 计数）——正是 pgr `paf-pangenome.md` "`--end-trim` 推迟"待办的直接可对照实现。
- `trimovp` 意图算法（`trimovp.rs` 整文件注释）：按 target 分组 → 按 query 分组 → 按
  `target_start` 二分插入排序 → 贪心保留与上一条不重叠且更长的记录。可作 pgr per-interval
  修剪的语义起点，但该命令在 1.1.0 未落地，需自行验证其正确性。

### 4.5 通用过滤（`filter.rs`）—— pgr 缺失子命令的参数形态

pgr `paf` 缺通用过滤。wgatools 三档阈值可参考：`-b` block 大小（target 跨度）、`-q` query
长度（contig 场景）、`-a` pair 累计比对大小（all-to-all，对应 pgr `query`/`to-bed` 的
`--min-tree-coverage` 语义）。

### 4.6 IO 工程技巧（`utils.rs`）

- **magic-byte 嗅探 + 自动解压**（`utils.rs:37-169`）：读文件头 6 字节识别 gz(`1f 8b 08`)/
  bz2(`42 5a 68`)/xz(`fd 37 7a 58 5a 00`)，`get_input_reader` 按需包 `MultiGzDecoder` 等；
  `get_output_writer`（`utils.rs:181`）按扩展名自动压缩（zlib level 6）。pgr 的
  `libs/io.rs` 已支持 gz/BGZF，可对照补齐 bz2/xz 的 magic 判定。
- **空 stdin 守卫**（`utils.rs:172 stdin_reader`）：用 `atty::is` 检测 stdin 非 tty，为空即报错，
  避免管道空输入时的阻塞/异常——pgr 的 `stdin` 约定可借鉴。
- **输出文件重写保护**（`utils.rs:231 check_outfile`）：缺 `-r/--rewrite` 时已存在输出文件报错，
  防止误覆盖。

### 4.7 并行聚合范式（`stat.rs`/`pafcov.rs`/`validate.rs`）

三处都采用 `reader.records().par_bridge().try_fold(init, acc|..).try_reduce(...)`：把**逐条
解析 + 并行聚合 + 错误传播**封装成一条链。pgr 大文件统计（`paf stat`）若需并行化可照此模式；
但注意 `par_bridge` 下的聚合归并需自定义 reducer（如 pafcov 的 `HashMap` 逐元素相加，
`pafcov.rs:46`）。

## 5. 结论

价值中等——pgr `paf` 已覆盖核心（图/索引/查询），且 pgr 的 validate 容错与 coverage
恒定段输出甚至优于 wgatools。wgatools 的增量参考是：**`filter` 三档阈值**、
**`parse_cigar_to_trim` 头尾 indel 计数（`--end-trim`）**、**`cs_to_cigar` 反向解析**、
**csv 制表符解析**与 **magic-byte IO**。按 `todo.md` §2 的优先级（`--min-tree-coverage`、
`--end-trim`）推进时再细读 §4.2/§4.4 对应文件即可，暂不立项。
