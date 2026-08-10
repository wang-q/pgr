# pgr fa / 2bit 命令族代码审核记录（2026-08-05）

> 原 `audit-fa-fq-2bit.md`。`fq` 部分已拆分至独立记录 `audit-fq.md`，
> 本文件现仅覆盖 fa / 2bit。

对 `pgr fa`（18 个子命令）与 `pgr 2bit`（5 个子命令）命令族及相关库文件
（`libs/fmt/fa`、`libs/fmt/twobit`、`libs/translate`、`libs/io`、
`libs/ds/range`、`libs/nt`、`libs/loc`、`libs/fasta/{chunk,stat,dedup,filter}`）
和全部测试/文档进行审核。`pgr fq` 命令族已拆分为独立记录 `audit-fq.md`。
以下仅保留有借鉴意义的结论；验证过程已精简。

## 与外部参考实现的语义一致性核对

2bit 家族对照 UCSC kent-tools 的 2bit 二进制格式规范逐字节核对，语义一致。
有意差异（均已记录）：
- 读取端同时支持 version 0（index 偏移 4 字节）与 version 1（index 偏移 8
  字节）两种布局；写入端产出 version 1。
- 打包时任何非 A/C/G/T 碱基（含 IUPAC 歧义码、U）按 N 处理并记录为
  hard-mask 块；软屏蔽小写 A/C/G/T 记录为 mask 块。`size --no-ns` 仅扣减
  N-block，与 UCSC "排除 hard-masked N 位置"语义一致。
- 畸形输入：kent-tools 部分输入直接崩溃/告警，pgr 统一友好出错/跳过
  （Zero-Panic）。

## 排除的疑点（安全不变量，经核验无需修复）

- `Range::from_str` 手写字节扫描器（不 panic），无 `:` 时回退 `chr` 为首个空白
  token；`range` 命令对含 `:` 的输入用 `rg.is_valid()` 兜底。
- `from_dna` 打包：`packed.len()` 恒等于 `ceil(len/4)`，与读取端
  `dna_size.div_ceil(4)` 一致；`bit_offset` 归零即推入整字节。
- `TwoBitWriter::write` 的 `record_size` 与实际写出的字节数一致；v1 偏移 8 字节，
  index 起算 `16 + Σ(1+N+8)` 正确。
- `read_sequence` 的 `cached` 记录缓存与 `get_sequence_len`/`get_sequence_blocks`
  的 seek 无位置串扰：所有读取方法均在使用前显式 seek，缓存命中时不再 seek。
- `translate` 反向链坐标：ORF 的 `end*3` 不超过该框碱基数，`dna_len - end*3`
  不会下溢，坐标恒 ≥ 1。

## 已知限制（有意保留）

- `2bit range` 的全序列反链请求 `seq_name(-)`（无坐标）不被解析：`Range::from_str`
  在无 `:` 时回退 `chr` 为整个 token，`strand` 为空，`has_sequence("seq_name(-)")`
  为假 → 告警跳过。此用法未文档化，且不会静默返回错误数据（仅跳过）。
- `2bit range` 在既无位置参数 `ranges` 也无 `--rgfile` 时静默输出空结果（与
  `fa range` 行为一致）。

## 记录项（未改，低风险 / 待决策）

- `fa split name`：`sanitize_filename` 把 `/`、`\`、`(`、`)`、`:` 替换为 `_`，
  两个不同名称可能清洗到同一文件名（如 `a/b` 与 `a_b`），会静默合并到同一文件。
- `fa one` 未找到序列时返回错误，但输出文件已先用 `File::create` 打开，会留下
  空文件（`-o`）。属一般的"出错留空文件"行为，非数据损坏，`-o` 与输入重叠已被
  `ensure_outfile_distinct` 前置拦截。
- `fa replace` 的 `read_replace_tsv` 不跳过 `#` 注释行，而 `read_names` 会跳过：
  若 TSV 带 `#` 表头会被误当作 key。命令文档未承诺注释支持，低风险，未改。
- `fa masked --gap` 的 `is_n` 会把 IUPAC 歧义码一并计为 N，帮助文本字面略窄于
  实现。与代码库"歧义码即 N"的统一口径一致，属术语精度问题。
- `fa split about -c N -o stdout`：`about` 模式按 `max_files`（默认 999）轮转并
  `break`，在 stdout 上仅输出前 `max_files` 分块容量内的记录。默认 `-c` 为
  `usize::MAX` 不受影响；流式 stdout 与按大小分块本就矛盾，极端组合。
- `read_u32_vec` 在畸形文件给出超大 `count` 时会尝试分配巨大缓冲（可能 OOM
  abort）。属全局既有模式与畸形输入鲁棒性范畴。

## 修复的缺陷（根因模式）

### Zero-Panic / 越界 / 溢出（fa）

- **`six-frame` 短序列 panic**：`&dna[frame..]` 在序列短于 frame 时 slice 越界；
  `dna_len - frame` 在 frame 超出序列长度时下溢。修复：`dna.get(frame..).unwrap_
  or(&[])`；`orfs.is_empty()` 提前返回。
- **`window` 窗口长度过大导致 usize 溢出**：`start + len` 在 `len` 接近
  usize::MAX 时溢出。修复：`start.saturating_add(len)`。
- **`mask_sequence` 对多字节 UTF-8 序列 panic**：以 `&str` 切片按字节偏移操作，
  遇多字节字符触发 char 边界 panic。修复：重写为直接操作 `&[u8]` 字节。

### 数据安全（`-o` 同输入保护，fa/2bit）

- **流式命令允许 `-o` 覆盖输入文件**（先打开输出截断、后读输入，静默清空数据）。
  修复：`count`/`size`/`n50`/`masked`/`one`/`some`/`order`/`rc`/`replace`/
  `range`/`window`/`dedup`（含 `--dups-file`）/`filter`/`to-2bit`/`gz` 及
  全部 5 个 `2bit` 子命令（含 `some` 的 list、
  `range` 的 rgfile）统一加入 `ensure_outfile_distinct`。`six-frame` 初漏，后补。
- **`split`/`window --chunk-records` 输出为目录/派生文件，可覆盖输入**。修复：
  在 `gen_fh`/`create_writer` 生成输出路径时用 `same_path` 与输入比对，命中即
  `bail!`（目录输出无法复用单文件的 `ensure_outfile_distinct`）。
- **`range` 的 `-o` 可覆盖 `.loc` 侧车索引**：先开 writer 后在 `open_indexed`
  读取 `infile.loc`；`-o` 名为 `infile.loc` 会先截断索引，随后 `loc_is_fresh` 因
  新 mtime 判定"新鲜"不重建，`load_loc` 读空索引，且 `.loc` 被永久损坏。修复：
  将 `infile.loc` 一并加入保护列表。
- **`gz` 的 `-o stdout` 创建字面文件 `stdout`** 而非写标准输出（未走 `io::writer`）。
  修复：`outfile == "stdout"` 时改向真实标准输出写 BGZF 流，并跳过无法生成的
  `.gzi` 索引。
- **`gz` 对目录输入会先建残缺输出文件**（`File::open` 对目录成功、错误延迟到首读，
  此时输出已创建）。修复：压缩分支打开输出前对非 `stdin` 输入做 `is_dir()` 前置
  检查。

### 输入校验 / 静默错误（fa）

- **`fa read_names` 未跳过 `#` 注释行**（`some`/`order`/`range` 等），修复：在
  `read_names` 跳过空行与 `#` 开头行。
- **`fa gz --reindex` 静默忽略 `-o`**（reindex 分支直接对 `infile` 建索引）。
  修复：reindex 分支检测到 `-o` 即报错（reindex 输出位置固定为输入旁）。

### 统计 / 算法正确性（fa）

- **`n50` 的 Nx 边界条件用 `>` 而非 `>=`**：累计长度恰好等于 goal 时顺延到下一
  片段。修复：改 `cumul_size >= *goal`。
- **`fa size --no-ns` 对 `-`/`*`/数字等非 IUPAC 字符计数为"有效碱基"**。修复：
  `--no-ns` 排除 N + IUPAC 歧义码 + Invalid，与"仅有效碱基"语义一致。
- **`rc` / `range` 负链对非 IUPAC 字符报错而非保留**：noodles `Sequence::
  complement()` 遇 `-`、`*` 报 invalid base，与文档 "Non-IUPAC preserved" 不符，
  且正/负链输出不一致。修复：改用 `NT_COMP` 表逐字节补全（反向遍历，未知字节
  NT_COMP=255 原样保留）。实测 `ACGT-*NR-` → `-YN*-ACGT`。

### 文档一致性（一次性小修，已精简）

`fa some/order/mask` 的 gzip 输出示例误导（`io::writer` 写端不压缩，压缩由专门
`fa gz` 子命令负责）→ 示例改 `-o output.fa` 并仅指明输入可为 gzipped；
`fa --no-ns`、`fa six-frame` 帮助文本与行为对齐。

## 结论

`fa`（18 子命令）与 `2bit`（5 子命令）命令族合计修复 23 处缺陷（Zero-Panic 3、
数据安全 8、输入校验/静默错误 2、文档一致性 4、统计正确性 2、算法正确性 2），
关键修复均附回归测试与文档澄清。经多轮纵深复审收敛，未再发现新问题。
