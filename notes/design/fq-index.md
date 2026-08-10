# pgr fq index：按 read name 的随机访问（设计稿）

> 状态：**一期已实现（2026-08）**——单文件 `pgr fq range`（API 与 `fa range`
> 对齐，含 name 归一化与交错 `#n` 消歧）；双端 S2 待二期。
>
> 配套：[seq-reader.md](seq-reader.md)（FAFQ 读取与 BGZF 基础设施）、
> [anchr-trim-replace.md](anchr-trim-replace.md)（fq 命令组现状，含 trim-qual）。

## 0. 动机

FASTA 侧已有按 name 的 `.loc` 索引（`src/libs/loc.rs`）：文本格式
`name\t明文偏移\trecord_size`，`fa range` 自动构建/复用，BGZF 文件经
`CachedBgzfReader` 按明文偏移 seek（依赖 `.gzi`）。用户希望 FASTQ 获得同样的
能力：**从超大 FASTQ 按 read name 随机提取子集**，而不是每次全文件流式过滤。

## 1. 核心问题：API 是否与 `fa range` 完全对齐

`fa range` 的形态：`pgr fa range in.fa "chr1:1-1000" "chr2(-):2000-3000"`，
自动建 `.loc`，支持 `-r`（范围文件）/`-c`（LRU 缓存）/`-u`（强制重建索引）/
`-o`。语义：`name`（整条）或 `name:start-end`（子段），`(-)` 负链反向互补。

两个方案：

| 方案 | 形态 | 收益 | 成本 |
|---|---|---|---|
| **A（推荐）：完全对齐** | `pgr fq range in.fq "read1:10-200" "read2"`，参数与 `fa range` 一致 | API 与 FA 一模一样，脚本/心智统一；长读（ONT/PacBio）可取 read 内部片段 | 子段提取需同时切 seq+qual；`.loc` 之外无额外索引信息 |
| B：最小 | `pgr fq index` + `pgr fq fetch`，只按 name 整条提取 | 实现最小 | API 与 FA 不一致，未来对齐要再改 |

**关键洞察：索引格式可以零成本对齐**。`.loc` 存的是"记录级"的
（明文偏移, record_size），FASTQ 的 4 行记录整体即一个 record；取子段时
`fetch_record` 先取整条再切分（`fa range` 正是 fetch 整条后 `slice_record`）。
所以方案 A 不需要在索引里多存任何字段。

**推荐 A**，理由：
1. 索引与查询架构直接复用 `loc.rs` 骨架，A 与 B 的差异只在"取到整条后要不要
   支持 `name:start-end` 切片"——增量极小。
2. 用户明确看重"API 一模一样"；对齐后 FA/FQ 脚本可互替。
3. read 子段对长读数据有真实用途（截取 read 内部区域），不是空泛灵活性。

## 2. `.loc` 格式与构建（FASTQ 版）

格式不变：`name\t明文偏移\trecord_size`（每 read 一行）。构建差异只在
**记录边界判定**：

- FASTA：遇 `>` 行首即新记录。
- FASTQ：4 行结构（`@name` / seq / `+` / qual），且质量行也可能以 `@` 开头，
  必须按结构逐条解析（复用 `SeqReader` 的解析逻辑，另累计明文 offset 与
  记录字节数），不能只看行首。
- `+` 行内容（空或重复 name）不影响：按行累计字节数即可。
- 折行 FASTQ（序列/质量多行）需在解析时处理，与 `SeqReader` 一致。

明文 / 普通 gzip / BGZF 三种输入：
- 明文、gzip：`GzReader` 流式解压扫描，`.loc` 存明文偏移。
- BGZF：同样流式扫明文偏移（构建侧）；查询侧 `CachedBgzfReader` + `.gzi`
  按明文偏移定位块——与 `fa range` 完全同路。

## 3. CLI 设计（方案 A）

```
pgr fq range in.fq "read1:10-200" "read2" ...

Options（与 fa range 一致）：
  -r, --rgfile    从文件读 name/range 列表
  -c, --cache     LRU 缓存容量（默认 1）
  -u, --update    强制重建 .loc 索引
  -o, --outfile   输出文件（默认 stdout）
```

语义：
- `name` → 输出整条 4 行记录（保留 `+` 行原内容）；
- `name:start-end` → 同时切 seq 与 qual，输出 4 行记录（`+` 行输出原内容或
  名称与切片一致？——见 §4）；
- `(-)` 负链：**不支持**（FASTQ 的 read 无基因组方向坐标语义，反向互补对
  read 切片无意义；若未来长读场景需要再补）。

## 4. 差异点与待定决策

| 主题 | FA 现状 | FQ 处理 | 状态 |
|---|---|---|---|
| 索引格式 | `name\t偏移\tsize` | 完全一致 | 定稿 |
| 子段语法 | `chr1:1-1000` | 复用 `ds::Range`，`read1:10-200` | 方案 A 定稿 |
| 负链 `(-)` | 反向互补 | 不支持（切片时忽略 strand） | 已定 |
| 子段输出 `+` 行 | — | 输出单个 `+`（与 trim-qual 一致） | 已定 |
| 重复 name | IndexMap 后者覆盖 | 归一化 key + `#n` 消歧（不丢数据） | 已定 |
| 索引重建 | mtime 判断（`loc_is_fresh`） | 同 | 定稿 |
| 缓存 | LRU<FastaRecord> | LRU<原始记录字节 Vec<u8>> | 定稿 |
| 普通 gzip | 不支持（需 BGZF） | 明确报错"only plain text and BGZF" | 已定 |

## 5. 双端与交错（重点，易错区）

### 5.1 命名模式与索引冲突

FASTQ 双端数据的 name 有三种常见模式：

| 模式 | R1 记录名 | R2 记录名 | 单文件 `.loc` 的 name key |
|---|---|---|---|
| CASAVA 老式（`/1` `/2`） | `read1/1` | `read1/2` | 不同 key，无冲突 |
| CASAVA 1.8+（description 区分） | `read1 1:N:...` | `read1 2:N:...` | 都是 `read1`（name 取首段） |
| 无后缀同名 | `read1` | `read1` | 都是 `read1` |

- **分离文件（R1.fq / R2.fq）**：各自建 `.loc`，单文件内 name 唯一，无冲突。
- **交错文件**：相邻 pair 的 name 规范化后相同 → 同一文件内同名多条，
  `IndexMap` 后者覆盖 → **静默丢一条**。必须消歧。

### 5.2 用户视角的核心问题

"取 read1，那 read1 的另一端（R2）怎么办？"——三个层面：

1. **分离文件**：用户要么在 R1.fq、R2.fq 各跑一次（相同 name 列表），要么
   命令层同时接受两个文件、一次取 pair；
2. **name 归一化**：`read1/1` / `read1/2` 应归一为 pair name `read1`
   （strip `/1` `/2`）；CASAVA 1.8+ 的 name 首段已一致，无需额外处理；
3. **交错文件**：一次查询应同时返回 pair 的两条（保持交错顺序）。

### 5.3 API 方案

| 方案 | CLI 形态 | 说明 |
|---|---|---|
| **S1 单文件（fa 对齐，最小）** | `fq range R1.fq "read1/1"` | 一次一个文件；配对由用户自行在 R2.fq 重复；交错文件需消歧 |
| **S2 双端感知** | `fq range R1.fq R2.fq "read1" -o r1.fq --outfile-2 r2.fq`；交错输入单文件 → 交错输出 | pair name 归一化；对齐 `fq trim-qual` 的双端输出形态 |

推荐 **S2**（符合"取 pair"直觉），分两期落地：

- **一期**：单文件 + 交错消歧。`fq range in.fq "read1"` 只作用于一个文件；
  双端用户用脚本对两个文件跑同一条命令。
- **二期**：双端感知。两个输入文件 + `-o/--outfile-2` 分离输出；交错输入
  单文件输出保持交错。

### 5.4 交错文件索引消歧

同 key 多条时，索引 key 追加序号：`read1#0`、`read1#1`（保持出现顺序）。
查询 `read1` 时返回精确 key + 全部 `#n` 变体（`read1`、`read1#1`）。效果：

- 不丢数据；
- 交错文件按名提取返回 pair 两条、顺序保持；
- 分离文件无影响（key 唯一，`#n` 不出现）。

`#n` 是内部消歧后缀（第一条不带 `#`，后续为 `#1`、`#2`...），用户查询时
无需写。CASAVA 1.8+ 的 `1:N:`/`2:N:` 不专门
识别（name 首段已一致，消歧兜底）。

### 5.5 决策点

1. ~~一期做 S1 还是 S2~~ → 已按 S1 实现（一期），S2 待二期。
2. ~~`#n` 消歧是否可接受~~ → 已按 `#n` 实现并验证。
3. 双端输出（二期）是否严格对齐 `fq trim-qual` 的 `-o/--outfile-2` 形态？

## 6. 测试计划（实现阶段）

- 单元：FASTQ 记录边界扫描（普通/折行/质量行含 `@`）、offset/size 累计、
  子段切分 seq+qual、pair name 归一化（strip `/1` `/2`）、交错同名消歧。
- 集成 `tests/cli_fq_range.rs`：明文/gzip/BGZF 建索引与提取、`name` 与
  `name:start-end`、`-r`/`-u`、不存在 name 报错、重复 name、索引自动重建、
  交错文件 pair 提取（两条同时返回）、双端分离提取（S2 时）。
- 吞吐 sanity：50 MB FASTQ 建索引耗时 + 按名提取小集合的耗时（记录到笔记，
  不建 criterion 基准）。

## 7. 二期（未做）

双端感知 S2：两个输入文件 + `-o/--outfile-2` 分离输出；交错输入单文件
输出保持交错。索引层已就绪（每文件独立 `.loc`、name 归一化），增量在命令层。

## 8. 一期实现记录（2026-08）

- 代码：`src/libs/loc.rs`（`normalize_pair_name`/`create_fq_loc`/
  `open_fq_indexed`/`query_fq_locs`）、`src/cmd_pgr/fq/range.rs`（参数与
  `fa range` 一致：infile + ranges + `-r`/`-c`/`-u`/`-o`）。
- 行为：4 行结构扫描建 `.loc`；name 归一化（strip `/1` `/2`），同 key 多条
  追加 `#n`；查询返回精确 key + 全部 `#n` 变体（交错/合并文件 pair 两条
  同时返回、保持顺序）；`name:start-end` 同时切 seq 与 qual，`+` 行输出
  单个 `+`；BGZF 复用 `.gzi`；普通 gzip 明确报错。
- 验证：lib 单元 4 个 + 集成 9 个（明文/子段/BGZF/交错 pair/`/1` `/2`
  归一化/缺失名警告/输出同名拒绝/FASTA 输入报错/普通 gzip 报错）；
  clippy clean；全量 lib 688 + 集成套件全绿。
- 吞吐 sanity（release，103 MB / 50 万条 100 bp 明文 FASTQ）：
  - 建索引 + 提取 1 条：0.43 s；`.loc` 12 MB（50 万行 × ~24 B）；
  - 索引已存在时 `-r` 提取 100 条（分散位置）：0.14 s。

---

*参考来源: [fa range](../../src/cmd_pgr/fa/range.rs) | [loc.rs](../../src/libs/loc.rs) | [seq-reader.md](seq-reader.md)*
