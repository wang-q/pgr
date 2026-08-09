# pgr fq index：按 read name 的随机访问（设计稿）

> 状态：**设计稿（未实现）**。核心问题是"是否与 `fa range` 的 API 完全对齐
> （支持按 read 取子段）"，待用户确认后定稿。
>
> 配套：[seq-reader.md](seq-reader.md)（FAFQ 读取与 BGZF 基础设施）、
> [fq-trim-q.md](fq-trim-q.md)（fq 命令组现状）。

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
| 负链 `(-)` | 反向互补 | 不支持 | 待用户确认 |
| 子段输出 `+` 行 | — | 建议保持原 `+` 行内容（切段后质量对齐 seq） | 待定 |
| 重复 name | IndexMap 后者覆盖 | 同（默认覆盖，必要时警告） | 待定 |
| 索引重建 | mtime 判断（`loc_is_fresh`） | 同 | 定稿 |
| 缓存 | LRU<FastaRecord> | LRU<FqRecord（4 行）> | 定稿 |

## 5. 测试计划（实现阶段）

- 单元：FASTQ 记录边界扫描（普通/折行/质量行含 `@`）、offset/size 累计、
  子段切分 seq+qual。
- 集成 `tests/cli_fq_range.rs`：明文/gzip/BGZF 建索引与提取、`name` 与
  `name:start-end`、`-r`/`-u`、不存在 name 报错、重复 name、索引自动重建。
- 吞吐 sanity：50 MB FASTQ 建索引耗时 + 按名提取小集合的耗时（记录到笔记，
  不建 criterion 基准）。

## 6. 实施步骤（定稿后）

1. `loc.rs` 泛化：抽取"记录扫描建索引"骨架，新增 FASTQ 版（4 行结构解析）；
2. `fq range` 命令（对齐 `fa range` 参数）；
3. 测试 + sanity；4. 更新 [fq-trim-q.md](fq-trim-q.md) 或本文档为已实现。

---

*参考来源: [fa range](../../src/cmd_pgr/fa/range.rs) | [loc.rs](../../src/libs/loc.rs) | [seq-reader.md](seq-reader.md)*
