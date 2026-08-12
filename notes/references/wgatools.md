# wgatools：Rust PAF/MAF 工具库（源码分析）

> 2026-08-12 整理，纯源码分析（`wgatools-1.1.0/`，Rust，BSD-2）。wgatools
> 是**全基因组比对格式（PAF/MAF）工具库**：解析 + 12+ 个子命令（stat/
> filter/index/pafcov/pileup/trimovp/dotplot/mafextra/rename/chunk/
> tview/validate/caller 等）。**与 pgr 的关系**：pgr `paf` 模块是自研 PAF
> 隐式图（index/query/to-gfa/to-vcf/stat/validate），wgatools 是最接近的
> 外部参照，可补 pgr 缺失的 PAF 通用操作。

## 1. 工具对照（pgr paf vs wgatools）

| wgatools | 作用 | pgr 现状 |
|---|---|---|
| `stat` | PAF 统计 | ✅ `paf stat` |
| `validate` | PAF 校验/修复 | ✅ `paf validate`（2026-08-06） |
| `index` | PAF 索引 | ✅ `paf index`（区间树 + `.paf.idx`） |
| `pafcov` | PAF → 覆盖度 | ✅ `pgr paf coverage`（2026-08-12 新增，cg:Z 扫描线） |
| `filter` | 按长度/身份过滤 PAF | ❌ pgr 无 PAF 过滤子命令 |
| `trimovp` | 修剪 overlap 末端 | ❌ pgr 无（`paf to-fas`/query 未做 per-interval 修剪，`todo.md` §2） |
| `rename` / `chunk` / `tview` / `mafextra` / `pileup` / `caller` | 命名/分块/查看/MAF 提取/堆叠/变异 | 部分与 `paf to-fas/to-maf/to-vcf`、`plot dot` 重叠 |

## 2. 借鉴点

- **`pafcov` 的 CIGAR → 覆盖度扫描线**（`tools/pafcov.rs` + `parser/cigar.rs`
  `update_cov_vec`）：PAF 带 `cg:Z` 时直接累积覆盖度向量——已落地为
  `pgr paf coverage`（2026-08-12，`libs/paf/cov.rs` 扫描线 + 恒定深度段
  合并），无需再经 SAM 中转。
- **`trimovp`**：overlap 末端修剪正是 `paf-pangenome.md` §"`--end-trim`
  推迟（需 per-interval 修剪 CIGAR）"的待办——wgatools 有现成语义可对照。
- **`filter`**：pgr `paf` 模块缺通用过滤子命令，`paf query`/`to-bed` 的
  过滤维度（`--min-tree-coverage` 等，`todo.md` §2）可参考 wgatools 的
  参数形态。
- **解析器容错**：wgatools 的 PAFReader 处理 tag/CIGAR 的边界情况丰富，
  pgr `paf::parser`/`cigar` 可作对照测试源（不引依赖，仅参考）。

> 结论：价值中等——pgr `paf` 已覆盖核心（图/索引/查询），wgatools 的
> 增量是 `filter`/`trimovp`/`pafcov` 三个通用操作与解析器对照。按
> `todo.md` §2 的优先级（`--min-tree-coverage`、`--end-trim`）推进时
> 再细读对应文件即可，暂不立项。
