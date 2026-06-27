# pgr 泛基因组第一步 — 实施记录

本文档记录第一步（V1）的实现过程和最终状态。路线决策见 [[paf-route.md]]，
代码级参考见 [[paf-implementation.md]]，图构建设计见 [[graph-design.md]]。

---

## 1. 完成后状态

三个命令全部实现，98 个测试通过：

```
pgr maf to-paf          ✅ MAF → PAF 转换
pgr paf index           ✅ 区间树索引（支持多文件合并 build_multi）
pgr paf query <region>  ✅ 坐标投影 + --transitive BFS
```

对应 impg 的能力栈：

```
索引层 ✅  |  查询层 ✅  |  图构建层 ← 见 graph-design.md |  应用层 ← 远期
```

## 2. 关键决策记录

| 决策 | 结论 | 理由 |
|------|------|------|
| 索引数据结构 | `IndexMap<String, u32>` + `HashMap<u32, Coitree>` | 复用 pgr 现有模式 |
| 区间树节点 | `PafMetadata`（target/query 坐标 + CIGAR） | 只存 u32 ID，不存序列名 |
| CIGAR 存储 | `Vec<CigarOp>` 全量存入节点 | MAF 转换产物规模可控 |
| PAF 解析 | 手写 tab 分割 | pgr 已有 csv 依赖但手写更轻 |
| 多文件索引 | `build_multi` 直接合并 | 不需要 impg 的 ForestMap/RwLock/ImpgIndex |
| 不输出 FAS/MAF | 只输出 PAF | Fas/MAF 是 MSA 产物，不是坐标投影产物 |

## 3. 源码结构

```
src/libs/paf/
├── cigar.rs        # CigarOp bit-packing + identity + stats
├── record.rs       # PafRecord（String 字段 + tags）
├── writer.rs       # PAF 行格式化
├── parser.rs       # 纯文本 PAF 解析
├── index.rs        # PafIndex + PafMetadata + SortedRanges
└── persist.rs      # .paf.idx 磁盘持久化

src/cmd_pgr/paf/
├── index.rs        # pgr paf index（含 -o 输出）
├── query.rs        # pgr paf query（--transitive）
└── mod.rs

src/cmd_pgr/maf/
└── to_paf.rs       # pgr maf to-paf
```

## 4. 测试覆盖

| 套件 | 测试数 |
|------|:------:|
| `paf::cigar` | 28 |
| `paf::parser` | 10 |
| `paf::index` | 16 |
| `paf::persist` | 12 |
| `cli_maf` | 5 |
| `cli_paf` | 23 |
| **总计** | **98** |

## 5. 变更日志

| 日期 | 内容 |
|------|------|
| 2026-06-27 | 三个命令全部实现（65 tests） |
| 2026-06-27 | 方向 A：持久化（12 tests） |
| 2026-06-27 | 方向 B：真实数据验证（3 tests） |
| 2026-06-27 | 覆盖率提升（5 tests） |
| 2026-06-27 | 方向 C：多文件索引（build_multi） |
| 2026-06-27 | 查询层打磨：--min-identity/--min-output-len/--merge-distance/--subset-sequence-list/--bed/--paf |
| 2026-06-27 | BED 成为默认输出，删除 TSV |
| 2026-06-28 | 文档整合：精简为实施记录 |
| 2026-06-28 | BED/TSV 删除，只输出 PAF；覆盖率补充（+4 tests） |
