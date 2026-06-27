# pgr 泛基因组第一步 — 实施记录

本文档记录第一步（V1）的实现过程和最终状态。路线决策见 [[paf-route.md]]，
代码级参考见 [[paf-implementation.md]]，图构建设计见 [[graph-design.md]]。

---

## 1. 完成后状态

三个命令全部实现，104 个测试通过：

```
pgr maf to-paf          ✅ MAF → PAF 转换
pgr paf index           ✅ 区间树索引（支持多文件合并 build_multi）
pgr paf query <region>  ✅ 坐标投影 + --transitive BFS + -o bed/paf + -b 批查
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
| 不输出 FAS/MAF | 只输出 PAF/BED | Fas/MAF 是 MSA 产物，不是坐标投影产物 |
| 默认输出 | PAF（非 impg 的 BED 默认） | 既有测试断言 PAF；PAF 含完整 CIGAR 更直接；BED 是轻量坐标用 `-o bed` 切换。见 [[graph-design.md]] §3.1 |
| 批查 | `-b regions.bed`（BED3 格式） | 对齐 impg `-b`，单/批 region 互斥 |

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
| `cli_paf` | 29 |
| **总计** | **104** |

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
| 2026-06-28 | 决策修订：BED 删除有误（impg 默认即 BED），待恢复为默认输出；见 [[graph-design.md]] §3 |
| 2026-06-28 | V1 实现：`-o bed` 可选 + `-b regions.bed` 批查（+6 tests，共 29）。**默认输出保持 PAF**（非 impg 的 BED 默认），理由：既有 23 测试断言 PAF、PAF 含完整 CIGAR 更直接、BED 是轻量坐标产物用 `-o bed` 显式切换。见 [[graph-design.md]] §3.1 |
| 2026-06-28 | V1 后处理过滤：`--min-degree N`（per-region distinct query 数）+ `--min-chain-length N`（per-query_id 累加对齐长度），+4 tests 共 33。`--end-trim` 推迟（需 per-interval CIGAR 修剪，与区间投影模型不兼容，待 V2）。见 [[paf-implementation.md]] §8、[[paf-route.md]] §4.4 |
| 2026-06-28 | **V4a 粗全局 GFA 物化（路线 B: seqwish DSU）**：新增 `pgr paf graph -f refs.fa --min-var-len 100`，输出 GFA v1.0（S/L/P）。`src/libs/paf/graph.rs` 470 行引擎（CIGAR 切分 → 段对 → DSU 传递闭包 → 节点序列 → 路径+novel 段补全 → 边派生）+ `src/cmd_pgr/paf/graph.rs` CLI 包装，5 单元 + 7 集成测试（共 40）。简化项：无 disk-backed interval tree / SparseBitVec / lock-free DSU；路径方向恒 `+`（反向已翻转坐标）；rGFA SN/SO/SR tag 暂缺（待兼容性需要再加）。见 [[graph-design.md]] §4.3.1、[[seqwish.md]] §6.2 |
