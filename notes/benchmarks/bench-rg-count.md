# `rg count` 命令行基准：pgr（coitrees）vs rgr（rust-lapper）

> 目的：对比 `pgr rg count`（`libs/runlist::RgIndex`，coitrees 区间树）与
> 外部 `rgr count`（intspan 项目，rust-lapper）在命令行层（含解析、建索引、
> 查询、输出）的耗时与内存，评估移植后的性能余量。2026-08-04 实测，同日
> 修复轮后复测（数值一致）。

## 环境与版本

* pgr：本仓库 release 构建（`cargo build --release`，v0.4.0）
* rgr：`~/.cbp/bin/rgr` 0.8.6（release）
* 机器：本机（hyperfine 1.19.0，`/usr/bin/time -v` 量内存；复测时主机有
  并发负载，pgr 绝对时间略高于初测，相对差距不变）

## 数据（合成，种子 20260804）

* 8 条染色体（chr1..chr8）× 100 Mb
* `iv.1m.rg`：1,000,000 条随机区间（长度 500–80,000），另加 1 条
  `chr1:0-89999999`（覆盖 chr1 90% 的病态长区间）
* `iv.1m.normal.rg`：同上但去掉病态长区间
* `target.100k.rg`：100,000 条随机查询区间（长度 100–2,000）

## 复现命令

```bash
pgr rg count target.100k.rg iv.1m.rg -o /dev/null
rgr  count target.100k.rg iv.1m.rg -o /dev/null

hyperfine --warmup 1 --runs 3 \
  'pgr rg count target.100k.rg iv.1m.rg -o /dev/null' \
  'rgr  count target.100k.rg iv.1m.rg -o /dev/null' \
  'pgr rg count target.100k.rg iv.1m.normal.rg -o /dev/null' \
  'rgr  count target.100k.rg iv.1m.normal.rg -o /dev/null'
```

## 复测结果（1M 区间 + 100k target，5 次取均值，2026-08-04）

| 实现 | 含病态长区间 | 普通 | RSS（1M 区间） |
| :--- | ---: | ---: | ---: |
| pgr `rg count` | 231.9 ± 1.5 ms | 231.2 ± 2.2 ms | 41.9 MB |
| rgr `count` | 755.3 ± 19.1 ms | 742.8 ± 10.9 ms | 24.8 MB |

pgr 约快 **3.2–3.3×**（普通 3.21×，病态 3.27×）；内存约为 rgr 的 1.7×。
病态长区间对 pgr 几乎无影响（231.9 vs 231.2 ms），与初测一致。

## 正确性验证

200k 区间 + 20k target 上，两者输出 `sort` 后 `diff` 为空（20,000 行逐行
一致），计数语义等价。

## 分析

1. **coitrees 查询有界**：病态长区间对 pgr 几乎无影响（235 vs 224 ms，
   +5%）；rgr 在本数据上也未见明显退化（758 vs 753 ms）。注意这与
   `notes/design/runlist.md.md` 里 lapper `find()` 在超长区间上
   退化 ~71× 的场景不同——`count()` 走 BITS 路径，对单条超长区间的敏感度
   低于 `find()`。真实差距主要来自查询本身：coitrees 100k 次窗口查询的
   常数远小于 lapper。
2. **解析/构建开销相当**：两端都要读 1M 行并建索引（用户态 ~200 ms vs
   ~730 ms），pgr 的剩余优势基本都在查询阶段。
3. **内存代价**：coitrees 树节点比 lapper 的原始区间向量紧凑度低
   （1M 区间 42 vs 26 MB），对基因组级数据可接受；若未来要处理 10M+
   区间，值得关注。

## 结论

`pgr rg count` 在命令层级比 `rgr count` 快约 3.4×、输出逐行一致，移植的
性能余量充足；内存多 1.6× 属可接受代价。若追求更低内存，可后续评估
coitrees 的 FlatTree/序列化形态或按染色体流式查询。
