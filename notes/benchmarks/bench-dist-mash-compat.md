# `dist mash` 与 Mash 2.3 字节级兼容对照

> 目的：验证 `pgr dist mash`（bottom-k MinHash）与参考工具 Mash（Ondov et
> al. 2016，本地 `Mash-master` + 系统 mash 2.3）在相同 k / sketch size 下
> 的距离输出完全一致。日期：2026-08-08。

## 算法核对结论（Mash-master 源码）

- canonical k-mer：正链与反向互补**字节比较取小**（`memcmp`），不是哈希取小；
- 哈希：`MurmurHash3_x64_128(kmer_bytes, seed=42)` 取前 64 位；
- 过滤：k-mer 含非 ACGT 直接跳过（整 k-mer 不参与）；
- bottom-k：全局保留 sketchSize（默认 1000）个最小唯一哈希；
- 距离：合并排序 → shared / denom（denom clamp 到 sketchSize）→
  Jaccard → `-ln(2J/(1+J))/k` clamp [0,1]。

**关键坑 1（Jaccard）**：Mash 的 Jaccard = `compareSketches` 的匹配数 /
sketchSize（合并两个排序 sketch、最多 sketchSize 步的匹配数），不是
标准集合 Jaccard。Mash-master/test 的 genome1×genome2 完整交集 581 个
哈希，但 Mash 报告 456/1000。`dist mash` 严格按 Mash 语义实现，
与 `mash dist` 的 shared/denom 完全一致（20/20 对，见对照 2）。

**关键坑 2（Containment，2026-08-08 发现并修复）**：`dist mash` 原实现
的 containment 用了与 Jaccard 相同的 merge-common（`common / a.len()`），
当集合满 sketchSize 时数值恰好等于 Jaccard——这不是标准 containment。
系统对照（5 株 × 20 对）显示它比标准值（完整集合交集 / query 集合大小，
Mash `within` 语义）系统性低估 0.12–0.18（相对偏差约 25–35%）。已修复：
containment 改为完整集合交集 / 第一个集合大小；修复后 20/20 对与标准值
完全一致（Δ=0.0000），Jaccard/距离不受影响。

**关键坑 3（Undersized sketch 的 denom，2026-08-08 审计修复）**：Mash
`compareSketches` 的 Jaccard 分母是 merge 遍历的 `denom`（一方提前耗尽时
补上剩余未遍历哈希，clamp 到 sketchSize），不是固定 sketchSize。原实现
除以 sketchSize，导致 sketch 未满（短 contig/质粒、高 k、小序列）时距离
错误——两个相同 46-hash sketch（k=15/s=1000）Mash 报 46/46、距离 0，
pgr 报 46/1000、距离 0.1621。已修复：分母与 `union` 输出均改为 Mash 的
denom 语义；修复后小 sketch 与 Mash 一致（2/2 → 0、0/92 → 1），满 sketch
不受影响。

**内存（2026-08-08 审计修复）**：草图构建原为全量物化（每条序列的完整
哈希 Vec + 缓冲、`--merge` 累积整个文件），内存 O(基因组长度)（4.6 Mb
单条 record 峰值 RSS 约 90 MB，100 Mb 基因组约 1 GB）。已改为滚动窗口
流式（O(k) 窗口）+ 增量 bottom-k 累积器（O(sketch_size)）：4.6 Mb 基因组
merge 模式峰值 RSS 降到 15.8 MB，且不随序列长度增长。流式实现与旧全量
逻辑逐哈希对照一致（单元测试 `test_for_each_mash_hash_matches_reference`）。

## 对照 1：Mash-master/test 三对（k=21, s=1000）

| 对 | mash dist | pgr dist mash | shared |
|---|---:|---:|---:|
| genome1 × genome2 | 0.0222766 | 0.0223 | 456/1000 |
| genome1 × genome3 | 0.0000000 | 0.0000 | 1000/1000 |
| genome2 × genome3 | 0.0222766 | 0.0223 | 456/1000 |

## 对照 2：5 株 E. coli 真实基因组 20 对（k=21, s=1000，整文件 `--merge`）

数据：`/tmp/pgr_cohort/data/{mg,sa,se,e2,cf}.fa.gz`（MG1655 / Sakai / SE11 /
E2348 / CFT073）。pgr 侧用 `--merge`（与 `mash dist` 的整文件语义一致）；
不加 `--merge` 时 pgr 是逐 contig 比较，行数会不同，注意区分。

| 对 | mash dist | pgr dist mash | diff |
|---|---:|---:|---:|
| mg×sa | 0.0227107 | 0.0227 | <1e-4 |
| mg×se | 0.0162749 | 0.0163 | <1e-4 |
| mg×e2 | 0.0304804 | 0.0305 | <1e-4 |
| mg×cf | 0.0312752 | 0.0313 | <1e-4 |
| sa×se | 0.0240597 | 0.0241 | <1e-4 |
| sa×e2 | 0.0328290 | 0.0328 | <1e-4 |
| sa×cf | 0.0350415 | 0.0350 | <1e-4 |
| se×e2 | 0.0339141 | 0.0339 | <1e-4 |
| se×cf | 0.0356220 | 0.0356 | <1e-4 |
| e2×cf | 0.0174761 | 0.0175 | <1e-4 |

全部 20 对一致，差异仅为 pgr 输出 4 位小数的舍入。结论：`dist mash`
与 Mash 字节级兼容，可作 Mash 的直接替代（同参数）。

## 性能

`benches/dist_sketch_benchmark.rs`（criterion，本机）：
`dist mash` 草图构建约 33 MiB/s（k=21 单线程）；对比 `dist mini` 约
122 MiB/s、`dist frac` 约 113 MiB/s（语义不同，仅量级参考）。
