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
`dist mash` 草图构建约 54 MiB/s（k=21 单线程、流式 bottom-k + 预筛，
2026-08-08 优化后；优化前 24 MiB/s）；对比 `dist mini` 约 104 MiB/s
（窗口采样，哈希次数少 5 倍）、`dist frac` 约 113 MiB/s（语义不同，
仅量级参考）。Mash C++ 单线程约 66 MiB/s（实测估算）。

### 并行扩展（2026-08-08，5 株 E. coli × 4 query，3 次取最快）

sketch 加载改为文件级并行后（`load_entries` rayon，三个草图命令共用）：

| -p | pgr dist mash（前） | pgr dist mash（后） | mash |
|---:|---:|---:|---:|
| 1 | 1.34 s | 0.57 s | 0.45 s |
| 4 | 1.34 s | 0.23 s | 0.18 s |
| 8 | 1.34 s | 0.24 s | 0.19 s |

单对（2×4.6 Mb）：pgr 0.24 s vs mash 0.18 s（1.3×）。剩余差距来自
纯 Rust murmur3 与 FASTA 读取（Mash 为 C++ 优化），未换哈希实现以保持
字节级兼容。`dist frac` 同样受益（0.70→0.31 s @ -p4）、`dist mini`
（0.53→0.30 s）。

### `dist frac` containment 同 k 对照（2026-08-08）

澄清此前"~10% Jensen 偏差"（`dist-cohort-validation.md`）：那是 frachash
k=21 vs 全 k=40 真值的 **k 不匹配假象**。同 k=21 全 canonical k-mer 集合
真值下，5 株 × 10 对：

| 对 | 真值 containment | frac (scale=1000) | 偏差 |
|---|---:|---:|---:|
| mg×sa | 0.6690 | 0.6591 | −1.5% |
| mg×se | 0.7474 | 0.7347 | −1.7% |
| sa×se | 0.5884 | 0.5840 | −0.7% |
| e2×mg | 0.4976 | 0.5012 | +0.7% |
| e2×sa | 0.5168 | 0.5049 | −2.3% |
| e2×se | 0.5021 | 0.5004 | −0.3% |
| cf×mg | 0.4819 | 0.4861 | +0.9% |
| cf×sa | 0.4875 | 0.4823 | −1.1% |
| cf×se | 0.4844 | 0.4795 | −1.0% |
| cf×e2 | 0.6680 | 0.6772 | +1.4% |

偏差正负对称、幅度 ~1–2%（scale=1000 时采样 SE≈0.7%），属采样方差而非
系统性偏差。Hera 2023 校正因子 (1−(1−s)^|A|) 对 |A|≥10⁵ ≈ 1，大肠杆菌
场景无校正效果；仅极短序列（|A|<~100）+ 大 scale 才有意义。`dist frac`
containment/ANI 保持现状。
