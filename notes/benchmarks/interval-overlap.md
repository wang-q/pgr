# 区间结构三方案基准：IntSpan vs coitrees vs rust-lapper

> 目的：回答"判断某个点是否落在基因组区间内"以及"枚举与窗口重叠的区间"
> 时，三种候选结构的取舍。数据生成镜像 rust-lapper 官方基准
> （`rust-lapper-master/benches/lapper_benchmark.rs`）：100 Mb 染色体、
> 区间长度 500..80 kb，另加一条覆盖 90% 染色体的超长区间作为病态用例。

## 三种结构（语义不同）

* `libs::ds::IntSpan`：合并后的 runlist 集合——只回答"点是否被覆盖"
  （`contains` 二分），丢弃区间身份；构造即做区间合并。
* `coitrees::BasicCOITree`：COITree 区间索引（PAF/pbit 已在用），枚举所有
  重叠区间，最坏情况有保证。
* `rust-lapper`：按起点排序 + 二分 + 最长区间补偿；`find`/`count`（BITS），
  外部 intspan 项目的 spanr coverage / rgr count 在用。

## 环境与执行

- `criterion 0.5.1`，`--measurement-time 2 --warm-up-time 1`，样本 100
- `rust-lapper 1.3.0`（dev-dependency，仅基准用）、`coitrees 0.4.0`
- 随机数据固定种子（`benches/interval_overlap_benchmark.rs` 内 SEED）

```bash
cargo bench --offline --bench interval_overlap_benchmark
```

## 结果（median，n = 1k / 10k / 100k 区间）

### 构造

| 方案 | 1k | 10k | 100k |
| --- | ---: | ---: | ---: |
| intspan add_pair（合并） | 0.047 ms | 1.066 ms | 2.417 ms |
| coitrees new | 0.028 ms | 0.323 ms | 6.304 ms |
| lapper new | 0.016 ms | 0.501 ms | 6.677 ms |

### 点成员查询（n 个随机点；用户核心场景）

| 方案 | 1k | 10k | 100k |
| --- | ---: | ---: | ---: |
| intspan contains（普通/超长） | 0.010 / 0.008 ms | 0.084 / 0.050 ms | 0.197 / 0.196 ms |
| coitree point query（普通/超长） | 0.017 / 0.034 ms | 0.694 / 0.902 ms | 27.975 / 31.905 ms |
| lapper count（普通/超长） | 0.015 / 0.014 ms | 0.803 / 0.804 ms | 11.980 / 12.098 ms |

### 区间重叠枚举（n 个 2 kb 查询窗口）

| 方案 | 1k | 10k | 100k |
| --- | ---: | ---: | ---: |
| coitree query（普通/超长） | 0.017 / 0.035 ms | 0.699 / 0.905 ms | 27.914 / 32.296 ms |
| lapper find（普通/超长） | 0.008 / 0.200 ms | 0.622 / 15.952 ms | 22.164 / **1570.356 ms** |

## 结论

1. **纯点成员查询用 IntSpan**：100k 时比 coitree 快 ~140x、比 lapper 快
   ~60x（合并后 span 数远小于区间数，`contains` 为二分）；代价是只有
   "是否覆盖"而没有区间身份。超长区间对它无影响。
2. **重叠枚举（需要区间身份）**：普通数据 lapper 最快（~1.2x），但
   rust-lapper README 自己警告的病态场景确实发生——一条覆盖 90% 染色体的
   区间使 lapper find 在 100k 时从 22 ms 恶化到 **1.57 s**（~71x）；
   coitrees 查询时间几乎不受影响。
3. **构造**：三者量级相当；100k 时 intspan 反而最快（add_pair 顺带合并）。

建议：pgr 的点成员 / 掩码路径维持 IntSpan；需要枚举重叠区间的生产代码
（PAF/pbit 索引）维持 coitrees（最坏情况有保证）；若某数据集确认不存在
超长区间，可考虑用 rust-lapper 换查询速度。
