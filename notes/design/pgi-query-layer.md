# PgiQuery 抽象层与 FastGA 顺序算法：解锁路径（设计稿 + 基准计划）

> 2026-08-05 起草。背景：`pgr align pgi` 的 merge 走 `PgiQuery`
> 抽象（resident 索引 / mmap 按需解码），FastGA 的 `new_merge_thread`
> 是"顺序指针 + 预计算 LBYTE + vlcp 表"模型。抽象层把顺序访问降级成
> 二分，导致完整 LCP（§3.6 慢 2.1×）与变长种子（§7.3）等 FastGA 算法
> 无法直接利用。本文给出诊断、候选路径与基准设计。
>
> **2026-08-05 定稿（方案 A 落地）**：§4 方案 A（归并式 merge）已实现并
> 验证——PSL 逐字节一致、merge 耗时不倒退（resident 反快 ~28%）、内存不
> 涨；见 §8 落地记录。
>
> **2026-08-05 复测**：resident 18.6 vs 27.2 ms（反快 ~32%）、流式
> 90-100 ms（§8.2）；`sd search pgi` 两用例由 wave trim 修复
> （commit `cd66774`）；`test_invalid_op_panics` 核对为普通 `panic!`
> （非 `debug_assert!`）。当前全量 `cargo test`（debug + release）
> 1300 个全部通过。

## 1. 问题诊断

### 1.1 现状：`PgiQuery` 抽象

`src/libs/pgi/mod.rs::PgiQuery` 暴露 10 个方法（方案 A 新增
`entry_lower_bound_ge`，§4/§8），两个实现 + 一个流式参考：

| 实现 | 存储 | `entry_range` | `entry_kmer` | `entry_next` |
|---|---|---|---|---|
| `PgiIndex`（resident） | `Vec<PgiEntry>`（24 B：kmer u128 + pos_start/freq u32×2） | `partition_point` 二分 | 数组直接读 | `i+1` |
| `PgiMmap`（query 侧） | mmap 字节流（`RecordLayout`：kmer 按需字节 + 位置按需字节） | `lower_bound` 二分（逐字节比较） | `unpack_kmer` 组装 u128 | `group_end` 顺序推进 |
| `PgiStream`（ref 侧） | 1 MiB 缓冲批解码 `Vec<(PgiEntry, Vec<u64>)>` | —（批内 Vec） | 直接读 | 批内 `i+1` |

`emit_entry_hits` 对**每个 a 条目**：`window(len)` → 定位窗口 → 扫窗口找
最大共享前缀 m → 再 `window(m)` 定位 → 发射。LCP 传播（`prev_kmer`）把
窗口起点提前到 `max(min_shared, lcp(prev, cur))`。方案 A 后窗口定位由
`MergeCursor` 批内顺序推进（大跳变回退二分），替代每条目 `entry_range`
二分（§4/§8）；本文 §1.2-1.4 的"降级"分析针对方案 A 前的二分模型。

### 1.2 对照：FastGA `new_merge_thread`（FastGA.c:610）的顺序模型

1. **cpre 分组**：T1（a 侧）按完整前缀分组，每组把 T2（b 侧）条目载入
   连续 `cache`（`kbyte` 字节/条），组内 `rcur/rend/suf1` **指针推进**；
2. **LBYTE 字段**：k-mer 表条目内嵌与前一条目的共享 lcp（排序后稳定），
   读取 O(1)，不即时计算；
3. **vlcp[plen] 表**：`plen ∈ [12, KMER]` 档位，维护"与当前 a 条目共享
   ≥ plen 的 b 侧窗口起点"，O(1) 摊销更新；`low = vlcp[plen]` 直接取窗口；
4. **频率过滤**：`kfreq = FREQ*kbyte` 用指针 `top = low + kfreq` 推进，
   无二分。

### 1.3 抽象层把什么降级了

| FastGA | pgr（PgiQuery） | 代价 |
|---|---|---|
| 指针顺序推进（O(1)/步） | 每条目二分（O(log n)） | 完整 LCP 尝试 1 慢 5×（543 vs 107 ms） |
| LBYTE 预计算 | `shared_prefix` 即时 XOR+clz | 可忽略（快指令），但无法支撑 vlcp 表 |
| vlcp[plen] 档位窗口 | `prev_kmer` 单档传播 | 完整 LCP 尝试 2 慢 2.1×（220-229 vs 107 ms） |
| cpre 组批处理 | rayon 批内独立二分 | 批边界游标无法复用 |

§3.6 结论"完整 LCP 正确但性能不可行"本质是**抽象层没有提供顺序访问
原语**，不是算法本身慢。

### 1.4 被卡住的算法清单

1. **完整 LCP（vlcp 表 / `.pgi` v3 LBYTE）**——§3.6 已失败、§7.6 暂缓；
2. **变长种子（完整 adaptamer，种子长度 > k）**——§7.3，merge 阶段需要
   顺序延展匹配（连续 kmer 延伸），当前每条目独立二分无法表达；
3. **merge 批处理的窗口复用**——FastGA cpre 组内 O(1) 窗口继承，pgr 批
   内仍逐条目二分。

## 2. 现状量化（基准基线，写入文档供对比）

> 基线为方案 A 前二分模型；当前数值见 §8.2（流式 merge 90-100 ms、
> resident 18.6 ms）。

| 指标 | 当前值（2026-08-05 实测） |
|---|---|
| merge 耗时（流式，E. coli 3 对） | ~90-100 ms（含 a 侧流式读；方案 A，§8） |
| merge 种子数 | 1,121,308（min_shared=12） |
| 完整 LCP（narrow_prefix 正确版） | 220-229 ms（慢 2.1×，§3.6） |
| 逐碱基递增（二分模拟） | 543 ms（慢 5×，§3.6） |
| 峰值 RSS（align 全流程） | 209-210 MB |

**拆分探针（阶段 0，已完成）**：在 `emit_entry_hits` 内分别计时
（a）`entry_range` 二分定位、（b）窗口扫描找 m、（c）`shared_prefix`
计算、（d）发射 positions。用 `examples/merge_mem_bench.rs` 的
split-profile 变体实测（mg1655 × nissle，resident，8 线程 release）：

| 分量 | 累计 CPU | 占比 | 说明 |
|---|---:|---:|---|
| (a) entry_range 二分 | 263.8 ms | **76.9%** | 2,287,166 次调用，0.87 次/条目 |
| (b) 窗口扫描找 m | 37.8 ms | 11.0% | 1,524,811 条目/次扫描 |
| (d) 发射 positions | 41.3 ms | 12.0% | — |
| (c) shared_prefix | — | 可忽略 | 快指令，纳入 (b) 计时 |

**结论**：主战场是 (a) `entry_range` 二分（76.9%），而非窗口扫描——
方案 A（顺序推进替代二分）是正确方向；方案 B/LBYTE 单独无收益（共享
前缀计算占比可忽略）。

## 3. 目标与成功标准

- **目标 1**：merge 顺序化（不动 `.pgi` 格式），窗口定位从"每条目二分"
  变为"批内顺序推进 + 大跳变自适应二分"——**✅ 已达成**（方案 A，§8）；
- **目标 2**：`.pgi` v3 内嵌 LBYTE（构建时算相邻条目 lcp），为 vlcp 表
  与变长种子铺路——**不做**（拆分探针显示 shared_prefix 占比可忽略；
  完整 LCP 已拒绝，§3.6/§8.4）；
- **目标 3**：完整 vlcp 表语义（FastGA `new_merge_thread`）性能 **≥
  简化版**（~90-100 ms），消除 §3.6 的 2.1×——**不做**（§3.6 结论：
  正确版仍慢 2.1×，性能不可行；简化版为最终形态）；
- **目标 4**：变长种子（§7.3）在顺序访问原语上可落地——**未立项**
  （顺序原语已就绪；**2026-08-12 实验否定，不做**——见
  `pgi-align.md` §7.3.1：短端/长端位置命中增量被链化 + 波扩展吸收，
  端到端 PSL 覆盖无收益）；
- **硬性成功标准**：每一步输出与当前版 **PSL 逐字节一致**（种子级允许
  FastGA 语义下的 0~-50 差异，被链化吸收）；chainnet 覆盖不降；
  merge 耗时不倒退（±5% 噪声内）；内存不涨。

## 4. 候选方案

### 方案 A：归并式 merge（顺序窗口推进，不动格式）——✅ 已落地

**机制**：a 条目按 kmer 全序单调，b 侧窗口下界/上界随之单调（不要求 a
与 b 全等，只要求窗口边界不倒退）。维护"b 游标"，批内每个 a 条目从上次
位置**顺序推进**到新窗口边界，替代 `entry_range` 二分。大跳变（a 步进
大、b 稀疏）时回退二分（自适应：顺序推进超过 `MAX_SEQ_SCAN=64` 组即转
二分）。

**改动面**：`PgiQuery` 新增 `entry_lower_bound_ge(key, from)` 顺序下界
原语（`PgiIndex` 数组顺序扫描 + 超限回退二分；`PgiMmap` 逐组
`group_end` 推进 + 超限回退 `lower_bound`）；`align.rs::emit_entry_hits`
新增 `MergeCursor`（f0/f1/first 三字段 ~24 B）批内维护 floor 窗口边界，
窄窗口从 floor 窗口内顺序推进；批边界（`first`）重置为二分定位。见 §8。

**成本/风险**：低（无格式改动）；风险是 b 侧稀疏时顺序推进比二分慢
（用 `MAX_SEQ_SCAN` 自适应缓解）；并行化下批边界二分开销可忽略（每批
一次）。

**收益**：merge 窗口定位从 O(m log n) → O(m + n)；这是解锁 vlcp 表的
前置（vlcp 需要窗口起点可顺序维护）。实测 resident 反快 ~28%（§8）。

### 方案 B：`.pgi` v3 内嵌 LBYTE（预计算相邻条目 lcp）

**机制**：构建时排序后对相邻条目算 lcp（O(n) `shared_prefix`），每条目
存 1 字节 LBYTE（`PgiEntry` 24→25 B；mmap `RecordLayout` +1 B）。merge
读 `entry_lcp(i)` O(1)，替代即时 `shared_prefix`。

**注意**：LBYTE 是"b 条目与其排序前驱的 lcp"，**不是**"b 条目与当前
a 条目的 lcp"——单独用省不了 m 扫描里的共享前缀计算；它的价值是支撑
vlcp 表的窗口继承（方案 C）与 a 侧传播（a 侧流也可带 LBYTE）。

**成本/风险**：格式版本升级（v2→v3，需兼容读 v2 或强制重建）；构建
O(n) 计算；mmap 布局调整。格式演进独立于性能验证，可与方案 A 并行做。

**状态（2026-08-05 更新）**：不做（§8.4 定稿）——拆分探针显示
shared_prefix 占比可忽略（§2），单独做无收益；完整 vlcp 表已拒绝
（§3.6），LBYTE 失去唯一立项理由。

### 方案 C：merge 侧物化 FastGA cache（b 侧中间表示）

**机制**：把 b 侧索引（resident/mmap）按 cpre 组解码成连续缓存段
（kmer 字节 + LBYTE + freq + positions 偏移），merge 在缓存上做 FastGA
式指针推进与 vlcp 表。统一 `PgiIndex`/`PgiMmap` 为一种表示。

**成本/风险**：内存（b 侧全部物化：E. coli ~1.1M 组 × ~12-24 B ≈
13-26 MB；人类规模需分块/cpre 粒度物化）；失去 mmap 按需读取优势的
部分场景。这是"完全解锁"路径，收益最大、工作量最大（~300-500 行 +
基准）。

### 方案 D：保持简化版（基线）

E. coli 级 LCP 零收益（§3.3 复测：种子差 50、时间/PSL/覆盖全一致）。
如果人类规模（§7.2）验证 LCP 仍无收益且变长种子不立项，则不需要任何
抽象层改动——**阶段 0 的基准应首先回答这个问题**。

## 5. 推荐路径（分阶段，每阶段可独立审核/回退）

```
阶段 0：拆分探针 + 人类规模 lcp 分布预研  ✅（拆分探针已完成）
  ├─ 验证点：merge 内部分布；相邻条目 lcp > 12 的比例（人类重复区）
  └─ 决策门：收益不出现 → 停（方案 D 定稿）
阶段 1：方案 A 归并式 merge（不动格式）  ✅（2026-08-05 落地，§8）
  ├─ 验证点：merge 时间、PSL 逐字节一致、内存
  └─ 决策门：时间不降或退化 → 回退，转阶段 3
阶段 2：方案 B `.pgi` v3 LBYTE（格式演进，可与 1 并行）
  ├─ 验证点：格式兼容策略、构建时间、读 LBYTE 路径
  └─ 决策门：v2/v3 兼容成本 > 收益 → 仅构建侧新增
阶段 3：方案 C / 完整 vlcp 表
  ├─ 验证点：完整 LCP ≥ 简化版性能（§3.6 目标）、语义等价
  └─ 决策门：性能达标 → 变长种子（§7.3）立项
阶段 4：变长种子（§7.3）——在顺序原语上实现 >k adaptamer
```

**2026-08-05 更新**：阶段 0/1 完成；阶段 2-3 经 §3.6/§8.4 定稿**不做**
（完整 LCP 性能不可行、LBYTE 无独立收益），阶段 4 未立项。路径剩余部分
仅在变长种子立项或人类规模数据（§7.2）显示新收益时重启。

**2026-08-12 更新**：阶段 4（变长种子）已实验否定——E. coli + 酵母实测
（`pgi-align.md` §7.3.1）：12–40 bp 短匹配的位置命中增量（31–35%）在端到端
PSL 覆盖上仅 +0.009–0.15%（链化 + 波扩展已吸收），>40 bp 长匹配全部落在
40-mer 已命中位置、对发现同源零增量。路径剩余仅在人类规模数据（§7.2）
显示新收益时重启。

## 6. 基准测试设计

### 6.1 微基准（`examples/merge_mem_bench.rs` 扩展）

- 输入：`tests/genome/{mg1655,nissle1917,ec958}.fa.gz` 的 `.pgi`；
- 变体：当前 / no-LCP / skip-scan / 方案 A / 方案 B（LBYTE 读）；
- 指标：merge 耗时（5 次取中位数）、种子数、峰值 RSS；
- 拆分探针：entry_range 定位 / m 扫描 / shared_prefix / 发射 各自占比。

### 6.2 端到端（对齐 §2.2 口径）

- 命令：`pgr align pgi <ref> <query> --parallel 8`（8 线程 release）；
- 指标：merge/chain_tubes/extend 阶段耗时、PSL 块数、chainnet syntenic
  覆盖（`pgr psl to-chain` + `pgr pl chainnet --syn`）；
- 一致性：与基线 PSL **逐字节 diff**（种子级差异容忍 0~-50，见 §3.3）。

### 6.3 人类规模预研（§7.2 数据到达后）

- 统计相邻条目 lcp 分布（`entry_lcp` / `shared_prefix` 采样），确认
  `lcp > 12` 的比例——决定 LCP/顺序化是否有收益；
- 记录 merge 时间与内存，与 E. coli 外推对比。

## 7. 决策点与开放问题

1. **值得做吗**：E. coli 零收益 → 收益只在人类规模（lcp 分布）或变长
   种子立项后。阶段 0 先回答——**已部分回答**：方案 A 在 E. coli 规模
   即有收益（resident 反快 ~28%），完整 LCP 仍无收益（§3.6 拒绝）；
2. **格式兼容**：`.pgi` v2/v3 双读还是强制重建？（v2 是 2026-08-05 刚
   稳定的格式，破坏性升级需谨慎）——**已解决**：LBYTE 不做，格式保持 v2
   （§8.4）；
3. **并行化**：rayon 批处理与顺序游标的结合方式（批内顺序、批间二分）；
   ——**已解决**：方案 A 的 `MergeCursor` 每批一个，批边界重置为二分（§8.1）；
4. **物化粒度**：方案 C 全量 vs cpre 分块（内存/IO 权衡）；
   ——**未立项**：方案 C 不做（§8.4）；
5. **与 dist/to-hv 的耦合**：`dist pgi`、`to_hv`、`count_unique` 也用
   `PgiQuery`，顺序化改动不应破坏它们的等价性（有现成逐字节测试）。
   ——**已验证**：方案 A 后 `sequential_merge_matches_binary_reference`
   与全量测试（1300）通过，dist/to-hv 逐字节测试未受影响。

## 8. 方案 A 落地记录（2026-08-05）

### 8.1 改动

| 文件 | 改动 |
|---|---|
| `src/libs/pgi/mod.rs` | `PgiQuery` 新增 `entry_lower_bound_ge(key, from)`；`MAX_SEQ_SCAN=64` |
| `src/libs/pgi/mmap.rs` | `PgiMmap` 实现 `entry_lower_bound_ge`（逐组 `group_end` 推进 + 超限回退 `lower_bound`） |
| `src/libs/pgi/align.rs` | `emit_entry_hits` 引入 `MergeCursor`（f0/f1/first），floor 窗口批内顺序推进、窄窗口子窗口内推进；新增等价测试 `sequential_merge_matches_binary_reference` |
| `examples/merge_mem_bench.rs` | 新增 `Probe` 拆分探针（range/scan/emit/shared_prefix 计时） |

语义不变：`entry_lower_bound_ge` 返回的下界与 `entry_range` 完全一致
（`from` 只是顺序推进提示，正确性不依赖它）；`MergeCursor` 每批一个
（~24 B），无内存增长。

### 8.2 实测（mg1655 × nissle1917，k=40 syncmer 8/5，8 线程 release）

`examples/merge_mem_bench.rs`（resident 双索引）：

| 变体 | merge 耗时 | 种子数 |
|---|---:|---:|
| **方案 A（lib）** | **18.6 ms** | 1,121,308 |
| 二分 + lcp（基线） | 27.2 ms | 1,121,308 |
| 二分 no-lcp | 27.6 ms | 1,121,358 |

resident 反快 **~32%**（18.6 vs 27.2 ms），种子数与二分 lcp 版完全一致。
拆分探针（§2）：range 76.3% / scan 11.0% / emit 12.7%。

流式路径（生产 `align pgi`，`RUST_LOG=debug` 探针，3 次）：

| 指标 | 方案 A | 基线（§2） |
|---|---|---:|
| merge 耗时 | 90 / 97 / 99 ms（nissle）；93-97 ms（sakai） | ~90-100 ms |
| 峰值 RSS（merge 后） | 118 MB | ~120 MB（同量级） |

无明显回归，峰值内存持平（`MergeCursor` 每批 24 B，可忽略）。

### 8.3 硬性成功标准核对

| 标准 | 结果 | 证据 |
|---|---|---|
| **PSL 逐字节一致** | ✅ | 端到端 `pgr align pgi` 输出 diff：mg1655×nissle（1299 块）、mg1655×sakai（738 块）与二分基线**逐字节一致**；单测 `sequential_merge_matches_binary_reference` 覆盖 resident + mmap 查询视图 |
| **merge 耗时不倒退** | ✅ | resident 18.6 vs 27.2 ms（反快 ~32%）；流式 90-100 ms |
| **内存不涨** | ✅ | 峰值 RSS 118 MB（持平）；游标 24 B/批 |

`cargo clippy -- -D warnings` clean；`cargo test --release --lib libs::pgi` 50 全过。

> 注：写作时的全量 `cargo test --release` 记录有 1 个 paf 失败
> （`test_invalid_op_panics`）与 2 个 `sd search pgi` 失败。**2026-08-05
> 复测后两者均已消除**：paf 用例经核对用普通 `panic!`（非 `debug_assert!`），
> release 下本就通过，原注表述有误；`sd search pgi` 两用例由 wave trim
> 修复（commit `cd66774`，见 design/pgi-align.md §3.5.7 勘误 7）。当前
> `cargo test`（debug）与 `cargo test --release` 均为 1300 全过。

### 8.4 决策与后续

- **方案 A 定稿**：顺序推进原语已落地，为变长种子（§7.3，阶段 4）铺路
  （vlcp 表已拒绝，见下）；
- **阶段 2（方案 B LBYTE）**：拆分探针显示 shared_prefix 占比可忽略、
  单独做无收益；且完整 vlcp 表已拒绝（§3.6），LBYTE 失去立项理由——
  **不做定稿**，格式保持 v2；
- **阶段 3（方案 C / vlcp）**：§3.6 结论"完整 LCP 正确但性能不可行
  （2.1×）"，不做；简化版（§3.3）为最终形态；
- **阶段 4（变长种子）**：**实验否定，不做**（2026-08-12，`pgi-align.md`
  §7.3.1）；顺序原语（`entry_lower_bound_ge` + `MergeCursor`）保留，若人类
  规模显示新场景可直接复用；
- **人类规模（§7.2）**：数据到达后复测 lcp 分布与顺序推进收益（高重复
  区相邻条目 lcp 可能 >12，窗口加速可能更明显）。
