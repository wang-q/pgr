# pgr align pgi：两基因组归并比对（设计定稿 + 决策记录）

> 定位：`.pgi` 的第一个比对消费者。输入两个已构建 .pgi，输出 PSL 块，
> 喂给 `pgr pl chainnet`（UCSC 链化由 pgr 承担，见 [[fastga.md]] §12.3 决策 3）。
>
> 状态：**2026-08-05 定稿**。管线完整、与 FastGA 对照稳定：chainnet 覆盖
> 持平、端到端反超 ~3.3×；08-05 完成 `PgiMmap` 复用与 lcp 传播两项内部
> 优化（PSL 输出逐字节一致），并对照源码裁定 Gap_Improver（§3.1）与完整
> LCP（§3.6）不移植。
>
> 结构：§0 当前状态 → §1 设计 → §2 验证与基准 → §3 决策记录 →
> §4 已排除方向 → §5 勘误与基准方法 → §6 FastGA 功能差距 →
> §7 未来方向 → §8 相关文档。

## 0. 当前状态（2026-08-05 定稿）

### 0.1 一句话

`pgr align pgi` 完整实现"种子归并 → tube 链化（FastGA `align_contigs`）→
mid-line wave 扩展 → PSL"，与 FastGA 对照：chainnet 覆盖持平（差
0.0-0.015%）、阶段耗时持平（~0.42-0.70s vs ~0.7s）、峰值内存更低（
~161-176 MB vs 332 MB）、端到端反超 ~3.3×。

### 0.2 关键数字（详见 §2）

| 对（MG1655 vs） | pgr 覆盖 | FastGA 覆盖 | pgr 耗时 | pgr 峰值内存 |
|---|---:|---:|---:|---:|
| Sakai | 89.30% | 89.3% | ~0.42s | ~176 MB（FastGA 332 MB） |
| EC958 | 86.28% | 86.3% | ~0.70s | ~161 MB |
| Nissle | 85.21% | 85.30% | ~0.46s | ~164 MB |

### 0.3 已落地（截至 2026-08-05）

- **核心管线**（2026-08-02）：种子归并 → tube 链化 → wave 扩展 → PSL；
  ref 流式 + query mmap；负链 PSL 帧按 UCSC 约定。
- **对齐 FastGA**（2026-08-03）：`pgi build --mask`（soft mask 感知种子）、
  单输入自比对。
- **内部优化**（2026-08-05）：`dist pgi`/`stat`/`to-hv` 复用 `PgiMmap`
  （§3.2）；种子合并 lcp 连续传播（§3.3）——两者 PSL 输出均逐字节一致
  （lcp 种子级有 0~-50 差异，被链化吸收，见 §3.3）。
- **greedy 移除 + self 特判**（2026-08-05）：实测 greedy 链化 + 窗口扩展
  在 chainnet syntenic 覆盖上比 tube 低 3.4-3.7 pp（两对验证），且是
  唯一需要 `--min-span/--max-gap/--band/--merge-gap` 与 `--workflow` 的
  复杂流程；无序列场景改用 chain_tubes + tube 几何块（§3.5.5），greedy
  整套已删除（§3.5.7）。self 模式补上 FastGA 的对角线 0 限制：同 contig
  正链 wave/回溯路径不允许跨越 diag 0，跨 0 的 tube 整管跳过（§3.5.6）。

### 0.4 明确不做

| 项 | 原因 |
|---|---|
| Gap_Improver 移植 | pgr wave 已精确回溯，架构冗余（§3.1 完整论证） |
| 完整 LCP（vlcp 表 / `.pgi` v3 LBYTE） | 正确版实测仍慢 2.1×，性能不可行（§3.6） |
| select 表达式 | 低价值便利，`fa range` 可替代（§3.4） |
| trace points / `.1aln` 紧凑存储 | 人类规模才需要（§6） |
| ALNchain | UCSC chain/net 更标准（§6） |
| GDB / scaffold / 完整 GIX 分片 | pgr 2bit + `.pgi` 已覆盖（§6） |

### 0.5 待办

- 人类规模（~3 Gb）验证——待数据（§7.2）。

## 1. 设计

### 1.1 范围

**做**：

1. 两个排序 .pgi 流的归并 → 种子命中（plen 最大共享前缀 + 频率过滤）；
2. anti-diagonal 空间链化 → 链（greedy 贪心 / tube 两种语义）；
3. 链扩展：banded 仿射 gap 局部比对（greedy）或 Myers wave（tube），无
   序列输入时每条链输出一个 PSL 块；
4. ref 流式 + query mmap 读取（E. coli 规模起不再整体载入内存）；
5. `pgr align pgi` CLI + 集成测试 + E. coli 三株系验证（§2）。

**不做**：见 §0.4（Gap_Improver、select、紧凑存储、ALNchain、GDB/GIX 分片）
与 §7（未来方向）。

### 1.2 数据流

```
ref.pgi（流式）+ query.pgi（mmap）
  │ 1. 归并（plen 最大共享前缀 + 频率过滤）→ 种子命中
  ▼
hit = (key, a_contig, a_pos, a_strand, b_contig, b_pos, b_strand)
  │ 2. 方向解析 + diag/anti 坐标变换
  ▼
fwd/rev 两个空间：(pos_a, pos_b')，diag = pos_a − pos_b'，anti = pos_a + pos_b'
  │ 3. 按 (contig_a, contig_b, 方向) 链化（greedy / tube）
  ▼
chain = 区间 + 对角线带（间距/带宽容忍，跨度过滤）
  │ 4. 扩展（banded 仿射 gap / mid-line wave）+ PSL 输出
  ▼
out.psl → pgr psl to-chain → pgr pl chainnet（现有字节级验证的链化主场）
```

### 1.3 关键语义

#### 1.3.1 种子方向（pgi 双链条目）

`.pgi` 每条位置带 strand 标记（0=正向 k-mer，1=RC k-mer）。两索引归并时
key 相等按 (a_strand, b_strand) 解析方向：

| (a, b) strand | 含义 | 变换 |
|---|---|---|
| (0,0) / (1,1) | 正链命中（两侧实际窗口相等） | `pos_b' = pos_b` |
| (0,1) / (1,0) | 负链命中（b 窗口 = RC(a 窗口)） | `pos_b' = b_len − k − pos_b` |

负链 `pos_b'` 是 RC(b) 空间坐标；PSL 输出时按 §1.3.4 还原。

#### 1.3.2 种子选择（FastGA `new_merge_thread` 语义，2026-08-02 移植）

1. **最大共享前缀（plen）**：a 条目只对其在 b 中的最长匹配发种子，仅共享
   plen 碱基的范围参与配对；短的部分匹配只有在它是该条目最长匹配时存活。
2. **扩展范围频率过滤**：`freq` 作用于 plen 处的出现数（`occ < freq` 才
   保留），而非固定窗口的条目数。
3. **canonical 去重**：`pgi build` 每个位置同时存 fwd/RC 两个 key，a 侧只
   保留 `kmer <= rc(kmer)` 的 canonical 条目（物理命中只发一次）。
4. **floor = 12**：默认 `min-shared` 为 FastGA 的 plen 下限 12。

> 实现注记（2026-08-05）：merge 默认以相邻条目 lcp 起步窗口（FastGA
> `vlcp` 传播的"匹配延续"近似，§3.3），会跳过 `[min_shared, lcp)` 的
> 更短独立匹配（种子级差异 0~-50，PSL 输出不受影响，见 §3.3/§3.6）。

#### 1.3.3 链化

**tube（唯一流程，2026-08-05 起）**：FastGA `align_contigs` 的忠实移植。
种子按对角线分桶（宽 64）→ 相邻桶对按 anti 归并（排序键 (diag 桶,
anti)，§3.5.3 修过顺序 bug）→ tube 维护 anti 覆盖与对角线范围，种子
anti 间隔超 `CHAIN_BREAK`（2000 bp，FastGA 内部值）断开、覆盖达
`CHAIN_MIN`（85 bp，单轴口径 = FastGA 170 anti）触发。tube 扩展用
mid-line wave（BUCK_ANTI=128 滑动），每个 tube 独立 `alast`（并行化
替代 FastGA 的逐对桶共享）+ 输出端 `dedupe_contained`（0.95 阈值，
§3.5.3 修过误删）。链化/扩展参数全部写死 FastGA 常量（BUCK=64、
BREAK=2000、MIN_COV=85、BUCK_ANTI=128、TUBE_MIN_LEN=50、
TUBE_MIN_RATE=0.35），CLI 仅暴露 `-f` 与 `--min-shared`。

**无序列输入**（`.pgi` 对 `.pgi`、不带 `--ref-seq/--query-seq`）：
`chain_tubes` 不需要序列，走 tube 几何块——每个 tube 按其种子跨度
（`a_start/a_end/b_start/b_end`）输出一个单块 PSL（§3.5.7）。有序列时
tube 用 mid-line wave 输出带真实身份率的多块 PSL。

**为什么删掉 greedy**（2026-08-05 实测，见 §3.5.5/§3.5.7）：greedy 用
精确 k-mer 种子 + 窗口扩展跳过低分区间，chainnet syntenic 覆盖比 tube
低 3.4-3.7 pp（Nissle 81.85% vs 85.29%、EC958 82.56% vs 86.32%，8 线程
release）；greedy 的身份率略高（97.73% vs 96.60%）但那是"挑容易的比对"，
对 ChainNet 目的（覆盖）不合算；tube 反而更快（4.9s vs 6.3s，含建索引）。
- **tube**：种子按对角线分桶（宽 64）→ 相邻桶对按 anti 归并（排序键
  (diag 桶, anti)，§3.5.3 修过顺序 bug）→ tube 维护 anti 覆盖与对角线范围，
  种子 anti 间隔超 `CHAIN_BREAK`（2000 bp，FastGA 内部值）断开、覆盖达
  `CHAIN_MIN`（85 bp，单轴口径 = FastGA 170 anti）触发。tube 扩展用
  mid-line wave（BUCK_ANTI=128 滑动），每个 tube 独立 `alast`（并行化
  替代 FastGA 的逐对桶共享）+ 输出端 `dedupe_contained`（0.95 阈值，
  §3.5.3 修过误删）。

#### 1.3.4 PSL 输出（UCSC 约定，负链是坑）

- q = query、t = ref（`pgr align pgi <ref> <query>`）；
- 每条链（greedy）/ 每窗（扩展）= 一个块；
- **负链**：qStart/qEnd 必须正链帧、内部 qStarts 必须 RC 帧（与 `psl chain`
  的 `calc_block_score` 一致）——§3.5.3 记录过整类 '-' 块被静默丢弃的 bug；
- match/mismatch 计数来自扩展（v1 链块为 0，由 `psl chain` 重算）。

#### 1.3.5 模糊碱基（N）与种子发射（2026-08-03 确认）

`pgi build`（`build.rs::collect_one_contig`）对 N/简并碱基是两层处理：

1. **s-mer 哈希：N 当作 A**（`let sb = if code == 4 { 0 } else { code }`，
   注释 "N treated as A, matching pgr dist"）。**syncmer 的位置选择完全
   不受 N 影响**——N 区域照常选出 closed syncmer 位置，也不会产生假种子。
2. **k-mer key：窗口内出现 N 即失效**（`kx/kxr/kvalid` 清零，发射条件
   `kvalid >= k`）。因此每个 N 会让它**前后约 k 个碱基（k=40 时 ±40 bp）
   范围内的种子全部缺失**（该范围内每个 syncmer 的 k-mer 窗口都含 N）。

影响评估（用于判断"清洗把 IUPAC 转 N"是否安全）：

*   种子空洞是局部的：链化用 `max-gap=1000` / `merge-gap=5000` 桥接，
    ±40 bp 的空洞远低于阈值，链照常跨过、区间照常合并。
*   只要 N 密度远低于 1/k（实际清洗后 repbase 约 1.4 万个 N、每条序列
    平均不到 1 个，tncentral/dfam 更少），种子密度损失可忽略；实测
    repbase 清洗后 e-align 输出正常（50,784 bp，与预期一致）。
*   边界：**长 N 区段**（scaffold 内的大段 N）会割出超过 merge-gap 的
    种子空洞，敏感度才真正下降——那是输入质量问题，不是清洗引入的。
*   反向安全：N 当 A 参与哈希只影响"选不选这个位置"，k-mer key 失效
    保证不会发射错误序列内容的种子，无假阳性。

### 1.4 CLI

```
pgr align pgi <ref> <query> -o out.psl
  [--freq 10] [--min-shared N] [--parallel 8] [--keep-index]
  [-k 40] [--smer 8] [--window 5]
  [--ref-seq ref.fa|2bit] [--query-seq query.fa|2bit]
```

- 输入 `<ref>/<query>` 可以是基因组（FASTA/.gz/2bit）或 `.pgi`，任意混用：
  基因组自动建索引（同目录同名 `.pgi` 存在则复用，否则临时目录，可选
  `--keep-index` 保留），序列本身兼作扩展输入；`.pgi` 直接使用，扩展序列
  由 `--ref-seq/--query-seq` 提供并做 contig 校验（防索引-序列不一致的
  静默错误）；
- 两侧参数必须一致（复用 `dist pgi` 校验）；`-k/--smer/--window` 仅序列
  输入生效（`.pgi` 读索引头），显式传入与复用缓存冲突时报错；
- 无序列输入（`.pgi` 对 `.pgi`）：每个 tube 一个几何 PSL 块（种子跨度）；
  有序列：mid-line wave 输出带真实身份率的多块 PSL；
- query 索引（或现场建的临时索引）必须是真实文件（mmap 不支持
  stdin/gzip）。

## 2. 验证与基准

### 2.1 测试基线

**912 测试全过**（2026-08-02 基线；后续 mmap 测试 +4）。回归覆盖：方向解析
（fwd/rev）、频率过滤、链化边界（band/gap/span）、tube anti 序、dedupe
0.95 保留延伸块、最大前缀/扩展范围过滤、负链 RC 帧、多 contig（RC contig
+ 2% 突变）、mmap 与全量读入等价性；集成测试 `tests/cli_align_pgi.rs`
（identical / RC / mutation / tube + 序列直入 / 复用 / 混用 / 校验拒绝）。

> **2026-08-05 更新**：`align pgi` 22 个 CLI 测试、`libs::pgi` 55 个单测
> 全通过，`cargo test` 全量 1255 通过（audit 后新增 crafted 索引/负链帧/
> 数据安全等回归，见 [[../audit/audit-pgi-align.md]]）。

### 2.2 当前基准（2026-08-05 复测，真实数据，8 线程，release）

| 对（MG1655 vs） | pgr chainnet 覆盖 | PSL 记录 | 比对耗时（预建索引） | 端到端（含建索引×2） | 峰值内存（预建索引） | FastGA 覆盖/耗时 |
|---|---:|---:|---:|---:|---:|---:|
| Sakai | 89.30%（582 块） | 738 | ~0.42 s | ~1.23 s | ~176 MB | 89.3% / ~0.7 s |
| EC958 | 86.28%（811 块） | 815 | ~0.70 s | ~1.3 s | ~161 MB | 86.3% / ~0.7 s |
| Nissle | 85.21%（820 块） | 1299 | ~0.46 s | ~1.23 s | ~164 MB | 85.30% / ~0.7 s |

> 口径：chainnet 覆盖 = `pgr psl to-chain` + `pgr pl chainnet --syn` 后
> 目标（mg1655）被 syntenic 块覆盖的碱基比例；pooled PSL identity =
> `(matches + rep_matches) / block_len`（block_len 含错配与插入）。
> 覆盖/块数随 wave trim 更新（2026-08-05，§3.5.7）。

阶段分布（预建索引 + `--ref-seq/--query-seq`，`RUST_LOG=debug` 探针，
2026-08-05 数据）：Sakai merge 93-97 / chain_tubes 74-77 / extend
236-242 ms；Nissle merge 90-99 / chain_tubes 70-72 / extend 284-286 ms；
EC958 merge 115 / chain_tubes 89 / extend 450 ms。merge 命中 1,121,308
（nissle）/ 1,181,074（sakai）/ 1,090,053（ec958）。resident 微基准
（`examples/merge_mem_bench.rs`）merge 18.6 ms（方案 A，见
[[pgi-query-layer.md]] §8.2）。

### 2.3 端到端管线验证

`pgr align pgi` → `pgr psl to-chain` → `pgr pl chainnet --syn` 全链路，与
FastGA 驱动版本对比 syntenic MAF。下表为早期管线快照（2026-08-02）：

| 输入对 | pgr 管线 | FastGA 管线 |
|---|---:|---:|
| MG1655 vs Sakai | 87.7%（392 块） | 89.3%（506 块） |
| MG1655 vs Nissle | 82.9%（541 块） | 85.3%（711 块） |

> **当前 chainnet 覆盖以 §2.2 为准**（Sakai 89.30% / Nissle 85.21%）——
> 87.7%/82.9% 与 89.30%/85.21% 的差距来自后续 merge-gap、种子选择与
> 负链 PSL 修复（§3.5.2/§3.5.3）。**2026-08-05 复测**：端到端（建索引 ×2
> + 比对）pgr 1.23 s vs FastGA 4.04 s，反超 ~3.3×，见
> [[../benchmarks/bench-pgi-align-vs-fastga.md]]。角色约定：`pgr align pgi
> <ref> <query>` 的 PSL 是 q=query/t=ref，FastGA 输出相反，喂 chainnet 前
> 需 `pgr psl swap`。

### 2.4 10 株 cohort 两两验证（45 对）

扩展块身份率矩阵（2026-08-05 复测；行×列 = ref×query；pooled PSL
identity = `(matches+rep)/block_len`，口径见 §2.2；45 对完整数据见
[[../benchmarks/dist-cohort-validation.md]]）：

- 分布 95.4-98.8%（均值 96.8%），与亲缘关系一致（最近 nissle–cft073
  98.8%、e24377a–se11 98.4%；最远 ec042–sakai 95.4%）；合并块后身份率
  低 ~0.2-0.5%；与 FastGA 同口径基本持平（mg1655×sakai：pgr 97.30% vs
  FastGA 97.28%）。
- 块数随 wave trim 更新（2026-08-05，§3.5.7）：低身份末端被裁，块数
  下降（如 mg1655–sakai 738 条记录）。

## 3. 决策记录（2026-08-05 重写重点）

> 除核心算法外，本模块还经历了几次"要不要做 X"的抉择。以下按决策顺序记录：
> 动机 → 尝试 → 实测 → 结论。判断原则：**对照源码语义；输出等价或可验证；
> 不因表面指标（如块数减少）牺牲结构正确性**。

### 3.1 Gap_Improver：尝试 → 放弃

**动机**：indel 复杂区覆盖缺口 ~0.7 kb（Nissle，噪声级）。Myers 的
`Gap_Improver`（align.c:6714）是 wave 输出后对 gap 区的精修，曾被视为
"indel 复杂区的答案"。

**尝试 A：仿射分数 DP（自选参数）**。盒检测 + banded 全局仿射重比对
（match 5 / mismatch -4 / gap open -8 / extend -6，即 pgr `AlignmentParams`
默认），接入 tube 扩展输出端。实测（MG1655 vs Sakai，`--workflow tube`）：

| 指标 | 无 | 仿射版 |
|---|---:|---:|
| chainnet syntenic 块 | 1204 | 912（-24%） |
| syntenic 覆盖 | 8,295,891 | 8,298,411（+2.5 kb） |
| `psl to-chain` 条数 | 23,106 | 13,964（-40%） |
| PSL identity | 0.9757 | 0.9741（-0.2%，mismatches +9k） |

表面看"块数大减、覆盖不降"似乎有益，但 **identity 下降是警示**：仿射参数
（gap open -8 远大于 mismatch -4）有强烈动机把**真实短 indel 折叠成
mismatch**。受控测试（1 个真实 indel + 1 个真实错配）显示长度差 1 的场景
保留 gap，但真实数据 mismatches +9k 说明结构失真确实在发生。

**尝试 B：unit-cost 盒内重算（FastGA 语义）**。复用 wave 的 `dandc_nd`
做盒内 unit-cost 最优路径，输出与不接入**逐字节一致**（无操作）。

**对照源码**（align.c:6714-7140）：FastGA 盒内是 Myers 最远到达点（F/G/H
数组）的 **unit-cost 稀疏 DP**（`n += 1`，mismatch 与 gap 同价，**无 gap
open/extend 参数**）。它有效的根源是 FastGA 的 wave 输出为 **tspace=100
采样的不完整 trace**，盒内 DP 负责补全；pgr 的 wave 用 `dandc_nd` 全量
精确回溯，无采样缺口——**源码语义下 Gap_Improver 对 pgr 必然无操作**。

**结论：不移植。** 两个尝试给出互补证据：

1. 源码语义（unit-cost）= 逐字节无操作 → 机制前提不成立（pgr wave 已精确）；
2. 自选参数（仿射）= 有副作用（identity -0.2%，结构失真）+ 无源码依据 →
   不可接受。块数减少迎合了 axtChain 评分对 gap 的厌恶（gap 惩罚远大于
   mismatch），**不代表生物学正确**。

**教训**：

* 块数减少本身不是成功标准；标准是覆盖不损失 + 结构不失真；
* 参数必须对照源码，不能自选；
* "看起来在改善"的启发式要复核其真实代价（与 §4 门控类启发式同款教训）。

### 3.2 dist/stat/to-hv 复用 PgiMmap：落地

**动机**：align 已用 `PgiMmap`（FastGA GIX 模型，内存 -23%），
`dist pgi`/`stat`/`to-hv` 仍全量 `PgiIndex::read`。

**做法**：`dist_between`/`index_to_hv` 泛型化为 `&impl PgiQuery`，三个命令
改用 `PgiMmap::open`；新增 `count_unique`（entry 组遍历）。新增
mmap/resident 等价性测试（dist 指标、hv 投影、unique 计数）。

**结论：落地**（2026-08-05）。输出与全量读入逐字节一致（测试固化），
内存 -23% 预期。这是 FastGA GIX 读取模型的直接对应，无自选参数。

### 3.3 种子合并 lcp 连续传播：落地

**动机**：FastGA `new_merge_thread`（FastGA.c:610）用 `vlcp[plen]` 表 +
相邻条目 LBYTE 做 O(1) 摊销的 lcp 传播；pgr `emit_entry_hits` 每个条目从
`min_shared` 窗口从头扫描。

**做法**：`emit_entry_hits` 新增 `prev_kmer` 参数，扫描窗口从
`max(min_shared, lcp(prev, cur))` 起步（窄窗口为空则回退 `min_shared`），
批内顺序维护前驱。

**结论：落地**（2026-08-05）。流式路径（唯一采用路径）merge 107 vs
109 ms（no-lcp，噪声级，性能中性）；PSL 输出与 no-lcp 逐字节一致
（734 块 / 4,948,358 aligned，3 对命令级复核，见 §3.6）。注意：它实际
采用 FastGA `vlcp` 的"匹配延续"近似——跳过 `[min_shared, lcp)` 的更短
匹配（种子级差异 0~-50，被链化吸收）；语义向 FastGA 对齐，非严格
"最长匹配"（§1.3.2 注记）。纯内存路径整体比流式慢（不采用），无附加
动作。是 §7.3 变长种子的语义前置。

> **tube 语境复测（2026-08-05，greedy 移除后）**：当时的评估基于 greedy
> （exact k-mer，min_shared=k），lcp ≤ k 使窗口无跳过空间，LCP 完全无
> 操作。tube 唯一流程后重测（mg1655 vs nissle，min_shared=12，8 线程
> release，各 3 次）：种子差 50（no-lcp 1,121,358 vs lcp 1,121,308，
> 正是跳过 `[12, lcp)` 短匹配的量级），但 **merge 耗时（~105 vs ~102
> ms）、PSL 输出（逐字节一致）、chainnet 覆盖（均 85.286%）全部无差异**。
> 原因：syncmer 稀疏采样下相邻条目 lcp 通常 < 12（start 经常回退
> floor），且跳过的短种子不形成独立链。**结论未变**：E. coli 级 LCP 零
> 收益；保留理由是语义对齐（tube 种子流 = FastGA adaptamer 口径）+ §7.3
> 变长种子前置，无害。完整 LCP（§3.6）2.1× 慢属 PgiQuery 抽象层，与
> workflow 无关，不做结论不变。**人类规模（§7.2）需复测**：高重复基因组
> 相邻条目 lcp 可能 > 12，窗口加速才可能显现。

### 3.4 select 表达式：不做（用处不大）

**动机**：FastGA `select.c` 支持"只比对选定 contig/区间"。

**评估**：对 pgr 用处不大，理由：

- pgr 的 `align pgi` 面向全基因组 pairwise 比对，contig/区间级定向场景可用
  `fa range` 截取 + 子索引间接实现，效果等同，只多一步中间文件；
- FastGA 的 select 主要为 ALNview 交互查看设计（`interpret_point` /
  `interpret_range` 服务于 dot plot 焦点定位），CLI 批量场景少用；
- 语法绑定 GDB 的 scaffold/contig 两级模型（`@`/`.` 符号），pgr 是平铺
  contig，移植需改语义；
- 对比对质量/性能无任何贡献。

**结论：不做**。避免为便利功能引入解析器 + CLI 集成 + 维护成本。

### 3.5 开发历史里程碑（2026-08-02 压缩保留）

> 早期逐条记录（v1→v3、5.x 编号）已精简，仅保留里程碑与结论；完整细节见
> git 历史与 [[../audit/audit-pgi-align.md]]。

#### 3.5.1 链化与扩展引擎演进

- **v1 链块 → v2 banded SW 扩展 → v3 分窗扩展**（16 kb 窗口 + 2 kb 重叠
  滑动）：自比对主链身份率精确 1.0000000（5.30M/0），跨株身份率 97.3%
  （pooled PSL identity；同口径 FastGA 97.28%）。
- **性能关键点**：banded DP 按带限列迭代（16000 列里只有 65 列有效，
  246× 浪费）；窗口级 rayon 负载均衡（自比对主链 332 窗口是单线程长尾 →
  摊平后自比对 37.8s→0.84s、跨株 2.0s→0.66s）。
- **仿射 gap**（open -8 + extend -6）：块数 -25%，身份率不变，仅 indel
  结构更干净。
- **Myers wavefront 移植**（`src/libs/pgi/wave.rs`，FastGA align.c；2026-08-05
  从 `libs/alignment/` 迁入，与 pgi 链化/扩展同属 FastGA 移植）：
  - 单独接入 banded 路径**不可用**（块数 3-4×、覆盖 71.6%/32.7%）——
    wave 依赖 tube 的锚定上下文；
  - 按 FastGA `align_contigs` 语义接入 tube：643 块/88.2%/425 MB vs
    FastGA 701/89.3%/~0.7s，质量逼近、内存 -3×；tube 并行 + anti 上限
    40 kb + 链排序 u128/rayon 后调用数 58k→~900，端到端 ~0.7s。

#### 3.5.2 种子语义演进

- **5.9 盲目部分匹配（无 plen 最大选择）：不可用**——min-shared 12 产生
  53015 块/0.9496 假阳性爆炸；部分种子收益依赖 FastGA 的 tube/wave 机制，
  不能直接移植到贪心链化。
- **5.29 plen 最大选择 + 扩展范围过滤 + canonical 去重落地**（§1.3.2）：
  种子减半、内存不升反降，三对覆盖 +0.05~0.24%（Sakai 89.26→89.31% 等）；
  tube 默认 floor 定为 12（12-19 bp 锚点补 indel 复杂区，与 5.9 的
  "无最大选择"机制不同）。

#### 3.5.3 bug 修复要点

| 条目 | 影响 | 修复 |
|---|---|---|
| tube 合并顺序 | Sakai 缺失 55 kb | 排序键 (diag 桶, anti)，按 anti 归并 |
| dedupe 误删 | 丢 3.1 kb 真实覆盖 | 双轴重叠阈值 0.80→0.95 |
| 大 tube 同源门控 | 误杀整管（7.3 kb，99% 身份） | 根因修复后直接移除（教训：门控类启发式要复核） |
| 负链 PSL 帧 | 所有 '-' 块被 chainnet 丢弃（Nissle 0.32% 差距主因） | qStart/qEnd 正链帧、qStarts RC 帧（UCSC 约定） |
| self 跨 diag 0（2026-08-05） | self 模式可能输出"同坐标自己 vs 自己"的假比对块 | wave 对角线限制（minp/maxp）+ banded 回溯路径（§3.5.6） |

#### 3.5.4 内存与性能优化结论

- 峰值内存：0.96 GB（双索引全量）→ 预建索引路径 **~161-176 MB**
  （EC958/Sakai；ref 流式 `PgiStream` + query mmap `PgiMmap` + positions
  位域打包 + 种子减半），端到端（含临时索引构建）~209-223 MB；均低于
  FastGA 实测 332 MB。
- 耗时：37.8s → 比对 **~0.42-0.70 s**（预建索引，Sakai/EC958/Nissle）；
  阶段分布（2026-08-05，Sakai）merge 93-97 / chain_tubes 74-77 / extend
  236-242 ms；端到端（含建索引×2）**~1.23 s**（§2.2）。
- 结构上限记录：`SeedHit` contig u16（build 守卫 >65535 报错）、
  `pack_position` cid 上限 2^20、`PgiEntry` 24 B/条——人类规模复核见 §7.2。
- 命令迁移（2026-08-03）：`pgr pgi align` → 顶层 `pgr align pgi`（`pgr pgi`
  收敛为纯索引管理 build/stat/to-hv）；基因组输入自动建索引（sibling 复用 +
  mtime 失效 + 参数一致性校验），`.pgi` 输入配 `--ref-seq/--query-seq`
  做 contig 校验。

#### 3.5.5 默认 workflow 转 tube：greedy vs tube 实测（2026-08-05）

**动机**：一直保留"greedy 默认、tube 可选"，但从未在同一代码上对比两者
的最终产出。2026-08-05 实测（当前代码，8 线程 release，MG1655 为 ref，
含自动建索引）：

| 指标 | greedy（原默认） | tube（现默认） |
|---|---:|---:|
| Nissle chainnet syntenic 覆盖 | 81.85%（769 块） | **85.21%（820 块）** |
| EC958 chainnet syntenic 覆盖 | 82.56%（724 块） | **86.28%（811 块）** |
| Nissle PSL 记录 | 1,574（全单块） | 1,299（多块带 gap） |
| Nissle PSL identity | 97.73% | 96.25% |
| 种子 | exact k-mer（min-shared=k） | partial（floor 12） |
| 端到端（含建索引×2） | 6.3 s | **1.23 s** |

**greedy 为什么覆盖低**：

1. 种子是精确 k-mer（§5.9 证明部分匹配在贪心链化下假阳性爆炸），indel
   复杂区（Nissle 每 ~300 bp 一个 indel）链直接断；
2. 16 kb 窗口 banded SW 在低分窗口返回 None，链回退几何块（identity 0）
   被 chainnet 过滤——正是 bench 笔记里"Sakai 差距来自分歧区的 wave 式
   补齐"所指。greedy identity 更高是"挑容易的比对"，对 ChainNet 覆盖
   不合算。

**结论**：默认与推荐流程全面转 tube（覆盖 +3.4-3.7 pp、更快）。无序列
场景（`.pgi` 对 `.pgi`）本就是一点小例外，为此保留整套 greedy
（链化 + 窗口扩展 + 5 个 CLI 参数）不合算——`chain_tubes` 不需要序列，
无序列时直接输出 tube 几何块即可。**greedy 整套已删除**（§3.5.7）。

#### 3.5.6 self 模式对角线 0 限制（2026-08-05）

**动机**：self 模式（单输入自比对）是本工具重要功能；对照 FastGA 发现
pgr 缺 `align_contigs` 的 self 分支（FastGA.c:3220-3240）：同 contig
正向自比对时，wave 对角线不允许跨越 0（精确自同线）——tube 全正则
`minp=1`、全负则 `maxp=-1`，跨 0 的 tube 整管跳过。

**实现**（对照源码语义）：

1. `wave.rs::local_alignment` 加 `selfie` 参数：`forward_wave_mid` 扩展
   时用 `minp/maxp` 硬夹对角线（FastGA `forward_wave` 的 `low>=minp`/
   `hgh<=maxp` 分支），反向 wave 的边界镜像换算（k' = m-n-k）；
2. 回溯路径：FastGA 盒内 DP 带宽 = 端点 diag 差（align.c Gap_Improver
   `Diag=|Fdag-d|+1`），pgr 的 `dandc_nd` 无带可能跨 0——新增
   `banded_edit_ops`（unit-cost 带限 DP + 回溯），self 模式下带 =
   **tube 原始带 ∪ 端点 diag** ∩ 单侧限制（端点会因拷贝边界漂移，如
   tandem repeat 处 diag 从 -400 漂到 -381，只锁端点带会丢失匹配
   对角线）；
3. `extend_tube` 只在 `a_contig==b_contig && strand==0` 时启用 selfie
   （FastGA `SELF && ctg1==ctg2 && !comp`）；跨 contig 的 repeat 拷贝
   不受限。

**验证**：wave 单测（正/负对角线路径不跨 0、跨 0 tube 返回 None）、
align 单测（tandem repeat self 输出块无 diag 0）、真实数据 mg1655 self
437 条记录 / 2,861 个 block（2026-08-05 复测；wave trim 后块数下降，
pooled 身份率 0.9803），**零个 diag 0 块**。附带修复
`align_to_psl_ext`（非流式）漏调 `drop_self_hits` 的不一致。

#### 3.5.7 greedy 流程删除（2026-08-05）

**动机**：上一轮把默认 workflow 转 tube 后，greedy 的唯一保留理由是
"无序列回退路径"。用户指出：无序列场景（`.pgi` 对 `.pgi` 不带
`--ref-seq/--query-seq`）很容易判断，为这一点例外保留整套复杂流程
（贪心链化 + 中间同源验证 + 相邻链合并 + 16 kb 窗口 banded 扩展 + 5 个
CLI 参数）不值得。

**做法**：

1. `Tube` 结构增加种子跨度字段（`a_start/a_end/b_start/b_end`，
   `chain_tubes`/`tubes_for_group` 归并时累计，语义同原 greedy 链的
   种子跨度）；
2. 新增 `tube_to_psl`（原 `chain_to_psl` 逻辑，字段换 tube）；无序列的
   `align_to_psl`/`align_to_psl_streaming` 改为 `chain_tubes` + 几何块；
3. 删除：`Chain`/`ChainCursor`/`chain_hits`/`merge_adjacent_chains`/
   `middle_is_homologous(_range)`/`push_chain`/`chain_to_psl`/
   `extend_chain`/`chain_windows`/`extend_window`/`WindowJob`/
   `EXTEND_WINDOW/STEP`/`GREEDY_MIDDLE_MIN_GAP`/`MIDDLE_*_CAP`/
   `SeqPair`、`Workflow` 枚举、`AlignParams` 的
   `min_span/max_gap/band/merge_gap/workflow` 字段；
4. CLI：`align pgi` 与 `rept e-align` 删除
   `--workflow/--min-span/--max-gap/--band/--merge-gap` 参数；
   `AlignParams` 只剩 `freq` + `min_shared`；
5. 测试：7 个链化单测 + 窗口扩展测试删除，`psl_block_coordinates` 改为
   tube 版；CLI 无序列测试直接走 tube 几何块。

**收益**：`align.rs` 从 ~2,800 行降到 ~1,800 行；CLI 更简单；无序列能力
保留（tube 几何块与 greedy 几何块同为种子跨度单块 PSL）。`rept e-align`
管线同步清理（RM 配方只用 `-f`/`--min-shared`，见 repeat-masking.md）。

**遗留**：`libs/alignment/banded.rs`（`align_banded_local`）的唯一调用者
`middle_is_homologous` 随 greedy 删除，成为孤儿，2026-08-05 已删除（无
消费者；git 历史可恢复）。

**副作用与修复（2026-08-05 定位，2026-08-05 修复，§5.1 勘误 7）**：本条
删除动作（commit `1a75965`）后 `tests/cli_sd.rs` 的两个倒位重复测试失败
（`command_sd_search_pgi_inverted_repeat` /
`command_sd_search_pgi_close_inverted_repeat`，断言 PSL 块 M 列应为 0，
实际 18/24/15）。git bisect 确认引入 commit = `1a75965`（`a2bfabd` PASS →
`1a75965` FAIL）。根因不是链化（tube 本身干净）：pgr 移植 FastGA
`forward_wave` 时 trim 判断只保留 `PATH_AVE`（60 列 ≥ 42 匹配），丢掉了
`TABLE`/`SCORE` 打分表检查（align.c `set_table`，TRIM_LEN=15，最近 30 列
比分"前缀正"才更新 trim 点），wave 端点因此爬进非同源侧翼 ~35-40 bp
（~50-72% 身份的噪声列）。修复：`wave.rs` 补 `TrimSpec`（复刻 `set_table`，
mscore/dscore 由参考碱基偏差 + 默认 ALIGN_RATE=0.3 推导），trim 更新加
`tip_ok` 检查；合成数据输出与 FastGA `-psl` 逐字节一致（1200 bp、M=0），
32768 项打分表与 C 版逐项一致。附带：mmap 重构 commit `2189b90` 本身是
broken 中间态（`to_hv.rs` 引用不存在的 `count_unique`，编译不过），被
下一个 commit `a2bfabd` 修复。

### 3.6 完整 LCP（vlcp 表 / `.pgi` v3 LBYTE）：尝试 → 不做

**动机**：FastGA `new_merge_thread`（FastGA.c:610）用 `vlcp[plen]` 表 +
相邻条目 LBYTE 做 O(1) 摊销的 lcp 传播；pgr 简化版（§3.3）仍逐条目算
`shared_prefix`。是否值得移植完整形态？

**尝试 1：逐碱基递增（二分模拟）**：`window(m+1)` 空即最大 m——merge
543 ms vs 简化版 107 ms（慢 5×）：每步 `window(m+1)` 都是 PgiQuery 二分
（O(log n)），FastGA 的 O(1) 靠 C 指针推进，PgiQuery 抽象无此优势；且
忽略频率时丢种子（715 vs 734 块——GIX 构建侧滤高频，pgr 在 merge 侧滤）。

**尝试 2：`narrow_prefix` 顺序接口（正确版）**：给 `PgiQuery` 加顺序收缩
方法（PgiIndex 二分 / PgiMmap 顺序字节），实现"逐碱基递增 + 频率回退"。
修复 3 个 bug：回退下限 `start`→`min_shared`、窗口回退后 `m` 重置、
PgiMmap 前缀掩码（`pack_kmer` 高位对齐，原 mask 取低位致 k%4≠0 时比较
错）。结果：语义完全等价（60 测试全过，含逐条目对比），但 merge
220-229 ms vs 107 ms（仍慢 2.1×）——顺序字节收缩在 Rust/PgiQuery 下
比 u128 位操作慢。

**性能分解（`examples/merge_mem_bench.rs`，两索引全量载入）**：

| 变体 | merge 耗时 | seed hits |
|---|---:|---:|
| 流式（命令路径，ref 流式 + query mmap） | 107 ms | 1,181,074 |
| 内存版（lib，lcp） | 24 ms | 1,181,074 |
| 内存版（no lcp，窗口从 min_shared） | 22 ms | 1,181,077 |
| 内存版（lcp + 跳过扫描） | 69 ms | 1,755,870 |

发现：

- **IO 占流式 merge ~78%**（内存计算仅 24 ms）；流式（107）比全量加载
  （~155 ms）更快——ref 流式 + query mmap 架构本身最优；
- **内存版 lcp 慢 2.3-2.4×**（4 对 E. coli，回退开销 ~8.5% 条目）；跳过
  扫描反而更慢（扫描是净收益：精确 max 缩小 occ 窗口/种子数）。纯内存
  路径整体比流式慢，不采用，故无需"禁用"动作；
- **种子语义影响**：lcp 跳过 `[min_shared, lcp)` 更短匹配，种子差异
  0~-50（mg–sakai -3、mg–nissle -50、mg–ec958 0、sakai–nissle -38），
  PSL 输出 3 对命令级逐字节一致（被链化吸收）。

**结论：不做**（2026-08-05 定稿）。完整 LCP 语义可做、性能不可行（正确版
仍慢 2.1×）；简化版（§3.3）是最终形态。`.pgi` v3 存 LBYTE 属独立格式
演进（非性能优化），暂缓（§7.6）。`narrow_prefix` 实验代码已全部回退。

## 4. 已排除方向（避免重试）

| 方向 | 结论 | 原因 |
|---|---|---|
| adaptamer 部分种子（盲目，无最大选择） | 不可用 | 弱种子大量假阳性：min-shared 12 → 53015 块/0.9496 |
| wave 单独接入 banded 路径 | 不可用 | 无 indel 偏好 + 贪心延伸；wave 依赖 tube 锚定上下文（§3.5.1） |
| 大 tube 同源门控（多对角线滑窗） | 已移除 | 采样漏真实对角线误杀整管；根因修复后收益消失 |
| 中心对角线滑窗身份率门控 | 不可用 | 覆盖 88.2%→70.9% |
| 种子覆盖密度 / 邻近门控 | 不可用/无效 | 无法干净区分生产性 tube |
| CHAIN_BREAK 调小（300/100） | 更差 | tube 碎片化，质量 87.3% |
| tube 默认 min-shared 30 | 更差 | 部分匹配噪声未被抑制 |
| pgi 解析并行化 | 无益 | 瓶颈是磁盘读取而非解析 |
| wave D&C 回溯每次调用全跑 | 已改 | FastGA 的 `dandc_nd` 是死代码 |
| **Gap_Improver 仿射版（自选参数）** | **已移除** | 参数非源码语义、identity -0.2% 结构失真（§3.1 完整论证） |

## 5. 勘误与基准方法

### 5.1 关键勘误清单

1. **FastGA 内存 "~0 MB" → 实测 332 MB**（2026-08-02 直接实测）：FAtoGDB
   ~7 MB、GIXmake ~160 MB、FastGA 比对主进程 332 MB；源码核查全项目无
   mmap（GIX 流式 read + GDB EXTERNAL 文件态）。pgr mmap 版 ~161-176 MB
   （预建索引路径，§2.2）已低于 FastGA。
2. **Nissle 基线曾误判无效 → 有效**（文件逐字节相同）；其 0.32% 差距主因
   是负链 PSL 坐标帧 bug（§3.5.3），非 chainnet 过滤。
3. **indel 复杂区不能用 naive 偏移身份判断**：每 ~300 bp 一个 indel，同偏移
   身份 ~25% 但共享 40-mer ≈ 99%（8990/10258）——种子不缺，缺口在扩展端
   （§7.1）。
4. **`is_minimal` 是 canonical 方向判断**，不是噪声抑制；真正的噪声抑制是
   plen 最大选择 + 扩展范围过滤。
5. **tube `CHAIN_BREAK` 口径**：FastGA `-s 1000` 内部翻倍为 anti 空间
   2000；修正后块数 +1.5%、三对覆盖 ±0.1% 内。
6. **验证数据统一**：早期 mmap 验证用合成随机序列 → tests/genome 真实数据
   （数字以 §2.2 为准）。
7. **wave trim 打分表丢失导致端点爬进非同源侧翼**（2026-08-05）：pgr 移植
   `forward_wave` 时只留 `PATH_AVE` 匹配数检查，丢掉 FastGA `set_table` 的
   最近 30 列"前缀正"比分检查（TRIM_LEN=15），倒位重复的互惠块端点各爬进
   侧翼 ~35-40 bp，M 列 18/24/15（§3.5.7）。修复后与 FastGA `-psl` 输出
   逐字节一致（§3.5.7）。

### 5.2 基准方法

- 数据：`tests/genome/{mg1655,sakai,nissle1917,ec958}.fa.gz`（另有
  cft073/e2348_69/e24377a/ec042/se11 等 cohort 株）；
- 命令：`pgr align pgi <ref> <query> --ref-seq --query-seq`（8 线程默认，
  2026-08-05 起 tube 是唯一流程，无 `--workflow`），release；
  `/usr/bin/time -v` + `RUST_LOG=debug` 阶段探针（merge/chain_tubes/
  extend + VmHWM）；
- 覆盖：`pgr psl to-chain` → `pgr pl chainnet --syn` syntenic 覆盖；
- 端到端（含建索引 ×2）见 [[benchmarks/bench-pgi-align-vs-fastga.md]]；
- 10 株 cohort 验证见 [[benchmarks/dist-cohort-validation.md]]（引用
  §2.4 的身份率矩阵）。

## 6. FastGA 功能差距

> 对照 [[fastga.md]]（参考笔记）与 `pgr align pgi` 现状。状态口径：已落地 /
> 不做（有结论） / 暂缓（未来可能）。

### 6.1 状态总览

| 功能 | 状态 |
|---|---|
| soft mask 感知的种子发现（`-M`） | **已落地**（`pgr pgi build --mask`） |
| 自比对模式（`FastGA A`） | **已落地**（`pgr align pgi` 单输入 self） |
| PAF `cs:Z` 输出（`-pafs/S`） | **已落地**（`pgr maf to-paf`） |
| Gap_Improver | **不做**（pgr wave 已精确，§3.1） |
| 完整 LCP（vlcp 表 / `.pgi` v3 LBYTE） | **不做**（正确版性能不可行，§3.6） |
| select 表达式 | **不做**（低价值，§3.4） |
| 多 mask union | 暂缓 |
| `-S` 对称 adaptamer | 暂缓（专门场景） |
| trace points / `.1aln` 紧凑存储 | 不做 |
| ALNchain | 不做 |
| GDB / scaffold / 完整 GIX 分片 | 不做 |

### 6.2 说明

- **Gap_Improver**：FastGA 盒内 unit-cost DP 针对 tspace 采样 trace 的补全；
  pgr wave 精确回溯，无需补全（§3.1）。
- **select 表达式**：`fa range` + 子索引可间接实现（§3.4）。
- **多 mask union / `-S` 对称 adaptamer**：专门场景，暂缓（§7）。
- **trace points / `.1aln`、ALNchain、GDB/GIX 分片**：pgr 用 PSL/MAF +
  UCSC chain/net + 2bit/`.pgi` 替代，不做（格式对比见 fastga.md §9/§10）。

## 7. 未来方向

> 2026-08-05 定稿。已落地部分稳定（§0-§3），剩余按价值排序。

### 7.1 indel 复杂区覆盖缺口（~0.7 kb，Nissle）

该区域每 ~300 bp 一个 indel，但共享 40-mer ≈ 99%（§5.1 勘误 3）：**种子
不缺，缺口在扩展端**。Gap_Improver 已裁定不做（§3.1）；候选方向：

* wave 扩展的边界处理（banded 窗口对低分区间跳过 → 参考 FastGA wave 的
  桥接方式，但保持 unit-cost 语义）；
* 或交由下游 chainnet 容忍（当前缺口噪声级，对 syntenic 覆盖无影响）。

### 7.2 人类规模（~3 Gb）验证

字段上限逐项核对：

| 字段 | 上限 | 人类规模评估 |
|---|---|---|
| `pos_start` u32 | 单 contig < 4.3 Gb | 人类最大染色体 ~250 Mb，OK |
| `SeedHit` contig u16 | >65535 contig 报错（build 守卫） | 人类 contig 数远小于此，OK |
| `PgiEntry` 24 B/条 | kmer u128 + pos_start/freq u32×2 | 人类 ~3 Gb 种子量级需实测 |
| positions 位域 | cid 上限 2^20（debug_assert） | OK |
| 总记录数 | `payloads.len() <= u32::MAX` 守卫 | 人类规模需实测是否接近 |

**执行**（待 GRCh38/CHM13 数据）：`pgi build` + `align pgi` 自比对/两染色体
比对，记录构建时间/峰值内存/磁盘/mmap 行为；核对 `.pgi` v2 打包位宽
（`ceil(k/4)` + position 按需字节）在长 contig 下无溢出。与 FastGA（-T8）
对照耗时/内存/覆盖。

### 7.3 完整 adaptamer：变长种子（>k）

FastGA 种子长度范围本身是 [12, KMER=40]，更长匹配靠链化/对齐覆盖。pgr 的
plen ∈ [12, k] 已等价；未做的是"种子 (start, len) 直接携带 >40 变长"——
benchmark 明确非当前优先级（Sakai 剩余差距来自分歧区的 wave 补齐，见
[[../benchmarks/bench-pgi-align-vs-fastga.md]]）。前置 lcp 机制已落地（§3.3）。

### 7.4 对称 adaptamer（`-S`）：不做

FastGA `-S`（FastGA.c:2340，README:199-207 已文档化）：`P1->maxp >
P2->maxp` 时交换 T1/T2，双向种子合并；用两方 adaptamer，结果与 A/B 顺序
基本无关，但通常发现 B 中更多重复比对。README 明确 **synteny 场景不建议，
仅重复结构分析时用**——pgr 的 `align pgi` 是 synteny 用途，且 `sd cross`
当前单向够用，故**不做**（与 §3.4 select 同类判断：专门场景 + 无消费者）。
若将来做对称跨基因组重复检测，按 FastGA.c:2340 语义实现（较小基因组做
种子侧，双向合并）。

**pbit 存储场景实测（2026-08-05，MG1655 ref × Sakai query）**：考虑
`align → PAF → pbit`（CIGAR delta 压缩，见 [[pbit.md]]）时，-S 理论上
可能通过"更多 query 覆盖（含重复）"改善压缩。实测：

- query 覆盖：单向 4,628,115 bp vs 对称 4,671,797 bp（**+0.9%**，重复区
  为主，符合 README）；
- pbit 归档（同名 contig 修正后）：单向 3,009,469 B vs 对称 3,009,470 B
  （**+1 字节，无收益**）。

原因：E. coli 重复区极少（~0.5%），-S 多覆盖的片段要么被 `min_match_len`
过滤，要么在 LZ-diff 下已近最优，CIGAR delta 边际收益≈0。

**重复遮蔽边界（2026-08-05 用户指出）**：以上实测用**未遮蔽** E. coli。
-S 的收益来源就是"更多重复比对"，而真实流程通常带重复遮蔽（`pgi build
--mask` / FastGA `-M` 滤掉 masked 种子）——遮蔽后 -S 的额外比对被移除，
覆盖差异趋零，归档更无差异。**结论（不做）在遮蔽场景下更稳健**。唯一例外
是"故意不遮蔽重复区做 delta 引用"的 pbit 流程，但细菌规模已实测无收益；
真核（重复区占比大）是唯一可能显著场景，需 §7.2 人类数据验证（且应连同
遮蔽与否两个版本一起测）。

### 7.5 多 mask union

FastGA `.1ano` 可叠加多个 soft mask；pgr `fa mask` 单 runlist。低优先便利，
暂缓。

### 7.6 完整 LCP（vlcp 表 / `.pgi` v3 LBYTE）

**现状**：简化版 lcp 已落地（§3.3）——`emit_entry_hits` 用相邻条目
`lcp(prev, cur)` 起步窗口，但**窗口内仍逐条目算 `shared_prefix` 找最大
plen**（O(窗口内条目数)）。

**FastGA 的完整机制**（`new_merge_thread`，FastGA.c:610）：

- GIX 每条记录带 **LBYTE**（与排序流中相邻记录的共享前缀，打包 1 字节，
  `.ktab` 记录 = kmer 后缀 7 B + mask 1 B + **lcp 1 B**）；
- merge 维护 **`vlcp[plen]` 表**：记录"共享前缀 = plen"时对方范围的起点；
  新条目从前一条目的 LBYTE 处继续，**逐碱基递增 plen**，每步范围收缩
  O(1)（指针推进），直到范围空或 plen = KMER——**O(1) 摊销**，不做窗口内
  逐条目比较。

**实现这个东西到底有什么好处**（2026-08-05 分析）：

1. **性能——E. coli 上收益小，人类规模才可能显著（推测，未实测）**：
   merge 瓶颈是磁盘 IO（简化版 A/B 实测仅 107 vs 109 ms），CPU 扫描不是
   主导；3 Gb 级索引 b 侧窗口内条目多（重复区动辄数百），逐条目
   `shared_prefix` 的 CPU 成本才显著，vlcp 的 O(1) 摊销才体现价值。
   **在 §7.2 人类数据到位前，性能收益无法证实**。
2. **`.pgi` v3 与 GIX 记录格式对齐**：补上 lcp 字段后，读取侧
   （PgiMmap/PgiStream）与 FastGA `Kmer_Stream` 形态一致，merge 无需现算
   `shared_prefix`；这是"实现形态对齐"的架构价值，独立于性能。
3. **变长种子（§7.3）的高效前置**：vlcp 表是 [12,40] 长度自适应种子选择
   的引擎（§3.3 结论）；完整版让"每个条目找最长共享前缀"从 O(窗口) 降为
   O(1) 摊销——若未来做变长种子/更长共线检测，这是现成机制。

**实现路径**（若重新评估，分两步；当前不做）：

1. **merge 侧**：`emit_entry_hits` 从 `start = max(min_shared, lcp(prev,
   cur))` 起步后改为**逐碱基递增收缩**（`window(m+1)` 为空即最大 m），
   替代窗口内逐条目 `shared_prefix`——语义等价（window 非空 ⟺ max m ≥ len），
   可先行 A/B；
2. **格式侧（若 1 有收益）**：`.pgi` v3 打包 per-entry LBYTE（排序后相邻
   entry 的 lcp，最后一个为 0），`PgiQuery` 增 `entry_lcp(i)`，merge 直接
   读字段。

**验证**：输出与简化版逐字节一致（等价性测试）；`libs::pgi` + cli_align_pgi
全绿；A/B merge 耗时（E. coli 预期小；人类规模待 §7.2 数据）。

**结论：不做**（2026-08-05 定稿，完整尝试与论证见 §3.6）。完整 LCP
语义可做、性能不可行（正确版仍慢 2.1×）；简化版（§3.3）是最终形态
（流式中性、PSL 零影响、FastGA 语义对齐）。`.pgi` v3 存 LBYTE 属独立
格式演进，若做需 §7.2 人类数据支撑。

## 8. 相关文档

- 索引格式与消费者规划：[[pbit.md]]（多参考节 + .pgi 距离消费者层级）
- FastGA 管线与简化移植评估：[[fastga.md]] §11/§12
- 泛基因组场景：[[ecoli-cohort.md]]、[[paf-pangenome.md]]
