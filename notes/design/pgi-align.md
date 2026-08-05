# pgr align pgi：两基因组归并比对（设计定稿 + 开发记录）

> 定位：`.pgi` 的第一个比对消费者。输入两个已构建 .pgi，输出 PSL 块，
> 喂给 `pgr pl chainnet`（UCSC 链化由 pgr 承担，见 [[fastga.md]] §12.3 决策 3）。
> 状态：2026-08-04 更新。质量（chainnet 覆盖）与 FastGA 持平（差
> 0.0-0.015%）、阶段耗时持平（~0.8s vs ~0.7s）、峰值内存更低（224 vs
> 332 MB，query 索引 mmap 零拷贝）；2026-08-04 端到端复测反超 FastGA
> ~2.3×（见 [[../benchmarks/bench-pgi-align-vs-fastga.md]]）。
>
> 结构：§0 当前状态 → §1 设计 → §2 验证与基准 → §3 开发历史 →
> §4 已排除方向 → §5 勘误与基准方法 → §6 FastGA 功能差距 →
> §7 相关文档。

## 0. 当前状态

**一句话**：`pgr align pgi` 已完整实现"种子归并 → 链化（greedy/tube）→
扩展（banded 仿射 gap / Myers wave）→ PSL"，MG1655 vs Sakai/EC958/Nissle
三项基准的 chainnet 覆盖与 FastGA 持平（差 0.0-0.015%）、耗时持平
（~0.8s vs ~0.7s）、峰值内存更低（224 vs 332 MB，Sakai）。

**实现构成**（2026-08-02 快照）：

- 索引读取：ref 走 `PgiStream` 流式分块；query 走 `PgiMmap` mmap 零拷贝
  （FastGA GIX 模型）；merge 经 `PgiQuery` trait 统一两种视图；
- 种子语义：最大共享前缀（plen）+ 扩展范围频率过滤 + canonical 去重
  （FastGA `new_merge_thread` 语义，§1.3.2）；
- 链化：greedy（默认）/ tube（`--workflow tube`，FastGA `align_contigs`）；
- 扩展：banded 仿射 gap（greedy 默认）与 mid-line Myers wave（tube）；
- PSL：负链 qStart/qEnd 用正链帧、内部 qStarts 用 RC 帧（UCSC 约定）。

**默认参数**（对齐 FastGA：-f 10 / -c 85 / -s 1000，本实现为未加倍空间）：

| 参数 | 默认 | 说明 |
|---|---|---|
| `-k` / `--smer` / `--window` | 40 / 8 / 5 | 与 FastGA GIX (12,8) 同源 |
| `-f` / `--freq` | 10 | 任一侧频率超限即跳过 |
| `-c` / `--min-span` | 85 | FastGA CHAIN_MIN |
| `-s` / `--max-gap` | 1000 | FastGA CHAIN_BREAK/2 |
| `--band` | 128 | 对角线带半宽 |
| `--merge-gap` | 5000 | 相邻共线链合并阈值（IS 元件断链） |
| `--min-shared` | tube=12 / greedy=k | tube 取 FastGA plen 下限 |
| `--workflow` | greedy | tube 需要 `--ref-seq/--query-seq` |
| `--parallel` | 8 | 专用 rayon 池（FastGA `-T` 默认） |

> 口径说明：FastGA 的 `-c`/`-s` 在内部**翻倍**为 anti 空间值
> （`CHAIN_MIN=170`、`CHAIN_BREAK=2000`，源码注释 "2x in anti-diagonal
> space"）。`-c/-s/--band` 只作用于 greedy 链化（单轴语义）；tube 路径用
> FastGA 语义的硬编码常量：`MIN_COV=85`（单轴口径 = FastGA 170 anti）、
> `BREAK=2000`（anti，= `-s 1000` 翻倍）、`BUCK_ANTI=128`。

**当前基准**（2026-08-02，tests/genome 真实数据，8 线程，release，详见 §2.2）：

| 对（MG1655 vs） | pgr 覆盖 | FastGA 覆盖 | pgr 耗时 | FastGA 耗时 | pgr 峰值内存 | FastGA 峰值内存 |
|---|---:|---:|---:|---:|---:|---:|
| Sakai | 89.33% | 89.3% | 0.77s | ~0.7s | **224 MB** | 332 MB |
| EC958 | 86.38% | 86.3% | 0.81s | ~0.7s | **205 MB** | — |
| Nissle | 85.28% | 85.30% | 0.65s | ~0.7s | **207 MB** | — |

**剩余工作**：

见 §7 未来方向（indel 复杂区扩展端、人类规模验证、`dist pgi` / `stat` /
`to-hv` 复用 `PgiMmap`、完整 adaptamer 等，按价值排序）。

**2026-08-03 新增**（对齐 FastGA，详见 §6）：

- `pgr pgi build --mask`：跳过 soft-mask 区（FASTA 小写 / 2bit mask_blocks
  转 N），抑制重复/低复杂度区种子（FastGA `-M`）；
- `pgr align pgi` 单输入自比对：`drop_self_hits` 过滤完全自身命中（同
  contig 同位置同方向，FastGA 跳过 diag=0 语义），保留内部重复与反向重复；
  MG1655 自比对 689 块、无全长主链、无完全自身子块。

**重要勘误索引**（详情见 §5.1）：

- FastGA 内存 "~0 MB" → 实测 332 MB（§5.1）；
- Nissle 基线曾误判无效 → 有效；其 0.32% 差距曾归因 chainnet 过滤 →
  实为负链 PSL 坐标 bug（§3.3）；
- `is_minimal` 曾读作噪声抑制 → canonical 方向判断；
- 早期验证曾用合成随机序列 → 统一为 tests/genome 真实数据。

## 1. 设计

### 1.1 范围（当前做/不做）

**做**：

1. 两个排序 .pgi 流的归并 → 种子命中（plen 最大共享前缀 + 频率过滤）；
2. anti-diagonal 空间链化 → 链（greedy 贪心 / tube 两种语义）；
3. 链扩展：banded 仿射 gap 局部比对（greedy）或 Myers wave（tube），无
   序列输入时每条链输出一个 PSL 块；
4. ref 流式 + query mmap 读取（E. coli 规模起不再整体载入内存）；
5. `pgr align pgi` CLI + 集成测试 + E. coli 三株系验证（§2）。

**不做（未来工作）**：

- lcp 连续传播的完整 adaptamer（固定 k + plen 最大选择已实现；最小种子
  选择未做，见 §4）；
- ~~pbit 内嵌索引段消费~~（已按决策 A 放弃，见 [[pbit.md]]）；
- 人类规模验证（见 §7）。

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
4. **floor = 12**：tube 默认 `min-shared` 为 FastGA 的 plen 下限 12；
   greedy 默认 = k（精确匹配）。

#### 1.3.3 链化

- **greedy（默认）**：按 (contig_a, contig_b, 方向) 分组，命中按
  (diag, pos_a) 排序后贪心延伸（`|diag−均线| ≤ band` 且 Δpos ≤ max_gap），
  跨度 ≥ min_span 才保留；`--merge-gap` 合并相邻共线链（IS 元件等对角线
  平移断链）。
- **tube**：种子按对角线分桶（宽 64）→ 相邻桶对按 anti 归并（排序键
  (diag 桶, anti)，§3.3 修过顺序 bug）→ tube 维护 anti 覆盖与对角线范围，
  种子 anti 间隔超 `CHAIN_BREAK`（2000 bp，FastGA 内部值）断开、覆盖达
  `CHAIN_MIN`（85 bp，单轴口径 = FastGA 170 anti）触发。tube 扩展用
  mid-line wave（BUCK_ANTI=128 滑动），每个 tube 独立 `alast`（并行化
  替代 FastGA 的逐对桶共享）+ 输出端 `dedupe_contained`（0.95 阈值，
  §3.3 修过误删）。

#### 1.3.4 PSL 输出（UCSC 约定，负链是坑）

- q = query、t = ref（`pgr align pgi <ref> <query>`）；
- 每条链（greedy）/ 每窗（扩展）= 一个块；
- **负链**：qStart/qEnd 必须正链帧、内部 qStarts 必须 RC 帧（与 `psl chain`
  的 `calc_block_score` 一致）——§3.3 记录过整类 '-' 块被静默丢弃的 bug；
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
  [--freq 10] [--min-span 85] [--max-gap 1000] [--band 128] [--merge-gap 5000]
  [--min-shared N] [--workflow greedy|tube] [--parallel 8] [--keep-index]
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
- 无序列：每条链一个 PSL 块；有序列：链细化成带真实身份率的块（16 kb
  窗口 + 2 kb 重叠滑动）；
- query 索引（或现场建的临时索引）必须是真实文件（mmap 不支持
  stdin/gzip）。

## 2. 验证与基准

### 2.1 测试基线

**912 测试全过**（2026-08-02 基线；后续 mmap 测试 +4）。
回归覆盖：方向解析（fwd/rev）、频率过滤、链化边界（band/gap/span）、
tube anti 序、dedupe 0.95 保留延伸块、最大前缀/扩展范围过滤、负链 RC 帧、
多 contig（RC contig + 2% 突变）、mmap 与全量读入等价性；集成测试
`tests/cli_align_pgi.rs`（identical / RC / mutation / tube + 序列直入 /
复用 / 混用 / 校验拒绝）。

> **2026-08-05 更新**：`align pgi` 22 个 CLI 测试、`libs::pgi` 55 个单测
> 全通过，`cargo test` 全量 1255 通过（audit 后新增 crafted 索引/负链帧/
> 数据安全等回归，见 [[../audit/audit-pgi-align.md]]）。

### 2.2 当前基准（2026-08-02，真实数据，8 线程，release）

| 对（MG1655 vs） | pgr chainnet 覆盖 | 块数 | pgr 耗时 | FastGA 覆盖/耗时 | pgr 峰值内存 | FastGA 峰值内存 |
|---|---:|---:|---:|---:|---:|---:|
| Sakai | 89.33% | 691 | 0.77s | 89.3% / ~0.7s | **224 MB** | 332 MB |
| EC958 | 86.38% | 756 | 0.81s | 86.3% / ~0.7s | **205 MB** | — |
| Nissle | 85.28% | 1213 | 0.65s | 85.30% / ~0.7s | **207 MB** | — |

> 注：块数/覆盖按当前默认 `min-shared=12`（tube）实测；早期记录的
> 588/794/793 块对应 k/2=20 时代（§3.2）。三对 PSL 与全量读入版
> **逐字节一致**（mmap 改动验证）。BREAK=1000→2000 对齐 FastGA 后
> 实测（§5.1 勘误 5）：块数 +1.5%，syntenic 覆盖 Sakai +0.02%、
> EC958 -0.09%、Nissle ±0.00%（噪声级），耗时/内存持平。
> FastGA 内存实测见 §5.1。

阶段分布（Sakai，`RUST_LOG=debug` 探针）：merge **198 ms / 171 MB**、
chain_tubes **237 ms / 218 MB**、extend **278 ms / 229 MB**，墙钟 0.77 s
（另含索引流式读取、序列加载、PSL 写盘）；merge 命中 2,471,561。

### 2.3 端到端管线验证（2026-08-02）

`pgr align pgi` → `pgr psl to-chain` → `pgr pl chainnet --syn` 全链路，与
FastGA 驱动版本对比 syntenic MAF：

| 输入对 | 指标 | pgr 管线 | FastGA 管线 |
|---|---|---:|---:|
| MG1655 vs Sakai | syntenic 覆盖 | **87.7%**（392 块） | 89.3%（506 块） |
| MG1655 vs Nissle | syntenic 覆盖 | **82.9%**（541 块） | 85.3%（711 块） |

结论：管线端到端可用（0.4s），块结构比 FastGA 更平滑（392 vs 506 /
541 vs 711）；覆盖差 1.6-2.4% 来自分歧区（FastGA 的 wave 能桥接 banded
窗口跳过的低分区间）。角色约定：`pgr align pgi <ref> <query>` 的 PSL 是
q=query/t=ref，FastGA 输出相反，喂 chainnet 前需 `pgr psl swap`。

> **2026-08-04 复测**：端到端（建索引 ×2 + 比对）pgr 1.67 s vs FastGA
> 3.86 s，反超 ~2.3×（初测 1.08× 持平），见
> [[../benchmarks/bench-pgi-align-vs-fastga.md]]。本表为早期管线快照，
> 当前 chainnet 覆盖以 §2.2 为准（Sakai 89.33% / Nissle 85.28%）——
> 87.7%/82.9% 与 89.33%/85.28% 的差距来自后续 merge-gap、种子选择与
> 负链 PSL 修复（§3.2/§3.3）。

### 2.4 10 株 cohort 两两验证（45 对）

扩展块身份率矩阵（初测 2026-08-02；行×列 = ref×query，块数为扩展块数，
默认参数含 `--merge-gap 5000`）：

- 分布 97.0-99.6%，与亲缘关系一致（e24377a/se11/ec2011c_3493 聚类
  99.1%+、nissle–cft073 99.6%）；合并块后身份率低 ~0.2-0.5%；整体比
  FastGA 高 ~0.5%（banded 局部取精确核心）。
- **2026-08-04 复测**（tube 链排序键 + syncmer 去重修复后）：块数整体下降
  （如 mg1655–sakai 862 → 791），身份率分布几乎不变（0.9702–0.9960）。
- 完整 45 对矩阵与复测数据见 [[../benchmarks/dist-cohort-validation.md]]。

## 3. 开发历史（里程碑；日期均为 2026-08-02）

> 早期逐条记录（v1→v3、5.x 编号）已随功能稳定而精简，仅保留里程碑与结论；
> 完整细节见 git 历史与 [[../audit/audit-pgi-align.md]]。

### 3.1 链化与扩展引擎演进

- **v1 链块 → v2 banded SW 扩展 → v3 分窗扩展**（16 kb 窗口 + 2 kb 重叠
  滑动）：自比对主链身份率精确 1.0000000（5.30M/0），跨株身份率 98.42%
  （FastGA 97.83%）。
- **性能关键点**：banded DP 按带限列迭代（16000 列里只有 65 列有效，
  246× 浪费）；窗口级 rayon 负载均衡（自比对主链 332 窗口是单线程长尾 →
  摊平后自比对 37.8s→0.84s、跨株 2.0s→0.66s）。
- **仿射 gap**（open -8 + extend -6）：块数 -25%，身份率不变，仅 indel
  结构更干净。
- **Myers wavefront 移植**（`src/libs/alignment/wave.rs`，FastGA align.c）：
  - 单独接入 banded 路径**不可用**（块数 3-4×、覆盖 71.6%/32.7%）——
    wave 依赖 tube 的锚定上下文；
  - 按 FastGA `align_contigs` 语义接入 tube：643 块/88.2%/425 MB vs
    FastGA 701/89.3%/~0.7s，质量逼近、内存 -3×；tube 并行 + anti 上限
    40 kb + 链排序 u128/rayon 后调用数 58k→~900，端到端 ~0.7s。

### 3.2 种子语义演进

- **5.9 盲目部分匹配（无 plen 最大选择）：不可用**——min-shared 12 产生
  53015 块/0.9496 假阳性爆炸；部分种子收益依赖 FastGA 的 tube/wave 机制，
  不能直接移植到贪心链化。
- **5.29 plen 最大选择 + 扩展范围过滤 + canonical 去重落地**（§1.3.2）：
  种子减半、内存不升反降，三对覆盖 +0.05~0.24%（Sakai 89.26→89.31% 等）；
  tube 默认 floor 定为 12（12-19 bp 锚点补 indel 复杂区，与 5.9 的
  "无最大选择"机制不同）。

### 3.3 bug 修复要点

| 条目 | 影响 | 修复 |
|---|---|---|
| tube 合并顺序 | Sakai 缺失 55 kb | 排序键 (diag 桶, anti)，按 anti 归并 |
| dedupe 误删 | 丢 3.1 kb 真实覆盖 | 双轴重叠阈值 0.80→0.95 |
| 大 tube 同源门控 | 误杀整管（7.3 kb，99% 身份） | 根因修复后直接移除（教训：门控类启发式要复核） |
| 负链 PSL 帧 | 所有 '-' 块被 chainnet 丢弃（Nissle 0.32% 差距主因） | qStart/qEnd 正链帧、qStarts RC 帧（UCSC 约定） |

### 3.4 内存与性能优化结论

- 峰值内存：0.96 GB（双索引全量）→ **224/206/210 MB**（Sakai/EC958/Nissle；
  ref 流式 `PgiStream` + query mmap `PgiMmap` + positions 位域打包 + 种子
  减半），已低于 FastGA 实测 332 MB。
- 耗时：37.8s → **~0.7-0.8s**；阶段分布（Sakai）merge 198 / chain_tubes
  244 / extend 273 ms。
- 结构上限记录：`SeedHit` contig u16（build 守卫 >65535 报错）、
  `pack_position` cid 上限 2^20、`PgiEntry` 24 B/条——人类规模复核见 §7。
- 命令迁移（2026-08-03）：`pgr pgi align` → 顶层 `pgr align pgi`（`pgr pgi`
  收敛为纯索引管理 build/stat/to-hv）；基因组输入自动建索引（sibling 复用 +
  mtime 失效 + 参数一致性校验），`.pgi` 输入配 `--ref-seq/--query-seq`
  做 contig 校验。

## 4. 已排除方向（避免重试）

| 方向 | 结论 | 原因 |
|---|---|---|
| adaptamer 部分种子（盲目，无最大选择） | 不可用 | 弱种子大量假阳性：min-shared 12 → 53015 块/0.9496 |
| wave 单独接入 banded 路径 | 不可用 | 无 indel 偏好 + 贪心延伸；wave 依赖 tube 锚定上下文（§3.1） |
| 大 tube 同源门控（多对角线滑窗） | 已移除 | 采样漏真实对角线误杀整管；根因修复后收益消失 |
| 中心对角线滑窗身份率门控 | 不可用 | 覆盖 88.2%→70.9% |
| 种子覆盖密度 / 邻近门控 | 不可用/无效 | 无法干净区分生产性 tube |
| CHAIN_BREAK 调小（300/100） | 更差 | tube 碎片化，质量 87.3% |
| tube 默认 min-shared 30 | 更差 | 部分匹配噪声未被抑制 |
| pgi 解析并行化 | 无益 | 瓶颈是磁盘读取而非解析 |
| wave D&C 回溯每次调用全跑 | 已改 | FastGA 的 `dandc_nd` 是死代码 |

## 5. 勘误与基准方法

### 5.1 关键勘误清单

1. **FastGA 内存 "~0 MB" → 实测 332 MB**（2026-08-02 直接实测）：FAtoGDB
   ~7 MB、GIXmake ~160 MB、FastGA 比对主进程 332 MB；源码核查全项目无
   mmap（GIX 流式 read + GDB EXTERNAL 文件态）。pgr mmap 版 224 MB 已低于
   FastGA。
2. **Nissle 基线曾误判无效 → 有效**（文件逐字节相同）；其 0.32% 差距主因
   是负链 PSL 坐标帧 bug（§3.3），非 chainnet 过滤。
3. **indel 复杂区不能用 naive 偏移身份判断**：每 ~300 bp 一个 indel，同偏移
   身份 ~25% 但共享 40-mer ≈ 99%（8990/10258）——种子不缺，缺口在扩展端
   （§7）。
4. **`is_minimal` 是 canonical 方向判断**，不是噪声抑制；真正的噪声抑制是
   plen 最大选择 + 扩展范围过滤。
5. **tube `CHAIN_BREAK` 口径**：FastGA `-s 1000` 内部翻倍为 anti 空间
   2000；修正后块数 +1.5%、三对覆盖 ±0.1% 内。
6. **验证数据统一**：早期 mmap 验证用合成随机序列 → tests/genome 真实数据
   （数字以 §2.2 为准）。

### 5.2 基准方法

- 数据：`tests/genome/{mg1655,sakai,nissle1917,ec958}.fa.gz`（另有
  cft073/e2348_69/e24377a/ec042/se11 等 cohort 株）；
- 命令：`pgr align pgi <ref> <query> --ref-seq --query-seq --workflow tube`
  （8 线程默认），release；`/usr/bin/time -v` + `RUST_LOG=debug` 阶段探针
  （merge/chain_tubes/extend + VmHWM）；
- 覆盖：`pgr psl to-chain` → `pgr pl chainnet --syn` syntenic 覆盖；
- 端到端（含建索引 ×2）见 [[benchmarks/bench-pgi-align-vs-fastga.md]]；
- 10 株 cohort 验证见 [[benchmarks/dist-cohort-validation.md]]（引用
  §2.4 的身份率矩阵）。

## 6. FastGA 功能差距

> 对照 [[fastga.md]]（参考笔记）与 `pgr align pgi` 现状，记录 FastGA 中相对
> 重要、pgr 尚未实现的功能；已落地的三项见 §0，此处仅作状态对照。
> 日期：2026-08-03。

### 6.1 状态总览

| 功能 | 状态 |
|---|---|
| soft mask 感知的种子发现（`-M`） | **已落地**（`pgr pgi build --mask`） |
| 自比对模式（`FastGA A`） | **已落地**（`pgr align pgi` 单输入 self） |
| PAF `cs:Z` 输出（`-pafs/S`） | **已落地**（`pgr maf to-paf`） |
| select 表达式 | 暂缓 |
| Gap_Improver | 暂缓 |
| 多 mask union | 暂缓 |
| `-S` 对称 adaptamer | 暂缓（专门场景） |
| trace points / `.1aln` 紧凑存储 | 不做 |
| ALNchain | 不做 |
| GDB / scaffold / 完整 GIX 分片 | 不做 |

### 6.2 当前差距（未完成）

| FastGA 功能 | pgr 现状 | 状态 | 说明 |
|---|---|---|---|
| **select 表达式**（只比对选定 contig/区间） | 无（需 fa range + 子索引间接实现） | 暂缓 | 低优先级 |
| **Gap_Improver**（wave 后 gap 区二次精修） | banded 仿射 gap 已覆盖；wave 路径无等价物 | 暂缓 | 质量微调，收益不确定 |
| **多 mask union**（.1ano 可叠加） | `fa mask` 单 runlist | 暂缓 | 低优先级 |
| **`-S` 对称 adaptamer**（双输入 A vs B 双向种子；FastGA 未文档化选项，仅 V1.5 源码支持） | pgr 双输入同样单向（canonical 半方向单发，`A B` ≠ `B A`；对称需双方向合并，未实现） | 暂缓 | 专门场景（对称的跨基因组重复/结构分析）才有价值；`sd` 的 cross 目前单向够用 |
| **trace points / ONEcode `.1aln` 紧凑存储** | PSL/MAF + BGZF | 不做 | 人类规模才需要，见 fastga.md §10 |
| **ALNchain（.1aln 链化）** | UCSC chain/net（更标准） | 不做 | — |
| **GDB 格式 / scaffold 语义 / 完整 GIX 分片** | pgr 2bit + `.pgi` | 不做 | 格式对比见 fastga.md §9 |

## 7. 未来方向

> 2026-08-05 复核。已落地部分稳定（§0-§5），剩余工作按价值排序：

1. **indel 复杂区覆盖缺口（~0.7 kb，Nissle，噪声级）——扩展端优先**。
   该区域每 ~300 bp 一个 indel，但共享 40-mer ≈ 99%（§5.1 勘误 3）：
   **种子不缺，缺在扩展端**（banded 窗口跳过低分区间，FastGA 的 wave 能
   桥接）。候选：FastGA 式 tube/wave 桥接、Gap_Improver（§6.2 暂缓）。
   benchmark 已明确 lcp/adaptamer 变长种子不是当前优先级
   （[[../benchmarks/bench-pgi-align-vs-fastga.md]]）。
2. **人类规模（~3 Gb）验证**：`pos_start` u32（单 contig > 4.3 Gb 不支持）、
   `SeedHit` contig u16（>65535 已守卫）、`PgiEntry` 24 B/条 等字段上限
   需按规模复核；`.pgi` 存储与 mmap 读取在人类规模的 IO 行为未实测。
3. **`dist pgi` / `stat` / `to-hv` 复用 `PgiMmap`**：目前仍全量
   `PgiIndex::read`，可复用 mmap 进一步降内存（align 已落地该路径）。
4. **完整 adaptamer（lcp 连续传播，变长种子 >k）**：当前种子长度上限 = k
   （40）；FastGA 靠排序流相邻条目 lcp 扩展到任意长度。仅在需要更高敏感度/
   更长共线检测时值得做（benchmark 结论见第 1 条）。
5. **对称 adaptamer（`-S`）**：双方向种子合并；专门场景（对称的跨基因组
   重复/结构分析），`sd` cross 目前单向够用（§6.2）。
6. **select 表达式 / 多 mask union**：低优先级便利功能（§6.2 暂缓项）。

## 8. 相关文档

- 索引格式与消费者规划：[[pbit.md]]（多参考节 + .pgi 距离消费者层级）
- FastGA 管线与简化移植评估：[[fastga.md]] §11/§12
- 泛基因组场景：[[ecoli-cohort.md]]、[[paf-pangenome.md]]
