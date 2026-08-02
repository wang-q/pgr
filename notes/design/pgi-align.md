# pgr pgi align：两索引归并比对（设计定稿 + 开发记录）

> 定位：`.pgi` 的第一个比对消费者。输入两个已构建 .pgi，输出 PSL 块，
> 喂给 `pgr pl chainnet`（UCSC 链化由 pgr 承担，见 [[fastga.md]] §12.3 决策 3）。
> 状态：2026-08-02 定稿。质量（chainnet 覆盖）与 FastGA 持平、速度持平、
> 峰值内存低于 FastGA（query 索引 mmap 零拷贝）。
>
> 结构：§0 当前状态 → §1 设计 → §2 验证与基准 → §3 开发历史 →
> §4 已排除方向 → §5 勘误与基准方法 → §6 相关文档。

## 0. 当前状态

**一句话**：`pgr pgi align` 已完整实现"种子归并 → 链化（greedy/tube）→
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

1. indel 复杂区覆盖 ~0.7 kb（Nissle，噪声级）——需 adaptamer 最小种子
   **选择**或 FastGA 式 tube/wave 桥接，见 §4；
2. `dist pgi` / `stat` / `to-hv` 仍全量 `PgiIndex::read`，可复用 `PgiMmap`
   进一步降内存；
3. 人类规模（~3 Gb）未验证：`pos_start` u32、`SeedHit` contig u16、
   `PgiEntry` 24 B/条 等字段上限需按规模复核（§3.4 有上限记录）。

**重要勘误索引**（详情见 §5.2）：

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
5. `pgr pgi align` CLI + 集成测试 + E. coli 三株系验证（§2）。

**不做（未来工作）**：

- lcp 连续传播的完整 adaptamer（固定 k + plen 最大选择已实现；最小种子
  选择未做，见 §4）；
- ~~pbit 内嵌索引段消费~~（已按决策 A 放弃，见 [[pbit.md]]）；
- 人类规模验证（见 §0 剩余工作 3）。

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

- q = query、t = ref（`pgr pgi align <ref> <query>`）；
- 每条链（greedy）/ 每窗（扩展）= 一个块；
- **负链**：qStart/qEnd 必须正链帧、内部 qStarts 必须 RC 帧（与 `psl chain`
  的 `calc_block_score` 一致）——§3.3 记录过整类 '-' 块被静默丢弃的 bug；
- match/mismatch 计数来自扩展（v1 链块为 0，由 `psl chain` 重算）。

### 1.4 CLI

```
pgr pgi align <ref.pgi> <query.pgi> -o out.psl
  [--freq 10] [--min-span 85] [--max-gap 1000] [--band 128] [--merge-gap 5000]
  [--min-shared N] [--workflow greedy|tube] [--parallel 8]
  [--ref-seq ref.fa|2bit] [--query-seq query.fa|2bit]
```

- 两侧索引参数必须一致（复用 `dist pgi` 校验）；
- 无序列文件：每条链一个 PSL 块；有 `--ref-seq/--query-seq`：链细化成
  带真实身份率的块（16 kb 窗口 + 2 kb 重叠滑动）；
- query 索引必须是真实文件（mmap 不支持 stdin/gzip）。

## 2. 验证与基准

### 2.1 测试基线

**912 测试全过**（2026-08-02；5.34 记录 908，+4 为 §3.6 的 mmap 测试）。
回归覆盖：方向解析（fwd/rev）、频率过滤、链化边界（band/gap/span）、
tube anti 序、dedupe 0.95 保留延伸块、最大前缀/扩展范围过滤、负链 RC 帧、
多 contig（RC contig + 2% 突变）、mmap 与全量读入等价性；集成测试
`tests/cli_pgi_align.rs`（identical / RC / mutation / tube）。

### 2.2 当前基准（2026-08-02，真实数据，8 线程，release）

| 对（MG1655 vs） | pgr chainnet 覆盖 | 块数 | pgr 耗时 | FastGA 覆盖/耗时 | pgr 峰值内存 | FastGA 峰值内存 |
|---|---:|---:|---:|---:|---:|---:|
| Sakai | 89.33% | 691 | 0.77s | 89.3% / ~0.7s | **224 MB** | 332 MB |
| EC958 | 86.38% | 756 | 0.81s | 86.3% / ~0.7s | **205 MB** | — |
| Nissle | 85.28% | 1213 | 0.65s | 85.30% / ~0.7s | **207 MB** | — |

> 注：块数/覆盖按当前默认 `min-shared=12`（tube）实测；早期记录的
> 588/794/793 块对应 k/2=20 时代（§3.2）。三对 PSL 与全量读入版
> **逐字节一致**（mmap 改动验证）。BREAK=1000→2000 对齐 FastGA 后
> 实测（§5.2 勘误 7）：块数 +1.5%，syntenic 覆盖 Sakai +0.02%、
> EC958 -0.09%、Nissle ±0.00%（噪声级），耗时/内存持平。
> FastGA 内存实测见 §5.1。

阶段分布（Sakai，`RUST_LOG=debug` 探针）：merge **198 ms / 171 MB**、
chain_tubes **237 ms / 218 MB**、extend **278 ms / 229 MB**，墙钟 0.77 s
（另含索引流式读取、序列加载、PSL 写盘）；merge 命中 2,471,561。

### 2.3 端到端管线验证（2026-08-02）

`pgr pgi align` → `pgr psl to-chain` → `pgr pl chainnet --syn` 全链路，与
FastGA 驱动版本对比 syntenic MAF：

| 输入对 | 指标 | pgr 管线 | FastGA 管线 |
|---|---|---:|---:|
| MG1655 vs Sakai | syntenic 覆盖 | **87.7%**（392 块） | 89.3%（506 块） |
| MG1655 vs Nissle | syntenic 覆盖 | **82.9%**（541 块） | 85.3%（711 块） |

结论：管线端到端可用（0.4s），块结构比 FastGA 更平滑（392 vs 506 /
541 vs 711）；覆盖差 1.6-2.4% 来自分歧区（FastGA 的 wave 能桥接 banded
窗口跳过的低分区间）。角色约定：`pgr pgi align <ref> <query>` 的 PSL 是
q=query/t=ref，FastGA 输出相反，喂 chainnet 前需 `pgr psl swap`。

### 2.4 10 株 cohort 两两验证（45 对）

扩展块身份率矩阵（行×列 = ref×query；块数为扩展块数；默认参数含
`--merge-gap 5000`）：

| pair | 块数 | 身份率 |
|---|---:|---:|
| mg1655–sakai | 862 | 0.9834 |
| mg1655–nissle1917 | 1378 | 0.9746 |
| mg1655–cft073 | 884 | 0.9739 |
| mg1655–e2348_69 | 880 | 0.9740 |
| mg1655–e24377a | 957 | 0.9859 |
| mg1655–ec042 | 934 | 0.9778 |
| mg1655–ec2011c_3493 | 903 | 0.9863 |
| mg1655–ec958 | 958 | 0.9747 |
| mg1655–se11 | 853 | 0.9870 |
| sakai–nissle1917 | 1592 | 0.9736 |
| sakai–cft073 | 1230 | 0.9700 |
| sakai–e2348_69 | 1138 | 0.9718 |
| sakai–e24377a | 1227 | 0.9807 |
| sakai–ec042 | 1240 | 0.9724 |
| sakai–ec2011c_3493 | 1261 | 0.9783 |
| sakai–ec958 | 1122 | 0.9712 |
| sakai–se11 | 953 | 0.9772 |
| nissle1917–cft073 | 1835 | 0.9962 |
| nissle1917–e2348_69 | 1291 | 0.9870 |
| nissle1917–e24377a | 1547 | 0.9736 |
| nissle1917–ec042 | 1621 | 0.9734 |
| nissle1917–ec2011c_3493 | 1608 | 0.9721 |
| nissle1917–ec958 | 1596 | 0.9880 |
| nissle1917–se11 | 1170 | 0.9735 |
| cft073–e2348_69 | 1040 | 0.9864 |
| cft073–e24377a | 1022 | 0.9733 |
| cft073–ec042 | 1121 | 0.9720 |
| cft073–ec2011c_3493 | 1102 | 0.9714 |
| cft073–ec958 | 1151 | 0.9867 |
| cft073–se11 | 816 | 0.9730 |
| e2348_69–e24377a | 981 | 0.9731 |
| e2348_69–ec042 | 1070 | 0.9727 |
| e2348_69–ec2011c_3493 | 1021 | 0.9714 |
| e2348_69–ec958 | 882 | 0.9847 |
| e2348_69–se11 | 868 | 0.9716 |
| e24377a–ec042 | 1200 | 0.9765 |
| e24377a–ec2011c_3493 | 1142 | 0.9911 |
| e24377a–ec958 | 1115 | 0.9726 |
| e24377a–se11 | 1040 | 0.9913 |
| ec042–ec2011c_3493 | 1201 | 0.9761 |
| ec042–ec958 | 1131 | 0.9728 |
| ec042–se11 | 877 | 0.9756 |
| ec2011c_3493–ec958 | 1211 | 0.9734 |
| ec2011c_3493–se11 | 829 | 0.9908 |
| ec958–se11 | 942 | 0.9727 |

全部 45 对 ~60s 完成。分布 97.0-99.6%，与亲缘关系一致（e24377a/se11/
ec2011c_3493 聚类 99.1%+、nissle–cft073 99.6%）。合并块后身份率比无合并
低 ~0.2-0.5%（分歧间隙计入计数，更接近真实）；身份率基于种子链发现的
块（偏保守区段），且整体比 FastGA 高 ~0.5%（banded 局部取精确核心）。

## 3. 开发历史（按主题；日期均为 2026-08-02）

> 历史条目保留原 "5.x" 编号以便对照旧引用。勘误以 §5.2 清单为准，
> 早期条目里的中间结论已被更正，阅读时注意。

### 3.1 链化与扩展引擎演进

- **v1 链块（最小闭环）**：两条目线性归并 + 贪心链化 + 每链一个 PSL 块。
  首个实测：MG1655 自比对主链 1 块覆盖 4,641,650/4,641,652 bp
  （99.9999%），其余为 rRNA/IS 等真实重复（负链 186 块）；vs Sakai
  1019 块（span 覆盖 95.7%）vs FastGA 701 块（99.7%）——真实并集覆盖
  Sakai 75.8% vs 78.2%、Nissle 77.3% vs 77.3%（span 求和有重复计数）。
- **5.2 v2 banded SW 扩展**：1015/1019 块被细化，身份率 98.41%
  （FastGA 97.83%）；并行化 20.5s→2.0s（32 核）。>30 kb 链当时回退为块。
- **5.3 株系验证（v2）**：nissle 97.62%（FastGA 97.09%）、sakai 98.41%
  （97.83%）、ec958 97.60%（FastGA 无基线）。
- **5.4 v3 分窗扩展**：16 kb 窗口 + 2 kb 重叠沿链对角线滑动；自比对主链
  331 窗口身份率精确 1.0000000（5.30M/0），整体 99.93%；Sakai 1093 块
  全部扩展，身份率 98.42%（FastGA 97.83%）。
- **5.6 性能**：banded DP 按带限列迭代（65 列而非 16000 列，246× 浪费）；
  窗口级负载均衡（自比对 37.8s→0.84s、跨株 2.0s→0.66s）。端到端
  （建索引 ×2 + 扩展）1.32s vs FastGA 1.22s（1.08×）。
- **5.7 链合并（--merge-gap）**：IS 插入使对角线平移超 band 断链。
  Sakai 1019→718 块 / 覆盖 4.44→4.51 Mb / 最大 58→83 kb；Nissle
  1634→1259 块 / 4.46→4.54 Mb。块数 -23~30%、覆盖 +1.6~1.8%。
- **5.8 仿射 gap**（M/I/D，open -8 + extend -6）：Sakai 块数 4710→3512、
  Nissle 7451→5609（-25%）；身份率不变，仅 indel 结构更干净。
- **5.11 Myers wavefront 独立移植（负结果）**：wave 核心单元测试通过，
  但单独接入 banded 路径更差（块数 3-4×，覆盖 71.6% / 32.7% 取决于锚点
  策略）——wave 依赖 tube 的锚定上下文；保留为独立实现
  （`src/libs/alignment/wave.rs`），按 FastGA 语义接入见 5.13。
- **5.12 tube 链化移植**（FastGA `align_contigs`）：对角桶（宽 64）+
  anti 归并 + tube（§1.3.3）。常量对齐 FastGA：`BREAK=2000`（anti，
  `-s 1000` 内部翻倍）、`MIN_COV=85`（单轴 = FastGA 170 anti）。
- **5.13 Myers wave + `Local_Alignment` 移植**：forward_wave_mid +
  Myers O(ND) D&C 回溯 + extend_tube（BUCK_ANTI=128；`alast` 每 tube
  独立，并行化替代 FastGA 的逐对桶共享，重叠由 `dedupe_contained` 兜底）。
  对照（tube）：
  banded 基线 862 块/87.7%/1.36s/1.38 GB → tube+wave 643 块/**88.2%**/
  8.7s/425 MB → FastGA 701 块/89.3%/~0.7s。质量逼近、内存 -3×，
  速度是短板（~4 万次 wave 调用，单次 ~0.18ms）。
- **5.14 调用数与 tube 结构对齐**：FastGA 全程仅 1062 次调用（tube 平均
  7.7 kb）；我们 ~4 万次（exact-40 种子把分歧岛桥接成巨型 tube，平均
  30 kb）。优化：tube 并行 + RC/互补预计算 + tube anti 上限 40 kb →
  8.7s→1.7-1.9s。
- **5.15 大 tube 同源门控（曾保留）**：>10 kb tube 做多对角线（9 条 ×
  64 bp 窗）滑窗身份检查，<50% 则跳过；Sakai 88.2%/1.52s、Nissle
  83.9%/1.55s，误杀率≈0。后被 5.22 移除（见 §3.3）。
- **5.18 排序优化**：链排序键打包 u128 + rayon `par_sort_unstable_by_key`
  （741→275 ms）；anti 序修复连锁使调用 58,370→883（FastGA 1062），
  顺序耗时 11.4s→2.7s；总 1.5s→1.0s。
- **5.23 速度与 FastGA 持平**：正链 tube 的 q 复制改 `Cow` 零拷贝
  （extend 470→179 ms）+ tube 形成并行（chain 312→198 ms）→
  **~0.7s**（Sakai 0.69-0.79s）。

### 3.2 种子语义演进

- **5.9 盲目部分匹配：负结果**。无 plen 最大选择的部分种子全部劣于精确
  匹配（MG1655 vs Sakai）：

  | min_shared | records | blocks | identity |
  |---|---|---:|---:|
  | 40（精确） | 862 | 3512 | 0.9834 |
  | 30 | 1130 | 4657 | 0.9827 |
  | 25 | 1302 | 6072 | 0.9810 |
  | 20 | 1614 | 8680 | 0.9781 |
  | 12 | 3375 | 53015 | 0.9496 |

  结论：部分种子的收益依赖 FastGA 的 tube/wave 种子质量机制，不能直接
  移植到贪心链化。
- **5.20 tube 默认 k/2=20**：Sakai 89.1%（+0.2%）、Nissle 84.7%（+0.3%），
  hit +20%；min-shared 30 反而更差。
- **5.29 种子选择移植**（§1.3.2 四项）：Sakai 247 万种子、89.26→89.31%；
  EC958 227 万、86.17→86.36%；Nissle 233 万、84.74→84.98%；内存不升反降
  （种子减半）；tube 默认 floor 定为 12（配合最大选择，12-19 bp 锚点补
  indel 复杂区，与 5.9 的"无最大选择"机制不同）。

### 3.3 bug 修复清单

| 条目 | 症状 | 根因 | 修复 | 回归测试 |
|---|---|---|---|---|
| 5.16 tube 合并顺序 | 高 anti 种子先处理抬高 ahgh，稠密种子 cov 不累计，整管丢弃（Sakai 缺失 55 kb） | 桶内按 (diag,a_pos) 排序，对角线漂移时 anti 序与 a_pos 序相反 | 排序键 (diag 桶, anti)，按 anti 归并 | `tube_merge_uses_anti_order_when_diagonal_drifts` |
| 5.17 dedupe 误删 | 同 tube 连续调用重叠 ~87% 的延伸块被删，丢 3.1 kb 真实覆盖 | 双轴 80% 重叠阈值过低 | 阈值 0.95 | `dedupe_keeps_blocks_that_extend_earlier_ones` |
| 5.22 门控误杀 | Sakai 最大缺失区 7.3 kb（99% 身份、4160 个共享 40-mer）整管被跳过 | 对角线采样步长 ~13 bp，真实 diag 落采样点之间（差 2 bp）时滑窗身份率 ~25% | 直接移除（anti 修复后调用数已 58k→~900，门控只剩误杀） | — |
| 5.30 负链 PSL 帧 | 所有 '-' 块被 `psl chain` 静默丢弃（MAF 中 '-' 块为 0） | qStart/qEnd 用 RC 帧、内部 qStarts 用正链帧，与 UCSC 约定相反，`calc_block_score` 得大额负分 | qStart/qEnd 正链帧、qStarts RC 帧 | `extend_chain_rc_query` |

5.16 修复效果：Sakai 88.2→88.9%、Nissle 83.9→84.4%，缺失 55→25 kb；
5.17：Sakai 覆盖 4,124,601→4,127,886 bp，缺失 24.9→21.7 kb；
5.22：Sakai 89.1→89.3%，缺失 15→7.7 kb（教训：门控类启发式在根因修复
后要及时复核，否则从"省时"变"漏块"）；
5.30：Sakai 89.31→89.33%、EC958 86.36→86.38%、Nissle 84.98→**85.28%**
——Nissle 的"0.32% 差距"绝大部分就是这个 bug。

### 3.4 内存优化时间线（Sakai 峰值，8 线程，另有标注除外）

| 条目 | 变更 | 峰值 |
|---|---|---:|
| 基线 | 双索引全量载入 + 并行 dandc | ~0.96 GB（32T） |
| 5.24 | `drop(hits)` + q `Cow` | 8T 960→639 / 32T ~825 MB |
| 5.27 | `align_to_psl_ext` 按值 + `mem::take` 释放 entries/positions | 875→639 MB（32T） |
| 5.28 | positions u64 位域打包（`pos 32\|cid 20\|strand 1`，12→8 B/条） | 639→607 MB（32T 表） |
| 5.29 | 种子选择（种子减半） | 607→586 MB（32T） |
| 5.31 | ref 索引流式 `PgiStream` | Sakai 604→398 / EC958 512→374 / Nissle 464→378 MB |
| 5.32 | `SeedHit` 24→16 B + 链排序基数化（u128 键 + u32 序数组） | ~381 / 321 / 310 MB（32T） |
| 5.33 | 共享 RC 预计算 + wave 预留收敛（4096→256×width）+ `--parallel 8` 专用池 | 296 / 281 / 284 MB |
| 5.35 | query mmap 零拷贝 `PgiMmap` | **224 / 206 / 210 MB** |

> 5.35 实测对照（同一输入、旧版全量读入）：Sakai 289→224、EC958
> 272→206、Nissle 277→210 MB（-23%）；阶段峰值（Sakai）：merge
> 171 MB、chain_tubes 218 MB、extend 229 MB。FastGA 同输入实测 332 MB
> （§5.1）——pgr 峰值已低于 FastGA。

结构上限记录：`SeedHit` contig u16（5.34 起 `build_from_seqs` 守卫
>65535 报错）；`pack_position` cid 上限 2^20（debug_assert）；
`PgiEntry` 24 B/条（kmer u128 + pos_start/freq u32×2）是 positions 位域
之后的下一优化对象。

### 3.5 性能优化时间线（Sakai，8 线程）

| 条目 | 变更 | 耗时 |
|---|---|---:|
| v3 分窗后 | 窗口级负载均衡（§3.1） | ~0.97s |
| 5.18 | 链排序 u128 键 + rayon（741→275 ms）+ 调用数骤降 | ~1.0s |
| 5.19 | pgi 批量解析（1 MB 分块、按 rec_size 对齐；加载 0.7→0.5s） | ~1.0-1.2s |
| 5.21 | merge 并行（4096 条分块，139→61 ms） | ~1.1s |
| 5.23 | q `Cow` + tube 并行 | ~0.7s |
| 5.34 | 阶段探针实测：merge 193 / chain_tubes 236 / extend 245 ms | ~0.89s |
| 当前 | 探针：merge 198 / chain_tubes 244 / extend 273 ms | 0.78s |

> 5.21 同时试过 pgi 解析并行化：**无益**（瓶颈是磁盘读取而非解析，已回退，
> 见 §4）。与 FastGA 的剩余差距是 wave 每调用成本（~0.47ms，Rust vs C）。

### 3.6 结构与读取演进（近期）

- **5.19 pgi 批量解析**：`PgiIndex::read` 由逐记录 `read_exact`（trait 对象
  虚拟分发）改为 1 MB 分块 + 切片解析，加载 0.7→0.5s。
- **5.31 `PgiStream`**：ref 侧流式分块、按条目批量产出、条目不跨批；
  `merge_seed_hits_from_stream` 用 rayon `par_bridge`；修正 5.24
  "惰性加载需新依赖或 unsafe"的结论（纯 std 流式即可）。
- **5.34 多 contig 支持**：`SeedHit` contig u16 + 守卫；3 contig 回归
  （c1 正 20000/20000、c2 负 15000/15000 RC 正确、c3 正 9858/10000 含
  2% 突变）；阶段耗时探针。
- **5.35 `PgiMmap`（mmap 零拷贝）**：query 侧记录驻留映射页，条目经
  packed k-mer 字节二分定位、位置按需解码；`PgiQuery` trait 统一
  resident/mmap 视图。修两个语义坑：prefix 哨兵键 `hi=2^(2k)` 需 clamp
  （打包后高位移出）；`entry_range` 返回记录区间、组内记录不能当独立
  条目（新增 `entry_next` 按组推进）。912 测试（新增 roundtrip/区间/
  截断/merge 等价性 4 项）。
- **5.36 FastGA 内存勘误**：见 §5.1。

## 4. 已排除方向（避免重试）

| 方向 | 结论 | 原因 / 出处 |
|---|---|---|
| adaptamer 部分种子（盲目，无最大选择） | 不可用 | 弱种子大量假阳性，min-shared 12 → 53015 块/0.9496（§3.2 5.9） |
| wave 单独接入 banded 路径 | 不可用 | unit-cost 波前无 indel 偏好 + 贪心延伸：块数 3-4×，覆盖 71.6% / 32.7%（锚点策略）；wave 依赖 tube 锚定上下文（§3.1 5.11） |
| 大 tube 同源门控（多对角线滑窗） | 已移除 | 采样漏真实对角线（差 2 bp 时窗口身份 ~25%）→ 误杀整管；根因修复后收益消失（§3.3 5.22） |
| 中心对角线滑窗身份率门控 | 不可用 | 覆盖 88.2%→70.9%，误杀薄保守区/偏移对角线（§3.1 5.14） |
| 种子覆盖密度门控 | 不可用 | 零块与生产性 tube 分布重叠，无法干净区分（5.14） |
| 种子邻近门控（amid ±300bp） | 无效 | 失败调用的 amid 本来就在种子附近（5.14） |
| CHAIN_BREAK 调小（300/100） | 更差 | tube 碎片化、重叠调用更多，质量 87.3%（5.14） |
| tube 默认 min-shared 30 | 更差 | 部分匹配噪声未被抑制（§3.2 5.20） |
| pgi 解析并行化 | 无益 | 瓶颈是磁盘读取而非解析（§3.5 5.21） |
| wave 死代码路径（D&C 回溯每次调用全跑） | 已改 | FastGA 的 `dandc_nd` 是死代码；我们曾每次调用跑完整 D&C（5.13，优化后非每次） |

## 5. 勘误与基准方法

### 5.1 FastGA 内存实测与 "~0 MB" 勘误（2026-08-02）

此前各节对照表把 FastGA 峰值内存记为 "~0 MB"，理由是"FastGA 用 mmap 保持
低 RSS"。直接实测（MG1655 vs Sakai，-T8，`/usr/bin/time -v`）推翻该结论：

| 阶段 | 进程 | 峰值 RSS |
|---|---|---:|
| FAtoGDB（FASTA → 2bit GDB） | 子进程 | ~7 MB |
| GIXmake -T8（k-mer 排序建索引） | 子进程 | ~160 MB |
| FastGA -psl 比对（-T8，主流程） | 主进程 | **332 MB** |

源码核查（仓库内 FASTGA-main）：

- **全项目无任何 mmap 调用**（rg 无命中）。GIX 经 `Open_Kmer_Stream`
  （libfastk.c:785）用 `read()` + `STREAM_BLOCK` 缓冲**流式读取**；
  GDB 序列默认 `seqstate == EXTERNAL`（不整库载入，每线程 fopen 一个
  .bps 句柄按需 fseeko），`Load_Sequences` 整库载入被 `#ifdef LOAD_SEQS`
  编译开关关闭；
- 比对阶段真实内存大头：`(nelmax+1)*swide` 种子排序数组
  （`swide = 2*DBYTE + JCONT + 2` ≈ 10 B，种子 ~360 万 → ~36 MB）、
  GIX 流缓冲、N 线程对齐缓冲与 contig 序列缓存。

"~0 MB" 的来源：**无测量记录**，系 9930353 起从参考笔记 fastga.md §9.4
的"可 mmap 整库"能力描述错误外推（该描述本身也不准确：GDB.c 的 COMPRESSED
态是 Malloc + read 整库读入，`Load_Sequences` 无 mmap），5.13/5.18/5.23/
5.24/5.27 沿用未核。

修正后的对照（8 线程，仅比对进程）：FastGA **332 MB** vs pgr mmap 版
**224 MB** / 全量读入版 289 MB——**pgr 的 mmap 版峰值已低于 FastGA**。

### 5.2 结论勘误清单

1. **Nissle 基线（5.25）**：曾怀疑 `nissle1917.fa.gz` 被替换、FastGA 基线
   无效 → 文件逐字节相同、基线有效。过程中验证"对含 indel 的对齐不能用
   naive 偏移身份判断"（该区域每 ~300 bp 一个 indel，同偏移身份 ~25% 但
   共享 40-mer 8990/10258 ≈ 99% 相关）。
   另：5.25 进一步抽查最大缺失区（t 2058021-2062339，2.2 kb）确认对齐
   本身正确，只是 chainnet 在重复区过滤了该块。
2. **EC958 基线（5.26）**：`ec958.fa.gz` 当天改过，但 FastGA 输出与当前
   文件坐标一致，基线有效；我们 86.2% vs FastGA 86.3%（缺失 18.9 kb、
   多覆盖 12.6 kb，净差很小）。
3. **Nissle 0.32% 差距归因（5.27→5.30）**：曾归因 chainnet 对重复区单块
   过滤 + indel 复杂区 → 主因是负链 PSL 坐标帧 bug（§3.3），修复后差
   0.015%（~0.7 kb）。
4. **`is_minimal` 语义（5.20→5.29）**：曾读作噪声抑制 → canonical 方向
   判断；真正的噪声抑制是 plen 最大选择 + 扩展范围过滤。
5. **FastGA 内存（5.13→5.36）**：曾记 "~0 MB" → 实测 332 MB（§5.1）。
6. **验证数据（5.35→5.36）**：早期 mmap 验证曾用合成 2×2 Mb 随机序列 →
   统一为 tests/genome 真实数据重测（数字以 §2.2 为准）。
7. **tube `CHAIN_BREAK` 口径（5.12→本次核对）**：曾把 `BREAK=1000`
   （未加倍值）用于 anti 空间间隔比较；FastGA 的 `-s 1000` 在内部翻倍为
   `CHAIN_BREAK=2000`（anti）。已修正为 2000（`MIN_COV=85` 因 cov 用
   单轴投影，与 FastGA 的 170 anti 等价，不动）。修正后实测：块数
   +1.5%，三对 syntenic 覆盖 ±0.1% 内（Sakai +0.02%、EC958 -0.09%、
   Nissle ±0.00%），耗时/内存持平（§2.2 数字以修正后为准）。

### 5.3 基准方法

- 数据：`tests/genome/{mg1655,sakai,nissle1917,ec958}.fa.gz`（另有
  cft073/e2348_69/e24377a/ec042/se11 等 cohort 株）；
- 命令：`pgr pgi align <ref> <query> --ref-seq --query-seq
  --workflow tube`（8 线程默认），release；`/usr/bin/time -v` +
  `RUST_LOG=debug` 阶段探针（merge/chain_tubes/extend + VmHWM）；
- 覆盖：`pgr psl to-chain` → `pgr pl chainnet --syn` syntenic 覆盖；
- 端到端（含建索引 ×2）见 [[benchmarks/bench-pgi-align-vs-fastga.md]]；
- 10 株 cohort 验证见 [[benchmarks/dist-cohort-validation.md]]（引用
  §2.4 的身份率矩阵）。

## 6. 相关文档

- 索引格式与消费者规划：[[pbit.md]]（多参考节 + .pgi 距离消费者层级）
- FastGA 管线与简化移植评估：[[fastga.md]] §11/§12
- 泛基因组场景：[[ecoli-cohort.md]]、[[paf-pangenome.md]]
