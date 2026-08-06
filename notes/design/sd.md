# pgr sd：分段重复（SD）检测与分析管线

> 实现笔记。BISER 本身的分析见 [[../references/biser.md]]；本文档记录 pgr 侧
> SD 管线的实现、验证与设计决策。日期：2026-08-03。

## 1. 命令与数据流

`pgr sd` 共 7 个子命令（对应 BISER 流程，检测阶段采用外部比对路线）：

```
sd search（--engine pgi|lastz，默认 pgi）→ putative SD hits（PSL）
  → sd align（pgr pl chainnet 非 --syn → MAF → PAF）
    → sd cluster（PAF → 重复家族 FASTA）
      → sd decompose（家族 → elementary SD BED）
        → sd cover（elementary + PAF → core duplicon 标记）
sd cross（跨基因组 SD 映射，--engine pgi|lastz）
sd run（全流程串联）
```

## 2. 库结构（`libs/sd/`）

- `mod.rs`：SD 公共过滤 `psl_block_len` / `psl_identity` / `passes_sd_filters`
  （T2T-CHM13 标准：block_len ≥ min_len、identity = (match+rep)/block_len，含 insert 碱基）；
- `search_pgi.rs`：原生 pgi 引擎（调用 `pgr align pgi` 自比对或两基因组比对），
  `SearchPgiOptions` + `pgi_to_hits` + `search_pgi`；
- `search_lastz.rs`：外部 lastz 引擎（`lastz --self` → LAV → PSL），
  `SearchLastzOptions` + `lastz_to_hits` + `search_lastz`；
- `cluster.rs` / `decompose.rs` / `cover.rs`：BISER cluster/decompose/cover 的 pgr 实现
  （`decompose` 用 k=10 完整 k-mer 索引 + 多拷贝 plane-sweep，`cover` 用贪心 set cover）。

两引擎共用 `libs/lastz::build_common_args`（query-depth/LAV 格式/preset + 打分矩阵临时文件），
`pgr align lastz` 与 `sd search --engine lastz` 同一来源。

## 3. 实施状态与验证

- **search**：`--engine pgi|lastz`（默认 pgi）。lastz 引擎 MG1655 81 秒 / 264 条；
  pgi 引擎 MG1655 8.3 秒 / 242 条（≥1 kb、≥90% 身份）。
- **align**：`pgr pl chainnet`（非 --syn）+ MAF→PAF 合并，MG1655 0.2 秒 / 90 条 PAF。
- **cluster**：PAF 双 mate 绑定 + 同染色体重叠区间 union-find，按连通分量输出
  `cluster_N.fa`（头 `{species}#{chrom}{strand}#{start}#{end}`）；90 条 PAF → 25 个 cluster。
- **decompose**：k=10 完整 k-mer 索引，共享 k-mer（≥5 个）经 Dsu 归组、`MAX_GAP=50`
  链式合并、`MIN_LEN=100` 过滤，输出 `species\tchrom\tbegin\tend\tset_id\tlength\tscore\tstrand`
  （基因组坐标，`set_id` 每 cluster 从 1 起，`run` 合并时跨 cluster 重编号）。
- **cover**：elementary set 覆盖 SD hit，贪心 set-cover 选最小集合，追加 `CORE`/`non-core`；
  MG1655 23 个 elementary sets，22 个标 CORE。
- **run**：search → align → cluster → decompose → cover 串联，MG1655 全程 82 秒。
- **cross**：两基因组映射（`--engine pgi|lastz`，默认 pgi）；lastz 引擎 MG1655×Sakai
  2 分 5 秒 / 303 条跨基因组同源区段。

## 4. 设计决策（自 BISER 迁移设计稿提炼）

### 4.1 为什么用外部比对路线，而非自研 k-mer + plane-sweep

BISER 原生 `search`（ordered Jaccard + plane-sweep，k=14/w=16 winnowing）面向低至 70%
同一性的古老 SD，自研 k-mer 索引 + plane-sweep 实现成本高；pgr 采用 T2T-CHM13 SD 标准
（>1 kb、>90% 同一性），外部比对器（lastz / pgi）默认参数正好匹配，无需为低同源性调参。
最终选择 `--engine pgi|lastz` 双引擎，`sd search` 收敛为一条命令。

### 4.2 为什么 PAF 替代 BISER `.align`

BISER 内部用 14 列 `.align` 中间格式；pgr 以标准 PAF（12 列 + `cg:Z:`/`cs:Z:` + 推荐 tag）
为统一中间格式，下游 cluster/decompose/cover 直接消费，避免私有格式与双向转换。

### 4.3 为什么 `pgr pl chainnet` 而非 `pgr pl ucsc`

chain/net 精修必须走**非 `--syn`**（`--syn` 会经 netFilter 只保留共线性比对，丢掉伴随重排的
SD）。`pgr pl ucsc` 依赖外部 kent-tools 且有 Linux 崩溃风险，改用原生
`pgr pl chainnet`（与 UCSC 字节级一致）。

### 4.4 cluster 的端点 coloring 语义

对 hit 的四个端点做区间 coloring，等价于 union-find 找重叠 hit 的连通分量；pgr 用
`ds::Dsu`（原 `paf/graph/dsu.rs` 迁移至 `libs/ds/`）实现，PAF 双 mate 绑定后按同染色体
重叠区间合并。

### 4.5 decompose 的 k=10 完整索引 + 多拷贝 plane-sweep

对 cluster FASTA 建 10-mer 完整索引（不做 winnowing），plane-sweep 时每个节点维护
`mappings`（各染色体当前最右边界）以跟踪同一 elementary SD 的多拷贝；`diff=50` 老化、
`merge gap=500` 合并相邻集合。pgr 实现直接在 `libs/sd/decompose.rs` 内完成
（未拆 k-mer 索引与 plane-sweep 两个独立模块）。

### 4.6 为什么不用覆盖度路线（pgr-repeat.sh）

覆盖度 + TE 差集 + 自匹配路线（早期备选设计，未采用）可复用 `pgr fa window`/`align lastz` 等现有命令，
但缺少 BISER 的 error model 与 elementary SD 分解，只适合重复屏蔽而非 SD 结构分析；
SD 检测/分解最终采用 search → chain/net → cluster/decompose 路线。

### 4.7 App-Egaz 的教训

App-Egaz/linkr 的历史流程（lastz 自比对 + link 图聚类 + blastn 扩展）效果不如 BISER：中间步骤
多、无 error model、无 elementary 分解；仅"lastz 自比对作为种子"和"图聚类思想"被吸收
（对应 `sd search --engine lastz` 与 cluster 的连通分量）。

### 4.8 忽略项：BISER 的 MAX_EXTEND 边界扩展未移植

**BISER 原版**在 `save_sd()` 输出 hit 前对边界做 `MAX_EXTEND`（5000 bp）填充，
再由后续局部比对（seed-and-extend + chaining + refinement）覆盖真实边界
（见 [[references/biser.md]] §3.3）。**pgr 移植时忽略了这一步**：

- `sd search`（pgi/lastz 引擎）：只比对 + 过滤，边界原样，无扩展；
- `sd align`（chain/net 非 `--syn`）：只链化已有块，不扩展边界；
- `sd cluster`：按 PAF 坐标直接提取序列，区间即比对块区间；
- `sd decompose`：只有序列内 `MAX_GAP=50` bp 的共享 k-mer 片段合并，
  无法发现比对块之外的同源（cluster FASTA 里没有那些碱基）。

**影响（2026-08-06 实测修正）**：10 个 E. coli 基因组对比 pgi/lastz 引擎，
pgi hit 左边界比 lastz 短、中位仅 2–6 bp、约一半短 >3 bp，右边界一致
（0 bp 中位）；大边界差（如 nissle）是 hit 结构差异（pgi 拆分 vs lastz
合并），非边界问题。旧的"1–11 bp"说法基于种子锚定理论假设，实测量级为
几 bp——刚过 `--min-len` 阈值的拷贝被过滤的风险存在但边际。
**后续可选改进（2026-08-06 定稿：不做）**：实测与 lastz 比较后确认收益边际
（pgi 左边界仅短 2–6 bp、右边界一致），且 `freq=50/k=31` 灵敏度优化已解决
更实质的漏检问题——不再移植 MAX_EXTEND（详见 §4.9）。

### 4.9 pgi 引擎灵敏度边界（2026-08-06 实测）

**方法**：10 个 E. coli 基因组（tests/genome/）各跑 `sd search` pgi 与
lastz 引擎，对 lastz 检出的每个 hit 按 identity 分箱，统计 pgi 覆盖
<50% 的漏检率。

| identity（分歧） | 漏检率 |
|---|---:|
| ≥0.95（<5%） | 0% |
| 0.93–0.95（5–7%） | 2.4–2.8% |
| 0.90–0.93（7–10%） | 11.4–18.2% |

**结论**：

- **pgi 在分歧 <5% 时稳定检出**（漏检率 ~0）；5–7% 开始零星漏检
  （~2.5%）；**≥7% 分歧（identity ≤93%）漏检率显著上升（11–18%）**——
  与设计里"SD 身份下限 ~90–93%"吻合，这是 pgi 引擎的灵敏度边界；
- **独立漏检通道**：e2348_69 有 562 个 identity 96–100% 的 hit 被漏
  （高拷贝重复被 `freq=10` 种子频率过滤），与分歧无关；
- **对后续方向的启示**：比边界扩展（§4.8，几 bp）更实质的问题是灵敏度
  （高分歧拷贝漏检）；提高方向 = 放宽 `freq`（代价是重复种子噪声）或为
  pgi 补一个 low-sensitivity 补充引擎（如 `align rest` 式）。

**✅ 灵敏度优化已落地（2026-08-06）**：`sd search` 增加 pgi 参数透传
（`--freq/--kmer/--smer/--window`），实测确定默认 `freq=50, kmer=31`：

- `freq` 10→50：解决高拷贝重复漏检（e2348_69 漏检 562→0）；
- `kmer` 40→31：解决 90–93% 分歧漏检（sakai 4→0、e24377a 2→0）；
- 10 个基因组整体漏检率 **13.1% → 0.26%**（579/4413 → 11/4215），
  剩余 11 个为每基因组 1–2 对边缘个案（identity 0.91–0.95，疑似
  低复杂度结构，k=21/window=3 均无效，记为已知限制）。

**实现**：`SearchPgiOptions` 增加 freq/kmer/smer/window（`pgr align pgi`
透传），CLI 默认 freq=50/kmer=31；集成测试
`command_sd_search_pgi_freq_retains_high_copy_repeats`（15 拷贝重复，
freq 50 > freq 10）。

**总体评估（2026-08-06，调低参数的副作用核查）**：

- **无假阳性引入**：k=31 的 hit 全部 ≥0.90 identity（各基因组 identity
  <0.90 的均为 0 个），中位 0.97–1.0——31-mer 特异性足够，未引入噪声；
- **hit 数量 +35%**（3510→4740）主要来自 e2348_69（freq 恢复高拷贝
  重复，180→1072）与 sakai/mg1655（高分歧 SD 恢复），是真漏检的恢复
  而非噪声；
- **性能持平**：k31/freq50 与 k40/freq10 相当甚至略快（sakai 0.74 s vs
  0.92 s），峰值内存无差异（~200 MB）；
- **与 lastz 一致性**：遮蔽后互相漏检 pgi 3.2% / lastz 6.0%，未遮蔽
  时 pgi 相对 lastz 漏检从 13.1% 降至 0.26%；
- **代价**：raw PSL 块数增加 ~35%（下游 `sd align`/cluster 处理量），
  属真实恢复的可接受成本。

### 4.10 重复遮蔽后的 SD 检测（2026-08-06 实测）

**动机**：10 个 E. coli 基因组此前未做重复遮蔽，SD 检测被散在重复
（IS 元件）干扰。pgr 原生遮蔽链路（`rept e-kmer` 三库 → `fa mask --hard`）
已可用，跑完整流程验证。

**方法**：对每个基因组 e-kmer 三库（tncentral/repbase/dfam）→ runlist 并集
→ `fa mask --hard`（N 遮蔽）→ `sd search`（pgi 与 lastz 引擎）。

**结果**：

- 遮蔽量 ~1.2%（E. coli 散在重复以 IS 元件为主，**tncentral 库主导**；
  dfam/repbase 仅额外贡献 ~16 bp）；
- 遮蔽后 pgi 检出 2200 / lastz 2069（未遮蔽 3510 / 4215）——IS 相关
  高 identity 重复（identity 中位 100%）被排除；
- 遮蔽后两引擎**互相漏检率 pgi 3.2% / lastz 6.0%**（未遮蔽时 pgi 相对
  lastz 漏检 13.1%）——遮蔽让 pgi 与 lastz 高度收敛，pgr 原生 pgi 引擎
  可作为 lastz 的替代（无外部依赖）；
- 剩余互相漏检（个别基因组 10–14%，如 e24377a/ec2011c_3493/se11）为
  两引擎的 hit 结构差异，需逐对核查（已知限制）。

**结论**：标准 SD 流程建议 = **先重复遮蔽（`rept e-kmer` + `fa mask`），
再 `sd search --engine pgi`（默认 freq=50/k=31）**；pgi 引擎在遮蔽后与
lastz 等价。

**soft mask vs hard mask（2026-08-06 修正：mask 只是发现阶段过滤）**：

- **mask 仅用于 `sd search` 发现候选**；后续 `sd align`/`cluster`/`decompose`
  全部读**原始 genome**（cluster 用 loc 索引提取），被 mask 的空缺自动补回，
  **任何引擎都不需要 hard mask**。BISER 同样：soft-mask 输入 → 内部转
  hard-mask 发现 → `translate` 映射回原基因组坐标（[[biser.md]] §1.2）。
- `pgr align pgi` 索引构建**硬编码 `mask=true`**（soft-mask 感知种子），
  pgi 引擎对 `fa mask` 默认小写天然有效（sakai soft 353 ≈ hard 351，
  且小写保留序列）；
- lastz 引擎大小写不敏感，soft mask 只部分过滤（sakai soft 386 vs
  hard 346）——多出的候选由下游 chain/net 精修吸收，**无需 hard mask**。

**实证（sakai，lastz 引擎，`sd align` 后对比）**：soft mask 候选 align 出
211 条 PAF，hard mask 189 条；**soft 独有的 22 条全部是真实 SD**
（含 12679 bp/id=0.981、1296 bp/0.977、1266 bp/0.973、1389 bp/0.940 等），
hard mask 因 N 区打断 lastz 比对而全部漏检。**hard mask 不仅不必要，还会
误伤真实 SD——pgr 只用 soft mask（`fa mask` 默认小写）**。
