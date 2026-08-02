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
