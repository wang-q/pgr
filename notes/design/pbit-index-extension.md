# pbit 多参考与内嵌索引段扩展（设计草案）

> 定位：pbit 在**泛基因组增量场景**下的扩展设计——多参考 + 参考索引内嵌。
> 独立成文以控制 [[pbit.md]] 的篇幅；pbit 格式基线（v1001 布局、Reference Index、
> append 语义）见主文档。
> 本文档是设计草案，格式细节待实现时定稿，落地时 bump 格式版本（1001 → 1002）。

> **实现状态（2026-08-02）**：单参考场景的内嵌索引段已落地为 pbit **v1002**——
> Reference Records 之后可选内嵌整参考 `.pgi`（`pgr pbit create --index`，
> 默认 k=40/syncmer 8/5），Footer 40 字节（`idx_offset`/`idx_size`），
> `pgr pbit to-index` 提取后可直接喂 `pgr pgi align`（E. coli 验证：
> 提取索引与独立构建字节一致、比对结果一致）；append 样本时索引段保留。
> 开发早期不做旧版本兼容（pbit.md 既有立场）。多参考（ref_id）与每参考
> 独立索引段仍未实现，见 §5 待定项。

## 9. v1003 多参考定稿（2026-08-02）

延续"不做旧版本兼容"，为多参考增量场景定稿：

1. **RefGroupEntry 增加 `ref_id: u32`**：段所属参考基因组序号；
   `ref_group_id` 全局唯一（跨参考连续编号）。
2. **Reference Index 节**：group 条目之后追加 Reference Table：
   ```
   u32 ref_count
   for each ref: (ref_name, idx_offset u64, idx_size u64,
                  group_start u32, group_count u32)
   ```
   每个参考的 `.pgi` 索引段偏移存于此处；v1002 的 footer idx 字段删除，
   **Footer 回到 24 字节**（ref_index/delta/sample 三偏移）。
3. **Sample Index 不加 ref_id 字段**：样本段已带全局 `ref_group_id`，
   反查即可得 ref_id（简洁优先，不存冗余）。
4. **路由：用户指定**：`create -r` 可重复（一个参考一个文件）；
   `--name` TSV 第 4 列为参考名/序号（默认 0）；`-i` 模式全部路由参考 0。
   自动路由（按 contig 覆盖/k-mer 相似度选参考）留作后续。
5. **`append-ref`**：追加参考 = Reference Records 尾部追加
   [新参考 2bit 段 | 新参考 .pgi]，重写 Reference Index（旧条目不动 +
   新条目 + ref table）+ delta + sample + footer（复用 open_for_append
   的"截断重写"模式）。
6. **版本 1003**。

> **实现状态（2026-08-02）**：v1003 已全部落地——`create -r` 可重复（每参考
> 一个文件）、`--name` TSV 第 4 列路由（名/序号，默认 0）、`-i` 模式路由
> 参考 0、`to-index --ref <name|N>` 按参考提取、`append-ref` 追加参考
> （`--index` 内嵌新参考 .pgi，旧样本/delta/索引全部保留）。
> E. coli 验证：mg1655+sakai 双参考归档（2501 参考段、2 样本），两参考
> 内嵌索引与独立构建字节一致，样本路由正确且重建完全精确
> （4,641,652 / 5,594,605 bp）。集成测试覆盖多参考路由、跨参考同名
> contig、append-ref 保留语义。

## 1. 背景与架构决策

泛基因组增量工作流：少量基因组起步，逐个添加新基因组；参考（锚）被反复比对，
样本逐个加入。

- **最终产物是图**（GFA / pangenome graph）；pbit 与索引都是**中间产物**，
  可重建、可丢弃。
- 因此格式演化没有长期兼容负担（软件随时更新，版本一起升级即可）；"内嵌单文件"的
  原子管理与分发简洁性 > 独立文件的解耦性。
- **边界**：内嵌的索引量 = **参考数**（几十个锚），**不随样本数增长**（4 万样本
  各自的索引临时建、比对完即弃，否则 pbit 会塞进 4 万索引段）。
- 参考的 2bit 序列记录本就在 pbit 参考层（随机访问），索引段只额外携带"查询结构"
  （参考级 syncmer 稀疏 k-mer 表）。

## 2. 布局（扩展后）

```
┌───────────────────────────────────────┐
│ Header (固定 36 字节 + 版本 bump)      │
├───────────────────────────────────────┤
│ 参考1: [2bit 记录 | 参考1 索引段]      │  ← 成对追加，append 单元
│ 参考2: [2bit 记录 | 参考2 索引段]      │
│ ...                                   │
├───────────────────────────────────────┤
│ Delta Data (样本 delta，append 在末尾) │
├───────────────────────────────────────┤
│ Reference Index + Sample Index         │
├───────────────────────────────────────┤
│ Footer                                 │
└───────────────────────────────────────┘
```

**"参考 + 索引段"是 append 单元**：每添加一个参考，把 [2bit 记录 | 索引段] 作为一个
整体追加到文件尾部，旧段一律不动。局部性好——加载单个参考时 2bit 与索引相邻，
seek 距离短。

## 3. 格式扩展点

1. **Reference Index 每条目扩展**（现为 `name + segment_offset`）：
   `seq_offset + seq_size`（2bit 记录）+ `idx_offset + idx_size`（索引段，可为 0/空）+
   `has_idx` flag。追加参考 = 表里加一条，旧条目不动。
2. **索引段内容 = 独立索引格式（.pgi）的字节**：pgr 定义自己的单文件索引格式
   `.pgi`（pgr genome index），独立文件形态供比对命令消费，内嵌 pbit 时原样存入
   索引段（与 2bit 的"独立格式 + pbit 参考层复用"先例一致，共享读写函数）。
   功能需求见 §3.1（目标：支撑类似 FastGA 的种子发现 / 链扫描 / 局部扩展）。
3. **Header / Footer**：版本 bump；追加参考后重写 Reference Index 与 Footer
   （pbit 已有"读 Footer → patch → 重写"的 append 模式，加样本时 patch sample_count；
   加参考同理 patch reference 条目与索引表偏移）。
4. **Sample Index 扩展**：每条样本加 `ref_id` 字段，记录该样本对照哪个参考
   （多参考下样本归属的显式表达）。
5. **索引段可选**：仅做解压访问的参考可以不建索引段（`has_idx = 0`），旧工具读不到
   索引段也照常解压数据。

### 3.1 索引格式（.pgi）功能需求

索引构建的目的：为后续实现**类似 FastGA 的功能**（输入基因组 → 种子 → 链 → 局部
扩展 → PSL，接入 `pgr pl chainnet`）。因此格式必须覆盖 FastGA 管线的全部需求，
内容对齐 [[fastga.md]] §11：

| FastGA 阶段 | 索引要提供 | .pgi 必须包含 |
|-------------|-----------|---------------|
| 种子发现（两索引归并）| 字典序排序的 k-mer 流 + 位置 | 排序条目：k-mer 编码（2-bit）+ 位置（contig, pos）|
| 频率过滤（-f 10）| 每个 k-mer 出现次数 | 条目携带频率，或排序流上即时统计 |
| 种子扩展（adaptamer）| 相邻相同 k-mer 的 lcp | 条目携带 lcp，或排序流即时计算 |
| 链扫描 | diag/anti 坐标、对角线桶 | 由位置即时计算（diag=i−j, anti=i+j）；预编码是实现选择 |
| 双链 | 正反链都索引 | 每序列正反链条目（或 canonical 2-bit 合并）|
| 两流归并 | 流式迭代 + 定位（GoTo）| 排序数组 + 流式读取接口（可 mmap）|

**设计原则**：
- 内容 = GIX 的功能集（syncmer 稀疏 + 2-bit + 位置 + 频率 + lcp + 双链 + 排序流）；
- 物理组织 = pgr 自研单文件（排序数组 `Vec<(u64, u32)>` + 附列），可 mmap，
  不用 GIX 的代理 + `.ktab` 分片 ensemble（见 [[fastga.md]] §10）；
- 采样参数（k 大小、syncmer 密度）写入 Header，供读取端对齐；
- 独立命令生成（`pgr pgi build`）与内嵌 pbit 共享同一读写实现。

## 4. 追加语义

- `append` 参考：尾部写入 [2bit | 索引段] → 更新 Reference Index → 重写 Footer。
- `append` 样本：不变（尾部 delta + patch sample_count），与参考索引段互不影响。
- 样本索引：比对时临时构建（FastGA 式，见 [[fastga.md]] §12），用完即弃，不落盘。
- 比对流程：新样本 → 临时索引 → 与目标参考的内嵌索引段比对（或两索引归并）→
  结果（PAF）进 `.paf.idx`、压缩结果进 pbit delta → 样本索引删除。

## 5. 待定项

- 索引段二进制格式定稿（k 大小、syncmer 参数、条目编码、mmap 布局）；
- 版本 bump（1001 → 1002）与旧文件兼容策略（读取时按版本分支，旧版忽略索引段）；
- 多参考下"样本 vs 参考"的比对路由（选哪个参考由用户指定还是自动？）。

## 6. .pgi 的距离消费者（dist 层级）

`.pgi` 除了服务 FastGA 式比对，还派生两类**距离计算**，形成距离精度/成本层级：

| 命令/模式 | 方式 | 复杂度 | 场景 |
|-----------|------|--------|------|
| `pgr dist seq`（现有）| 序列 sketch（minimizer/syncmer 采样）| O(序列长) | 大规模粗筛（4 万去冗余的 Mash 近邻）|
| `pgr dist hv`（现有）| 序列 → hypervector → 比较 | O(序列长) + O(dim) | 更快粗筛（向量比较）|
| `pgr dist pgi`（新增）| 两 .pgi 排序 k-mer 流**归并** | O(\|K1\|+\|K2\|) | .pgi 已存在时的**精确**距离（确定性、零采样方差）|
| `pgr pgi to-hv` + HV 比较（新增）| .pgi → hypervector（一次投影）→ 比较 | O(\|K\|) 生成 + O(dim) 比较 | .pgi 已存在时的超大规模两两距离（4 万级 KNN）|

### 6.1 `pgr dist pgi`（索引归并距离）

- 输入：两个 `.pgi`（或目录 + 通配符两两）；
- 算法：两排序 k-mer 流归并 → `total1/total2/inter/union` → **精确** Jaccard /
  containment / Mash 距离（输出与 `dist seq` 对齐）；
- **校验**：两索引采样参数（k、w、密度、canonical）必须一致（Header 比对），
  不一致报错——否则共享 k-mer 比没有距离语义；
- 无额外采样参数（索引 Header 自带）。

### 6.2 `pgr pgi to-hv`（索引 → hypervector）+ HV 比较

- **动机**：索引归并是 O(\|K\|)（排序流线性），对 4 万级两两仍贵；把 .pgi 的 k-mer
  哈希集**一次投影**到固定维度向量后，HV 比较是 O(dim)，便宜一个量级。
  采样已在索引构建时完成，HV 生成无需再读序列。
- `pgr pgi to-hv <idx.pgi> -o out.hv`：遍历索引 k-mer 集，折叠到固定向量
  （维度、哈希折叠方式写入 HV Header）。
- HV 比较：复用 `dist hv` 的向量距离核心（新入口接收 `.hv` 文件，或新增子命令），
  输出同样式（jaccard/containment/mash 从向量内积/距离推导）。
- **校验**：HV 的采样参数与维度必须一致才可比（Header 记录，比较时比对）。
- **与 `dist hv` 的关系**：`dist hv` 从序列采样生成 HV；`pgi to-hv` 从索引生成。
  两者 HV **参数不同则不可互比**（.pgi 是 FastGA 系参数，如 k=40、(12,8) syncmer；
  `dist hv` 默认 k=8/w=55），只在各自参数体系内比较。

### 6.3 层级定位

- 粗筛（海量 pair）：`dist seq` / `dist hv`（sketch，快）；
- 精筛（候选 pair）：FastGA 比对（比对器）；
- 精确距离（.pgi 已存在）：`dist pgi`（归并）；
- 超大规模粗筛（.pgi 已存在）：`pgi to-hv` + HV 比较（O(dim)）。

> **实现状态（2026-08-02）**：`dist pgi` 已实现并验证（45 对 cohort）；
> `pgi to-hv` + `dist hv <a.hv> <b.hv>` 已实现（`.hv` 直接比较）。
> **关键实证**（详见 [[../benchmarks/dist-cohort-validation.md]]）：
> `dist pgi` 是"采样集合的精确距离"（确定性，但与身份率 Spearman 0.54）；
> `dist hv` 已修复为**稀疏投影 + 余弦**（`.hv` v2：每 k-mer 更新 `--sparse`
> 个随机维度、头存 sparse/n_kmer、比较用余弦估计共享数）：与 `dist pgi`
> 的 mash Spearman 0.97、共享数平均误差 2.4%、45 对比较快 50×——
> "粗筛近似层"定位成立；`dist seq`（k=8）仍是与身份率最贴近的草图层
> （0.82 vs 0.51）。

## 7. CLI 设计（pgr pgi 命令族）

```
pgr pgi — Manages pgr genome index (.pgi) files

Commands:
  build   构建 .pgi（输入 .2bit 优先，FASTA fallback）
  stat    索引统计（k/syncmer 参数、条目数、文件大小）
  to-hv   投影为 hypervector（供 dist 快速比较）
  show    查看索引条目（调试用，可选）
```

> **明确不做**：`import-gix`（从 FastGA GIX 导入）——格式逆向成本高、收益低，
> 互通通过"调用 FastGA 时保留其 GIX"实现，pgr 不读 GIX 字节。连 TODO 也不留。

### 7.1 `pgr pgi build`

```
pgr pgi build <infile> -o out.pgi
  <infile>       FASTA（.fa/.fa.gz）或 .2bit（推荐，构建最快）
  -k, --kmer <40>        k-mer 大小
      --smer <12>        syncmer 长度
      --window <8>       窗口 s-mer 数（密度 2/(w+1)）
      --no-rev           只索引正链（默认正反链）
  -t, --threads <4>      线程
```

- 默认参数与 FastGA GIX 对齐（k=40、(12,8) syncmer），服务"类似 FastGA 功能"
  与后续互通；
- `to-hv` / `dist pgi` 的距离语义跟随索引 Header 参数，不额外指定。

### 7.2 配套入口

```
pgr dist pgi <idx1> <idx2> [--list files.txt]   # 索引归并距离（见 §6.1）
pgr pgi to-hv <in.pgi> -o out.hv [--dim 1024]   # 索引 → HV（见 §6.2）
pgr dist hv <out1.hv> <out2.hv>                 # HV 比较（复用现有向量距离核心）
pgr pbit create/append ... --index              # 内嵌索引段（触发方式待定，见 §8）
```

### 7.3 待定

- pbit 内嵌索引段的 CLI 触发方式（`--index` flag？自动？细粒度控制？）——暂不决定；
- `show` 是否纳入首版（调试期可能不需要）。

## 8. 相关文档

- pbit 格式基线：[pbit.md](pbit.md)（文件格式规范 v1001、决策记录）
- 序列索引选型：[fastga.md](../references/fastga.md) §10/§12（GIX 评估、索引选型）
- 泛基因组场景：[ecoli-cohort.md](../ecoli-cohort.md)、[paf-pangenome.md](../paf-pangenome.md)
