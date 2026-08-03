# pbit 设计笔记（含多参考扩展）

> **✅ 已决策（2026-08-03）**：作者按推荐确认全部 6 项开放项（详见
> [[pbit-decisions.md]] 文末拍板清单）：路由保持手动 + 多参考未指定时警告
> （已实现）、Sample Index 不加 ref_id、append-ref 不重路由、维持单参考
> 压缩模型、版本策略维持"不兼容 + bump"、HV sketch 内嵌暂缓。暂停解除，
> 可继续开发。
>
> **已决策（2026-08-02）**：~~内嵌索引~~ → **决策 A：索引不进 pbit**。
> 实测索引 ~92 MB vs 压缩归档 ~1.1 MB（79×），内嵌会让"压缩格式"失去意义；
> `.pgi` 保持独立临时工作对象（用时现建 0.3s 或旁路缓存），与 FastGA 的
> GIX"独立文件、用完即删"定位一致。HV sketch 内嵌（决策 B）暂缓，其算法
> 设计仍需思考。格式已按决策 A 落地为 **v1004**（无索引字段）。
>
> **已决策的开放项（2026-08-03 确认）**：
> 1. **样本 vs 参考的路由**：保持手动（TSV 第 4 列参考名/序号，默认参考 0）；
>    多参考且未指定时 `log::warn` 警告（已实现）。自动路由留待多样性 cohort
>    数据证明收益后再做。
> 2. **Sample Index 不加 ref_id**：样本段已带全局 `ref_group_id`，经 Reference
>    Table 反查即可。
> 3. **append-ref 不重路由**：只加参考、不改已有样本路由；"换锚"需求出现时
>    做显式 `re-anchor` 子命令。
> 4. **多参考压缩模型**：维持单参考/样本；格式已支持每段 ref_group_id，
>    未来做按 contig 路由无需格式变更。
> 5. **版本策略**：不做旧版本兼容，格式改动直接 bump 版本；长期归档需求
>    出现时用 `convert` 逃生舱。
> 6. **决策 B（HV sketch 内嵌）**：暂缓，触发条件 = 出现"无源 FASTA、仅归档、
>    需距离粗筛"的真实工作流（设计稿见 [[pbit-decisions.md]]）。

## 当前状态（v1004，2026-08-02）

pbit 为原生"2bit 参考 + delta 样本"群体基因组压缩格式（区别于 C++ AGC 的
`.agc`）。已实现：

- `pgr pbit create`（单/多参考，`-r` 可重复，TSV 第 4 列路由样本到参考）、
  `append`（追加样本）、`append-ref`（追加参考）、`stat` / `to-fa` /
  `some` / `range`（读取）；
- 多参考（每参考一个 2bit 段组 + Reference Table），样本路由到指定参考；
  E. coli 双参考归档验证（样本路由正确、重建精确）；
- **不内嵌索引**（决策 A）：`.pgi` 为独立临时工作对象，比对/距离时现建；
- 版本 1004，仅当前版本可读写（不做旧版本兼容）。

## 快速参考

| 子命令 | 分组 | 用途 | 关键参数 |
|--------|------|------|----------|
| `create` | build | 创建归档（单/多参考） | `-r ref.fa`（可重复）, `-i sample.fa` / `--name tsv`, `-o out.pbit` |
| `append` | build | 追加样本 | `in.pbit`, `-i sample.fa`, `-o out.pbit`（可选） |
| `append-ref` | build | 追加参考 | `in.pbit`, `-r ref.fa`, `-o out.pbit`（可选） |
| `to-fa` | transform | 提取所有样本为 FASTA | `in.pbit`, `-o out_dir/` |
| `some` | subset | 按样本名列表提取 | `in.pbit`, `sample_list.txt`, `-o out.fa` |
| `range` | subset | 按 contig/区间提取 | `in.pbit`, `chr1:1-1000`, `-o out.fa` |
| `stat` | info | 统计/列表 | `in.pbit`, `--samples` / `--refs` / `--contigs` |

样本名默认取输入 FASTA basename（`--name` TSV 可覆盖）。TSV 列：
`sample_name<TAB>fasta_path[<TAB>paf_path][<TAB>ref_name]`。

## 文件格式规范（v1004）

所有整数固定大小小端序（u32/u64），字符串为 u32 长度前缀 + UTF-8，不用
varint/null 终止。参考层直接复用标准 2bit 记录（`read_2bit_record` /
`write_2bit_record`，与 twobit.rs 共享代码，保留 N/mask blocks）。

### 文件结构总览

```
┌─────────────────────────────────────┐  ← offset 0
│ Header (固定 36 字节)               │
├─────────────────────────────────────┤
│ Reference Records                   │  ← 每段一个标准 2bit 记录（跨参考连续）
├─────────────────────────────────────┤  ← footer.ref_index_offset
│ Reference Index                     │  ← 参考段条目 + Reference Table
├─────────────────────────────────────┤  ← footer.delta_data_offset
│ Delta Data                          │  ← 每参考组的 delta 列表（flate2 压缩）
├─────────────────────────────────────┤  ← footer.sample_index_offset
│ Sample Index (collection)           │  ← flate2(序列化 samples/contigs/segments)
├─────────────────────────────────────┤
│ Footer (固定 24 字节)               │
└─────────────────────────────────────┘  ← EOF
```

### Header（36 字节）

```
0  4  magic              0x54494250 ('PBIT')
4  4  version            major*1000 + minor（当前 1004）
8  4  segment_size       分段大小（bp，如 4096）
12 4  kmer_len           LZ-diff 哈希 k-mer 长度（如 15）
16 4  min_match_len      LZ-diff 最小匹配长度（如 18）
20 4  ref_group_count    参考段总数（每段一个 group，跨参考）
24 4  sample_count       样本数
28 8  ref_records_offset Reference Records 起始偏移（通常 36）
```

### Reference Records（标准 2bit 记录，连续存储）

每段一个标准 2bit 记录：`dna_size + n_blocks + n_starts/n_sizes +
mask_block_count + mask_starts/mask_sizes + reserved + packed_dna`。
参考段不二次压缩（2bit 已压缩 4 倍，delta 层负责进一步压缩）。

### Reference Index（参考段条目 + Reference Table）

```
u32 ref_group_count
for each ref_group:
  str contig_name
  u32 ref_id            该段所属参考基因组序号
  u64 segment_offset    该 2bit 记录的文件偏移
u32 ref_count           Reference Table
for each ref:
  str ref_name          参考名（FASTA basename）
  u32 group_start       该参考首段的 group id
  u32 group_count       该参考段数
```

> `ref_group_id` 全局唯一（跨参考连续编号），样本段通过它反查所属参考，
> 故 Sample Index 不重复存 ref_id。Reference Table 只负责参考命名与段范围
> （索引不内嵌，决策 A）。

### Delta Data

```
u32 ref_group_count
for each ref_group:
  u32 delta_count
  for each delta:
    u8  is_rev_comp
    u32 raw_length      样本段长度（query 轴）
    u32 packed_size
    u8  encoding        0 = LZ-diff, 1 = CIGAR
    bytes packed_data   flate2(编码字节)
```

delta 为变长记录，不在文件存 per-delta 偏移表；`Decompressor::new` 顺序
扫描头部（10 字节/条）构建内存 `delta_offsets`，不解压数据。

### Sample Index（collection，flate2 整块压缩）

```
u32 sample_count
for each sample:
  str name
  u32 contig_count
  for each contig:
    str contig_name
    u32 segment_count
    for each segment（固定 16 字节）:
      u32 ref_group_id
      u32 delta_id
      u32 ref_start     参考段内相对起始（CIGAR 用；LZ-diff 填 0）
      u32 ref_end       参考段内相对结束（CIGAR 用；LZ-diff 填 0）
str cmd_line
```

### Footer（24 字节）

```
u64 ref_index_offset / u64 delta_data_offset / u64 sample_index_offset
```

> **设计要点**：无文件尾 magic（Header 已含）；全部固定大小字段（解析简单、
> 破坏可校验）；`ref_group_count`/`sample_count` 在 Header/Index/Delta/
> Sample 中重复——读取各 section 时就地校验一致性，避免损坏文件越界。

## 编码算法

### LZ-diff（默认路径）

样本段按 `segment_size` 分段，与参考段（整条 2bit 记录）做 LZ-diff
（k-mer 哈希索引找最长匹配，`kmer_len`/`min_match_len` 控制），差异
编码后 flate2 压缩。无 PAF 时所有样本段走此路径。

### PAF 驱动的 CIGAR delta（`--paf` 路径）

用 PAF（含 `cg:Z:` CIGAR，建议 `--eqx`）驱动压缩：样本段被 PAF alignment
完整覆盖时，把 CIGAR 按段切片存储（`ref_start/ref_end` 定位参考区间），
`packed_data = flate2(u32 op_count + [CigarOp; op_count] + u32 base_count +
[u8; base_count])`，CigarOp 为 `(op << 29) | len` 的 u32。

关键决策（详见旧版决策记录，已实现）：
- **段级回退**：最佳 alignment 未完整覆盖段、段跨多条 alignment 衔接、
  CIGAR target 投影跨参考段边界 → 整段回退 LZ-diff（不拆段、不合并）；
- **`M` 处理**：`maf to-paf` 源头已产出 `=/X`；对 minimap2 未用 `--eqx`
  的 `M` CIGAR 保留 `M → =/X` 拆分 fallback（压缩端读样本序列拆分）；
- **query-side 索引**：pbit 自建 `BasicCOITree<PafMetadata>`（不复用
  `PafIndex.reverse_trees`——其缺 `-` 链、元数据被交换、无公开查询接口）；
- **无/坏 CIGAR**：记录级错误跳过 + 回退 LZ-diff（log 警告）；整个 PAF
  不可用则报错终止（避免"以为生效实为全回退"）。

## 多参考（v1003/v1004 扩展）

泛基因组增量场景：少量基因组起步、逐个添加；参考（锚）被反复比对，样本
逐个加入。最终产物是图（GFA），pbit 与索引都是可重建的中间产物，故格式
演化不做长期兼容负担。

**索引定位（决策 A，2026-08-02）**：`.pgi` **不进 pbit**。实测单参考
归档 1.1 MiB、内嵌索引后 92.8 MiB（索引/压缩数据 = 79×）——内嵌会让
pbit 失去"压缩格式"的意义。参考索引在需要比对/距离时现建（~0.3s）或由
工作流在归档旁缓存为 `ref.pgi` 兄弟文件；与 FastGA"GIX 独立文件、用完即删"
定位一致。HV sketch 内嵌（决策 B）暂缓，其算法设计待后续思考。

**追加语义**：
- `append` 样本：尾部追加 delta → patch sample_count（v1001 起即有）；
- `append-ref` 参考：Reference Records 尾部追加新参考 2bit 段，重写
  Reference Index（旧条目不动 + 新条目 + Reference Table）+ delta + sample
  + Footer（截断重写模式，旧样本保留）；样本路由由压缩时的 `ref_group_id`
  决定，追加参考不影响已有样本。

**样本路由**：v1 为用户指定（TSV 第 4 列参考名/序号，默认参考 0）；自动
路由（k-mer 相似度/contig 覆盖）为开放项（见顶部）。

### .pgi 距离消费者层级（已实现并验证，索引独立于 pbit）

| 命令/模式 | 方式 | 复杂度 | 与身份率的 Spearman | 定位 |
|---|---|---:|---:|---|
| `dist seq`（k=8 syncmer） | 草图 | O(序列长) | 0.82（最优） | 大规模粗筛 |
| `dist pgi` | 两排序流归并 | O(\|K1\|+\|K2\|) | 0.54 | 已建索引时的确定性精确距离 |
| `dist hv`（.hv v2 稀疏） | 稀疏投影 + 余弦 | O(dim) | 0.51（与 `dist pgi` 的 ρ=0.97） | `dist pgi` 的 50× 快近似 |

> 注意：k=40 syncmer 集合受采样位置漂移影响，与真实身份率的排序相关只有
> ~0.5；`dist hv` 已修复为稀疏投影（每 k-mer 更新 `--sparse` 个维度），
> 详细数据见 `notes/benchmarks/dist-cohort-validation.md`。

## 参考资料

- 多参考扩展设计过程与 .pgi 消费者规划：旧 `pbit-index-extension.md` 已并入
  本文（该文件现为跳转 stub）；
- 距离消费者验证数据：`notes/benchmarks/dist-cohort-validation.md`；
- pgi 比对管线（`.pgi` 的消费者）：`notes/design/pgi-align.md`；
- AGC 算法参考：`notes/references/agc-cpp.md`。
