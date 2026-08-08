# pbit 设计笔记（含多参考扩展）

## 本文档结构

1. 设计基础（2026-08-09）：从目标推导的 7 条原则，所有决策的前提；其后
   为历史决策记录
2. 核心设计：编码模型——参考分段 + 样本三级编码 + 段级混合 + 主链/碎链
   + 无损保证（快速总览）
3. 当前状态：命令族、能力演进 v1004–v1010、行为增强
4. 快速参考：子命令表
5. 编码模型详解：LZ-diff / CIGAR / Identity / Raw 各路径
6. 大链与碎链：定义、旧现状、判定细节定稿、实现与验证
7. 文件格式规范：字节级布局（Header/Reference/Delta/Sample/PAF
   Recovery/Footer）
8. PAF 驱动编码的演进：#14 诊断 → v1010 的三条路线与约束解除
9. 专题：遮蔽（v1005）、统一序列访问 API（不做）、多参考
10. 附录：早期开放项决策稿（2026-08-03，选项对比与拍板记录）

---

## 设计基础（2026-08-09 明确）

> 本文档所有设计决策的**前提**。原则：决策必须能回溯到下面的基础；没有
> 依据的决策不进文档，新增决策先说明它由哪条基础推出。

1. **目标：样本尽量复用参考序列**——参考是压缩锚，样本只存"相对参考的
   差异"，不原样存储样本内容。无法归到参考的内容（Raw 2bit）只是兜底，
   不是目标路径。
   - **算法来源（AGC 策略）**：k-mer（minimizer）引导的参考段选择，把
     样本序列压缩到参考上——样本段先在参考中定位最佳参考段，再用
     LZ-diff 编码差异；delta 为空即直接指向参考（零载荷，对应 v1010
     Identity）。详见 [agc-cpp.md](../references/agc-cpp.md)。
2. **坐标系约定：参考固定为 target，样本为 query**——所有生成的 PAF 都
   以参考为基础、样本向上面对（链路固定 `align pgi ref sample`，方向
   不可反转；正反方向不等价，实测 809 vs 1237 条链）。
   - 理由 ①：用户理解——样本的差异以参考为参照物；
   - 理由 ②：PAF 查询 / mapping 必须共享同一参考坐标系（`pgr paf` 把
     query 区间投影到 target，跨样本同源查询的前提）。
   - 推论：所有样本的比对都在参考坐标下组织；查询、比较、pbit 编码
     都以参考坐标为统一基准。
3. **无损保证**：样本逐碱基还原（`to-fa`）；输入 PAF 可还原
   （`to-paf`，v1009 起）。
4. **PAF 存储基础**：PAF 只存坐标 + 匹配统计 + 标签（不存序列）；
   **主链重建行**的 `qlen/tlen/qname/tname` 不重复存——由归档 contig 表
   恢复（v1009；碎链行是例外，整行原样存，见第 5 条）。
   实测恢复区字段占比（00_3076 vs 00_3230，809 行真实 PAF）：
   **`cg:Z` + `cs:Z` ≈ 72%**（cs 46% + cg 37%），坐标（q+t）≈ 9%，
   `gi`/`bi` ≈ 3.4%，`ms` ≈ 1.3%，matches/block/mapq ≈ 1.5%。
   推论：PAF 存储/压缩优化的焦点在 **cg/cs 文本**；碎链行恢复区实测
   ~100 KB/样本（flate2 后），使 delta/gzip +9 pp（0.356 → 0.448），
   见 [benchmarks/bench-scale-and-pbit.md](../benchmarks/bench-scale-and-pbit.md) #14k。
5. **计算优先于存储（可计算的不存储）**：归档只存"无法从其他信息推导"
   的字段；能从归档内容重算的字段一律不存、导出时重算。大链行按此
   原则：`qname/qlen/qstart/qend/strand/tname/tlen/tstart/tend/matches/
   block/mapq/gi/bi/cg/cs` 全部由归档重建/重算，只存 `paf_id`（归组键）
   + `ms`（比对器打分，不可重算）。**碎链行是例外**：整行原样存
   （"存进去什么，出来什么"的最大保真 + 实现简化，2026-08-09 用户
   拍板），含可重算的 name/len——该例外与压缩算法（flate2 现状）一并
   维持，暂不优化。
6. **PAF 路由/载体定位**：pbit 归档是"参考 + 样本序列 + 比对信息
   （PAF）"的统一载体——主链整链编码+重建、碎链原样存行，完整比对
   信息都保存在归档里，承担 PAF 路由作用；外界程序通过 `to-paf`（未来
   可加直接查询）访问。pbit 的压缩主链判定**不替 PAF 查询做路由**：
   查询层的 primary/secondary 由查询语义自己决定（参考一个位点可合法
   对应多条链，多映射信息由碎链行完整保留）。
7. **匹配区间可访问性**：比对的**坐标与链结构**必须以可访问形式保留
   （外界程序可直接访问，用于 PAF mapping/路由）——碎链行原样可逐字
   还原、主链可重算还原。**CIGAR 是例外**：它占 PAF 字段大头（cg:Z
   ≈ 37%），允许压缩存储（如未来位打包候选），只需保证导出时重建为
   标准 `cg:Z` 即可。

**由基础推出的结构**（详见下文）：

- 参考按 `segment_size` 分段 → 每段一条 2bit 记录，是 LZ-diff 的索引基础；
- 样本段三级编码：CIGAR/Identity（能对到参考）→ LZ-diff（内容相似）→
  Raw 2bit（兜底）；
- 段级混合编码（v1007+）：一段内 PAF 覆盖部分 CIGAR、未覆盖部分
  LZ/Raw 补齐；
- 大链 = 能匹配上 Reference 的链（参与编码）；碎链 = 被其他大链覆盖的链
  （即使相似度更高，原样存行以便 `to-paf` 还原）。

---

以下为历史决策记录（前提见上"设计基础"）：

> **✅ 已决策（2026-08-03）**：作者按推荐确认全部 6 项开放项（详见
> §附录 文末拍板清单）：路由保持手动 + 多参考未指定时警告
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
>    需距离粗筛"的真实工作流（设计稿见 §附录）。
>    **（2026-08-09 变更为明确不做**——HV 评测未达预期，后续换其他形式，
>    见 `todo.md` §5）。

## 核心设计：编码模型（2026-08-09 统一理解）

pbit 是"2bit 参考 + delta 样本"的群体基因组压缩格式。编码模型如下：

**参考层**：参考基因组按 `segment_size`（默认 4096 bp）分段，每段一条
标准 2bit 记录（Reference Record），同时是 LZ-diff 的索引基础。

**样本层（三级编码，逐段决策）**：

1. **CIGAR / Identity**：样本段能对到参考（PAF 记录覆盖）→ 用 CIGAR
   表示差异（纯 `=` 段用 Identity 零载荷指向参考区间）；
2. **LZ-diff**：对不上但内容与某参考段相似 → 差异编码（k-mer 内容匹配）；
3. **Raw（2bit）**：都匹配不上 → 原文 2bit 存储，保证严格无损。

**PAF 驱动为唯一路径（2026-08-09 拍板）**：`pbit create`/`append` 强制
要求 PAF（`--paf` 或 `--name` TSV 第 3 列；CLI 校验已实现），无 PAF 的
独立压缩路径退役；PAF 记录必须带 `cg:Z`（推荐链路
`chainnet --t-name '' --q-name '' → maf to-paf` 天然满足）。LZ-diff /
Raw 仍是段级兜底（PAF 未覆盖/无匹配的段），不是独立入口。无 `cg:Z` 的
记录**跳过编码、原样存行还原**（决策 7 维持，2026-08-09 确认——不参与
CIGAR 编码，但行仍在碎链恢复区，`to-paf` 可还原）。

**从 MAF 建立 pbit（管道，2026-08-09 验证）**：chainnet 原生输出 MAF，
可经 `pgr maf to-paf in.maf | pgr pbit create -r ref.fa -i sample.fa
-p stdin -o out.pbit` 一行建立（`-p stdin` 从管道读 PAF，`maf to-paf`
默认输出 stdout）。pbit **不直接消费 MAF**——保持 PAF 为唯一输入语义
（强制 PAF 架构一致），MAF → PAF 转换由 `maf to-paf` 承担。

**段级混合（v1007+）**：同一段内 PAF 覆盖的部分用 CIGAR，未覆盖部分用
LZ-diff / Raw 补齐——不要求整段对齐、不整段回退。

**大链 / 碎链（2026-08-09 拍板）**：

- 主链/碎链都按 **sample（query）坐标** 定义——链在 sample 上的 query
  区间是判定主体；reference 是编码锚与坐标系，不是主/碎判定的主体。
- 大链 = 能匹配上 Reference 的比对链（**整链编码 + 重建**），判定不看长度；
- 碎链 = 被其他大链覆盖的链（即使相似度更高也归碎链，原样存行以便
  `to-paf` 还原；**碎链仍可编码其覆盖的段**（分类只决定"重建 vs
  存行"，不决定"能否编码段"，2026-08-09）。**PAF 逐条可还原是硬约束**
  （先保证还原，再优化压缩）。
  详见"大链与碎链：定义与判定"。

**无损保证**：`to-fa` 逐碱基还原样本；`to-paf`（v1009 起）还原输入 PAF。

## 当前状态（v1010，2026-08-09）

pbit 为原生"2bit 参考 + delta 样本"群体基因组压缩格式（区别于 C++ AGC 的
`.agc`）。已实现：

- 命令族：`create`（单/多参考，`-r` 可重复，TSV 第 4 列路由样本到参考）、
  `append`、`append-ref`、`stat` / `to-fa` / `some` / `range`（读取）、
  `to-paf`（v1009 起导出内嵌比对）；
- 多参考（每参考一个 2bit 段组 + Reference Table），样本路由到指定参考；
  E. coli 双参考归档验证（样本路由正确、重建精确）；
- **不内嵌索引**（决策 A）：`.pgi` 为独立临时工作对象，比对/距离时现建；
- 版本 1010，仅当前版本可读写（不做旧版本兼容）。

能力演进（v1004 → v1010）：
- v1005：soft mask 随样本存储（`ContigSegs.mask_blocks`，继承 2bit
  mask_blocks 语义），存进存出一致；
- v1006：严格无损——无参考匹配段 Raw 存储；LZ 兜底内容匹配化
  （canonical k-mer 倒排，无 PAF/同名也可 ~100% 无损）；
- v1007：CIGAR 任意参考区间 + 段级混合编码（PAF 驱动生效，
  54%→39%，见下文"PAF 驱动编码的演进"）；
- v1008：Raw 段改标准 2bit 记录（语义与参考层统一）；
- v1009：`to-paf` 无损还原输入 PAF（大链重建 + 碎链原样存行）；
- v1010：Identity 零载荷指向参考区间（纯 `=` 段，见下文对应章节）。

行为增强（2026-08-09，格式不变，仍为 v1010）：
- `create`/`append` 强制 PAF（无 PAF 独立路径退役，空 PAF 可禁用 CIGAR）；
- 主链/碎链按覆盖关系链级判定（`BIG_CHAIN_MIN_LEN` 退役）、碎链也可编码
  其覆盖的段、无 `cg:Z` 记录存行还原（见"大链与碎链"章节）。

## 快速参考

| 子命令 | 分组 | 用途 | 关键参数 |
|--------|------|------|----------|
| `create` | build | 创建归档（单/多参考） | `-r ref.fa`（可重复）, `-i sample.fa` / `--name tsv`, **`-p sample.paf`（必填）**, `-o out.pbit` |
| `append` | build | 追加样本 | `in.pbit`, `-i sample.fa`, **`-p sample.paf`（必填）**, `-o out.pbit`（可选） |
| `append-ref` | build | 追加参考 | `in.pbit`, `-r ref.fa`, `-o out.pbit`（可选） |
| `to-fa` | transform | 提取所有样本为 FASTA | `in.pbit`, `-o out_dir/` |
| `to-paf` | transform | 导出内嵌比对为 PAF（v1009） | `in.pbit`, `-s sample`（可选）, `-o out.paf` |
| `some` | subset | 按样本名列表提取 | `in.pbit`, `sample_list.txt`, `-o out.fa` |
| `range` | subset | 按 contig/区间提取 | `in.pbit`, `chr1:1-1000`, `-o out.fa` |
| `stat` | info | 统计/列表 | `in.pbit`, `--samples` / `--refs` / `--contigs` |

样本名默认取输入 FASTA basename（`--name` TSV 可覆盖）。TSV 列：
`sample_name<TAB>fasta_path[<TAB>paf_path][<TAB>ref_name]`。

## 编码模型详解

### LZ-diff（兜底路径）

样本段按 `segment_size` 分段，与参考段（整条 2bit 记录）做 LZ-diff
（k-mer 哈希索引找最长匹配，`kmer_len`/`min_match_len` 控制），差异
编码后 flate2 压缩。PAF 为强制输入；**空 PAF**（无记录）时所有样本段
走此路径（等价于无 CIGAR 输入）。

### PAF 驱动的 CIGAR delta（`--paf` 路径）

用 PAF（含 `cg:Z:` CIGAR，建议 `--eqx`）驱动压缩：样本段被 PAF 记录
覆盖的部分按段切片存 CIGAR（`ref_start/ref_end` 定位参考区间，v1007 起
为参考文件全局坐标、可跨参考段），同一段内未覆盖部分由 LZ-diff / Raw
补齐（**段级混合编码**，不要求整段对齐）。`packed_data =
flate2(u32 op_count + [CigarOp; op_count] + u32 base_count +
[u8; base_count])`，CigarOp 为 `(op << 29) | len` 的 u32。

**X/I 差异碱基存储（2026-08-05 核实）**：CIGAR 本身只存操作 + 长度，
`X`（mismatch）/`I`（插入）的碱基内容收集进 **`xi_bases`**（差异碱基
流），随 CIGAR 一起 flate2 压缩存储（`u32 base_count + [u8; base_count]`）；
解码 `apply_cigar` 按 X/I 出现顺序从 xi_bases 取回，`=` 段从参考取、
`D` 段跳过参考。**mismatch/插入碱基不丢失**：解码端校验
`X/I 消费数 == xi_bases 长度` 且 `参考消费 == ref 长度`，不一致即报错
（数据损坏时拒绝输出，而非静默丢碱基）。`=` 段依赖参考 2bit 记录完整。

关键决策（详见旧版决策记录，已实现）：
- **段级回退（v1006 及之前）**：最佳 alignment 未完整覆盖段、段跨多条
  alignment 衔接、CIGAR target 投影跨参考段边界 → 整段回退 LZ-diff
  （不拆段、不合并）。**v1007 起该路径退役**——改为段级混合编码
  （覆盖部分 CIGAR + 未覆盖部分 LZ/Raw 补齐，见上）；
- **`M` 处理**：`maf to-paf` 源头已产出 `=/X`；对 minimap2 未用 `--eqx`
  的 `M` CIGAR 保留 `M → =/X` 拆分 fallback（压缩端读样本序列拆分）；
- **query-side 索引**：pbit 自建 `BasicCOITree<PafMetadata>`（不复用
  `PafIndex.reverse_trees`——其缺 `-` 链、元数据被交换、无公开查询接口）；
- **无/坏 CIGAR**：无 `cg:Z` → 跳过编码、原样存行还原（决策 7，
  2026-08-09 确认）；malformed CIGAR → 记录级跳过（log 警告）；整个 PAF
  不可用则报错终止（避免"以为生效实为全回退"）。

### Identity 段存储：零开销指向参考区间（v1010，2026-08-09）

**背景**：AGC 中样本段与参考完全相同（LZ-diff delta 为空）时，段描述符
直接指向组内参考序列（`in_group_id = 0`），不产生任何 delta 载荷
（`segment.cpp::add` 的 `IMPROVED_LZ_ENCODING` 分支）。pbit 此前即使
纯 `=` 段也打包 CIGAR（`u32 op_count + op + u32 base_count`，再 flate2），
同一参考区间被多个样本重复引用时各存一份。

**实现（v1010）**：新增 `DeltaEncoding::Identity = 3`。CIGAR 切片经
`split_m_to_eqx` 后若全为 `=`（无 X/I/D），不打包 CIGAR，改为零载荷
delta：`packed_data` 为空，`ref_start/ref_end`（全局坐标）即样本段内容；
解码端直接 `read_ref_interval` 取参考区间（`-` 链照常 `rev_comp`）。
delta 按 `(encoding, is_rev_comp, raw_length)` 去重：同一参考组内所有
同向等长 Identity 段共享一个空载荷条目，区间差异由段描述符携带——
与 AGC "指向参考" 同构，而 pbit 因 delta 表按组存储，不需要特殊的
"id 0 = 参考"约定。

**一致性**：Identity 与 CIGAR 的 `=` 段语义相同（大小写不敏感比对、
soft-mask 由 `mask_blocks` 恢复），`to-fa`/`get_sample` 长度校验不变；
`pbit to-paf` 对 Identity 段合成单条 `=N` 参与链重建，`cg/cs/gi/bi`
与打包 CIGAR 完全一致（`cs:Z` 输出 `:N`）。

### Raw 段存储：为什么用 2bit 而不是 ASCII+flate2（2026-08-08，v1008）

**背景**：v1007 的段级混合编码中，PAF 未覆盖且 LZ-diff 也匹配不上的段
用 Raw 存储（无损兜底）。初版（v1006/v1007）Raw = `flate2(ASCII 碱基)`
（DEFLATE = LZ77 + Huffman）。用户提出疑问：为什么不用 2bit？

**实测（10 个真实 4 kb 段，0 N）**：

| 存储方式 | 大小（相对 ASCII） |
|---|---|
| ASCII 原始 | 100% |
| flate2(ASCII)（v1006/07 的 Raw） | 32.2% |
| 2bit 原始 | **25.0%** |
| 2bit + flate2 | 25.3% |

**结论**：
1. **2bit 比 flate2(ASCII) 小 7 pp**——flate2 压不到 DNA 的理论熵限
   （4 符号 = 2 bit/碱基），Huffman 码长限制 + DEFLATE 开销让它停在
   ~32%；2bit 精确 25%，且**再压几乎无空间**（25.3%，DEFLATE 头尾
   开销反超收益）。
2. **v1008 起 Raw 段 = 标准 2bit 记录**（与参考段同构：
   `write_2bit_record` / `read_2bit_record`，N 用 2bit 的 N-block 结构
   表示，样本 soft-mask 仍由 `mask_blocks` 单独存储）。不再 flate2。
3. 差异小的原因是 Raw 段在 LZ 兜底后占比极小（实测 13/1947 段）；
   但 2bit 实现几乎零成本（复用参考段的编解码路径），且语义上
   "Raw = 参考同构的 2bit" 更自洽。

## 大链与碎链：定义与判定（2026-08-09 澄清）

> 状态：**定义与判定细节已定稿并实现（2026-08-09）**，真实数据验证见
> 本节末"实现与验证"。本节先记录 v1009 的旧现状（含失误），再给出定稿。

### 定义（用户拍板，2026-08-09）

**主链/碎链都指 sample（query）轴上的概念**：判定主体是链在 sample 上的
query 区间（覆盖、重叠、贪心选择都按 query 坐标）；reference 是编码锚
与坐标系，不参与主/碎判定（reference 坐标只用于 CIGAR 的参考区间定位）。

- **大链（主链）**：**能匹配上 Reference 的比对链**——**整链参与编码**
  （段级 CIGAR/Identity），`to-paf` 时按 `paf_id` 归组重建回原记录。
  判定**不看长度**（v1009 的 10 kb 阈值是失误，见下）。
- **碎链**：**被其他大链覆盖的链**——即使其相似度更高，只要落在主链
  覆盖范围内，就是碎链；**原样存行**（PAF 恢复区），`to-paf` 原样还原。
  **碎链仍可编码其覆盖的段**（2026-08-09 确认）——主/碎分类只
  决定"重建 vs 存行"，不决定"能否编码段"：碎链覆盖的段照常可用它的
  CIGAR 编码（压缩样本），只是 `to-paf` 输出存的行而非重建。相似度
  **不参与**大/碎的主体判定。

语义上等同于 primary/secondary。**分类本质（2026-08-09 讨论明确）**：
每条 PAF 记录（链）二选一——主链 = 整链编码 + 重建；碎链 = 原样存行。
两种都在归档里，**PAF 逐条可还原是硬约束**（"存进去什么，出来什么"）：
先保证还原，再优化压缩；压缩率是主链集合内的次级目标。因此**不能跨记录
合并编码**——合并后无法还原出原来的多条记录。

**主链选择标准（2026-08-09 讨论明确）**：同一 sample 区间与 reference
多次匹配时，选择哪一条作为压缩主链，**不纯粹由相似度决定**（相似度高 →
CIGAR 短，只是次级收益）。"连续"是单条 PAF 记录的天然属性（记录即连续
对齐单位），**不设独立判定**（详见判定细节第 3 条）；相似度降级为
**同覆盖下的次级 tiebreak**。

### v1009 旧实现现状（2026-08-09 已被新判定取代，保留作历史记录）

`compressor.rs` v1009 时**不是**按覆盖关系判定，而是：

1. **候选主链按长度过滤**：`BIG_CHAIN_MIN_LEN = 10_000`（span ≥ 10 kb
   才算候选，`try_encode_segment_cigar`）；span < 10 kb 的记录不参与
   CIGAR 编码，段直接走 LZ/Raw 兜底。
2. **best 主链按"覆盖度 → 相似度"选**：对每个 4 kb 段，在候选主链中选
   覆盖度最大者，覆盖度相同再比 `gap_compressed_identity`——**相似度
   参与了主链选择**，与"碎链即使相似度更高也算碎链"相反。
3. **恢复区分类也按长度**：`append_sample_with_paf` 里
   `is_complete_big = span ≥ 10 kb && 所有跨段都被编码` → 重建（ms 表）；
   其余（span < 10 kb，或被更大重叠链抢走部分段的大链）→ 原样存行。
4. **覆盖关系只在段级隐式存在**："嵌套小链不许抢大链"靠 10 kb 过滤
   实现，并没有链级的主/碎消解步骤。

实测影响（`benchmarks/bench-scale-and-pbit.md` #14k，00_3076 vs 00_3230）：809 条
PAF 记录中按长度分类，碎链行原样存储约 ~100 KB/样本（大头 cg/cs 文本），
PAF 恢复区使 delta/gzip 0.356 → 0.448（+9 pp）。

### 与定义的差距（失误点）

1. **10 kb 阈值是武断的**：能匹配参考但 span < 10 kb 的链被排除在编码
   之外（回退 LZ/Raw），不能匹配参考的判定被长度替代。
2. **相似度参与主链选择**：与"碎链即使相似度更高也算碎链"矛盾。
3. **碎链定义被长度化**：v1009 的"碎链"实际是"小链 + 不完整大链"，
   不是"被主链覆盖的链"。

### 判定细节（定稿，2026-08-09）

> 2026-08-09 澄清：此前讨论的"覆盖轴"指 **以谁为坐标**（reference vs
> query），不是链与链之间的重叠；编码以 reference 为坐标。本定稿按
> 该理解 + 既定原则（还原优先、参考是压缩锚、相似度仅作次级 tiebreak）
> 写成，由代码效果验证。

1. **还原优先**：每条链要么"整链编码 + 重建"（主链），要么"原样存行"
   （碎链）——PAF 逐条可还原是硬约束，压缩率在此前提下优化；不能跨
   记录合并编码。"存行"不排斥编码：碎链的行原样存（供还原），其覆盖
   的段仍可用它的 CIGAR 编码（压缩样本），两者互不影响。
2. **主链选择 = 链级贪心（无长度阈值）**：
   - 对每个样本，取全部有 `cg:Z` 的记录（无 `cg:Z` 者不参与编码，行
     仍保留以便还原）；
   - 排序键：**query 覆盖段数**（`ceil((qe - qs) / segment_size)`，即
     该链能让多少个 4 K 段走 CIGAR）降序 → **相似度**（gi）降序 →
     **输入序**（record_id 升序，可复现）；
   - 贪心依次取链：query 区间与**任何已选主链重叠**（相交即重叠）→
     碎链；与所有已选主链不相交 → 主链。
   - 完全包含、部分重叠、完全重叠统一由此规则覆盖：
     覆盖段数多者先选为主，与其重叠者（即使相似度更高）归碎——与
     "被覆盖的链 = 碎链"一致；完全重叠时 span 无区分度，落到相似度
     tiebreak。
3. **"连续"不设独立判定（2026-08-09 确认）**：PAF 记录（链）本身就是
   比对器保证的连续对齐单位；链内部 gap 由 CIGAR 的 `I`/`D` op 表达，
   段级切片照常处理——**不需要断链判定**。主链 = 整条 PAF 记录（天然
   连续）。
4. **段级编码（主链 + 碎链都可编码）**：每个段在覆盖它的**所有有
   `cg:Z` 的链**（主链 + 碎链）中选覆盖度最大者，用其 CIGAR/Identity
   编码；任何链都不覆盖的段 → LZ/Raw 兜底（语义不变）。主链的跨段若
   被覆盖度更大的碎链抢走 → 该主链不完整，按第 5 条降级存行。
5. **恢复区**：主链所有跨段都被它自己编码 → 整链重建（只存 ms 表）；
   主链有跨段未编码（被碎链抢走或切片失败）→ 降级存行；碎链 + 无
   `cg:Z` 记录 → 原样存行（碎链即使编码了段也存行）。
6. **`BIG_CHAIN_MIN_LEN` 退役**：主链判定不再依赖长度，常量删除。

### 实现与验证（2026-08-09）

- 实现：`PafQueryIndex` 保留无 `cg:Z` 记录的原始行（存行还原）；
  `Compressor::select_main_chains` 链级贪心 + 段级编码放开到所有有
  `cg:Z` 的链（主 + 碎）；恢复区按主链完整与否分类。`BIG_CHAIN_MIN_LEN`
  已删除。
- 真实数据（00_3076 vs 00_3230，809 行 PAF）：`to-paf` 核心字段
  （qname/坐标/strand/matches/block/cg/ms）**100% 逐字段一致**；编码分布
  Cigar 1231 + Identity 15 + LzDiff 688 + Raw 13，与旧实现（10 kb 阈值）
  一致——判定改造未改变压缩结果（该数据中小链均嵌套于大链覆盖内）。

## 文件格式规范（v1010）

> 本规范描述 v1010 现状；各版本格式变更见对应章节：遮蔽（v1005）、
> CIGAR 任意参考区间（v1007）、Raw 段 2bit 化（v1008）、PAF 恢复区 +
> `paf_data_offset`（v1009）、Identity（v1010）。

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
├─────────────────────────────────────┤  ← footer.paf_data_offset
│ PAF Recovery (v1009)                │  ← flate2(PafRecovery：大链 ms 表 + 碎链行)
├─────────────────────────────────────┤
│ Footer (固定 32 字节)               │
└─────────────────────────────────────┘  ← EOF
```

### Header（36 字节）

```
0  4  magic              0x54494250 ('PBIT')
4  4  version            major*1000 + minor（当前 1010）
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
    u8  encoding        0 = LZ-diff, 1 = CIGAR, 2 = Raw, 3 = Identity
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
    for each segment（固定 24 字节）:
      u32 ref_group_id
      u32 delta_id
      u32 ref_start     参考文件全局坐标（v1007 起，CIGAR/Identity 用；LZ-diff/Raw 填 0）
      u32 ref_end       同上（参考文件全局坐标）
      u32 q_start       样本 contig 内偏移（v1007 起）
      u32 paf_id        源 PAF 记录 id（v1009 起；LZ-diff/Raw 填 u32::MAX）
str cmd_line
```

### Footer（32 字节）

```
u64 ref_index_offset / u64 delta_data_offset / u64 sample_index_offset /
u64 paf_data_offset（v1009 起）
```

> **设计要点**：无文件尾 magic（Header 已含）；全部固定大小字段（解析简单、
> 破坏可校验）；`ref_group_count`/`sample_count` 在 Header/Index/Delta/
> Sample 中重复——读取各 section 时就地校验一致性，避免损坏文件越界。

## PAF 驱动编码的演进（#14 诊断 → v1010）

> 本节的路线与约束来自 `genome-nn-query.md` §8.5（原"pbit CIGAR 编码重构
> 建议"），因属 pbit 实现/设计，已收拢至此统一维护。

**起点（#14 诊断）**：初版 CIGAR 路径要求**段相位对齐**（单条 PAF 记录
全覆盖 4096 bp 段），真实基因组的 indel 即破坏、段大小调参无效；
LZ-diff 兜底又要求样本 contig 与参考**同名**，跨组装样本无法走通。
三条可选路线（按"改动小 → 收益大"排序）：

1. **LZ 兜底内容匹配化**（→ **v1006，2026-08-08 已实现**）：把 LZ-diff
   的"按 contig 名找参考段"改成"按内容找相似参考段"（canonical k-mer
   倒排索引 + `best_ref_group`）。不改归档格式、不依赖长链对齐，任意
   组装命名可用；压缩率低于 CIGAR 但**立即可用**（真实近缘样本 delta =
   gzip-9 的 53%，~100% 无损，`benchmarks/bench-scale-and-pbit.md` #14f）。
2. **跨相位 CIGAR 编码**（→ **v1007，2026-08-08 已实现**）：delta 引用
   "任意参考区间"而非固定段（`SegmentDesc.ref_start/ref_end` 改参考文件
   全局坐标 + 新增 `q_start`），按链/段混合编码（CIGAR 覆盖部分 +
   Raw/LZ 补齐）。压缩率最高（真实 98.6% 对 delta/gzip 54%→39%）。
   后续：v1009 `to-paf` 无损还原输入 PAF（大链重建 + 碎链原样存行，
   delta/gzip 0.448）；v1010 Identity 零载荷（纯 `=` 段）。
3. **pgi 长链链化**（**明确不做，2026-08-09 用户裁定**）：minimap2 式
   chaining 产出长链、满足"单记录全覆盖段"的路线**不实现**——项目优势
   是引入 UCSC chainnet 经典链化管线（效果最好），自研 chain 始终不如
   它；链化依赖由 chainnet 承担。对重排多的基因组链化仍会失败，收益
   有限。（`psl to-paf` 无 cg:Z；链级 `cg:Z` 生产者同样**明确不做**，
   推荐链路 `chainnet → maf to-paf` 自带 cg:Z。）

**约束现状**（随版本解除）：
- contig 同名约束 → v1006 内容匹配 + v1007 PAF 驱动，已解除；
- 段相位对齐约束 → v1007 任意参考区间，已解除；
- 长链全覆盖约束 → v1009 起不再要求单记录全覆盖段；**2026-08-09 主链/
  碎链按覆盖关系判定（链级贪心），长度阈值（`BIG_CHAIN_MIN_LEN`）退役**，
  碎链也可编码其覆盖的段（见"大链与碎链"章节）。

## 遮蔽（soft mask）处理：设计决策（2026-08-05）

> 状态：**已实现（v1005，2026-08-05）**。`collection.rs`（格式 + 版本
> 1005）、`compressor.rs`（两条编码路径提取 mask）、`decompressor.rs`
> （`get_sample` 应用 mask 还原小写）、docs/pbit.md、回归测试（lib
> `test_get_sample_roundtrip_soft_mask` + cli 两处断言更新）均已落地。

### 现状与问题

**参考遮蔽**：完整保留。参考 FASTA 小写（soft mask）经
`write_2bit_record(do_mask=true)` 存入 2bit `mask_blocks`；查询/读序列时
（`read_2bit_record` 默认 `no_mask=false`）还原小写。**但 delta 编码读
参考时 `no_mask=true`（全大写）**——mask_blocks 只服务读取/查询，不参与
压缩编码。

**样本（query）遮蔽**：半保留、**不对称**（当前实现的副作用，非设计）。
`read_fasta` 原样读入小写；CIGAR 编码时：

- `M` 段用 `eq_ignore_ascii_case` 比较（大小写不敏感）→ `=` / `X`；
- `X`/`I` 差异段的样本碱基**原样**进 `xi_bases`（小写保留，解码时原样
  还原）；
- `=` 匹配段从参考解码（大写）——样本遮蔽信息丢失。

结果：同一样本存进 pbit 再取出，**匹配段大写、差异段小写，遮蔽状态取决于
该段是否与参考匹配，完全随机**。这不是可用语义。

**LZ-diff 路径（对照）**：**不保留**样本遮蔽——样本与参考编码端即转 2bit
（`encode_base` 大小写不敏感），解码走 `decode_base` 恒输出大写
A/C/G/T/N；`test_lowercase_input`（segment.rs）与
`test_encode_decode_roundtrip_lowercase`（lz_diff.rs）已固化"小写输入 →
大写输出"。方案 B 下 LZ-diff 编码路径同样不用改（2bit 层处理大写），
样本遮蔽统一在解码端应用 mask_blocks 还原。

### 方案对比

| 方案 | 做法 | 代价 | 适用 |
|---|---|---|---|
| **A：样本统一大写** | 编码前样本序列 `make_ascii_uppercase`（含 `xi_bases`），遮蔽不进 delta | 小 | **否决**——与 2bit 血统相悖，存进存出不一致 |
| **B：样本存 mask blocks** | `Collection` 每样本/contig 存 mask_blocks（同 2bit 语义），解码还原小写 | 格式加字段 + 每样本存储 | **决策**——继承 2bit 遮蔽存储 |
| **C：遮蔽对齐（工作流层）** | 参考与样本用**同一套重复注释**遮蔽（`pgi build --mask` / FastGA `-M` / `fa mask`），pbit 不感知遮蔽 | 无代码，流程约定 | 与 A/B 配合，保证上游比对与存储语义一致 |

### 决策：方案 B + C（修正，2026-08-05）

**pbit 的格式语义：继承 2bit 的遮蔽存储**——pbit 从 2bit 演化，2bit 的
`mask_blocks` 本就是遮蔽（小写区间）的标准存储；参考已如此（2bit 记录
保留 mask_blocks），**样本同样应存储遮蔽，存进存出保留小写**（存小写、
取小写），不给用户"存遮蔽、取无遮蔽"的意外。

1. **样本存储遮蔽（方案 B）**：`Collection` 为每个样本/contig 增加
   `mask_blocks`（小写区间，格式与 2bit 的 mask_blocks 一致：starts +
   sizes）。编码时 `read_fasta` **保留小写**：序列转大写做 delta（2bit
   层语义），小写区间提取为 mask_blocks 随样本存储；解码时 delta 还原
   大写序列后应用样本 mask_blocks 还原小写——**存进存出一致**。
2. **参考遮蔽不变**：参考仍用 2bit `mask_blocks`（现状），读取时还原
   小写。
3. **遮蔽是工作流层概念（方案 C）**：上游比对（`align pgi` / FastGA
   `--mask`）用与参考一致的遮蔽生成 PAF；pbit 存储层只负责"遮蔽随样本
   存进存出"，比对/压缩不因遮蔽改变编码路径。
4. 方案 A（样本丢弃遮蔽、统一大写）**否决**——与 pbit 的 2bit 血统
   相悖（2bit 存遮蔽），存进存出不一致会让用户疑惑。

### 实施步骤（2026-08-05 已全部落地）

1. **`collection.rs`**：`ContigSegs` 增加 `mask_blocks: Vec<(u32, u32)>`
   （小写区间，0-based，与 2bit mask_blocks 同语义）；`Collection` 序列化
   写入/读取该字段（格式变更，版本 1004 → 1005，仅新版本可读写）。
2. **`compressor.rs`**：`read_fasta` 读取样本时**保留小写**，提取每个
   contig 的小写区间为 mask_blocks；序列转大写后进入 delta 编码（CIGAR
   delta 与 LZ-diff 路径不变，2bit 层本就处理大写）；mask_blocks 写入
   `Collection`。
3. **解码**：`to-fa` / `some` / `range` 还原大写序列后应用样本
   mask_blocks（转小写）——存进存出一致（小写保留）。参考仍按现有
   `read_2bit_record` 的 mask_blocks 还原。
4. **文档**：`docs/pbit.md` 注明"pbit 样本/参考均存储遮蔽（小写区间，
   同 2bit mask_blocks 语义），存进存出一致"。
5. **测试**：
   - 新回归：小写样本 FASTA → create → `to-fa` 解码 = 小写还原（遮蔽
     保留、存进存出一致）；含 X/I 差异段与 LZ-diff 段混合样本；
   - 现有 `cli_pbit*` 全量回归（大写输入行为不变）。

### 与 `-S`（对称 adaptamer）的关联

遮蔽对齐后（方案 C），`-S` 多找到的"更多重复比对"在 pbit 场景被遮蔽过滤，
且 E. coli 实测（未遮蔽）已无归档收益（覆盖 +0.9%、归档 +1 字节，见
[[pgi-align.md]] §7.4）——**`-S` 对 pbit 无帮助的结论在遮蔽流程下更稳健**。

## 统一序列访问 API（内部 pbit、暴露 twoBit）：评估 → 不做（2026-08-05）

> 状态：**不做**（通盘评估后否决；替代方案见下）。与遮蔽方案（v1005）
> 正交，不影响遮蔽实现。

**方案**：定义统一 `SequenceSource` trait（contigs / read_range / blocks /
has，twoBit 风格），底层 `TwoBitSource`（包装 TwoBitFile）+ `PbitSource`
（pbit 归档参考段，2bit 记录读取），现有 twoBit 消费者（psl chain、twobit
命令族、net/chain、pgi build 等）迁移到 trait——实现"内部 pbit、对外暴露
twoBit 风格接口"。

**动机（用户提出）**：pbit 归档作为群体基因组的**最终分发格式**时，消费方
可能只有 .pbit，需要能像 twoBit 一样访问参考序列。

**评估（通盘，含自我质疑）**：

1. **"分发格式" ≠ "消费者必须直接读归档"**：pbit 参考段本就是标准 2bit
   记录，导出一个标准 twoBit/FASTA 是一行命令的量级，且外部工具也能消费；
   统一 API 省掉的只是"先导出"一步，价值密度低，却要承担 10+ 消费者
   trait 化的改造与回归风险。
2. **私有格式做"最终分发"有生态锁定风险**：pgr 在 fastga.md 中批评过
   FastGA 的 GDB"生态锁定"（私有、外部无法消费）；pbit 私有分发同样如此。
   分发物应尽量标准（twoBit/FASTA）或提供标准导出，而非让消费方必须用
   pgr 内部 API。统一 API 强化了"只能在 pgr 里消费"的锁定。
3. **统一抽象掩盖 pbit 本质**：twoBit 是独立随机访问的序列文件，pbit 是
   压缩容器（参考 + delta 样本）；让归档"伪装成序列文件"丢失"在读归档"
   的语义，长期是混淆而非整洁。
4. **违反简洁原则**：为"分发后直接随机查询归档参考"这一未被真实工作流
   证实的场景做大型抽象，且样本段（delta）本就不在 twoBit 语义内（统一
   接口实际只覆盖参考段，收益面窄）。

**结论：不做统一 API 层**。若"pbit 作为最终分发格式"是真实方向，正确的事
是让分发物可被标准工具消费，而非让 pgr 内部 API 适配私有归档。

**替代方案（若分发需求出现）**：

- `pbit` → 参考 twoBit / 样本 FASTA **导出命令**（复用现有
  `read_2bit_record` / `to-fa`），外部工具可消费；
- 个别消费者需要直接读归档参考时，做 `PbitRefReader` 轻量适配（单实现，
  不 trait 化），按需接入。

**与遮蔽方案的关系**：遮蔽（v1005）是格式层变更，统一访问 API 是访问层
构想，两者正交；本否决不影响遮蔽方案实施。

## 多参考（v1003/v1004 扩展）

泛基因组增量场景：少量基因组起步、逐个添加；参考（锚）被反复比对，样本
逐个加入。最终产物是图（GFA），pbit 与索引都是可重建的中间产物，故格式
演化不做长期兼容负担。

**索引定位（决策 A，2026-08-02）**：`.pgi` **不进 pbit**。实测单参考
归档 1.1 MiB、内嵌索引后 92.8 MiB（索引/压缩数据 = 79×）——内嵌会让
pbit 失去"压缩格式"的意义。参考索引在需要比对/距离时现建（~0.3s）或由
工作流在归档旁缓存为 `ref.pgi` 兄弟文件；与 FastGA"GIX 独立文件、用完即删"
定位一致。HV sketch 内嵌（决策 B）暂缓，其算法设计待后续思考。
（**2026-08-09 变更为明确不做**：HV 评测未达预期，后续换其他形式。）

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
> 详细数据见 `notes/design/hv.md`。

## 附录：早期开放项决策稿（2026-08-03，决策已全部落地）

> 目的：把 [[pbit.md]] 顶部暂停的 5 个开放项 + 决策 B（HV sketch 内嵌）
> 收敛为"选项 + 证据 + 推荐"，作者逐项确认即可解锁继续开发。
> 日期：2026-08-03。所有证据来自 `src/libs/pbit/`、`src/cmd_pgr/pbit/`
> 与 `notes/design/hv.md` 的当前实现。

> **✅ 已确认（2026-08-03）**：作者按推荐全部采纳（1A 2A 3A 4A 5A 6A）。
> 1A 的"多参考未指定路由时警告"已实现于 `cmd_pgr/pbit/mod.rs`
> `resolve_ref_id`。本文保留选项对比供未来重新评估。

> **状态核对（2026-08-09）**：6 项决策全部仍有效且已落地（详见
> [[pbit.md]] 顶部"已决策"摘要）；唯一变更是**决策 B（HV sketch 内嵌）
> 由"维持暂缓"改为"明确不做"**（HV 评测未达预期，后续换其他形式，见
> `todo.md` §4 与 `pbit.md`）。本文档定位 = 早期决策过程记录（选项对比 +
> 证据 + 拍板清单），权威状态以 `pbit.md` 为准。

### 背景事实（影响全部决策）

1. **路由只影响压缩率，不影响正确性**：任何参考都能压缩任何样本（LZ-diff
   找最长匹配），路由只决定"锚"选得好不好 → 压缩率差异。
2. **pbit 是中间产物**：泛基因组流程中归档、索引、图都可重建，格式演化
   不做长期兼容负担（已确认）。
3. **样本段已带全局 `ref_group_id`**（Sample Index 每条 segment 16 字节
   固定字段），跨参考路由在**格式上已经可表达**，只是压缩器按"单参考"
   实现。

---

### 开放项 1：样本 vs 参考的路由

**现状**：用户指定（TSV 第 4 列参考名/序号）；默认参考 0；
`-i` 模式固定参考 0（`cmd_pgr/pbit/create.rs:107`、`append.rs:75`
经 `resolve_ref_id` 解析）。

**选项**：

| 选项 | 含义 | 优点 | 缺点 |
|---|---|---|---|
| A. 保留手动（推荐） | TSV 第 4 列指定；多参考且未指定时**警告**而非静默用参考 0 | 零行为变更；用户有完全控制 | 多参考时漏填第 4 列会静默选错锚 |
| B. 自动路由 | 按 k-mer 草图相似度（`dist seq` k=8 syncmer 风格，O(序列长)）自动选最相似参考 | 消除人为错误；无 TSV 依赖 | 每样本多一次参考草图标定成本；与手动的结果可能不一致，回放困难 |
| C. 手动 + `--auto-route` 开关 | A 为主，未指定时可用开关走 B | 兼容现有流程，行为显式 | CLI 多一个参数 |

**证据**：`dist seq`（k=8 syncmer 草图）是当前与身份率最贴近的粗筛层
（Spearman 0.82，见 design/hv.md 证据附录），成本 O(序列长)；
大肠杆菌任意两株 ≥90% 身份，锚的选择对压缩率影响很小。自动路由的真正
价值在**跨物种/多样性 cohort**（此时参考间差异大）。

**推荐**：**选项 A 起步**——只补一个"多参考 + 未指定路由 → 警告"的
改进；`--auto-route`（选项 C）留到真实多样性 cohort 出现、且能用
压缩率数据证明收益时再实现（避免推测性开发）。

### 开放项 2：Sample Index 是否加 ref_id

**现状**：不加。segment 存 `ref_group_id`（全局唯一，跨参考连续编号），
经 Reference Table（`group_start`/`group_count`）反查 ref_id。

**分析**：
- 加 ref_id = 每条 segment 16→20 字节（+25% 索引体积）+ 与 Reference
  Index 重复数据、存在不一致风险。
- 反查成本：参考数通常个位数，线性/二分扫 Reference Table 可忽略。
- "反查可得"的简洁决策已被 v1004 落地，无已知消费方需要热路径反查。

**推荐**：**保持不加**。若未来出现每段热路径反查（不预期），在
Reference Table 上补二分查找即可，无需格式变更。

### 开放项 3：`append-ref` 的语义

**现状**：只加参考、不改已有样本路由（`append_ref.rs`：追加 2bit 段 +
重写 Reference Index，旧样本段与 delta 原样保留）。

**选项**：

| 选项 | 含义 | 优点 | 缺点 |
|---|---|---|---|
| A. 保持"不重路由"（推荐） | append-ref 只追加，样本锚不变 | O(参考段数) 轻量；语义清晰 | 换到更好的锚需手动重建 |
| B. 自动重路由/重压缩 | append-ref 后把已有样本对全参考重压缩 | 压缩率随参考集优化 | 一次 append-ref 变 O(全部样本重压缩)；隐式巨量耗时，违反最小惊讶 |

**证据**：增量流程里样本通常在其最佳参考已存在后才加入（先 `create -r ref`
再逐批 `append`），"换锚"是例外场景而非常规。重压缩 = 重新读全部样本 +
重跑 LZ-diff，在 4 万级 cohort 上是小时级操作，不应藏在 `append-ref` 里。

**推荐**：**选项 A**。若"换锚"需求出现，做成显式独立子命令
（如未来的 `pbit re-anchor`），不在 append-ref 中隐式触发。

### 开放项 4：多参考压缩模型

**现状**：样本整体路由到一个参考（`Compressor::set_cur_ref_id` 每样本
设置；`append_sample` 的 contig→ref_group_ids 查找限定在当前参考的
段范围）。**格式已支持每段 ref_group_id**，跨参考在格式层无阻碍。

**选项**：

| 选项 | 含义 | 压缩率影响 | 成本 |
|---|---|---|---|
| A. 单参考/样本（推荐，现状） | 样本全部 contig 用一个锚 | 多样 cohort 中"染色体 vs 质粒"等跨参考混合场景吃亏 | 无 |
| B. 按 contig 路由 | 每个 contig 独立选参考 | 中等 | 压缩器改造：逐 contig 做方向检测 + 参考段查找 |
| C. 按段/混合参考 | 段内跨参考匹配 | 理论最高 | LZ-diff 跨参考索引、PAF 投影歧义、复杂度大增 |

**证据**：压缩器以"contig 级方向检测 + 单参考段集合"组织（`compressor.rs`
`append_sample`），B 需要把方向检测与段查找按 contig×参考重做；C 会破坏
"段→单参考"的简单映射，且收益需实测证明。

**推荐**：**选项 A**，但把"格式已支持每段 ref_group_id"写进设计笔记
（v1005 若做 B 无需格式变更，只是压缩器改造）。触发条件：出现真实
跨参考混合 cohort 且压缩率差距可量化。

### 开放项 5：版本策略

**现状**：已确认不做旧版本兼容，格式改动直接 bump 版本（当前 1004）。

**推荐**：**维持**。补充一条逃生舱：若未来出现"长期归档"需求，用
`convert`（读旧版写新版）显式迁移，而不是让读取端做多版本兼容。
这是决策而非开放项，建议在 pbit.md 中把该条从"待决策"移到"已决策"。

---

### 决策 B：HV sketch 内嵌（算法设计笔记）

**目标**：让 pbit 归档自带样本间距离能力（`dist hv` 风格），不依赖
`.pgi`（决策 A 已排除内嵌索引）。

**现状与证据**：
- `.hv` 目前只能从 `.pgi` 投影得到（`pgi to-hv`，`libs/pgi/to_hv.rs`），
  消费端 `dist hv a.hv b.hv` 直接比较；
- 修复后的稀疏投影 `.hv` 是 `dist pgi` 的 50× 快、排序 ρ=0.97 的近似层；
  但 dim=1024 的旧稠密投影在大规模集合上饱和退化（Spearman −0.05，不可用，
  见 design/hv.md 证据附录）；
- `dist seq`（k=8 syncmer 草图）不依赖任何索引，直接从 FASTA 计算，已是
  最优粗筛层（Spearman 0.82）。

**核心判断**：`dist seq` 已经从"源 FASTA 在时"的粗筛需求；内嵌 sketch
只对"源 FASTA 已删除、只剩归档"的场景有价值。这是窄场景，支撑"暂缓"。

**若实现（作为设计稿，不立即做）**：

1. 在 `create`/`append` 压缩样本时，直接从样本序列计算 sketch，
   **不经过 pgi**：k=8 syncmer 采样（复用 `libs/syncmer`）→ 每 k-mer
   更新 `--sparse` 个随机维度（复用 `libs/hv::hash_hv_sparse`）→ 稀疏
   HV 存入 Sample Index（flate2 块内，v1005 格式 bump）。
2. 规模估算（需实测校准）：dim=4096 × i32 = 16 KB/样本（稠密全量）；
   稀疏存储只存被更新维度 ≈ (k-mer 数 × sparse) 条。4 万样本量级约
   0.6–2 GB，占归档总量比例需 cohort 数据实测。
3. 消费者：新增 `pgr pbit dist`（归档内样本两两距离/近邻），复用
   `libs/hv::calc_distances` 的余弦实现。
4. 对比基线：实现后须与 `dist seq`（FASTA 侧）在 10 株 cohort 上对比
   Spearman，证明内嵌值不值得那部分体积。

**推荐**：**维持暂缓**。触发条件 = 出现"无源 FASTA、仅归档、需距离
粗筛"的真实工作流；届时按上述设计稿实现并先跑 4 的对比验证。

---

### 拍板清单（逐项回复即可）

1. 路由：A（手动 + 多参考警告）/ C（加 `--auto-route`）？
2. Sample Index 加 ref_id：否（推荐）/ 是？
3. append-ref：不重路由（推荐）/ 需显式 re-anchor？
4. 多参考压缩：单参考/样本（推荐）/ 按 contig 路由？
5. 版本策略：维持"不兼容 + bump"（推荐）？
6. 决策 B（HV sketch 内嵌）：维持暂缓（推荐）/ 按设计稿实现？

## 参考资料

- 多参考扩展设计过程与 .pgi 消费者规划：早期设计稿已并入本文（详见
  "多参考"章节与"附录：早期开放项决策稿"）；
- 距离消费者验证数据：`notes/design/hv.md`；
- pgi 比对管线（`.pgi` 的消费者）：`notes/design/pgi-align.md`；
- AGC 算法参考：`notes/references/agc-cpp.md`。
