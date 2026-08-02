# FastGA 功能差距与落地计划

> 对照 [[fastga.md]]（参考笔记）与 `pgr align pgi` 现状，列出 FastGA 中相对
> 重要、pgr 尚欠缺的功能，按优先级给出落地计划。日期：2026-08-03。

## 1. 差距清单

| FastGA 功能 | pgr 现状 | 重要性 | 建议 |
|---|---|---:|---|
| **soft mask 感知的种子发现**（`-M`/.1ano，mask 区不产生种子） | `pgi build` 把小写当正常碱基（codes 表映射 a/c/g/t）；`fa mask`、2bit mask_blocks 读取已存在 | 高 | **P1 做**：建索引跳过 soft-mask 区，抑制重复/低复杂度区假阳性 |
| **自比对模式**（`FastGA A`，单输入检测内部重复） | `pgr align pgi` 需两个输入；`sd search` 用外部 lastz `--self` | 中高 | **P1 做**：`pgr align pgi` 支持单输入 self 语义 |
| **PAF `cs:Z` 字符串输出**（`-pafs/S`） | PAF 已有 `cg:Z`，无 `cs:Z` | 中低 | **P2 可选**：CIGAR 可逆压缩标签，下游工具兼容 |
| **select 表达式**（只比对选定 contig/区间） | 无（需 fa range + 子索引间接实现） | 低 | 可选，暂缓 |
| **Gap_Improver**（wave 后 gap 区二次精修） | banded 仿射 gap 已覆盖；wave 路径无等价物 | 低 | 暂缓（质量微调，收益不确定） |
| **多 mask union**（.1ano 可叠加） | `fa mask` 单 runlist | 低 | 暂缓 |
| **trace points / ONEcode `.1aln` 紧凑存储** | PSL/MAF + BGZF | 低（人类规模才需要） | **不做**（规模不匹配，见 fastga.md §10） |
| **`-S` 对称 adaptamer** | `sd` 管线已覆盖重复分析 | 低 | **不做**（更慢，专门场景） |
| **ALNchain（.1aln 链化）** | UCSC chain/net（更标准） | — | **不做** |
| **GDB 格式 / scaffold 语义 / 完整 GIX 分片** | pgr 2bit + `.pgi` | — | **不做**（格式对比见 fastga.md §9） |

## 2. P1：soft mask 感知的索引构建（已实施 2026-08-03）

**目标**：`pgr pgi build` 支持跳过 soft-mask 区域（与 FastGA `-M` 语义一致），
使重复/低复杂度区不产生种子。

**方案**：

- 新增 `--mask` 选项（默认关闭，避免行为变化；FastGA 同样默认忽略 mask）：
  - FASTA 输入：小写碱基视为 masked；
  - 2bit 输入：读 `mask_blocks`（`fmt/twobit.rs` 已有读取与软 mask 应用）。
- 实现位置：`libs/pgi/build.rs`——masked 碱基按 N（code 4）处理，k-mer 窗口
  含 masked 即跳过（等价于 FastGA"采样点落在 mask 内不产生种子"，且更严格）。
  `read_fasta` 保留大小写；`read_2bit` 需把 mask_blocks 随序列返回。

**验证**：
1. 单元测试：含重复区的序列，`--mask` 建索引后 masked 区 k-mer 不存在；
   无 mask 时行为与现状逐字节一致（不回归）；
2. 集成测试：`pgr pgi build --mask` 小写 FASTA 与 2bit（mask_blocks）等价；
3. 真实数据：对含 IS/重复区的 E. coli 建 mask 索引，`pgr align pgi` 的
   重复区假阳性链减少、主链覆盖不降。

**结果（2026-08-03）**：`pgi build --mask` 已落地（`build_from_path` 增加
`mask` 参数；FASTA 小写与 2bit mask_blocks 统一转 N 跳过）。单元测试
`mask_skips_lowercase_fasta_kmers` / `mask_skips_2bit_mask_blocks` 验证
mask 后 k-mer 是未 mask 的子集；真实数据（MG1655，100 kb 区间小写）positions
下降 ~2.2%，与区间比例一致；CLI 集成测试
`command_pgi_build_mask_fasta_2bit_equivalent` 验证小写 FASTA 与 2bit
mask_blocks 在 `--mask` 下产出相同索引。923 测试全过。

## 3. P1：自比对模式（已实施 2026-08-03）

**目标**：`pgr align pgi A -o out.psl`（单输入，ref=query）检测基因组内部重复
与单倍型间同源（FastGA `FastGA A` 语义）。

**方案**：

- CLI：`query` 变为可选；仅给一个输入时进入 self 模式（同一索引两侧）；
- merge：self 时 canonical 去重防同一物理命中重复发射；
- 链化/扩展：跳过与自身完全重合的对角线段（diag=0 且区间相同的命中，
  FastGA self 语义），保留不同位置的重复（diag≠0 与反向重复）。

**验证**：
1. 单元测试：含串联重复的序列自比对，输出含重复块且不含 trivially 自身块；
2. ~~集成测试：与 `sd search` 的重复区结果对照~~（**已取消 2026-08-03**：
   `sd search` 引擎随后切换为 pgi（替代 lastz），lastz 对照失去意义）；
3. E. coli 自比对：主链外出现 rRNA/IS 等真实重复块（对照 pgi-align.md
   §3.1 v1 自比对 745 块、负链 186 块的已知结构）。

**结果（2026-08-03）**：`pgr align pgi` 支持单输入 self（query 可选；merge
后 `drop_self_hits` 过滤完全自身命中 = 同 contig 同位置同方向，与 FastGA
跳过 diag=0 一致）。单元测试 `drop_self_hits_filters_exact_identity` +
集成测试 `command_align_pgi_self`（串联重复无完全自身块）；MG1655 自比对
689 块、无全长主链（双输入对照 1021 块含主链）、无任何完全自身子块，重复块
（rRNA 等）保留。923 测试全过。

## 4. P2：PAF `cs:Z` 输出（已实施 2026-08-03）

**目标**：PAF 输出追加 `cs:Z`（minimap2/FastGA `-pafs` 规范：`:N` 匹配游程、
`*` 错配、`+seq` 插入、`-seq` 删除），与现有 `cg:Z` 并存。

**验证**：从 CIGAR 生成 `cs:Z` 的单元测试（含 `=`/`X`/`I`/`D` 混合），
roundtrip 可逆；`pgr paf` 下游（to-gfa/to-vcf）不受影响。

**结果（2026-08-03）**：`libs/paf/cigar.rs` 新增 `cs_from_alignment`
（FastGA `-pafs` 风格：`:N` 匹配游程、`*<ref><qry>` 错配、`+<qry>` 插入、
`-<ref>` 删除，可逆），`maf_block_to_paf`（MAF→PAF，`pgr maf to-paf` /
`pgr sd align` 共用）追加 `cs:Z:` tag，与 `cg:Z` 并存。单元测试
`cs_from_alignment_mixed_ops` / `cs_from_alignment_indels_and_length_check`
（= / X / I / D 混合、长度校验）+ `command_maf_to_paf_basic` 断言 cs:Z；
端到端 `cs:Z::8+A+A`。925 测试全过。

> 三个实施项合计 930 测试全过（2026-08-03 复核；sd 交叉验证测试随
> `sd search` 引擎切换方案取消，见 §3 验证点 2）。

## 5. 不做项与理由

- **trace points / `.1aln`**：ONEcode 私有格式 + 紧凑存储收益在人类规模才显著，
  pgr 的 PSL/MAF + BGZF 已覆盖当前规模（fastga.md §10.3）；
- **`-S` 对称模式**：重复分析由 `sd` 管线承担，对称模式更慢；
- **ALNchain / GDB / scaffold / 完整 GIX 分片**：pgr 有更标准的 UCSC
  chain/net 与 2bit/`.pgi` 等价物，格式对比见 fastga.md §9。

## 6. 相关文档

- 参考：[[fastga.md]]（FastGA 源码分析）
- 开发：[[pgi-align.md]]（比对管线设计 + 开发记录）
- 格式：[[../references/gfa.md]] 等（paf 生态）
