# pgr 近期待办

> 依据 `project-understanding.md`、各设计笔记与审计记录整理。
> 功能层基本齐备，近期大头是验证与数据驱动的扩展。
> 已完成条目只留一行结论，细节见链接文档。

## 0. 会话交接（2026-08-10：BBTools 替换逐命令核对与修复）

**已完成（trim.era.sh 8 步 M0-M8 全部移植，测试 1594 全绿，fmt/clippy
clean，代码已提交）**：

- `fq sample`/`trim-adapter`/`norm` 与 BBTools 39.38 逐字节一致复核完成：
  修复 changequality、qtrim 空 read 边界、`--ref` 可选、norm minq=6；
  reformat/bbduk 参数全景盘查（唯一缺口 = ihist，用户决定放着）；
  fairy/khmer 源码分析不改变 norm 精确/近似决策；bbnorm 深度分箱暂不做。
  细节见 `design/anchr-trim-replace.md` §4.0/§4.8–4.11、
  `references/fairy.md`、`references/khmer.md`。
- `pgr kmer` 用 Lambda 真实数据验证（table/hist/gsize/gc/qhist/qcheck/
  profile 全跑通；修复 gsize peak 估计：全局众数 → CallPeaks 主峰，
  peak=56 与 BBTools 一致、genome_size 47786 bp）。细节见
  `design/kmer.md` §10.8。
- `fq split` 多参数核对完成：stdout 输出 R1 与 repair golden 逐字节一致、
  尾记录（无 --outfile-single）warning + 丢弃均有测试；确认 `pgr::writer`
  不做 `.gz` 自动压缩（设计如此，压缩由 shell 管道负责）、`fq interleave`
  默认重命名 reads（roundtrip 逐字节不适用）。
- bbduk 第一梯队功能补齐（2026-08-10）：`fq trim-adapter` 新增 qtrim
  r/l/rl/w、polymer（poly-A/G/C + filter）、maq/mbq/maxnrate/mcb/mlf/
  maxlength、GC 过滤、forcetrim、kmask（N/lc/fully-covered）——14 组
  Lambda 真实数据黑盒对照 39.38/40.01 全部逐字节一致。细节见
  `design/anchr-trim-replace.md` §4.12。
- bbduk 两次调用拆分（2026-08-11）：`fq trim-adapter` 拆为
  `fq clean`（kmer 修剪 + 质量/组成过滤，对应 bbduk 第一次调用）+
  `fq filter`（kmer 污染过滤，对应第二次调用）；参数名统一 pgr 长名
  风格并在帮助标注 `(bbduk: 原名)`；`fq trim-qual`（sickle 语义）保留。

**挂账/待决**：

1. **norm 精确 vs 近似定稿**（§4.8 未定）：pgr 走精确表 + 外部桶（1TB
   答案）；bbnorm bits=16 近似表结果依赖 -Xmx。差异 = 定义差异不是 bug。
2. **anchr 模板替换**（用户自己处理，命令已齐）：把 trim.era.sh 的
   bbtools 调用换成 pgr 命令 + 管道串联（原语路线，pgr 不内置 pl trim）；
   merge.era.sh 的 bbduk 纯 qtrim → `pgr fq trim-adapter ... --no-ktrim
   --no-tbo --no-tpe --max-ns=-1 --force-trim-mod 0 --trim-quality <qual> --minlen <len>`。
   trim.era.sh 的 bbduk 两次调用 → `pgr fq clean`（trim）+ `pgr fq filter`
   （filter）。
3. **ihist**（2_insert_size.era.sh 的 reformat ihist）：SAM→insert size
   直方图，pgr 无 SAM 命令，用户决定放着。
4. **bbnorm 深度分箱**：暂不做（§4.9）。
5. `fq split` 多参数核对已完成（见上）；`fq clump` 多参数 golden 验证
   **不做**（体积控制，见 `design/anchr-trim-replace.md` §4.4 M1 注）。
   `kmer` 系列已用 Lambda 真实数据验证完成（见上）。

---

**历史会话（已完成，一行结论，细节见各设计笔记）**：

- **2026-08-09 `pgr kmer` 七子命令全量交付**（table/profile/hist/gc/
  qhist/qcheck/gsize，含 GenomeScope 2.0 原生迁移 `--model`、plot
  heat/spectra 拆分、`.pkt`/`.pkp`/`.hist` 三格式定稿），测试 1544 全绿
  → `design/kmer.md` §10。
- **fq 系列（2026-08）**：`trim-qual`（sickle 替代）、`range`（FASTQ
  `.loc` 索引一期 + 双端二期）、BGZF 写侧基准、FAFQ reader 笔记 →
  `design/fq-trim-qual.md`、`design/fq-index.md`、`design/seq-reader.md`。

## 1. 等数据/场景到位再启动

- [ ] 稀疏 s=1 完整 45 对 cohort 复测：等 10 株数据（5 株 10 对已验 ρ=0.988）。
- [ ] 4 万 E. coli cohort 端到端：核心步骤就绪，等真实 cohort 数据
      （来源：`ecoli-cohort.md`）。
- [ ] 人类规模（GRCh38/CHM13）验证：`.pgi` 字段上限、内存/耗时与 FastGA
      对照（来源：`design/pgi-align.md` §7.2）。
- [ ] pbit 自动路由：等多样性 cohort 数据证明收益（来源：`design/pbit.md`）。
- [ ] `--sym` 场景开关：先量化方向偏差再实现（来源：`design/pgi-align.md` §7.4.1）。
- [ ] 完整 adaptamer（变长种子 >k）：前置 lcp 已落地，只差立项
      （来源：`design/pgi-query-layer.md`）。
- [ ] `dist mash` 序列级并行：等单文件多 contig 大规模场景（文件级并行
      已覆盖多文件，见 `design/hv.md` 性能节）。
- [ ] 物种内聚类选参考 + PBit 归档（核心用例，UI 待讨论）：场景工作流见
      `design/genome-nn-query.md` §5；缺口 = `dist` 输出与 Necom 格式对齐、
      参考挑选、pbit 自动路由（`design/pbit.md` 决策点 1 触发场景）；
      HV 最近邻/SQLite/ANN 调研结论见 §6 与 `design/hv.md`（证据附录）、
      `benchmarks/bench-scale-and-pbit.md` #10b。
- [ ] 标定/检索剩余：E. coli NR 全量（15,574）实跑（2,088 已实测，
      `benchmarks/bench-scale-and-pbit.md` #8b）；pbit **多参考**样本验证待做
      （高分歧已在 #14g 实测：E. albertii ANI≈90%，delta/gzip 78–81%）；
      §7.4 #10/#9/#14/#19 状态见 `genome-nn-query.md`。

## 2. 低风险审计记录项（可顺手修）

- [ ] PAF `query_length`/`target_length` 恒 0：**待决策**——改 `.paf.idx`
      持久化 src_size，影响索引格式兼容性（来源：`audit/audit-paf.md`）。
- [ ] `syncmer.rs` 重复发射同一位置：**暂缓**——消费方已去重，收益小风险高
      （来源：`audit/audit-rept-sd.md`）。

## 3. 技术债（有空再议）

- [ ] `fas` 模块职责过重（20 子命令）考虑拆分（来源：`genome-nn-query.md` §8.6）。

## 4. 明确不做（避免重复立项）

- Gap_Improver、完整 LCP、`.1aln`、trace points、ALNchain、GDB/GIX 分片
  （`design/pgi-align.md` §6）；多 mask union（§7.5）；`-S` 对称 adaptamer
  （§7.4.1）；hybrid 逻辑留 `cmd_pgr/`（commit `d5281bc`，有意为之）。
- sd 边界扩展（BISER MAX_EXTEND）：已评估不做（`design/sd.md` §4.8/§4.9）。
- `dist mash` murmur3 SIMD：已评估不做（1.3× 差距可接受，风险大于收益）。
- 命令树 dispatch 宏简化：已评估不做（用户裁定——宏抹消对注册/分发代码
  的显式理解；`tests/cli_consistency.rs` 已核对注册一致性）。
- **pbit HV sketch 内嵌（决策 B）：明确不做**（2026-08-09 用户裁定——HV
  评测未达预期，后续换其他形式；原暂缓触发条件 = 无源 FASTA 归档需距离
  粗筛，见 `design/pbit.md` 决策 B）。
- **链级 `cg:Z` 生产者：明确不做**（2026-08-09 用户裁定——推荐链路
  `chainnet → maf to-paf` 自带 cg:Z；`psl to-paf` 无 cg:Z 记录走
  "跳过编码、存行还原"）。
- **碎链 cg 位打包：暂缓/不做**（2026-08-09 用户明确——别纠结碎链，
  整行压缩量可接受，见 `design/pbit.md` §大链与碎链）。
- **pgi 长链链化（pbit 路线 3）：明确不做**（2026-08-09 用户裁定——项目
  优势 = 引入 UCSC chainnet 经典链化管线，自研 chain 效果始终不如它；
  链化依赖由 chainnet 承担，见 `design/pbit.md` §PAF 驱动编码的演进）。
- **gzip 并行解压 / zlib-ng / libdeflate：明确不做**（2026-08-09 用户裁定——
  程序常被 shell 包裹并行执行，pgr 侧 `fa` 保持单线程；inflate 内部已是
  zlib-rs AVX2，见 `benchmarks/bench-profile-hotspots.md` 场景 1）。

## 5. 待实现 / 待决策（2026-08-09 文档扫描补充）

**已完成（一行结论，细节见链接）**：

- **spanr cover 名字截断**（2026-08-09 代码核对）：pgr 内建 runlist 区间
  操作替代外部 spanr；`rept/trf.rs` 带点 contig 名映射规避 →
  `design/repeat-masking.md`、`audit/audit-runlist-rg.md`。
- **wide 128-bit 化 / SIMD 梯队**（2026-08-09）：linalg（双累加器）/poa/hv
  改 128-bit；`paf::cigar` 分类掩码一次扫描（40 M 列 0.347 s，~37%）；
  `twobit::from_dna` 三级分类（pbit create 83→58 ms）；hv_benchmark 拆分
  现役/历史对照 → `design/simd-optimization.md` §6、
  `benchmarks/bench-simd-hv-jaccard.md`。
- **pbit 文件非确定性排查**（2026-08-09）：虚惊——collection 元数据嵌入
  完整命令行（含 `-o` 文件名），换输出名导致字节差；输出确定性已确认 →
  `design/pbit.md`。
- **用户文档改动清单 #22**（2026-08-09）：dist/pbit/align-pgi/pgi 文档
  已落地 → `design/genome-nn-query.md` §8.6。
- **FastK/Profex 原生迁移**（2026-08-09）：`rept s-kmer`/`e-kmer` 原生化，
  `--keep-index` 缓存升级单文件 `.pkt` → `design/kmer.md`。

待实现：

- [ ] **repeat masking：pgi 参数标定 + 真核验证**：CLI 透传已实现
      （`align pgi` 的 `-f/--min-shared/-k/--smer/--window`），但默认值
      未按 §2.5 调整（`--freq` 10 → 100、`--min-shared` 12 → 16 待验证）；
      真核（拟南芥/玉米等转座子丰富）与 RepeatMasker masked 输出对比
      recall（E. coli 无转座子无参考价值）；polyA/卫星低复杂度缺口由
      `rept trf` 兜底（来源：`design/repeat-masking.md` §2.4/§2.5）。
- [ ] **paf 查询层扩展（待实现）**：`--min-tree-coverage`（Caf Tree
      Coverage 过滤维度，查询时无法全图计算，作传递闭包后处理过滤）；
      `--end-trim` 推迟（需 per-interval 修剪 CIGAR，待序列输出引入时
      一并处理）（来源：`paf-pangenome.md` §Caf 过滤维度对照表）。
- [ ] **chain 算法待验证（低优先）**：KD-tree 已实现并用于 `psl chain`
      （`libs/ds/kdtree.rs`）；`best_crossover` 已接入 `fas_multiz` merge
      （`libs/ds/crossover.rs`）——两者的**真实数据验证**待做；KD-tree
      用于 PAF 链式化 / POA 排序仍待评估（PAF 当前未明确需要链式化）
      （来源：`chain-algorithms.md` §12.3）。
