# pgr 近期待办

> 依据 `project-understanding.md`、各设计笔记与审计记录整理。
> 功能层基本齐备，近期大头是验证与数据驱动的扩展。
> 已完成条目只留一行结论，细节见链接文档。

**挂账/待决**：

1. **norm 精确 vs 近似定稿**（§4.8 未定）：pgr 走精确表 + 外部桶（1TB
   答案）；bbnorm bits=16 近似表结果依赖 -Xmx。差异 = 定义差异不是 bug。
2. **anchr 模板替换**（用户自己处理，命令已齐）：trim.era.sh 的 bbtools
   调用换成 pgr 命令 + 管道串联（原语路线，pgr 不内置 pl trim）：
   clumpify dedupe → `fq clump --dedupe`；bbduk trim/filter → `fq clean`/
   `fq filter`；bbduk 纯 qtrim（merge）→ `fq trim-adapter --no-ktrim
   --no-tbo --no-tpe --max-ns=-1 --force-trim-mod 0 --trim-quality <qual>
   --minlen <len>`；reformat sample → `fq sample`；kmercountexact →
   `kmer hist`；repair → `fq split --repair`（含 8_spades/8_mr_spades 的
   hnsm 管道形态）；tadpole contig → `asm contig`（unitigs/2_insert_size，
   含 0_cleanup 文件名同步）；anchors bbwrap → `asm map` + `sam to-rg` +
   `rg coverage`；2_insert_size bbmap ×2 → `asm map --paired --max-reads
   {{opt.reads}}` + `sam ihist`（reformat ihist 替代，Picard 保留外部；
   完美匹配对足以估计插入长度，见 `asm-map.md` §2.6）。已知偏差：contig
   的 bubble 解析与 tadpole 有少量差异（`-Xmx` 相关），已接受确定性输出。
3. **bbnorm 深度分箱**：暂不做（§4.9）。
4. **kmer table 的 k 上限 64 vs anchr 2_fastk 的 k=81**（2026-08-11 记录，
   `design/kmer.md` §11）：anchr `2_fastk.tera.sh` 用 `FastK -t1
   -k<21|51|81>`，k=21/51 没问题，**k=81 超出 `pgr kmer table` 当前
   u128 上限 64**，属已知缺口。若替代 2_fastk 需要 k=81，给 `libs/kmer`
   扩表示（参考 FastK 字节打包或 u128 双字），当前未做；`asm
   contig`/`unitig` 走 tadpole 多字 Kmer，k=81 已可用，不受影响。

---

**历史会话（已完成，一行结论，细节见各设计笔记）**：

- **asm 命令族 + SAM 工具（2026-08-11/12 提交）**：`pgr asm`
  contig/unitig/map（含 `--min-count-seed`、`--links`/`--gfa`、
  `--paired`/`--max-reads`）、`pgr sam` ihist/to-rg（noodles-sam 0.81
  解析）；basecov 移出 map（SAM 派生）；map 流式分块 + 头对称；contig
  计数并行 + 排序快照（576→157 ms）；写出端手写（refname 全头字段与
  noodles 严格字符集冲突）→ `design/asm-map.md`、`design/fq-assemble.md`、
  `design/kmer.md` §11。
- **fq 纠错命令拆分（2026-08-11 提交）**：`fq ecc`→`ec-kmer`、
  `merge --ecco`→`ec-overlap`，golden 对照逐字节一致 →
  `design/anchr-merge-replace.md`。
- **trim 8 步 M0-M8 全部移植**（fq sample/clump/split/clean/filter/norm/
  trim-adapter/kmer hist），与 BBTools 39.38 逐字节一致，代码已提交 →
  `design/anchr-trim-replace.md`。
- **anchr merge 流程 7 步全部移植**（fq merge/ec-kmer/ec-overlap/extend/assemble +
  split --repair + s-filter），golden/统计对照完成；clumpify ecc 按计划
  跳过 → `design/anchr-merge-replace.md`、`design/fq-assemble.md`。
- **2026-08-09 `pgr kmer` 七子命令全量交付**（table/profile/hist/gc/
  qhist/qcheck/gsize，含 GenomeScope 2.0 原生迁移 `--model`、plot
  heat/spectra 拆分、`.pkt`/`.pkp`/`.hist` 三格式定稿），测试 1544 全绿
  → `design/kmer.md` §10。
- **fq 系列（2026-08）**：`trim-qual`（sickle 替代）、`range`（FASTQ
  `.loc` 索引一期 + 双端二期）、BGZF 写侧基准、FAFQ reader 笔记 →
  `design/anchr-trim-replace.md`（含 trim-qual）、`design/fq-index.md`、`design/seq-reader.md`。

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

- [x] PAF `query_length`/`target_length` 恒 0：**已修复**（2026-08-11，
      `.paf.idx` v5 持久化每序列长度，旧索引报错重建；
      `audit/audit-paf.md`）。
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
- `fq clump` 多参数 golden 验证：**不做**（体积控制，见
  `design/anchr-trim-replace.md` §4.4 M1 注）。
- **filterbytile.sh / 光学去重：不做**（2026-08-10/12）——flowcell tile
  坐标质量过滤（trim `--tile`，默认关）与光学去重同源，需坐标解析且无
  真实坐标数据可验证；`design/anchr-trim-replace.md` §3/§4 已记录。
  另注：8_spades/8_mr_spades 的 `repair.sh`（hnsm filter 管道）由
  `pgr fq split --repair` 覆盖，模板改写时验证 stdin/interleaved 形态。
- **asm contig 计数表 radix 化：已评估不做**（2026-08-11 基准：Lambda
  20k 下 radix 比 `cmp_bases` 比较排序慢，几十万 k-mer 规模不划算；
  数百万级再评估，见 `fq-assemble.md` §7）。

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
- [ ] **pgr asm unitig 真实数据验证**：已实现为独立命令（借鉴 BCALM
      graph3，顺序无关/无气泡；设计 `fq-assemble.md` §8）；待用 anchr
      `pe.cor.fa` 与 bcalm 输出对照 unitig 集合与连续性。
- [ ] **pgr asm unitig L: 边/GFA 真实对照**：`--links`/`--gfa` 已实现
      （LinkTigs 语义，方向规则见 `fq-assemble.md` §8 + 单测）；待与
      bcalm `LinkTigs` 输出（`.unitigs.fa` 的 `L:` 头）在真实数据上
      对照方向与边集合。
- [ ] **pgr asm map 真实数据验证**：已实现（`asm/map.rs`，设计
      `asm-map.md`）；待用真实 UT.fasta + reads 与 bbmap 输出对照
      mapped 比例与覆盖度（`sam to-rg` + `rg coverage`，本机 Java 配对
      读 gz 失败，黑盒对照暂缓）。
- [ ] **chain 算法待验证（低优先）**：KD-tree 已实现并用于 `psl chain`
      （`libs/ds/kdtree.rs`）；`best_crossover` 已接入 `fas_multiz` merge
      （`libs/ds/crossover.rs`）——两者的**真实数据验证**待做；KD-tree
      用于 PAF 链式化 / POA 排序仍待评估（PAF 当前未明确需要链式化）
      （来源：`chain-algorithms.md` §12.3）。
