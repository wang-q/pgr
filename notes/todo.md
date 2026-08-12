# pgr 近期待办

> 依据 `project-understanding.md`、各设计笔记与审计记录整理。
> 功能层基本齐备，近期大头是验证与数据驱动的扩展。
> 已完成条目只留一行结论，细节见链接文档。
> 按类型组织（已完成 / 待实现 / 挂账待决 / 待验证等数据 / 低风险审计 /
> 技术债 / 明确不做），不按会话轮次。

## 0. 会话交接（2026-08-12 晚，依赖清理 + 二进制瘦身 + 下一步规划）

> 会话交接材料，供下一次会话恢复上下文；读取后按用户指示清理。

**当前状态**：1741 测试通过，fmt/clippy 干净，工作树干净。
**最近提交**：`638fbf3`（依赖清理收尾：probminhash/isal-rs 删除）、
`6290d7d`（tera 移除 + regex 挪 dev-deps + pgi no-seq 警告）、`d167ac3`
（pgi A1 + 键错位修复）、`f5bb2e7`（unitig 级 contain 预过滤）、
`d9b478b`（paf coverage + OLC 修复 + 输出级去重）。

**本会话成果（依赖审计 + 瘦身）**：
- **release 17.0 → 12.62 MB**：bio（死依赖）、tera（plot 4 条渲染路径
  改纯 Rust 拼接/替换，8 个 golden 逐字节一致）、regex（生产代码 0 使用
  → dev-deps，A/B 实测 -0.48MB）、probminhash/isal-rs（零引用删除）。
- **依赖全量审计**：生产依赖全部真实使用、无死依赖；regex 无法彻底移出
  二进制（env_logger→env_filter 传递依赖，仍 ~0.61MB）；xlsx 输出链
  （rust_xlsxwriter+zip+zopfli+zlib-rs ~0.72MB）为最大第三方贡献者，属
  真实功能未动；`sqlite3-ext-vtab` 的 bundled 特性为 sqlite-vector-rs
  基准必需（保留）。
- **二进制构成（nm 符号级 9.40MB）**：pgr 3.68 / std+core+alloc 1.95 /
  异常+解栈元数据 0.91 / xlsx 链 0.72 / regex 家族 0.61 / rayon 0.37 /
  clap 0.26 / 压缩 C 0.18 / crossbeam 0.12 / anyhow 0.11 / serde_json
  0.10 / noodles 0.08 MB。

**关键裁定（必须遵守，自旧交接保留）**：
- k-mer 只保留一套，以 FASTK-master 为准（不用 tadpole `Vec<u64>`、不做
  定长对象键）——**已落地（M1–M5 完成）**；存储 = FastK 式连续打包；
  u128 仅剩算法中间量（pgi 构建滚动、qcheck 判定、FastGA 移植位运算）
- OLC 已落地（见 §1/§3）；v1 覆盖度 repeat breaking 待真实宏基因组数据调参
- 气泡处理不做（明确不做区）

**下一步（按优先级，用户已确认方向）**：
1) **paf 查询层 `--min-tree-coverage`**（`paf-pangenome.md` §Caf）：纯代码，
   查询时传递闭包后处理过滤。
2) **fas_multiz best_crossover 真实多基因组验证**（`chain-algorithms.md`
   §12.3）：合成测试已覆盖，真实数据待确认本地是否有合适多基因组比对。
3) 等数据：OLC 长读（本地无 HiFi/ONT）、4 万 E. coli cohort、人类规模
   pgi（GRCh38/CHM13）。pgi 结构性重写（收集时分桶）挂账，等真实大基因
   组场景证明必要性。

**参考源码（本地，gitignore 参考目录）**：`FASTK-master/`（长 k 第一参考）、
`FASTGA-main/`（pgi 参考，k-mer 与 FastK 同套）、`canu-2.3/`、`wgs-8.3rc2/`、
`metaMDBG-metaMDBG-1.4/`、`bcalm/`。

## 1. 已完成（一行结论，细节见链接）

- **依赖瘦身与审计**（2026-08-12）：release 17.0→12.62 MB。bio（死依赖）、
  tera（plot 模板改纯 Rust，8 个 golden 逐字节一致）、regex（生产 0 使用
  → dev-deps，-0.48MB 实测）、probminhash/isal-rs（零引用）；生产依赖全量
  审计无死依赖；`sqlite3-ext-vtab` bundled 特性为 sqlite-vector-rs 基准
  必需（保留）→ 详见 Cargo.toml 注释与 §0。
- **asm 命令族 + SAM 工具**（2026-08-11/12 提交）：`pgr asm`
  contig/unitig/map（含 `--min-count-seed`、`--links`/`--gfa`、
  `--paired`/`--max-reads`）、`pgr sam` ihist/to-rg（noodles-sam 0.81
  解析）；basecov 移出 map（SAM 派生）；map 流式分块 + 头对称；contig
  计数并行 + 排序快照（576→157 ms）；写出端手写（refname 全头字段与
  noodles 严格字符集冲突）→ `design/asm-map.md`、`design/fq-assemble.md`、
  `design/kmer.md` §11。
- **fq 纠错命令拆分**（2026-08-11 提交）：`fq ecc`→`ec-kmer`、
  `merge --ecco`→`ec-overlap`，golden 对照逐字节一致 →
  `design/anchr-merge-replace.md`。
- **trim 8 步 M0-M8 全部移植**（fq sample/clump/split/clean/filter/norm/
  trim-adapter/kmer hist），与 BBTools 39.38 逐字节一致 →
  `design/anchr-trim-replace.md`。
- **anchr merge 流程 7 步全部移植**（fq merge/ec-kmer/ec-overlap/extend/
  assemble + split --repair + s-filter），golden/统计对照完成；clumpify ecc
  按计划跳过 → `design/anchr-merge-replace.md`、`design/fq-assemble.md`。
- **`pgr kmer` 七子命令全量交付**（2026-08-09）：table/profile/hist/gc/
  qhist/qcheck/gsize，含 GenomeScope 2.0 原生迁移 `--model`、plot
  heat/spectra 拆分、`.pkt`/`.pkp`/`.hist` 三格式定稿，测试 1544 全绿
  → `design/kmer.md` §10。
- **fq 系列**（2026-08）：`trim-qual`（sickle 替代）、`range`（FASTQ
  `.loc` 索引一期 + 双端二期）、BGZF 写侧基准、FAFQ reader 笔记 →
  `design/anchr-trim-replace.md`（含 trim-qual）、`design/fq-index.md`、
  `design/seq-reader.md`。
- **spanr cover 名字截断**（2026-08-09）：pgr 内建 runlist 区间操作替代
  外部 spanr；`rept/trf.rs` 带点 contig 名映射规避 →
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
- **PAF `query_length`/`target_length` 恒 0**（2026-08-11）：已修复，
  `.paf.idx` v5 持久化每序列长度，旧索引报错重建 → `audit/audit-paf.md`。
- **k-mer 表示统一到 FastK + 长 k（k=81）落地**（2026-08-12）：
  `libs/kmer/key.rs` Kmer 字节键（FastK 字节序实测锁定）+ radix 泛化 +
  FastK golden 对照（k=21/51/81 逐条一致）；KmerTable/quality/norm/map/
  pgi/tadpole 全部迁移到打包字节（无 per-key 对象头），qcheck 查询侧
  转 Kmer；`.pkt`/`.pgi` 格式不变（字节与 FastK 一致，旧缓存直接兼容、
  不 bump）；k=81 table/hist/profile 端到端（2939 条与 FastK 一致）；
  → `design/kmer.md` §12。
- **kmer 计数性能修复**（2026-08-12）：count_mg1655 297 → 158 ms
  （双窗口滚动 + canonical 半长比较 + 收集直落字节 + `append` 移动
  缓冲 + emit 传引用）；单线程 347 ms vs FastK 481 ms（-28%）、
  8 线程 158 vs 188 ms（-16%）；基准拆解（radix 45 ms / 分组 13 ms）
  保留在 kmer_benchmark → `design/kmer.md` §12 基准段。
- **repeat masking 真核标定（2026-08-12）**：`rept e-align` 默认参数定稿
  f100/ms16（`--freq` 50→100、`--min-shared` 12→16 已落地）；S288c（去
  soft-mask）+ repbase 对 `pgr rept masker` 参考（RM 复刻）recall 67.6%、
  over-mask 0.029%；f10 明显低敏感、f100≈f500、k31≈k40≈k21、
  min-identity 0.60≡0.70；漏检 192 kb 以 <100 bp 真实 Ty LTR 碎片为主
  （种子敏感度边界，非参数问题），MITO 为参考侧 AT 富集 over-mask；
  MG1655 对照 1.29% 更接近 RM 1.06% → `design/repeat-masking.md` §2.5.1。
- **pgi 真实数据验证 + 优化**（2026-08-12）：E. coli 4 基因组
  （~20 Mb）build 817 ms（vs GIXmake 505 ms，1.62×）、峰值内存
  145 MB（迁移前 ~180 MB，-19%）；并行收集（+17%）、分组
  HashSet → 组内排序去重（-30%）、pack_kmer 大端拷贝（-47%）；
  align mg1655×sakai 738 PSL 与迁移前逐条一致（`.pgi` 兼容确认）
  → `benchmarks/bench-pgi-vs-gixmake.md`、`bench-pgi-align-vs-fastga.md`。
- **pgi build 打包字节键 + 关键 bug 修复（2026-08-12，A1）**：排序中间
  表示 u128 → FastK 打包字节（内存 -10%，wall 0.84 s vs GIXmake
  505 ms ≈ 1.66×，未达 1.2×）；collect 去重 HashSet → 位图（-30 ms）、
  partition 逐字节复制改 copy_from_slice（sort -12%）、group 键比较
  u64 双段（-5%）。**抓出并修复严重 bug**：打包键重构致 entry 键错位
  （`.pgi` 静默损坏），经基因组回验修复 + 回归测试
  `grouped_entries_match_positions`；进一步收敛需结构性改动
  （GIXmake 式收集时分桶）→ `bench-pgi-vs-gixmake.md`。

## 2. 待实现

- [ ] **paf 查询层扩展**：`--min-tree-coverage`（Caf Tree Coverage 过滤维度，
      查询时无法全图计算，作传递闭包后处理过滤）；`--end-trim` 推迟（需
      per-interval 修剪 CIGAR，待序列输出引入时一并处理）
      （来源：`paf-pangenome.md` §Caf 过滤维度对照表）。

## 3. 挂账 / 待决

- **norm 精确 vs 近似定稿**（`anchr-trim-replace.md` §4.8 未定）：pgr 走
  精确表 + 外部桶（1TB 答案）；bbnorm bits=16 近似表结果依赖 -Xmx。
  差异 = 定义差异不是 bug。
- **anchr 模板替换**（用户自己处理，命令已齐）：trim.era.sh 的 bbtools
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
- **多 k unitig 的 OLC 拼接**：**已实现**（2026-08-12，`pgr asm olc`
  等四命令，`design/olc.md`，承接 `references/canu.md` §8）。剩余待真实
  宏基因组数据验证：overlap 是否允许少量错配、不同 k unitig 冗余去重、
  repeat breaking 覆盖度证据阈值（Canu 6/15 的单元化版本）、列投票
  consensus（v1）。

## 4. 待验证 / 等数据或场景到位

- [x] 稀疏 s=1 完整 45 对 cohort 复测（2026-08-12，10 株本地数据已齐）：
  mash 排序 Spearman 0.9814 / Pearson 0.9969 / max |Δ| 0.0035
  （10 对子集为 0.988/0.0025，全 cohort 略降但一致）→ `design/hv.md`。
- [ ] 4 万 E. coli cohort 端到端：核心步骤就绪，等真实 cohort 数据
      （来源：`ecoli-cohort.md`）。
- [ ] 人类规模（GRCh38/CHM13）验证：`.pgi` 字段上限、内存/耗时与 FastGA
      对照（来源：`design/pgi-align.md` §7.2）；E. coli 多基因组
      （~20 Mb，8 contig）已实测（见 §1），人类规模仍待。
- [ ] pbit 自动路由：等多样性 cohort 数据证明收益（来源：`design/pbit.md`）。
- [x] `--sym` 场景开关：方向偏差已量化（2026-08-12，带 `--ref-seq/--query-seq`
  正确路径，A→B vs B→A 的 aligned 差异：mg1655×sakai +0.49%、
  mg1655×se11 −1.81%、sakai×se11 +0.51%、mg1655×cft073 +0.08%——
  多数 <1%，接近"噪声级"标准，全对场景大概率不需要 `--sym`，默认关
  保持不变。⚠️ 注意：**不带序列的 `align pgi` 路径输出的是未精化的
  geometric blocks**（帮助文档明示，每个种子链管一个粗 block），不能当
  精化比对用；本项首次量化误用了该路径，数字作废重测
  （来源：`design/pgi-align.md` §7.4.1）。
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
- [x] **pgr asm unitig 真实数据验证**（2026-08-12）：MG1655 1M 纠错 reads
      （`/home/wangq/data/anchr/mg1655/2_illumina/merge/pe.cor.fa.gz` 采样）
      × k=31，与 `/home/wangq/.cbp/bin/bcalm` 对照——**unitig 序列 100%
      一致**（2403/2403 canonical 归一后逐条相同）。→ `fq-assemble.md` §8.1。
- [x] **pgr asm unitig L: 边/GFA 真实对照**（2026-08-12）：无向边集 pgr
      3801 vs bcalm 3331，共同 2577（bcalm 的 77%）、pgr 多 1224/缺 754
      ——bcalm 边 = canonical 端点共享图（5019）的子集，过滤条件不在
      README/本地源码；尝试按 README 语义重写后边集未对齐且回归真实
      边，已回退。**结论：保持简化语义，逐边对齐需 bcalm 链接实现
      源码，暂不立项**（`fq-assemble.md` §8.1 已记录）。
- [x] **pgr asm map 真实数据验证（自一致性，2026-08-12）**：MG1655 参考 ×
      1M 纠错 reads（anchr `pe.cor.fa` 采样）——完美贴回 505,305/1M
      （50.5%，539,045 hits，1.5% 多映射）；**mapped reads 坐标抽查
      10/10 与参考区间逐碱基一致**；mapped 子集平均覆盖 38.4×。50%
      贴回率是设计使然（完美匹配拒绝任何带错 reads；anchr bwa 容错映射
      mosdepth 271× 佐证 reads 来自该参考）。⚠️ bbmap 黑盒对照仍暂缓
      （本机 Java 配对读 gz 失败）→ `asm-map.md` §4。
- [ ] **OLC 宏基因组/长读真实数据验证**：Lambda（Illumina 108bp）已验证
      （`design/olc.md` §12，含 40× 原始与 9× 纠错两种路径）；长读
      （HiFi/ONT）与宏基因组数据到位后调 v1 参数（overlap 错配容忍、
      repeat breaking 覆盖度阈值、多 k 反馈）。
- [x] chain 算法验证（2026-08-12）：**KD-tree 已由 UCSC 管线字节级验证
      覆盖**——`verify-ucsc-pipeline.sh:71` 跑 `pgr psl chain`（KD-tree
      链式化，`libs/chain/connect.rs`）并与 axtChain 逐字节一致（E. coli
      mg1655×sakai 主流程 + --syn + medium）。`best_crossover` 仍有
      fas_multiz 真实多基因组数据验证待做（合成测试已覆盖交叉合并语义，
      `tests/cli_fas_multiz.rs`）；KD-tree 用于 PAF 链式化 / POA 排序
      仍待评估（PAF 当前未明确需要链式化）（来源：`chain-algorithms.md`
      §12.3）。

## 5. 低风险审计记录项（可顺手修）

- [x] `pgr align pgi` 不带序列路径（geometric blocks，未精化）易被误当
      精化比对使用——已加运行时 `log::warn!` 显式警告（2026-08-12：
      "writing unrefined geometric blocks ... pass --ref-seq/--query-seq"），
      集成测试 `command_align_pgi_warns_without_extension_sequences` /
      `_no_warning_with_extension_sequences` 锚定。
- [ ] `syncmer.rs` 重复发射同一位置：**暂缓**——消费方已去重，收益小风险高
      （来源：`audit/audit-rept-sd.md`）。

## 6. 技术债（有空再议）

- [ ] `fas` 模块职责过重（20 子命令）考虑拆分（来源：`genome-nn-query.md` §8.6）。

## 7. 明确不做（避免重复立项）

- Gap_Improver、完整 LCP、`.1aln`、trace points、ALNchain、GDB/GIX 分片
  （`design/pgi-align.md` §6）；多 mask union（§7.5）；`-S` 对称 adaptamer
  （§7.4.1）；hybrid 逻辑留 `cmd_pgr/`（commit `d5281bc`，有意为之）。
- **bbnorm 深度分箱：暂不做**（`anchr-trim-replace.md` §4.9，未定、暂不实现）。
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
